import json
import numpy as np
import pandas as pd
import pytest
from ichnos.population import summarize_cells, SUMMARY_METHOD
from ichnos.snapshot_adapter import observation_from_summary, decode_summary
from test_snapshot_decoder import model, config, observation
from ichnos.snapshot_decoder import build_snapshot_calibration, decode_snapshot
from ichnos.posterior import fit_replicate_noise, decode_posterior


def summary_fixture(model):
    table, ref, ratio = model
    cfg = config(table, ref, ratio)
    cfg['acquisition'].update(nd_filter_green=0., nd_filter_red=0.)
    cfg['summary_method'] = SUMMARY_METHOD
    calibration = build_snapshot_calibration(table, **cfg)
    obs = observation(table)
    cells = []
    for i in range(4):
        g = obs['corrected_green'] * (i + 1)
        r = obs['image_ratio_red_green']
        cells.append(dict(session_id='session', sample_id='sample', specimen_id='target',
            biological_replicate_id='culture', condition_id='condition', timepoint=1, cell_id=i,
            acquisition_order=1, qc_pass=i < 3, corrected_mean_green=g,
            corrected_mean_red=g*r, ratio_red_green=r, sampling_time_hours=1.,measurement_time_hours=1.,
            exposure_ms_green=100.,exposure_ms_red=100.,nd_filter_green=0.,nd_filter_red=0.,objective='40X',
            acquisition_json=json.dumps(cfg['acquisition'])))
    frame = pd.DataFrame(cells)
    summary = summarize_cells(frame, data_kind='synthetic').iloc[0].to_dict()
    return table, calibration, frame, summary


def test_paired_green_and_ratio_use_same_qc_cells(model):
    table, calibration, frame, summary = summary_fixture(model)
    assert summary['corrected_green_median'] == frame.corrected_mean_green.iloc[1]
    obs = observation_from_summary(summary, calibration, source='cells hash')
    assert obs['corrected_green'] == summary['corrected_green_median']
    assert obs['image_ratio_red_green'] == summary['ratio_red_green_median']
    options = dict(ratio_tolerance=.01, green_tolerance=.01, green_floor=0)
    a = decode_summary(table, calibration, summary, source='cells hash', **options)
    b = decode_snapshot(table, calibration, obs, **options)
    assert a == b


@pytest.mark.parametrize('field,value', [('biological_replicate_id',None),('summary_status','insufficient_cells'),
    ('summary_method','mean'),('acquisition_json',None)])
def test_missing_or_incompatible_summary_rejected(model, field, value):
    _, calibration, _, summary = summary_fixture(model)
    summary[field] = value
    with pytest.raises(ValueError):
        observation_from_summary(summary, calibration, source='input')


def test_adapter_posterior_and_training_estimator_contract(model):
    table, calibration, _, summary = summary_fixture(model)
    obs = observation_from_summary(summary, calibration, source='input')
    rows = []
    for i, (dr,dg) in enumerate(((.98,.99),(1.,1.02),(1.03,.99),(1.01,1.03))):
        rows.append({**obs, 'specimen_id':f'train-{i}','biological_replicate_id':f'train-{i}',
            'condition_id':'train', 'dose_uM':25.,'time_hours':1.,
            'image_ratio_red_green':obs['image_ratio_red_green']*dr,
            'corrected_green':obs['corrected_green']*dg})
    noise = fit_replicate_noise(table,calibration,rows,green_floor=0)
    options = dict(prior={'kind':'uniform_grid'},max_mahalanobis_squared=25)
    result = decode_summary(table,calibration,summary,source='input',noise=noise,**options)
    direct = decode_posterior(table,calibration,noise,obs,**options)
    assert result == direct
    rows[0].pop('summary_method')
    incompatible = fit_replicate_noise(table,calibration,rows,green_floor=0)
    with pytest.raises(ValueError,match='summary_method'):
        decode_summary(table,calibration,summary,source='input',noise=incompatible,**options)
