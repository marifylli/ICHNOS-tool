from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from ichnos.calibration import fit_session_calibration, model_reference_id
from ichnos.decoder import DoseTable, build_dose_table
from ichnos.snapshot_decoder import (
    ObservableTable, build_observable_table, build_snapshot_calibration,
    decode_snapshot, load_observable_table, save_artifact,
)

ROOT=Path(__file__).resolve().parents[2]
ACQUISITION=dict(objective='40X',green_units='corrected_image_units',
                 extraction_green='green_component',extraction_red='red_component',
                 gain_setting='fixed synthetic setting',exposure_ms_green=100,exposure_ms_red=100)


@pytest.fixture(scope='module')
def model(tmp_path_factory):
    directory=tmp_path_factory.mktemp('snapshot-reference')/'reference'
    process=subprocess.run([sys.executable,str(ROOT/'scripts/run_protocol.py'),
        '--variant','ox','--dose','25','--dose-units','uM','--times-hours','.5','1','8',
        '--initialization','equilibrium','--out-dir',str(directory)],capture_output=True,text=True)
    assert process.returncode==0,process.stderr
    reference=json.loads((directory/'metadata.json').read_text())
    table=build_observable_table(variant='ox',doses_uM=[0,25,75],times_hours=[.5,1,8],initialization='equilibrium')
    import pandas as pd
    ratios=pd.read_csv(directory/'fluorescence.csv').ratio_red_green.to_numpy()
    calibration=fit_session_calibration(ratios,1.5*ratios,session_id='session',
        reference_id=model_reference_id(reference),source_kind='synthetic',source='separate synthetic reference')
    return table,reference,calibration


def config(table,reference,ratio):
    return dict(session_id='session',source_kind='synthetic',acquisition=deepcopy(ACQUISITION),
        low_reference=dict(specimen_id='low',biological_replicate_id='low-culture',source='synthetic low reference',
                           dose_uM=0,time_hours=.5,corrected_green=float(20+100*table.green[0,0])),
        high_reference=dict(specimen_id='high',biological_replicate_id='high-culture',source='synthetic high reference',
                            dose_uM=75,time_hours=8,corrected_green=float(20+100*table.green[-1,-1])),
        min_image_span=1.,min_model_span=.01,ratio_calibration=ratio.to_dict(),
        ratio_reference_metadata=reference,ratio_reference_specimen_ids=['ratio-reference'],
        ratio_reference_replicate_ids=['ratio-culture'])


def observation(table,i=1,j=1):
    return dict(session_id='session',specimen_id='target',biological_replicate_id='target-culture',
        data_kind='synthetic',acquisition=deepcopy(ACQUISITION),ratio_scale='image',
        image_ratio_red_green=float(1.5*table.ratio_table.ratios[i,j]),
        corrected_green=float(20+100*table.green[i,j]))


def controlled(model):
    real,reference,ratio=model
    table=ObservableTable(DoseTable([0,25,75],[.5,1,8],
        [[.2,.4,.6],[.2,.4,.6],[.3,.5,.7]],deepcopy(real.ratio_table.provenance)),
        [[1,5,3],[4,2,6],[7,2,9]])
    return table,build_snapshot_calibration(table,**config(table,reference,ratio))


def decode(table,calibration,snapshot=None,**kwargs):
    options=dict(ratio_tolerance=1e-6,green_tolerance=1e-6,green_floor=0.)
    options.update(kwargs)
    return decode_snapshot(table,calibration,snapshot or observation(table),**options)


def test_each_observable_is_ambiguous_but_pair_resolves_snapshot(model):
    table,calibration=controlled(model)
    result=decode(table,calibration)
    assert result['ratio_only_candidate_count']==2
    assert result['green_only_candidate_count']==2
    assert result['status']=='unique_grid_pair'
    assert result['dose_estimate_uM']==25
    assert result['elapsed_time_estimate_hours']==1
    assert result['posterior_computed'] is False
    assert result['experimentally_validated'] is False


def test_two_observables_can_remain_ambiguous(model):
    table,_=controlled(model)
    table.green[0,1]=2
    calibration=build_snapshot_calibration(table,**config(table,model[1],model[2]))
    result=decode(table,calibration)
    assert len(result['candidate_pairs'])==2
    assert result['status']=='ambiguous'
    assert result['dose_estimate_uM'] is None
    assert result['elapsed_time_estimate_hours'] is None


def test_reference_normalization_and_no_clipping(model):
    table,calibration=controlled(model)
    low=decode(table,calibration,observation(table,0,0))
    high=decode(table,calibration,observation(table,2,2))
    assert low['green_norm']==0
    assert high['green_norm']==1
    snapshot=observation(table)
    snapshot['corrected_green']=100.
    result=decode(table,calibration,snapshot)
    assert result['green_norm']<0
    assert result['green_norm_outside_reference_interval'] is True
    assert result['green_was_clipped'] is False


@pytest.mark.parametrize('field,value',[
    ('specimen_id','low'),('specimen_id','ratio-reference'),
    ('biological_replicate_id','high-culture'),('biological_replicate_id','ratio-culture'),
    ('session_id','wrong'),('data_kind','experimental'),('ratio_scale','model'),
    ('corrected_green',float('nan')),('image_ratio_red_green',-1),
])
def test_independence_scale_session_and_values(model,field,value):
    table,calibration=controlled(model)
    snapshot=observation(table)
    snapshot[field]=value
    with pytest.raises(ValueError): decode(table,calibration,snapshot)


def test_acquisition_change_and_modified_calibration(model):
    table,calibration=controlled(model)
    snapshot=observation(table)
    snapshot['acquisition']['exposure_ms_green']=200
    with pytest.raises(ValueError,match='acquisition'): decode(table,calibration,snapshot)
    modified=deepcopy(calibration)
    modified['image_green_anchors'][0]+=1
    with pytest.raises(ValueError,match='modified'): decode(table,modified)


@pytest.mark.parametrize('minimum,image_minimum',[(100.,1.),(.01,1e9)])
def test_nearly_flat_or_insufficient_reference_span(model,minimum,image_minimum):
    table,ref,ratio=model
    settings=config(table,ref,ratio)
    settings.update(min_model_span=minimum,min_image_span=image_minimum)
    with pytest.raises(ValueError,match='span'): build_snapshot_calibration(table,**settings)


def test_mismatched_reference_model_rejected(model):
    table,ref,ratio=model
    settings=config(table,ref,ratio)
    settings['ratio_reference_metadata']=deepcopy(ref)
    settings['ratio_reference_metadata']['fluorescence']['f']*=2
    with pytest.raises(ValueError,match='identity'): build_snapshot_calibration(table,**settings)


def test_abstains_for_floor_and_no_joint_match(model):
    table,calibration=controlled(model)
    assert decode(table,calibration,green_floor=1000)['status']=='below_detection_floor'
    snapshot=observation(table)
    snapshot['image_ratio_red_green']=.7*1.5
    assert decode(table,calibration,snapshot)['status']=='no_grid_pair_match'
    snapshot['image_ratio_red_green']=10
    assert decode(table,calibration,snapshot)['status']=='out_of_response_domain'


def test_forward_builder_matches_legacy_ratio_table(model):
    table,_,_=model
    legacy=build_dose_table(variant='ox',doses_uM=[0,25,75],times_hours=[.5,1,8],initialization='equilibrium')
    np.testing.assert_allclose(table.ratio_table.ratios,legacy.ratios,rtol=0,atol=0)


def test_fresh_model_simulation_recovers_known_snapshot(model):
    table,ref,ratio=model
    calibration=build_snapshot_calibration(table,**config(table,ref,ratio))
    independent=build_observable_table(variant='ox',doses_uM=[25,100],times_hours=[1],initialization='equilibrium')
    snapshot=observation(independent,0,0)
    result=decode(table,calibration,snapshot,ratio_tolerance=1e-7,green_tolerance=1e-7)
    assert result['status']=='unique_grid_pair'
    assert result['dose_estimate_uM']==25
    assert result['elapsed_time_estimate_hours']==1


def test_artifacts_and_cli_roundtrip_no_overwrite(model,tmp_path):
    table,ref,ratio=model
    save_artifact(table.to_dict(),tmp_path/'table.json')
    restored=load_observable_table(tmp_path/'table.json')
    np.testing.assert_equal(restored.green,table.green)
    settings=config(table,ref,ratio)
    save_artifact(settings,tmp_path/'config.json')
    base=[sys.executable,str(ROOT/'scripts/run_snapshot_decoder.py')]
    run=subprocess.run(base+['calibrate','--table',str(tmp_path/'table.json'),'--config',str(tmp_path/'config.json'),
        '--out',str(tmp_path/'calibration.json')],capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    save_artifact(observation(table),tmp_path/'observation.json')
    command=base+['decode','--table',str(tmp_path/'table.json'),'--calibration',str(tmp_path/'calibration.json'),
        '--observation',str(tmp_path/'observation.json'),'--ratio-tolerance','0.000001',
        '--green-tolerance','0.000001','--green-floor','0','--out',str(tmp_path/'decoded.json')]
    run=subprocess.run(command,capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    result=json.loads((tmp_path/'decoded.json').read_text())
    assert result['status']=='unique_grid_pair'
    assert len(result['observation_file_sha256'])==64
    assert subprocess.run(command,capture_output=True).returncode!=0
    modified=json.loads((tmp_path/'table.json').read_text())
    modified['green'][0][0]+=1
    (tmp_path/'table.json').write_text(json.dumps(modified))
    with pytest.raises(ValueError,match='modified'): load_observable_table(tmp_path/'table.json')


def test_image_verification_cli_adds_green_information(tmp_path):
    output=tmp_path/'verification'
    command=[sys.executable,str(ROOT/'scripts/verify_snapshot_decoder.py'),'--out-dir',str(output)]
    run=subprocess.run(command,capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    report=json.loads((output/'verification.json').read_text())
    assert report['passed'] is True
    assert report['ratio_only_candidate_count']>1
    assert report['joint_candidate_count']==1
    assert report['actual']==dict(dose_uM=400.,elapsed_time_hours=8.)
    assert report['posterior_computed'] is False
    with np.load(output/'target.npz') as target, np.load(output/'green-low.npz') as low:
        assert target['green'].max()>low['green'].max()
    assert subprocess.run(command,capture_output=True).returncode!=0


@pytest.mark.parametrize('change',[
    dict(ratio_reference_specimen_ids=[]),dict(ratio_reference_replicate_ids=['x','x']),
    dict(source_kind='experimental_reference'),dict(min_image_span=0),
])
def test_reference_contract_cannot_be_omitted(model,change):
    table,reference,ratio=model
    settings=config(table,reference,ratio)
    settings.update(change)
    with pytest.raises(ValueError):
        build_snapshot_calibration(table,**settings)
