import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd
from PIL import Image
import pytest
from tests.model.test_snapshot_decoder import model, config
from ichnos.snapshot_decoder import build_snapshot_calibration
from ichnos.posterior import fit_replicate_noise
from ichnos.population import SUMMARY_METHOD
from ichnos.artifacts import write_json, sha256_file
from ichnos.workflow import run_workflow


def workflow_inputs(model, tmp_path):
    table, reference, ratio = model
    cfg = config(table, reference, ratio)
    cfg.update(summary_method=SUMMARY_METHOD)
    cfg['acquisition'].update(extraction_green='scalar',extraction_red='scalar',
        nd_filter_green=0.,nd_filter_red=0.)
    yy,xx = np.mgrid[:128,:128]
    mask = np.zeros((128,128),bool)
    for cy,cx in ((32,32),(64,96),(96,32)):
        mask |= (yy-cy)**2+(xx-cx)**2 < 100
    g = 20+100*table.green[1,1]
    r = 1.5*table.ratio_table.ratios[1,1]
    for name,array in [('green',np.where(mask,g,0.)),('red',np.where(mask,g*r,0.)),
                       ('bf',np.where(mask,.2,.5))]:
        Image.fromarray(array.astype(np.float32)).save(tmp_path/f'{name}.tif')
    # References undergo the same segmentation/summary estimator as targets.
    # Mean image intensity is affected by the actual measured mask area.
    from ichnos_image.pipeline import ImageSet, process_image_set
    for anchor, value in ((cfg['low_reference'],20+100*table.green[0,0]),
                          (cfg['high_reference'],20+100*table.green[-1,-1])):
        plane = np.where(mask,value,0.).astype(np.float32).astype(float)
        bright = np.where(mask,.2,.5).astype(np.float32).astype(float)
        records = process_image_set(ImageSet(plane,plane,'reference',0,0,100,100,0,0,'40X',1,30,
            bright_field=bright),bleed_green_to_red=0.)
        anchor['corrected_green'] = float(np.median([x.corrected_mean_green for x in records]))
    calibration = build_snapshot_calibration(table, **cfg)
    for name, data in [('table',table.to_dict()),('calibration',calibration)]:
        write_json(tmp_path/f'{name}.json',data)
    ilow,ihigh = calibration['image_green_anchors']
    mlow,mhigh = calibration['model_green_anchors']
    measured_g = ilow+(table.green[1,1]-mlow)/(mhigh-mlow)*(ihigh-ilow)
    manifest = []
    for sample, warmup in [('accepted',30.),('rejected',0.)]:
        manifest.append(dict(session_id='session',sample_id=sample,specimen_id=sample,
            biological_replicate_id=sample+'-culture',condition_id='target',timepoint=0,
            acquisition_order=len(manifest),green_path='green.tif',red_path='red.tif',bright_field_path='bf.tif',
            exposure_ms_green=100.,exposure_ms_red=100.,nd_filter_green=0.,nd_filter_red=0.,
            objective='40X',burner_hours=1.,lamp_warmup_minutes=warmup,gain_setting='fixed synthetic setting'))
    pd.DataFrame(manifest).to_csv(tmp_path/'images.csv',index=False)
    training = []
    for i,(dr,dg) in enumerate(((.98,.99),(1,1.02),(1.03,.99),(1.01,1.03))):
        training.append(dict(session_id='session',specimen_id=f'train-{i}',biological_replicate_id=f'train-{i}',
            data_kind='synthetic',ratio_scale='image',acquisition=cfg['acquisition'],
            observation_unit='biological_replicate_summary',summary_method=SUMMARY_METHOD,
            image_ratio_red_green=float(r*dr),corrected_green=float(measured_g*dg),condition_id='training',dose_uM=25.,time_hours=1.))
    noise = fit_replicate_noise(table,calibration,training,green_floor=0)
    write_json(tmp_path/'noise.json',noise)
    write_json(tmp_path/'prior.json',{'kind':'uniform_grid'})
    configuration = dict(image_manifest='images.csv',bleed_default=0.,data_kind='synthetic',
        table='table.json',calibration='calibration.json',noise='noise.json',prior='prior.json',
        max_mahalanobis_squared=25.,image_options={'saturation_value':65535.})
    write_json(tmp_path/'config.json',configuration)
    return tmp_path/'config.json'


def test_real_cli_images_to_posterior_and_rejection_report(model,tmp_path):
    config_path = workflow_inputs(model,tmp_path)
    output = tmp_path/'run'
    cmd = [sys.executable,str(Path(__file__).resolve().parents[2]/'scripts/run_workflow.py'),
           '--config',str(config_path),'--out-dir',str(output)]
    result = subprocess.run(cmd,capture_output=True,text=True)
    assert result.returncode == 0,result.stderr
    report = json.loads((output/'qc-report.json').read_text())
    assert report['n_cells'] == 6 and report['n_cells_qc_pass'] == 3
    accepted,rejected = report['samples']
    assert accepted['status'] == 'posterior_available'
    assert rejected['status'] == 'rejected' and 'insufficient_cells' in rejected['reason']
    posterior = json.loads((output/accepted['result_file']).read_text())
    assert posterior['posterior_computed']
    assert np.sum(posterior['posterior_mass']) == pytest.approx(1.)
    manifest = json.loads((output/'run-manifest.json').read_text())
    assert manifest['outputs']['cells.csv'] == sha256_file(output/'cells.csv')
    assert 'source_tree_sha256' in manifest['code']
    assert (output/'report.html').exists()
    before = sha256_file(output/'qc-report.json')
    assert subprocess.run(cmd,capture_output=True).returncode != 0
    assert sha256_file(output/'qc-report.json') == before


def test_experimental_requires_enforced_focus_and_controls(model,tmp_path):
    path = workflow_inputs(model,tmp_path)
    cfg = json.loads(path.read_text());cfg['data_kind']='experimental';write_json(path,cfg)
    with pytest.raises(ValueError,match='focus'):
        run_workflow(path,tmp_path/'run')
    assert not (tmp_path/'run').exists()
