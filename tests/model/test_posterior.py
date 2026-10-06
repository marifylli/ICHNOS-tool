from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy.special import logsumexp
from ichnos.posterior import fit_replicate_noise, decode_posterior, prior_grid
from ichnos.snapshot_decoder import build_snapshot_calibration, prepare_snapshot, save_artifact
from tests.model.test_snapshot_decoder import model, controlled, observation, config

ROOT=Path(__file__).resolve().parents[2]


def rows(table,calibration):
    deviations=np.array([[-.03,-.02],[-.02,.01],[-.01,-.03],[.01,.03],[.02,-.01],[.03,.02]])
    result=[]
    for group,(i,j) in enumerate([(0,0),(2,2)]):
        for k,noise in enumerate(deviations):
            row=observation(table,i,j)
            row.update(observation_unit='biological_replicate_summary',condition_id=f'c{group}',
                dose_uM=float(table.ratio_table.doses_uM[i]),time_hours=float(table.ratio_table.times_hours[j]),
                specimen_id=f'train-{group}-{k}',biological_replicate_id=f'culture-{group}-{k}')
            row['image_ratio_red_green']*=float(np.exp(noise[0]))
            low,high=calibration['image_green_anchors']
            row['corrected_green']+=float(noise[1]*(high-low))
            result.append(row)
    return result


def setup(model):
    table,calibration=controlled(model)
    training=rows(table,calibration)
    noise=fit_replicate_noise(table,calibration,training,green_floor=0)
    target=observation(table)
    target['observation_unit']='biological_replicate_summary'
    return table,calibration,training,noise,target


def posterior(table,calibration,noise,target,**kwargs):
    options=dict(prior={'kind':'uniform_grid'},max_mahalanobis_squared=25.)
    options.update(kwargs)
    return decode_posterior(table,calibration,noise,target,**options)


def test_covariance_is_within_condition_and_preserves_correlation(model):
    table,calibration,training,noise,target=setup(model)
    latent=[]
    for row in training[:6]:
        prepared,_=prepare_snapshot(table,calibration,row,green_floor=0)
        latent.append([np.log(prepared['ratio_red_green']),prepared['green_norm']])
    expected=np.cov(np.array(latent).T,ddof=1)
    np.testing.assert_allclose(noise['covariance'],expected,rtol=1e-12,atol=1e-15)
    assert noise['residual_degrees_of_freedom']==10
    assert noise['covariance'][0][1]!=0
    assert noise['observation_unit']=='biological_replicate_summary'


def test_posterior_agrees_with_direct_gaussian_and_normalizes(model):
    table,calibration,training,noise,target=setup(model)
    result=posterior(table,calibration,noise,target)
    prepared,green=prepare_snapshot(table,calibration,target,green_floor=0)
    delta=np.stack([np.log(table.ratio_table.ratios)-np.log(prepared['ratio_red_green']),green-prepared['green_norm']],axis=-1)
    logp=-.5*np.einsum('...i,ij,...j->...',delta,np.linalg.inv(noise['covariance']),delta)
    expected=np.exp(logp-logsumexp(logp))
    np.testing.assert_allclose(result['posterior_mass'],expected,rtol=1e-10,atol=1e-14)
    assert np.sum(result['posterior_mass'])==pytest.approx(1)
    assert result['map_grid_pairs'][0]['dose_uM']==25
    assert result['map_grid_pairs'][0]['elapsed_time_hours']==1
    assert result['credible_mass_achieved']>=.95
    assert result['status']=='posterior_available'
    assert result['covariance_uncertainty_propagated'] is False


def test_density_prior_uses_irregular_grid_cell_sizes(model):
    table,_=controlled(model)
    actual=prior_grid(table,dict(kind='uniform_density',dose_bounds_uM=[0,75],time_bounds_hours=[0,8]))
    expected=np.outer([12.5,37.5,25],[.75,3.75,3.5])/(75*8)
    np.testing.assert_allclose(actual,expected)
    assert not np.allclose(actual,prior_grid(table,{'kind':'uniform_grid'}))


def test_credible_set_keeps_equal_modes_and_cutoff_ties(model):
    table,_=controlled(model)
    table.green[0,1]=2
    calibration=build_snapshot_calibration(table,**config(table,model[1],model[2]))
    noise=fit_replicate_noise(table,calibration,rows(table,calibration),green_floor=0)
    target=observation(table)
    target['observation_unit']='biological_replicate_summary'
    result=posterior(table,calibration,noise,target,credible_mass=.2)
    assert len(result['map_grid_pairs'])==2
    assert len(result['credible_set'])==2
    assert result['credible_mass_achieved']==pytest.approx(1)
    masses=np.zeros((3,3)); masses[0,1]=.2; masses[1,1]=.8
    result=posterior(table,calibration,noise,target,prior=dict(kind='explicit_mass',masses=masses.tolist()))
    assert result['status']=='prior_only'
    np.testing.assert_allclose(result['posterior_mass'],masses)


@pytest.mark.parametrize('change',['cell_unit','duplicate_culture','mixed_condition','few_replicates','zero_variance'])
def test_noise_rejects_pseudoreplication_and_invalid_fits(model,change):
    table,calibration,training,noise,target=setup(model)
    if change=='cell_unit': training[0]['observation_unit']='cell'
    elif change=='duplicate_culture': training[-1]['biological_replicate_id']=training[0]['biological_replicate_id']
    elif change=='mixed_condition': training[0]['dose_uM']=25
    elif change=='few_replicates': training=training[:2]
    elif change=='zero_variance':
        for row in training: row['corrected_green']=220.
    with pytest.raises(ValueError): fit_replicate_noise(table,calibration,training,green_floor=0)


def test_target_must_not_be_used_to_fit_noise(model):
    table,calibration,training,noise,target=setup(model)
    with pytest.raises(ValueError,match='overlaps noise-fitting'):
        posterior(table,calibration,noise,training[0])


def test_bad_model_fit_is_flagged_even_with_normalized_posterior(model):
    table,calibration,training,noise,target=setup(model)
    target['image_ratio_red_green']=20.
    result=posterior(table,calibration,noise,target)
    assert result['status']=='poor_model_fit'
    assert result['model_fit_check_passed'] is False
    assert np.sum(result['posterior_mass'])==pytest.approx(1)
    assert result['model_fit_check_is_a_calibrated_p_value'] is False


def test_negative_normalized_green_is_allowed_but_log_zero_ratio_is_not(model):
    table,calibration,training,noise,target=setup(model)
    target['corrected_green']=100.
    result=posterior(table,calibration,noise,target)
    assert result['green_norm']<0
    target['image_ratio_red_green']=0
    with pytest.raises(ValueError,match='positive'):
        posterior(table,calibration,noise,target)
    target['corrected_green']=0
    result=posterior(table,calibration,noise,target)
    assert result['status']=='below_detection_floor'
    assert result['posterior_computed'] is False


@pytest.mark.parametrize('prior',[
    {},dict(kind='explicit_mass',masses=[[1,2]]),dict(kind='explicit_mass',masses=np.zeros((3,3)).tolist()),
    dict(kind='uniform_density',dose_bounds_uM=[1,75],time_bounds_hours=[0,8]),
])
def test_prior_requires_valid_explicit_choice(model,prior):
    table,_=controlled(model)
    with pytest.raises(ValueError): prior_grid(table,prior)


def test_noise_artifact_and_cli_roundtrip(model,tmp_path):
    table,calibration,training,noise,target=setup(model)
    changed=deepcopy(noise); changed['covariance'][0][0]*=2
    with pytest.raises(ValueError,match='modified'): posterior(table,calibration,changed,target)
    for name,data in [('table',table.to_dict()),('calibration',calibration),('replicates',training),('target',target),('prior',{'kind':'uniform_grid'})]:
        save_artifact(data,tmp_path/f'{name}.json')
    base=[sys.executable,str(ROOT/'scripts/run_posterior.py')]
    common=['--table',str(tmp_path/'table.json'),'--calibration',str(tmp_path/'calibration.json')]
    run=subprocess.run(base+['fit-noise']+common+['--replicates',str(tmp_path/'replicates.json'),'--green-floor','0','--out',str(tmp_path/'noise.json')],capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    command=base+['decode']+common+['--noise',str(tmp_path/'noise.json'),'--observation',str(tmp_path/'target.json'),
        '--prior',str(tmp_path/'prior.json'),'--max-mahalanobis-squared','25','--out',str(tmp_path/'posterior.json')]
    run=subprocess.run(command,capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    result=json.loads((tmp_path/'posterior.json').read_text())
    assert result['posterior_computed'] is True
    assert len(result['observation_file_sha256'])==64
    assert subprocess.run(command,capture_output=True).returncode!=0


def test_synthetic_posterior_verifier_separates_training_and_holdout(tmp_path):
    snapshot=tmp_path/'snapshot'
    run=subprocess.run([sys.executable,str(ROOT/'scripts/verify_snapshot_decoder.py'),'--out-dir',str(snapshot)],capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    output=tmp_path/'posterior'
    command=[sys.executable,str(ROOT/'scripts/verify_posterior.py'),'--snapshot-dir',str(snapshot),
             '--out-dir',str(output),'--heldout-count','4']
    run=subprocess.run(command,capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    report=json.loads((output/'verification.json').read_text())
    assert report['computational_checks_passed'] is True
    assert report['n_training_cultures']==36
    assert report['n_heldout']==4
    assert 0<=report['heldout_inclusion_fraction']<=1
    training=json.loads((output/'replicates.json').read_text())
    heldout=json.loads((output/'heldout-observations.json').read_text())
    assert set(row['biological_replicate_id'] for row in training).isdisjoint(row['biological_replicate_id'] for row in heldout)
    assert (output/'posterior.png').is_file()
    assert subprocess.run(command,capture_output=True).returncode!=0
