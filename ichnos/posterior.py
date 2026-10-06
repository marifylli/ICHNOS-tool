"""Conditional grid posterior with correlated biological-replicate noise.

Observation coordinates are (log calibrated ratio, normalized green). This
plug-in Gaussian model does not propagate calibration/covariance uncertainty.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from .calibration import _nonempty
from .decoder import _digest, _vector
from .snapshot_decoder import _check_digest, _number, prepare_snapshot, validate_snapshot_calibration

TRANSFORM=['log_ratio_red_green','green_norm']
UNIT='biological_replicate_summary'


def _point(table,dose,time):
    dose=_number(dose,'training dose')
    time=_number(time,'training time')
    i=np.flatnonzero(np.isclose(table.ratio_table.doses_uM,dose,rtol=0,atol=1e-12))
    j=np.flatnonzero(np.isclose(table.ratio_table.times_hours,time,rtol=0,atol=1e-12))
    if len(i)!=1 or len(j)!=1:
        raise ValueError('training condition must identify one model grid point')
    return int(i[0]),int(j[0])


def fit_replicate_noise(table,calibration,replicates,*,green_floor,min_replicates_per_condition=3):
    """Pool within-condition sample covariance, one independent culture per row.

    Culture IDs cannot recur, including across times/conditions. Group means
    remove between-condition signal from variance estimation, not model bias
    from inference. At least three replicates per condition is a numerical
    minimum, not evidence of sufficient biological sample size.
    """
    validate_snapshot_calibration(table,calibration)
    floor=_number(green_floor,'green_floor')
    minimum=min_replicates_per_condition
    if isinstance(minimum,bool) or not isinstance(minimum,int) or minimum<3:
        raise ValueError('min_replicates_per_condition must be an integer >= 3')
    if not isinstance(replicates,list) or not replicates:
        raise ValueError('replicates must be a nonempty list of biological summaries')
    groups={}
    cultures=set()
    specimens=set()
    coordinates={}
    for row in replicates:
        if row.get('observation_unit')!=UNIT:
            raise ValueError('noise rows must be biological_replicate_summary, not cells or fields')
        condition=_nonempty(row.get('condition_id'),'condition_id')
        prepared,predicted_green=prepare_snapshot(table,calibration,row,green_floor=floor)
        if prepared['status']=='below_detection_floor':
            raise ValueError('noise replicate is below the declared detection floor')
        ratio=prepared['ratio_red_green']
        if ratio<=0:
            raise ValueError('log-ratio noise requires strictly positive ratios')
        culture=row['biological_replicate_id']
        specimen=row['specimen_id']
        if culture in cultures or specimen in specimens:
            raise ValueError('replicate/specimen IDs must be independent and unique across all rows')
        cultures.add(culture)
        specimens.add(specimen)
        point=_point(table,row['dose_uM'],row['time_hours'])
        if condition in coordinates and coordinates[condition]!=point:
            raise ValueError('a condition mixes different dose/time points')
        if condition not in coordinates and point in coordinates.values():
            raise ValueError('one dose/time condition must not be split into arbitrary groups')
        coordinates[condition]=point
        groups.setdefault(condition,[]).append([float(np.log(ratio)),prepared['green_norm']])
    residuals=[]
    summaries={}
    for condition,values in groups.items():
        if len(values)<minimum:
            raise ValueError('insufficient independent replicates in a condition')
        values=np.asarray(values)
        mean=values.mean(axis=0)
        residuals.extend(values-mean)
        i,j=coordinates[condition]
        predicted_ratio=table.ratio_table.ratios[i,j]
        if predicted_ratio<=0:
            raise ValueError('training model ratio must be strictly positive')
        model_mean=[float(np.log(predicted_ratio)),float(predicted_green[i,j])]
        summaries[condition]=dict(n_replicates=len(values),grid_index=[i,j],
                                 mean=mean.tolist(),model_mean=model_mean,
                                 mean_minus_model=(mean-model_mean).tolist())
    residuals=np.asarray(residuals)
    degrees=len(replicates)-len(groups)
    covariance=residuals.T@residuals/degrees
    if not np.isfinite(covariance).all():
        raise ValueError('replicate covariance is non-finite')
    try:
        np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError('replicate covariance must be positive definite; no automatic regularization') from exc
    if np.linalg.cond(covariance)>1e12:
        raise ValueError('replicate covariance is numerically ill-conditioned')
    data=dict(schema_version=1,kind='replicate_noise',
        table_artifact_id=table.to_dict()['artifact_id'],calibration_artifact_id=calibration['artifact_id'],
        transform=TRANSFORM,observation_unit=UNIT,
        method='pooled within-condition sample covariance; fixed across grid',
        n_replicates=len(replicates),n_conditions=len(groups),residual_degrees_of_freedom=degrees,
        min_replicates_per_condition=minimum,green_floor=floor,covariance=covariance.tolist(),
        condition_summaries=summaries,biological_replicate_ids=sorted(cultures),specimen_ids=sorted(specimens),
        source_rows=replicates,source_kind=calibration['source_kind'],
        calibration_uncertainty_propagated=False,covariance_uncertainty_propagated=False,
        experimentally_validated=False)
    data=json.loads(json.dumps(data,allow_nan=False))
    return {**data,'artifact_id':_digest(data)}


def validate_noise(table,calibration,noise):
    _check_digest(noise,'replicate_noise')
    expected=fit_replicate_noise(table,calibration,noise['source_rows'],green_floor=noise['green_floor'],
                                min_replicates_per_condition=noise['min_replicates_per_condition'])
    for key,value in expected.items():
        if key in {'artifact_id','covariance','condition_summaries'}:
            continue
        if noise.get(key)!=value:
            raise ValueError(f'inconsistent noise metadata: {key}')
    if not np.allclose(noise['covariance'],expected['covariance'],rtol=1e-10,atol=0):
        raise ValueError('covariance differs from the recorded replicate fit')
    if noise['condition_summaries'].keys()!=expected['condition_summaries'].keys():
        raise ValueError('condition summaries differ from replicate groups')
    for group,summary in expected['condition_summaries'].items():
        for key,value in summary.items():
            if not np.allclose(noise['condition_summaries'][group][key],value,rtol=1e-10,atol=1e-12):
                raise ValueError('condition summary differs from recorded replicate fit')


def _widths(axis,bounds,name):
    bounds=_vector(bounds,name)
    if len(bounds)!=2 or bounds[1]<=bounds[0] or bounds[0]>axis[0] or bounds[1]<axis[-1]:
        raise ValueError(f'{name} must be increasing and enclose all grid nodes')
    edges=np.concatenate([[bounds[0]],(axis[:-1]+axis[1:])/2,[bounds[1]]])
    widths=np.diff(edges)
    if (widths<=0).any() or not np.isfinite(widths).all():
        raise ValueError(f'{name} gives non-positive or non-finite grid cell widths')
    return widths


def prior_grid(table,prior):
    """Return explicitly chosen node masses; never assume uniform silently."""
    shape=table.ratio_table.ratios.shape
    kind=prior.get('kind')
    if kind=='uniform_grid':
        masses=np.ones(shape)
    elif kind=='uniform_density':
        dose_widths=_widths(table.ratio_table.doses_uM,prior['dose_bounds_uM'],'dose_bounds_uM')
        time_widths=_widths(table.ratio_table.times_hours,prior['time_bounds_hours'],'time_bounds_hours')
        masses=dose_widths[:,None]*time_widths[None,:]
    elif kind=='explicit_mass':
        raw=np.asarray(prior['masses'],dtype=object)
        if raw.shape!=shape:
            raise ValueError('prior masses must match dose/time table shape')
        masses=_vector(raw.ravel(),'prior masses').reshape(shape)
    else:
        raise ValueError('prior kind must be uniform_grid, uniform_density or explicit_mass')
    if not np.isfinite(masses).all() or masses.max()<=0:
        raise ValueError('prior must have finite positive total mass')
    masses=masses/masses.max()
    return masses/masses.sum()


def decode_posterior(table,calibration,noise,observation,*,prior,max_mahalanobis_squared,credible_mass=.95):
    """Conditional posterior for one culture summary; retains every grid mode."""
    validate_noise(table,calibration,noise)
    if observation.get('observation_unit')!=UNIT:
        raise ValueError('target must be a biological_replicate_summary')
    if observation.get('biological_replicate_id') in noise['biological_replicate_ids'] or observation.get('specimen_id') in noise['specimen_ids']:
        raise ValueError('target overlaps noise-fitting biological replicate/specimen IDs')
    threshold=_number(max_mahalanobis_squared,'max_mahalanobis_squared',positive=True)
    level=_number(credible_mass,'credible_mass',positive=True)
    if level>=1:
        raise ValueError('credible_mass must lie strictly between 0 and 1')
    masses=prior_grid(table,prior)
    prepared,predicted_green=prepare_snapshot(table,calibration,observation,green_floor=noise['green_floor'])
    retained=('table_artifact_id','calibration_artifact_id','observation','ratio_red_green','green_norm',
              'green_norm_outside_reference_interval','green_was_clipped','green_floor_image_units','reference_independence')
    result={key:prepared[key] for key in retained}
    result.update(schema_version=1,kind='snapshot_posterior',noise_artifact_id=noise['artifact_id'],
        prior=json.loads(json.dumps(prior,allow_nan=False)),prior_mass=masses.tolist(),
        transform=TRANSFORM,credible_mass_requested=level,max_mahalanobis_squared=threshold,
        posterior_computed=False,posterior_mass=None,credible_set=[],map_grid_pairs=[],
        experimentally_validated=False,calibration_uncertainty_propagated=False,
        covariance_uncertainty_propagated=False,model_parameter_uncertainty_propagated=False,
        scope='conditional discrete posterior with fixed fitted covariance and calibration',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if prepared['status']=='below_detection_floor':
        return {**result,'status':'below_detection_floor'}
    if prepared['ratio_red_green']<=0 or (table.ratio_table.ratios<=0).any():
        raise ValueError('log-ratio likelihood requires positive observation and every model ratio')
    z=np.array([np.log(prepared['ratio_red_green']),prepared['green_norm']])
    predictions=np.stack([np.log(table.ratio_table.ratios),predicted_green],axis=-1)
    differences=(z-predictions).reshape(-1,2)
    covariance=np.asarray(noise['covariance'])
    chol=np.linalg.cholesky(covariance)
    standardized=np.linalg.solve(chol,differences.T).T
    with np.errstate(over='ignore',invalid='ignore'):
        distances=np.sum(standardized**2,axis=1).reshape(masses.shape)
    if not np.isfinite(distances).all():
        raise ValueError('observation/model discrepancy exceeds numeric likelihood range')
    log_likelihood=-.5*(distances+2*np.log(np.diag(chol)).sum()+2*np.log(2*np.pi))
    support=masses>0
    log_prior=np.full(masses.shape,-np.inf)
    log_prior[support]=np.log(masses[support])
    log_posterior=log_likelihood+log_prior
    posterior=np.exp(log_posterior-logsumexp(log_posterior))
    posterior/=posterior.sum()
    flat_likelihood=bool(np.ptp(log_likelihood[support])<=1e-10)
    supported_min=float(distances[support].min())
    poor_fit=supported_min>threshold
    sorted_mass=np.sort(posterior.ravel())[::-1]
    cutoff_index=min(int(np.searchsorted(np.cumsum(sorted_mass),level)),len(sorted_mass)-1)
    cutoff=sorted_mass[cutoff_index]
    region=(posterior>=cutoff) | np.isclose(posterior,cutoff,rtol=1e-12,atol=0)
    map_mask=np.isclose(posterior,posterior.max(),rtol=1e-12,atol=0)
    def pairs(mask):
        return [dict(dose_uM=float(table.ratio_table.doses_uM[i]),
                     elapsed_time_hours=float(table.ratio_table.times_hours[j]),
                     posterior_mass=float(posterior[i,j])) for i,j in np.argwhere(mask)]
    result.update(posterior_computed=True,posterior_mass=posterior.tolist(),
        doses_uM=table.ratio_table.doses_uM.tolist(),times_hours=table.ratio_table.times_hours.tolist(),
        dose_marginal_mass=posterior.sum(axis=1).tolist(),time_marginal_mass=posterior.sum(axis=0).tolist(),
        credible_set=pairs(region),credible_mass_achieved=float(posterior[region].sum()),
        credible_set_definition='highest node-mass set; all cutoff ties retained; not a continuous HPD region',
        map_grid_pairs=pairs(map_mask),likelihood_uninformative_on_prior_support=flat_likelihood,
        minimum_mahalanobis_squared=supported_min,model_fit_check_passed=not poor_fit,
        model_fit_check_is_a_calibrated_p_value=False,
        status='poor_model_fit' if poor_fit else ('prior_only' if flat_likelihood else 'posterior_available'))
    return result
