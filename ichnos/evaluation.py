"""Synthetic robustness benchmark; off-grid cell inclusion is not continuous coverage."""
from pathlib import Path
import numpy as np
import pandas as pd

from .artifacts import new_output_directory, write_json, code_provenance, sha256_file
from .calibration import fit_session_calibration, model_reference_id
from .snapshot_decoder import build_observable_table, build_snapshot_calibration
from .posterior import fit_replicate_noise, decode_posterior
from .population import SUMMARY_METHOD


SCENARIOS = {
    'matched': (1., 0., 1.),
    'ratio_gain_minus_10pct': (.9, 0., 1.),
    'ratio_gain_plus_10pct': (1.1, 0., 1.),
    'green_offset_minus_10pct_span': (1., -.1, 1.),
    'green_offset_plus_10pct_span': (1., .1, 1.),
    'green_span_minus_10pct': (1., 0., .9),
    'green_span_plus_10pct': (1., 0., 1.1),
}


def synthetic_calibration(table):
    """Known synthetic image mapping; no experimental references are inferred."""
    provenance = table.ratio_table.provenance
    run = provenance['runs'][0]
    reference = dict(model_sbml_sha256=run['prepared_model_sbml_sha256'],
        parameter_profile=provenance['parameter_profile'],protocol={'variant':provenance['variant'],
        'dose':run['dose_uM'],'dose_units':'uM'},
        observation_times_hours=table.ratio_table.times_hours.tolist(),
        initialization=run['initialization'],solver_used=run['solver_used'],exposure=run['exposure'],
        fluorescence={'ratio_direction':'red/green','f':provenance['f'],'eps':provenance['eps']})
    ratios = table.ratio_table.ratios[0]
    ratio = fit_session_calibration(ratios,1.5*ratios,session_id='synthetic-evaluation',
        reference_id=model_reference_id(reference),source_kind='synthetic',
        source='known scale; reference metadata from recorded model run; not independent biological validation')
    time = float(table.ratio_table.times_hours[-1])
    return build_snapshot_calibration(table,session_id='synthetic-evaluation',source_kind='synthetic',
        acquisition=dict(objective='synthetic',green_units='arbitrary synthetic units',
            extraction_green='synthetic',extraction_red='synthetic',gain_setting='synthetic',
            exposure_ms_green=100,exposure_ms_red=100),
        low_reference=dict(specimen_id='low',biological_replicate_id='low',source='synthetic known anchor',
            dose_uM=0.,time_hours=time,corrected_green=10000.),
        high_reference=dict(specimen_id='high',biological_replicate_id='high',source='synthetic known anchor',
            dose_uM=float(table.ratio_table.doses_uM[-1]),time_hours=time,corrected_green=11000.),
        min_image_span=1.,min_model_span=1e-10,ratio_calibration=ratio.to_dict(),ratio_reference_metadata=reference,
        ratio_reference_specimen_ids=['ratio-ref'],ratio_reference_replicate_ids=['ratio-ref'],
        summary_method=SUMMARY_METHOD)


def synthetic_observation(calibration, coordinates, name):
    ilow,ihigh = calibration['image_green_anchors']
    return dict(session_id=calibration['session_id'],specimen_id=name,biological_replicate_id=name,
        data_kind='synthetic',ratio_scale='image',acquisition=calibration['acquisition'],
        observation_unit='biological_replicate_summary',summary_method=SUMMARY_METHOD,
        image_ratio_red_green=float(np.exp(coordinates[0])*calibration['ratio_calibration']['c_session']),
        corrected_green=float(ilow+coordinates[1]*(ihigh-ilow)),
        source='synthetic Gaussian culture-summary draw; not noisy image simulation')


def transformed_means(table, calibration):
    low, high = calibration['model_green_anchors']
    return np.stack([np.log(table.ratio_table.ratios),(table.green-low)/(high-low)],axis=-1)


def posterior_metrics(result, doses, times, truth_dose, truth_time, on_grid):
    """Keep all MAP ties; off-grid inclusion means the containing midpoint cell."""
    metrics = dict(posterior_computed=result['posterior_computed'],status=result['status'],
        truth_in_discrete_set=None,truth_cell_in_set=False,credible_area_fraction=None,
        map_dose_error_min=None,map_dose_error_max=None,map_time_error_min=None,map_time_error_max=None,
        credible_mass_achieved=None)
    if not result['posterior_computed']:
        return metrics
    maps = result['map_grid_pairs']
    dose_errors = [abs(p['dose_uM']-truth_dose) for p in maps]
    time_errors = [abs(p['elapsed_time_hours']-truth_time) for p in maps]
    # Nearest-node/Voronoi cells, bounded by the declared domain. This extends
    # discrete nodes into cells only for a separately labeled diagnostic.
    i,j = int(np.argmin(abs(doses-truth_dose))),int(np.argmin(abs(times-truth_time)))
    region = {(p['dose_uM'],p['elapsed_time_hours']) for p in result['credible_set']}
    cell_hit = (float(doses[i]),float(times[j])) in region
    dose_widths = np.diff(np.r_[doses[0],(doses[:-1]+doses[1:])/2,doses[-1]])
    time_widths = np.diff(np.r_[times[0],(times[:-1]+times[1:])/2,times[-1]])
    area = sum(dose_widths[a]*time_widths[b] for a,d in enumerate(doses) for b,t in enumerate(times) if (d,t) in region)
    metrics.update(truth_cell_in_set=cell_hit, truth_in_discrete_set=cell_hit if on_grid else None,
        credible_area_fraction=float(area/((doses[-1]-doses[0])*(times[-1]-times[0]))),
        map_dose_error_min=float(min(dose_errors)),map_dose_error_max=float(max(dose_errors)),
        map_time_error_min=float(min(time_errors)),map_time_error_max=float(max(time_errors)),
        credible_mass_achieved=result['credible_mass_achieved'])
    return metrics


def run_evaluation(out_dir, *, variants=('ox','er'), dose_nodes=9,time_nodes=17,
                   trials=12,noise_scales=(.5,1.,2.),seed=20261007,
                   scenarios=tuple(SCENARIOS),refine=False):
    if dose_nodes < 3 or time_nodes < 3 or trials < 1:
        raise ValueError('dose/time nodes >=3 and trials >=1 required')
    if not variants or not set(variants) <= {'ox','er'}:
        raise ValueError('variants must contain ox and/or er')
    if not noise_scales or any(not np.isfinite(x) or x <= 0 for x in noise_scales):
        raise ValueError('noise scales must be finite and positive')
    if not scenarios or not set(scenarios) <= set(SCENARIOS):
        raise ValueError('unknown calibration scenario')
    rng = np.random.default_rng(seed)
    base_covariance = np.array([[.02**2,-.35*.02*.035],[-.35*.02*.035,.035**2]])
    rows, observations = [], []
    prior = dict(kind='uniform_density',dose_bounds_uM=[0,800],time_bounds_hours=[0,8])
    with new_output_directory(out_dir) as out:
        for variant in variants:
            # Same target draws and truths across coarse/refined grids. Each
            # noise level is fitted independently with matching synthetic noise.
            doses,times = np.linspace(0,800,dose_nodes),np.linspace(0,8,time_nodes)
            table = build_observable_table(variant=variant,doses_uM=doses,times_hours=times,initialization='equilibrium')
            calibration = synthetic_calibration(table)
            means = transformed_means(table,calibration)
            tables = [('base',table,calibration)]
            if refine:
                finer = build_observable_table(variant=variant,doses_uM=np.linspace(0,800,2*dose_nodes-1),
                    times_hours=np.linspace(0,8,2*time_nodes-1),initialization='equilibrium')
                tables.append(('refined',finer,synthetic_calibration(finer)))
            for grid,forward,cal in tables:
                write_json(out/f'{variant}-{grid}-table.json',forward.to_dict())
                write_json(out/f'{variant}-{grid}-calibration.json',cal)
            truths = []
            for k in range(trials):
                i,j = int(rng.integers(dose_nodes)),int(rng.integers(time_nodes))
                truths.append(('on_grid',float(doses[i]),float(times[j]),means[i,j]))
            # Independent forward simulations, not interpolation from decoder table.
            off_doses = np.sort(rng.uniform(0,800,trials))
            off_times = np.sort(rng.uniform(0,8,trials))
            off = build_observable_table(variant=variant,doses_uM=off_doses,times_hours=off_times,initialization='equilibrium')
            off_means = transformed_means(off,calibration)
            write_json(out/f'{variant}-offgrid-truth-table.json',off.to_dict())
            for k in range(trials):
                truths.append(('off_grid',float(off_doses[k]),float(off_times[k]),off_means[k,k]))
            # Permute time association to avoid a built-in dose/time correlation.
            permutation = rng.permutation(trials)
            truths[trials:] = [('off_grid',float(off_doses[k]),float(off_times[permutation[k]]),
                               off_means[k,permutation[k]]) for k in range(trials)]
            for scale in noise_scales:
                covariance = base_covariance * scale**2
                training = []
                for group,(i,j) in enumerate(((0,0),(dose_nodes//2,time_nodes//2),(dose_nodes-1,time_nodes-1))):
                    for replicate in range(12):
                        row = synthetic_observation(calibration,rng.multivariate_normal(means[i,j],covariance),
                                                    f'train-{variant}-{scale}-{group}-{replicate}')
                        row.update(condition_id=f'group-{group}',dose_uM=float(doses[i]),time_hours=float(times[j]))
                        training.append(row)
                fitted = {}
                for grid,forward,cal in tables:
                    noise = fit_replicate_noise(forward,cal,training,green_floor=0.)
                    fitted[grid]=noise
                    write_json(out/f'{variant}-{grid}-noise-{scale:g}.json',noise)
                for k,(truth_kind,dose,time,mean) in enumerate(truths):
                    perturbation = rng.multivariate_normal(np.zeros(2),covariance)
                    for scenario in scenarios:
                        ratio_gain,green_offset,green_gain = SCENARIOS[scenario]
                        # Identical latent noise across calibration scenarios; fixed
                        # decoding calibration and training, with perturbed targets.
                        z = mean + perturbation
                        z = np.array([z[0]+np.log(ratio_gain),green_gain*z[1]+green_offset])
                        name = f'target-{variant}-{scale}-{k}-{scenario}'
                        target = synthetic_observation(calibration,z,name)
                        observations.append(dict(trial_id=name,truth_dose_uM=dose,truth_time_hours=time,
                            truth_kind=truth_kind,observation=target))
                        for grid,forward,cal in tables:
                            result = decode_posterior(forward,cal,fitted[grid],target,prior=prior,
                                max_mahalanobis_squared=25.,credible_mass=.95)
                            rows.append(dict(variant=variant,grid=grid,noise_scale=scale,scenario=scenario,
                                truth_kind=truth_kind,trial_id=name,truth_dose_uM=dose,truth_time_hours=time,
                                **posterior_metrics(result,forward.ratio_table.doses_uM,
                                    forward.ratio_table.times_hours,dose,time,truth_kind=='on_grid')))
                print(f'{variant}: noise scale {scale:g} complete',flush=True)
        frame = pd.DataFrame(rows)
        frame.to_csv(out/'trials.csv',index=False)
        summaries = []
        for key,group in frame.groupby(['variant','grid','noise_scale','scenario','truth_kind']):
            summary = dict(zip(['variant','grid','noise_scale','scenario','truth_kind'],key))
            summary.update(n_trials=len(group), n_posteriors=int(group.posterior_computed.sum()),
                truth_cell_inclusion_fraction=float(group.truth_cell_in_set.mean()),
                poor_fit_fraction=float((group.status=='poor_model_fit').mean()))
            for name in ('map_dose_error_min','map_dose_error_max','map_time_error_min','map_time_error_max','credible_area_fraction'):
                summary[name] = float(group[name].mean()) if group[name].notna().any() else None
            summary['discrete_set_inclusion_fraction'] = (float(group.truth_in_discrete_set.mean())
                if key[-1]=='on_grid' else None)
            summaries.append(summary)
        report = dict(schema_version=1,experimentally_validated=False,source_kind='synthetic',seed=seed,
            dose_nodes=dose_nodes,time_nodes=time_nodes,trials_per_truth_kind=trials,
            noise_scales=list(noise_scales),calibration_scenarios={k:SCENARIOS[k] for k in scenarios},
            refined_grid=refine,n_trials=len(rows),summaries=summaries,
            synthetic_noise_covariance=base_covariance.tolist(),prior=prior,
            offgrid_metric='inclusion of the nearest-node midpoint cell; NOT continuous 95% coverage',
            limitations=['summary-space Gaussian noise, not microscopy simulation',
                         'calibration perturbations test robustness; uncertainty is not propagated',
                         'small trial counts are diagnostic, not proof of nominal coverage',
                         'default profiles and equilibrium; no experimental accuracy established'])
        write_json(out/'evaluation.json',report)
        write_json(out/'observations.json',observations)
        write_json(out/'manifest.json',dict(code=code_provenance(),
            outputs={p.name:sha256_file(p) for p in sorted(out.iterdir()) if p.is_file()}))
    return report
