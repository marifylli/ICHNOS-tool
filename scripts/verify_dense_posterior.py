"""Dense synthetic posterior benchmark; conditional inference, not wet-lab validation."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from ichnos.posterior import decode_posterior, fit_replicate_noise, prior_grid
from ichnos.snapshot_decoder import build_observable_table, build_snapshot_calibration, save_artifact


def ratio_posterior(ratio, table, variance, prior):
    """Marginal Gaussian likelihood: green is omitted, not conditioned on."""
    masses = prior_grid(table, prior)
    log_weights = -.5 * (np.log(ratio) - np.log(table.ratio_table.ratios)) ** 2 / variance
    support = masses > 0
    log_weights[~support] = -np.inf
    log_weights[support] += np.log(masses[support])
    return np.exp(log_weights - logsumexp(log_weights))


def summarize_mass(mass, doses, times, truth_index, level=.95):
    """Retain cutoff/MAP ties; report error range rather than arbitrary tie-break."""
    ordered = np.sort(mass.ravel())[::-1]
    cutoff = ordered[min(np.searchsorted(np.cumsum(ordered), level), mass.size - 1)]
    region = (mass >= cutoff) | np.isclose(mass, cutoff, rtol=1e-12, atol=0)
    maps = np.argwhere(np.isclose(mass, mass.max(), rtol=1e-12, atol=0))
    i, j = truth_index
    dose_errors = np.abs(doses[maps[:, 0]] - doses[i])
    time_errors = np.abs(times[maps[:, 1]] - times[j])
    return dict(map_tie_count=len(maps), map_dose_abs_error_min_uM=float(dose_errors.min()),
        map_dose_abs_error_max_uM=float(dose_errors.max()),
        map_time_abs_error_min_hours=float(time_errors.min()),
        map_time_abs_error_max_hours=float(time_errors.max()),
        credible_node_count=int(region.sum()), credible_node_fraction=float(region.mean()),
        credible_mass_achieved=float(mass[region].sum()), truth_in_credible_set=bool(region[i, j]))


def run(snapshot_dir, out_dir, *, dose_nodes=17, time_nodes=33, repeats=3,
        seed=20261006, clearance_rate=.2):
    if dose_nodes < 3 or dose_nodes % 2 == 0 or time_nodes < 3 or repeats < 1:
        raise ValueError('use an odd dose_nodes >= 3, time_nodes >= 3 and repeats >= 1')
    if not np.isfinite(clearance_rate) or clearance_rate <= 0:
        raise ValueError('clearance_rate must be finite and positive')
    out = Path(out_dir)
    if out.exists():
        raise FileExistsError(f'output directory already exists: {out}')
    config = json.loads((Path(snapshot_dir) / 'calibration-config.json').read_text())
    if config['source_kind'] != 'synthetic':
        raise ValueError('benchmark requires synthetic reference configuration')
    doses, times = np.linspace(0, 800, dose_nodes), np.linspace(0, 8, time_nodes)
    table = build_observable_table(variant='ox', doses_uM=doses, times_hours=times, initialization='equilibrium')
    calibration = build_snapshot_calibration(table, **config)
    truth_clearance = build_observable_table(variant='ox', doses_uM=doses, times_hours=times,
        initialization='equilibrium', clearance_rate_per_hour=clearance_rate)
    ilow, ihigh = calibration['image_green_anchors']
    mlow, mhigh = calibration['model_green_anchors']
    ratio_scale = calibration['ratio_calibration']['c_session']
    covariance = np.array([[.02**2, -.35*.02*.035], [-.35*.02*.035, .035**2]])
    rng = np.random.default_rng(seed)

    def observation(source, i, j, name):
        mean = [np.log(source.ratio_table.ratios[i, j]), (source.green[i, j] - mlow)/(mhigh-mlow)]
        z = rng.multivariate_normal(mean, covariance)
        return dict(session_id=calibration['session_id'], specimen_id=name,
            biological_replicate_id=name+'-culture', data_kind='synthetic', ratio_scale='image',
            acquisition=calibration['acquisition'], observation_unit='biological_replicate_summary',
            image_ratio_red_green=float(np.exp(z[0])*ratio_scale),
            corrected_green=float(ilow+z[1]*(ihigh-ilow)),
            source='synthetic Gaussian summary draw; not experimental biological replicate')

    training = []
    for group, (i, j) in enumerate([(0, 0), (dose_nodes//2, time_nodes//2), (dose_nodes-1, time_nodes-1)]):
        for k in range(12):
            row = observation(table, i, j, f'train-{group}-{k}')
            row.update(condition_id=f'group-{group}', dose_uM=float(doses[i]), time_hours=float(times[j]))
            training.append(row)
    noise = fit_replicate_noise(table, calibration, training, green_floor=100.)
    prior = dict(kind='uniform_density', dose_bounds_uM=[0, 800], time_bounds_hours=[0, 8])
    out.mkdir(parents=True)
    for name, data in [('table', table.to_dict()), ('clearance-truth-table', truth_clearance.to_dict()),
                       ('calibration', calibration), ('noise', noise), ('prior', prior)]:
        save_artifact(data, out/f'{name}.json')
    rows, observation_rows, illustrations = [], [], {}
    for scenario, source in [('matched', table), ('clearance_misspecified', truth_clearance)]:
        for i in range(dose_nodes):
            for j in range(time_nodes):
                for repeat in range(repeats):
                    name = f'{scenario}-{i}-{j}-{repeat}'
                    target = observation(source, i, j, name)
                    result = decode_posterior(table, calibration, noise, target, prior=prior,
                                              max_mahalanobis_squared=25, credible_mass=.95)
                    observation_rows.append(dict(scenario=scenario, dose_uM=float(doses[i]),
                        elapsed_time_hours=float(times[j]), repeat=repeat, observation=target))
                    if not result['posterior_computed']:
                        raise ValueError('synthetic target below floor; benchmark must report rather than drop it')
                    joint = np.asarray(result['posterior_mass'])
                    ratio = ratio_posterior(target['image_ratio_red_green']/ratio_scale, table,
                                             noise['covariance'][0][0], prior)
                    for method, mass in [('ratio_only', ratio), ('ratio_and_green', joint)]:
                        metrics = summarize_mass(mass, doses, times, (i, j))
                        rows.append(dict(scenario=scenario, method=method, dose_uM=float(doses[i]),
                            elapsed_time_hours=float(times[j]), repeat=repeat,
                            poor_fit=result['status']=='poor_model_fit' if method=='ratio_and_green' else '',
                            **metrics))
                    if i == dose_nodes//2 and j == time_nodes-1 and repeat == 0:
                        illustrations[scenario] = dict(ratio_only=ratio.tolist(), ratio_and_green=joint.tolist(),
                            observation=target, joint_result=result)
            print(f'{scenario}: {i+1}/{dose_nodes} doses complete', flush=True)
    with (out/'trials.csv').open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    import pandas as pd
    frame = pd.DataFrame(rows)
    numeric = list(summarize_mass(np.ones((2,2))/4, np.array([0,1]), np.array([0,1]), (0,0)))
    frame.groupby(['scenario','method','dose_uM','elapsed_time_hours'])[numeric].mean().to_csv(out/'per-node.csv')
    summaries = []
    for (scenario, method), group in frame.groupby(['scenario','method']):
        summary = dict(scenario=scenario, method=method, n_trials=len(group))
        summary.update({name:float(group[name].mean()) for name in numeric})
        summary['poor_fit_fraction'] = float(group.poor_fit.mean()) if method=='ratio_and_green' else None
        summaries.append(summary)
    report = dict(source_kind='synthetic', experimentally_validated=False, seed=seed,
        dose_nodes=dose_nodes, time_nodes=time_nodes, repeats_per_node=repeats, prior=prior,
        known_covariance=covariance.tolist(), clearance_rate_per_hour=clearance_rate,
        decoding_model='default ox, constant exposure, equilibrium initialization', summaries=summaries,
        limitations='on-grid synthetic summary-space draws; fixed calibration and fitted covariance; no image-noise or experimental coverage validation; few repeats per node',
        truth_sampling='equal repeats at every node; aggregate inclusion is not prior-predictive coverage')
    save_artifact(report, out/'verification.json')
    save_artifact(observation_rows, out/'observations.json')
    save_artifact(illustrations, out/'illustrations.json')
    render(illustrations, doses, times, out)
    return report


def render(illustrations, doses, times, out):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    vmax = max(np.max(data[method]) for data in illustrations.values() for method in ('ratio_only','ratio_and_green'))
    for row, (scenario, data) in enumerate(illustrations.items()):
        for col, method in enumerate(('ratio_only','ratio_and_green')):
            axis = axes[row, col]
            mass = np.asarray(data[method])
            raster = axis.pcolormesh(times, doses, mass, shading='nearest', vmin=0, vmax=vmax, cmap='viridis')
            axis.plot([8], [400], 'rx', markersize=9, label='Generating dose/time')
            axis.set(xlabel='Elapsed time (h)', ylabel='Initial dose (uM)', title=f'{scenario.replace("_", " ")} — {method.replace("_", " ")}')
            axis.legend(loc='upper left', fontsize=8)
    fig.colorbar(raster, ax=axes, label='Conditional probability mass per grid node')
    fig.suptitle('Same synthetic target per row; shared prior and color scale\nSummary-space noise; not experimental validation')
    fig.savefig(out/'posterior-comparison.png', dpi=180)
    fig.savefig(out/'posterior-comparison.pdf')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot-dir', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--dose-nodes', type=int, default=17)
    parser.add_argument('--time-nodes', type=int, default=33)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--seed', type=int, default=20261006)
    parser.add_argument('--clearance-rate', type=float, default=.2)
    args = parser.parse_args()
    try:
        report = run(**vars(args))
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f'Error: {exc}\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
