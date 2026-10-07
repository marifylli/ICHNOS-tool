"""Traceable paired population medians -> snapshot compatibility/posterior."""
import json
import math
from .population import SUMMARY_METHOD
from .snapshot_decoder import decode_snapshot
from .posterior import decode_posterior


def observation_from_summary(summary, calibration, *, source):
    """Never infer biological IDs or acquisition metadata from reference artifacts."""
    if summary.get('summary_status') != 'ready':
        raise ValueError('insufficient_cells: summary is not ready')
    if summary.get('summary_method') != SUMMARY_METHOD:
        raise ValueError('unsupported summary_method')
    if calibration.get('summary_method') != SUMMARY_METHOD:
        raise ValueError('calibration must explicitly declare the same summary_method; rebuild from matching references')
    for name in ('session_id', 'specimen_id', 'biological_replicate_id'):
        if not isinstance(summary.get(name), str) or not summary[name].strip():
            raise ValueError(f'{name} is required; cells/fields are not independent cultures')
    acquisition = summary.get('acquisition_json')
    if not isinstance(acquisition, str):
        raise ValueError('acquisition_json is required')
    acquisition = json.loads(acquisition)
    for name in ('objective', 'exposure_ms_green', 'exposure_ms_red', 'nd_filter_green', 'nd_filter_red'):
        if name not in acquisition or acquisition[name] != summary.get(name):
            raise ValueError(f'acquisition_json disagrees with summary {name}')
    ratio, green = float(summary['ratio_red_green_median']), float(summary['corrected_green_median'])
    if not math.isfinite(ratio) or not math.isfinite(green) or ratio <= 0 or green < 0:
        raise ValueError('invalid paired summary: posterior requires positive ratio and non-negative green')
    if not isinstance(source, str) or not source.strip():
        raise ValueError('source is required')
    return dict(session_id=summary['session_id'], specimen_id=summary['specimen_id'],
        biological_replicate_id=summary['biological_replicate_id'], data_kind=summary['data_kind'],
        acquisition=acquisition, ratio_scale='image', image_ratio_red_green=ratio,
        corrected_green=green, observation_unit='biological_replicate_summary',
        summary_method=SUMMARY_METHOD, source=source,
        summary_selection=dict(sample_id=summary['sample_id'], timepoint=int(summary['timepoint']),
                               n_cells_used=int(summary['n_cells_used'])))


def decode_summary(table, calibration, summary, *, source, noise=None, prior=None,
                   max_mahalanobis_squared=None, ratio_tolerance=None, green_tolerance=None,
                   green_floor=None, credible_mass=.95):
    observation = observation_from_summary(summary, calibration, source=source)
    if noise is not None:
        if any(row.get('summary_method') != SUMMARY_METHOD for row in noise['source_rows']):
            raise ValueError('noise training must declare the same summary_method as target')
        if green_floor is not None and green_floor != noise['green_floor']:
            raise ValueError('population green floor must match noise fitting')
        return decode_posterior(table, calibration, noise, observation, prior=prior,
            max_mahalanobis_squared=max_mahalanobis_squared, credible_mass=credible_mass)
    return decode_snapshot(table, calibration, observation, ratio_tolerance=ratio_tolerance,
        green_tolerance=green_tolerance, green_floor=green_floor)
