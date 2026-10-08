"""Explicit, identity-checked manual review sensitivity; never changes baseline QC."""
import numpy as np
import pandas as pd
from .population import GROUP_COLUMNS, summarize_cells

IDENTITY = ['session_id', 'sample_id', 'timepoint', 'acquisition_order', 'cell_id']


def review_summaries(cells, review, *, data_kind, min_cells=3, green_floor=0.0):
    required = IDENTITY + ['review_reason']
    if not set(required).issubset(review.columns):
        raise ValueError(f'review requires columns: {required}')
    for frame in (cells, review):
        if frame[IDENTITY].isna().any().any() or frame.duplicated(IDENTITY).any():
            raise ValueError('missing or duplicate object identity')
    if not review.review_reason.map(lambda v: isinstance(v, str) and bool(v.strip())).all():
        raise ValueError('review_reason must be a nonempty string')
    joined = review[required].merge(cells[IDENTITY], on=IDENTITY, how='left',
                                    indicator=True, validate='one_to_one')
    if joined['_merge'].ne('both').any():
        raise ValueError('review contains identities absent from cell CSV')
    baseline = summarize_cells(cells, data_kind=data_kind, min_cells=min_cells,
                               green_floor=green_floor)
    marked = cells.merge(review[required], on=IDENTITY, how='left', validate='one_to_one')
    filtered = marked.copy()
    filtered['qc_pass'] = filtered['qc_pass'] & filtered.review_reason.isna()
    sensitivity = summarize_cells(filtered, data_kind=data_kind, min_cells=min_cells,
                                  green_floor=green_floor)
    metrics = ['n_cells_used', 'ratio_red_green_median', 'corrected_green_median', 'summary_status']
    comparison = baseline[GROUP_COLUMNS + metrics].merge(
        sensitivity[GROUP_COLUMNS + metrics], on=GROUP_COLUMNS,
        suffixes=('_baseline', '_without_flagged'), validate='one_to_one')
    for metric in ['ratio_red_green_median', 'corrected_green_median']:
        denominator = comparison[metric + '_baseline'].replace(0, np.nan)
        comparison[metric + '_change_percent'] = (
            100 * (comparison[metric + '_without_flagged'] - denominator) / denominator)
    return baseline, sensitivity, comparison, marked
