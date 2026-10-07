"""Isolate saturation-registration impact using one reprocessed cell table."""
import numpy as np
import pandas as pd


def saturation_impact(cells):
    required = {'sat_flag', 'sat_flag_legacy', 'qc_pass', 'qc_pass_legacy_saturation',
                'ratio_red_green', 'session_id', 'sample_id', 'timepoint'}
    if not required.issubset(cells.columns):
        raise ValueError('rerun raw images with the registered saturation masks before auditing')
    def median(frame, flag):
        values = pd.to_numeric(frame.loc[frame[flag], 'ratio_red_green'], errors='coerce')
        values = values[np.isfinite(values)]
        return float(values.median()) if len(values) else None
    rows = []
    for key, group in cells.groupby(['session_id', 'sample_id', 'timepoint'], dropna=False):
        rows.append(dict(session_id=str(key[0]), sample_id=None if pd.isna(key[1]) else str(key[1]),
            timepoint=int(key[2]), n_cells=len(group),
            newly_saturated=int((group.sat_flag & ~group.sat_flag_legacy).sum()),
            no_longer_saturated=int((~group.sat_flag & group.sat_flag_legacy).sum()),
            n_qc_before=int(group.qc_pass_legacy_saturation.sum()), n_qc_after=int(group.qc_pass.sum()),
            median_ratio_before=median(group, 'qc_pass_legacy_saturation'),
            median_ratio_after=median(group, 'qc_pass')))
    return dict(scope='saturation registration only; other processing held fixed',
                historical_results_reproduced=False, samples=rows)
