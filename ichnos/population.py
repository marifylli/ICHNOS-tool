"""QC-filtered sample summaries; cells are not biological replicates."""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import pandas as pd

from .calibration import SessionCalibration, calibrate_image_ratios


@dataclass(frozen=True)
class CalibrationBinding:
    calibration: SessionCalibration
    reference_id: str


GROUP_COLUMNS = ["session_id", "sample_id", "condition_id", "timepoint"]
ACQUISITION_COLUMNS = [
    "sampling_time_hours", "measurement_time_hours", "exposure_ms_green",
    "exposure_ms_red", "nd_filter_green", "nd_filter_red", "objective",
]
SUMMARY_COLUMNS = [
    *GROUP_COLUMNS, *ACQUISITION_COLUMNS,
    "n_cells_total", "n_cells_qc_pass", "n_cells_used", "n_images",
    "ratio_red_green_median", "ratio_red_green_q25", "ratio_red_green_q75",
    "calibrated_ratio_red_green_median", "calibrated_ratio_red_green_q25",
    "calibrated_ratio_red_green_q75", "c_session", "calibration_source_kind",
    "calibration_reference_id", "calibration_source", "data_kind",
]


def summarize_cells(
    cells: pd.DataFrame, *, data_kind: str,
    calibrations: dict[str, CalibrationBinding] | None = None,
    min_cells: int = 3, green_floor: float = 0.0,
) -> pd.DataFrame:
    """Median of cell ratios per sample/timepoint, with descriptive quartiles.

    A missing calibration stays missing. Supplying calibrations requires
    coverage of every session. Never combine sessions or biological samples.
    green_floor is in corrected image units, not a calibrated detection limit.
    """
    if data_kind not in {"synthetic", "experimental"}:
        raise ValueError("data_kind must be synthetic or experimental")
    if isinstance(min_cells, bool) or not isinstance(min_cells, int) or min_cells < 1:
        raise ValueError("min_cells must be a positive integer")
    if isinstance(green_floor, bool) or not math.isfinite(green_floor) or green_floor < 0:
        raise ValueError("green_floor must be finite and non-negative")
    required = {
        "session_id", "sample_id", "timepoint", "cell_id", "acquisition_order",
        "qc_pass", "corrected_mean_green", "corrected_mean_red", "ratio_red_green",
        *ACQUISITION_COLUMNS,
    }
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"missing cell columns: {sorted(missing)}")
    frame = cells.copy()
    if "condition_id" not in frame:
        frame["condition_id"] = None
    for name in ("session_id", "sample_id"):
        if not frame[name].map(lambda value: isinstance(value, str) and bool(value.strip())).all():
            raise ValueError(f"{name} is required for every cell; no inferred grouping")
    if not frame["qc_pass"].map(lambda value: isinstance(value, (bool, np.bool_))).all():
        raise ValueError("qc_pass must contain booleans, not truthy strings or numbers")
    identity = ["session_id", "sample_id", "timepoint", "acquisition_order", "cell_id"]
    if frame.duplicated(identity).any():
        raise ValueError("duplicate cell identity; refusing to count cells twice")
    if calibrations is not None:
        missing_sessions = set(frame.session_id) - set(calibrations)
        if missing_sessions:
            raise ValueError(f"missing calibration for sessions: {sorted(missing_sessions)}")
        for session in set(frame.session_id):
            binding = calibrations[session]
            calibration = binding.calibration
            if calibration.session_id != session or calibration.reference_id != binding.reference_id:
                raise ValueError("calibration session or reference_id mismatch")
            expected = "synthetic" if data_kind == "synthetic" else "experimental_reference"
            if calibration.source_kind != expected:
                raise ValueError("calibration source_kind does not match data_kind; synthetic calibration cannot calibrate experimental images")
            session_rows = frame[frame.session_id == session]
            for name in ("exposure_ms_green", "exposure_ms_red", "nd_filter_green", "nd_filter_red", "objective"):
                if session_rows[name].isna().any() or session_rows[name].nunique() != 1:
                    raise ValueError(f"calibration requires one known {name} setting per session")

    rows = []
    for key, group in frame.groupby(GROUP_COLUMNS, dropna=False, sort=False):
        row = dict(zip(GROUP_COLUMNS, key))
        for name in ACQUISITION_COLUMNS:
            if group[name].nunique(dropna=False) != 1:
                raise ValueError(f"inconsistent {name} within sample/timepoint {key}")
            value = group[name].iloc[0]
            row[name] = None if pd.isna(value) else value
        green = pd.to_numeric(group.corrected_mean_green, errors="raise").to_numpy(float)
        red = pd.to_numeric(group.corrected_mean_red, errors="raise").to_numpy(float)
        ratios = pd.to_numeric(group.ratio_red_green, errors="raise").to_numpy(float)
        eligible = (
            group.qc_pass.to_numpy(bool) & np.isfinite(green) & np.isfinite(red)
            & (green > green_floor) & (red >= 0) & np.isfinite(ratios) & (ratios >= 0)
        )
        if not np.allclose(ratios[eligible], red[eligible] / green[eligible], rtol=1e-8, atol=1e-12):
            raise ValueError("ratio_red_green differs from red/green channel intensities")
        used = ratios[eligible]
        row.update(
            n_cells_total=len(group), n_cells_qc_pass=int(group.qc_pass.sum()),
            n_cells_used=len(used), n_images=int(group.acquisition_order.nunique()),
            data_kind=data_kind,
        )
        for column in SUMMARY_COLUMNS:
            if column not in row:
                row[column] = None
        if len(used) >= min_cells:
            q25, median, q75 = np.quantile(used, [.25, .5, .75])
            row.update(ratio_red_green_median=median, ratio_red_green_q25=q25, ratio_red_green_q75=q75)
        if calibrations is not None:
            binding = calibrations[row["session_id"]]
            calibration = binding.calibration
            row.update(
                c_session=calibration.c_session,
                calibration_source_kind=calibration.source_kind,
                calibration_reference_id=calibration.reference_id,
                calibration_source=calibration.source,
            )
            if len(used) >= min_cells:
                corrected = calibrate_image_ratios(
                    used, calibration, session_id=row["session_id"],
                    reference_id=binding.reference_id,
                )
                q25, median, q75 = np.quantile(corrected, [.25, .5, .75])
                row.update(
                    calibrated_ratio_red_green_median=median,
                    calibrated_ratio_red_green_q25=q25,
                    calibrated_ratio_red_green_q75=q75,
                )
        rows.append(row)
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)
