"""Assemble per-cell records and write the pipeline CSV against ichnos_image.schema (Stage 8)."""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import pandas as pd

from .instrument import (
    FOCUS_SCORE_THRESHOLD,
    LAMP_WARMUP_THRESHOLD_MINUTES,
    REGISTRATION_SHIFT_THRESHOLD_PX,
)

from .extract import CellFeatures, compute_ratio
from .schema import CellRecord, CSV_COLUMNS, validate


def build_records(
    features: list[CellFeatures],
    *,
    session_id: str,
    timepoint: int,
    edge_flagged_ids: set[int],
    focus_score: float,
    registration_shift_px: float,
    exposure_ms_green: float,
    exposure_ms_red: float,
    nd_filter_green: float,
    nd_filter_red: float,
    objective: str,
    burner_hours: float,
    lamp_warmup_minutes: float,
    acquisition_order: int,
    focus_score_threshold: float = FOCUS_SCORE_THRESHOLD,
    registration_shift_threshold_px: float = REGISTRATION_SHIFT_THRESHOLD_PX,
    lamp_warmup_threshold_minutes: float = LAMP_WARMUP_THRESHOLD_MINUTES,
    sampling_time_hours: float | None = None,
    measurement_time_hours: float | None = None,
) -> list[CellRecord]:
    """Turn per-cell features + per-image acquisition/QC metadata into schema
    rows. QC flags here are image-level (focus_score, registration_shift_px)
    or set already on the feature (sat_flag) / passed in (edge_flagged_ids);
    qc_pass folds them into one column but never discards rows — dropping
    failed cells is a caller decision, not export's.

    Threshold defaults come from ichnos.config -- see that module for how
    each was calibrated and its caveats, rather than duplicating that here.
    """
    records = []
    lamp_flag = lamp_warmup_minutes < lamp_warmup_threshold_minutes
    for f in features:
        ratio = compute_ratio(f.corrected_mean_green, f.corrected_mean_red)
        edge_flag = f.cell_id in edge_flagged_ids
        qc_pass = (
            not f.sat_flag
            and not edge_flag
            and not lamp_flag
            and focus_score >= focus_score_threshold
            and registration_shift_px <= registration_shift_threshold_px
        )
        record = CellRecord(
            session_id=session_id,
            timepoint=timepoint,
            cell_id=f.cell_id,
            area_px=f.area_px,
            edge_flag=edge_flag,
            raw_mean_green=f.raw_mean_green,
            raw_mean_red=f.raw_mean_red,
            corrected_mean_green=f.corrected_mean_green,
            corrected_mean_red=f.corrected_mean_red,
            integrated_green=f.integrated_green,
            integrated_red=f.integrated_red,
            ratio_red_green=ratio,
            sat_flag=f.sat_flag,
            focus_score=focus_score,
            registration_shift_px=registration_shift_px,
            lamp_flag=lamp_flag,
            qc_pass=qc_pass,
            exposure_ms_green=exposure_ms_green,
            exposure_ms_red=exposure_ms_red,
            nd_filter_green=nd_filter_green,
            nd_filter_red=nd_filter_red,
            objective=objective,
            burner_hours=burner_hours,
            lamp_warmup_minutes=lamp_warmup_minutes,
            acquisition_order=acquisition_order,
            sampling_time_hours=sampling_time_hours,
            measurement_time_hours=measurement_time_hours,
        )
        validate(record)
        records.append(record)
    return records


def export_csv(records: list[CellRecord], out_path: str | Path, append: bool = False) -> Path:
    """Write records to CSV, columns ordered per ichnos_image.schema.CSV_COLUMNS.

    append=True adds rows to an existing file (e.g. accumulating timepoints
    across a time-lapse run) without repeating the header.
    """
    out_path = Path(out_path)
    df = pd.DataFrame([asdict(r) for r in records], columns=CSV_COLUMNS)
    write_header = not (append and out_path.exists())
    df.to_csv(out_path, mode="a" if append and out_path.exists() else "w", header=write_header, index=False)
    return out_path
