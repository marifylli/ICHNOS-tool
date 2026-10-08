"""Per-cell CSV schema for the ICHNOS pipeline's final dataset (Stage 8).

The one definition of CellRecord -- ichnos_image.export builds these and
writes them to CSV; anything else that reads that CSV (a future Stage 6
ratio/decoder package, analysis notebooks) should treat this as the schema,
not re-derive it from the CSV's header.

Field list reflects the current (Olympus, no ApoTome) plan: no
grid_setting/section_thickness, per-session instrument metadata instead
(session_id, exposure/ND per channel, burner_hours, lamp_warmup_minutes,
acquisition_order), plus the QC flags from Stage 7 (focus_score, sat_flag,
edge_flag, registration_shift_px, lamp_flag).
"""
from __future__ import annotations

from dataclasses import dataclass, fields
from typing import Optional

import math
from numbers import Real


def validate_elapsed_hours(value, *, field_name: str) -> float | None:
    """Optional elapsed time in hours from stress onset."""
    if value is None:
        return None

    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(
            f"{field_name} must be a finite non-negative number or None"
        )

    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError(
            f"{field_name} must be a finite non-negative number or None"
        )

    return value


def validate_optional_identifier(value, *, field_name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field_name} must be a non-empty string or None")
    return value

# gia pragmatikous xronous peiramatos

@dataclass
class CellRecord:
    # identity
    session_id: str
    timepoint: int
    cell_id: int

    # segmentation-derived
    area_px: int
    edge_flag: bool

    # raw per-channel intensities (pre-correction)
    raw_mean_green: float
    raw_mean_red: float

    # corrected per-channel intensities (background + optional flat-field + crosstalk; no photobleaching)
    corrected_mean_green: float
    corrected_mean_red: float
    integrated_green: float
    integrated_red: float

    # derived
    ratio_red_green: Optional[float]

    # QC (Stage 7)
    sat_flag: bool
    focus_score: float
    registration_shift_px: float
    lamp_flag: bool
    qc_pass: bool

    # acquisition / instrument metadata (per session, repeated per row for convenience)
    exposure_ms_green: float
    exposure_ms_red: float
    nd_filter_green: float
    nd_filter_red: float
    objective: str
    burner_hours: float
    lamp_warmup_minutes: float
    acquisition_order: int

    sampling_time_hours: float | None = None
    measurement_time_hours: float | None = None
    sample_id: str | None = None
    condition_id: str | None = None
    specimen_id: str | None = None
    biological_replicate_id: str | None = None
    acquisition_json: str | None = None
    sat_flag_legacy: bool | None = None
    qc_pass_legacy_saturation: bool | None = None
    focus_score_green: float | None = None
    focus_score_red: float | None = None
    contrast_to_noise_green: float | None = None
    contrast_to_noise_red: float | None = None
    focus_agreement: float | None = None
    focus_status: str = 'not_evaluated'
    focus_qc_mode: str = 'report'
    qc_reasons: str = ''


CSV_COLUMNS = [f.name for f in fields(CellRecord)]


def validate(record: CellRecord) -> None:
    """Cheap structural sanity checks; not a substitute for Stage 7 QC itself."""
    if record.area_px <= 0:
        raise ValueError(f"cell {record.cell_id}: non-positive area_px")
    if record.raw_mean_green < 0 or record.raw_mean_red < 0:
        raise ValueError(f"cell {record.cell_id}: negative raw intensity")

    for name in ("sample_id", "condition_id", "specimen_id", "biological_replicate_id"):
        validate_optional_identifier(getattr(record, name), field_name=name)

    for name in ("sampling_time_hours", "measurement_time_hours"):
        validate_elapsed_hours(
            getattr(record, name),
            field_name=name,
        )
