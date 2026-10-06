"""Bind one calibrated sample trajectory to a compatible dose table."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .calibration import model_reference_id, session_calibration_from_dict
from .decoder import decode_dose, load_dose_table


TIME_FIELDS = ("sampling_time_hours", "measurement_time_hours")


def check_reference_compatibility(table, reference):
    """Compare the model/protocol, not the reference dose with unknown dose."""
    provenance = table.provenance
    if provenance["variant"] != reference["protocol"]["variant"]:
        raise ValueError("variant differs from calibration reference")
    if provenance["parameter_profile"] != reference["parameter_profile"]:
        raise ValueError("parameter profile differs from calibration reference")
    runs = provenance["runs"]
    if not runs:
        raise ValueError("dose table has no protocol records")
    for run in runs:
        digest = run.get("prepared_model_sbml_sha256")
        if digest is None:
            raise ValueError("dose table lacks prepared-model identity; rebuild the table")
        if digest != reference["model_sbml_sha256"]:
            raise ValueError("prepared model differs from calibration reference")
        if run["solver_used"] != reference["solver_used"]:
            raise ValueError("solver differs from calibration reference")
        exposure = {key: value for key, value in run["exposure"].items() if key != "initial_dose"}
        expected = {key: value for key, value in reference["exposure"].items() if key != "initial_dose"}
        if exposure != expected:
            raise ValueError("exposure protocol differs from calibration reference")
        for key in ("method", "starting_state", "zero_stress_duration_hours", "criterion", "observables_checked"):
            if run["initialization"][key] != reference["initialization"][key]:
                raise ValueError(f"initialization {key} differs from calibration reference")
    for name in ("f", "eps"):
        if provenance[name] != reference["fluorescence"][name]:
            raise ValueError(f"{name} differs from calibration reference")
    if reference["fluorescence"]["ratio_direction"] != "red/green":
        raise ValueError("calibration reference uses a different ratio direction")


def decode_sample(
    *, samples_csv, summary_metadata, dose_table, calibration_reference,
    session_id: str, sample_id: str, time_field: str, ratio_tolerance,
    condition_id: str | None = None,
) -> dict:
    """Decode one explicitly selected trajectory without recalibrating it."""
    if time_field not in TIME_FIELDS:
        raise ValueError(f"time_field must be one of {TIME_FIELDS}")
    for name, value in (("session_id", session_id), ("sample_id", sample_id)):
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{name} must be explicit and non-empty")
    path = Path(samples_csv)
    metadata = json.loads(Path(summary_metadata).read_text())
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if metadata.get("samples_csv_sha256") is None:
        raise ValueError("summary lacks samples CSV identity; regenerate the summary")
    if metadata["samples_csv_sha256"] != digest:
        raise ValueError("samples CSV differs from summary metadata")
    if metadata.get("calibration_equation") != "calibrated_ratio = original_image_ratio / c_session" or metadata.get("f_applied_again") is not False:
        raise ValueError("unsupported summary calibration convention")
    data_kind = metadata.get("data_kind")
    if data_kind not in {"synthetic", "experimental"}:
        raise ValueError("summary data_kind is missing or invalid")
    artifact = metadata.get("calibrations", {}).get(session_id)
    if artifact is None:
        raise ValueError("no session calibration in summary metadata")
    calibration = session_calibration_from_dict(artifact)
    if calibration.session_id != session_id:
        raise ValueError("calibration session differs from selected session")
    expected_source = "synthetic" if data_kind == "synthetic" else "experimental_reference"
    if calibration.source_kind != expected_source:
        raise ValueError("calibration source does not match summary data kind")
    reference = json.loads(Path(calibration_reference).read_text())
    if model_reference_id(reference) != calibration.reference_id:
        raise ValueError("reference metadata does not match session calibration reference_id")
    table = load_dose_table(dose_table)
    check_reference_compatibility(table, reference)
    frame = pd.read_csv(path, dtype={"session_id": str, "sample_id": str, "condition_id": str})
    selected = frame[(frame.session_id == session_id) & (frame.sample_id == sample_id)]
    if condition_id is not None:
        selected = selected[selected.condition_id == condition_id]
    if selected.empty:
        raise ValueError("selected sample/session/condition has no rows")
    if selected.condition_id.nunique(dropna=False) != 1:
        raise ValueError("sample has several conditions; select condition_id explicitly")
    selected = selected.sort_values(time_field)
    for name, expected in (
        ("calibration_reference_id", calibration.reference_id),
        ("calibration_source_kind", calibration.source_kind),
        ("calibration_source", calibration.source), ("data_kind", data_kind),
    ):
        if not (selected[name] == expected).all():
            raise ValueError(f"sample rows have inconsistent {name}")
    if not np.allclose(selected.c_session.to_numpy(float), calibration.c_session, rtol=1e-12, atol=0):
        raise ValueError("sample coefficient differs from recorded session calibration")
    minimum = metadata.get("min_cells")
    if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 1:
        raise ValueError("invalid summary minimum cell count")
    used = selected.n_cells_used.to_numpy(float)
    if not np.isfinite(used).all() or (used < minimum).any() or (used != np.floor(used)).any():
        raise ValueError("selected sample has insufficient or invalid usable cell counts")
    raw = selected.ratio_red_green_median.to_numpy(float)
    corrected = selected.calibrated_ratio_red_green_median.to_numpy(float)
    if not np.isfinite(raw).all() or not np.isfinite(corrected).all():
        raise ValueError("selected sample has missing or non-finite ratio summaries")
    if not np.allclose(corrected, raw / calibration.c_session, rtol=1e-8, atol=1e-12):
        raise ValueError("calibrated median differs from original median / c_session")
    result = decode_dose(
        table, times_hours=selected[time_field].to_numpy(),
        calibrated_ratios=corrected, ratio_tolerance=ratio_tolerance,
    )
    result["sample_linkage"] = {
        "session_id": session_id, "sample_id": sample_id,
        "condition_id": (None if pd.isna(selected.condition_id.iloc[0]) else selected.condition_id.iloc[0]),
        "time_field": time_field, "timepoints": selected.timepoint.tolist(),
        "samples_csv_sha256": digest, "data_kind": data_kind,
        "calibration_reference_id": calibration.reference_id,
        "c_session": calibration.c_session, "calibration_applied_again": False,
        "n_cells_used": used.astype(int).tolist(),
        "summary_green_floor_corrected_image_units": metadata.get("green_floor_corrected_image_units"),
        "summary_green_floor_is_experimentally_calibrated": False,
        "original_ratio_medians": raw.tolist(),
        "experimental_instrument_validated": False,
        "reference_dose_is_not_assumed_to_be_sample_dose": True,
    }
    return result
