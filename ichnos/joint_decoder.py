"""Joint discrete dose/elapsed-time compatibility with known measurement spacing."""
from __future__ import annotations

import hashlib
from pathlib import Path
import numpy as np
from .decoder import DoseTable, _vector, decode_dose


def decode_joint_dose_time(
    table: DoseTable, *, relative_times_hours, elapsed_time_grid_hours,
    calibrated_ratios, ratio_tolerance, green_measurements=None, green_floor=None,
):
    """Estimate elapsed time at the first measurement, retaining every grid pair.

    Relative times start at zero; spacing is known. Exposure precedes the first
    measurement. Every shifted time must already exist in the response table.
    No interpolation, extrapolation, continuous identifiability or CI is implied.
    """
    relative = _vector(relative_times_hours, "relative_times_hours")
    if relative[0] != 0 or (np.diff(relative) <= 0).any():
        raise ValueError("relative times must start at zero and strictly increase")
    elapsed = _vector(elapsed_time_grid_hours, "elapsed_time_grid_hours")
    if (np.diff(elapsed) <= 0).any():
        raise ValueError("elapsed time grid must strictly increase")
    slices = []
    for offset in elapsed:
        shifted = relative + offset
        indices = []
        for time in shifted:
            matches = np.flatnonzero(np.isclose(table.times_hours, time, rtol=0, atol=1e-12))
            if len(matches) != 1:
                raise ValueError("each elapsed + relative time must match exactly one table time")
            indices.append(int(matches[0]))
        slices.append((float(offset), DoseTable(
            table.doses_uM, table.times_hours[indices], table.ratios[:, indices], table.provenance,
        )))
    candidates = []
    predictions = []
    checked = []
    for offset, view in slices:
        result = decode_dose(
            view, times_hours=view.times_hours, calibrated_ratios=calibrated_ratios,
            ratio_tolerance=ratio_tolerance, green_measurements=green_measurements,
            green_floor=green_floor,
        )
        checked.append(result)
        predictions.extend(view.ratios.tolist())
        for dose in result["candidate_doses_uM"]:
            candidates.append({"dose_uM": dose, "elapsed_time_at_first_measurement_hours": offset})
    first = checked[0]
    output = {
        "artifact_id": table.to_dict()["artifact_id"],
        "decoder_code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "relative_times_hours": relative.tolist(), "elapsed_time_grid_hours": elapsed.tolist(),
        "calibrated_ratios": first["calibrated_ratios"], "ratio_tolerance": first["ratio_tolerance"],
        "candidate_pairs": candidates, "dose_estimate_uM": None,
        "elapsed_time_estimate_hours": None, "time_origin": "first measurement",
        "discrete_grid_only": True, "time_estimated": False,
        "experimentally_validated": False, "tolerance_is_confidence_interval": False,
        "detection_floor_checked": first["detection_floor_checked"],
        "dose_grid_uM": table.doses_uM.tolist(),
    }
    if green_floor is not None:
        output.update(green_measurements=first["green_measurements"], green_floor=first["green_floor"])
    if first["status"] == "below_detection_floor":
        return {**output, "status": "below_detection_floor"}
    if len(candidates) == 1:
        return {**output, "status": "unique_grid_pair", "time_estimated": True,
                "dose_estimate_uM": candidates[0]["dose_uM"],
                "elapsed_time_estimate_hours": candidates[0]["elapsed_time_at_first_measurement_hours"]}
    if candidates:
        return {**output, "status": "ambiguous"}
    predicted = np.asarray(predictions)
    observed = np.asarray(first["calibrated_ratios"])
    tolerance = np.asarray(first["ratio_tolerance"])
    outside = ((observed < predicted.min(axis=0) - tolerance) |
               (observed > predicted.max(axis=0) + tolerance)).any()
    return {**output, "status": "out_of_response_domain" if outside else "no_grid_pair_match"}
