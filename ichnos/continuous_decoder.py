"""Piecewise-linear dose compatibility, with explicit interpolation checks."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .decoder import _digest, _vector, build_dose_table, decode_dose


FRACTIONS = np.array([.25, .5, .75])


def _tolerances(value, n_times, name):
    if np.asarray(value, dtype=object).ndim == 0:
        values = np.repeat(_vector([value], name, positive=True), n_times)
    else:
        values = _vector(value, name, positive=True)
    if values.shape != (n_times,):
        raise ValueError(f"{name} must be scalar or one value per time")
    return values


def _holdout_doses(table):
    return (table.doses_uM[:-1, None] + np.diff(table.doses_uM)[:, None] * FRACTIONS).ravel()


def validate_interpolation(table, *, allowed_error):
    """Check each segment at quarter, half and three-quarter grid doses.

    This is sampled computational verification, not a certified error bound
    over all doses and not experimental validation.
    """
    tolerance = _tolerances(allowed_error, len(table.times_hours), "allowed_error")
    provenance = table.provenance
    first = provenance["runs"][0]
    initialization = first["initialization"]
    method = initialization["method"]
    if method not in {"equilibrium_assumption", "finite_preincubation_assumption"}:
        raise ValueError("unsupported interpolation-check initialization")
    duration = initialization["zero_stress_duration_hours"]
    holdouts = _holdout_doses(table)
    all_doses = np.sort(np.concatenate([table.doses_uM, holdouts]))
    # The serialized criterion carries the original requested tolerance.
    criterion = initialization["criterion"]
    try:
        equilibrium_tolerance = float(criterion.rsplit("<=", 1)[1])
    except (ValueError, IndexError) as exc:
        raise ValueError("cannot recover equilibration criterion") from exc
    direct = build_dose_table(
        variant=provenance["variant"], profile_name=provenance["parameter_profile"]["name"],
        doses_uM=all_doses, times_hours=table.times_hours,
        initialization=("equilibrium" if method == "equilibrium_assumption" else "finite-preincubation"),
        preincubation_hours=(duration if method == "finite_preincubation_assumption" else None),
        equilibration_horizon_hours=duration,
        equilibration_tolerance=equilibrium_tolerance,
        clearance_rate_per_hour=first["exposure"]["clearance_rate_per_hour"],
    )
    if direct.provenance["parameter_profile"] != provenance["parameter_profile"]:
        raise ValueError("interpolation check profile differs from the table")
    expected_model = first.get("prepared_model_sbml_sha256")
    if expected_model is None:
        raise ValueError("table lacks model identity; rebuild it before interpolation checks")
    for run in direct.provenance["runs"]:
        if run["prepared_model_sbml_sha256"] != expected_model or run["solver_used"] != first["solver_used"]:
            raise ValueError("interpolation check forward model or solver differs from table")
    anchors = direct.ratios[np.searchsorted(all_doses, table.doses_uM)]
    if not np.allclose(anchors, table.ratios, rtol=1e-9, atol=1e-12):
        raise ValueError("re-simulated grid anchors differ from the stored table")
    actual = direct.ratios[np.searchsorted(all_doses, holdouts)]
    predicted = np.column_stack([
        np.interp(holdouts, table.doses_uM, table.ratios[:, index])
        for index in range(len(table.times_hours))
    ])
    errors = np.abs(predicted - actual).reshape(len(table.doses_uM) - 1, 3, -1).max(axis=1)
    return {
        "schema_version": 1, "table_artifact_id": table.to_dict()["artifact_id"],
        "source_kind": "model_generated", "experimentally_validated": False,
        "method": "piecewise-linear; direct simulation at fractions 0.25, 0.5, 0.75",
        "tested_doses_uM": holdouts.tolist(), "simulated_ratios": actual.tolist(),
        "interpolated_ratios": predicted.tolist(),
        "max_abs_error_per_segment_time": errors.tolist(),
        "allowed_error": tolerance.tolist(), "passed": bool((errors <= tolerance).all()),
        "sampled_check_only": True, "certified_uniform_error_bound": False,
        "direct_simulation_provenance": direct.provenance,
    }


def save_interpolation_validation(validation, path):
    payload = {**validation, "validation_id": _digest(validation)}
    with Path(path).open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def load_interpolation_validation(path):
    data = json.loads(Path(path).read_text())
    identifier = data.pop("validation_id")
    if _digest(data) != identifier:
        raise ValueError("modified interpolation validation artifact")
    return data


def _check_validation(table, validation):
    if validation.get("schema_version") != 1 or validation.get("table_artifact_id") != table.to_dict()["artifact_id"]:
        raise ValueError("interpolation validation belongs to a different table")
    if validation.get("source_kind") != "model_generated" or validation.get("experimentally_validated") is not False:
        raise ValueError("invalid interpolation validation provenance")
    holdouts = _holdout_doses(table)
    if not np.array_equal(np.asarray(validation["tested_doses_uM"]), holdouts):
        raise ValueError("interpolation validation must cover all prescribed holdout doses")
    actual = np.asarray(validation["simulated_ratios"], dtype=float)
    if actual.shape != (len(holdouts), len(table.times_hours)) or not np.isfinite(actual).all() or (actual < 0).any():
        raise ValueError("invalid direct-simulation validation ratios")
    predicted = np.column_stack([
        np.interp(holdouts, table.doses_uM, table.ratios[:, index])
        for index in range(len(table.times_hours))
    ])
    errors = np.abs(predicted - actual).reshape(len(table.doses_uM) - 1, 3, -1).max(axis=1)
    allowance = _tolerances(validation["allowed_error"], len(table.times_hours), "allowed_error")
    passed = bool((errors <= allowance).all())
    if validation.get("passed") is not passed or not np.allclose(errors, validation["max_abs_error_per_segment_time"], rtol=1e-12, atol=1e-15):
        raise ValueError("inconsistent interpolation check results")
    return allowance, passed


def decode_continuous_dose(
    table, *, times_hours, calibrated_ratios, ratio_tolerance,
    interpolation_validation, green_measurements=None, green_floor=None,
):
    """Intersect ratio compatibility constraints on every linear segment."""
    result = decode_dose(
        table, times_hours=times_hours, calibrated_ratios=calibrated_ratios,
        ratio_tolerance=ratio_tolerance,
        green_measurements=green_measurements, green_floor=green_floor,
    )
    result.pop("candidate_doses_uM")
    result.update(
        dose_estimate_uM=None, discrete_grid_only=False,
        method="piecewise_linear", compatible_dose_regions_uM=[],
        interpolation_checked=True, interpolation_check_is_sampled=True,
        interpolation_error_is_certified_bound=False,
    )
    allowance, passed = _check_validation(table, interpolation_validation)
    result["interpolation_error_allowance"] = allowance.tolist()
    result["interpolation_validation_id"] = _digest(interpolation_validation)
    if result["status"] == "below_detection_floor":
        return result
    if not passed:
        return {**result, "status": "interpolation_check_failed"}
    observations = np.asarray(calibrated_ratios, dtype=float)
    tolerance = np.asarray(result["ratio_tolerance"]) + allowance
    result["effective_ratio_tolerance"] = tolerance.tolist()
    segments = []
    for index in range(len(table.doses_uM) - 1):
        start = table.ratios[index]
        slope = table.ratios[index + 1] - start
        lower, upper = 0., 1.
        for prediction, delta, measured, margin in zip(start, slope, observations, tolerance):
            if delta == 0:
                if abs(prediction - measured) > margin:
                    lower, upper = 1., 0.
                    break
            else:
                ends = sorted(((measured - margin - prediction) / delta,
                               (measured + margin - prediction) / delta))
                lower, upper = max(lower, ends[0]), min(upper, ends[1])
        if lower > upper:
            continue
        left, right = table.doses_uM[index:index + 2]
        denominator = float(np.sum((slope / tolerance) ** 2))
        alpha = (float(np.sum(slope * (observations - start) / tolerance ** 2)) / denominator
                 if denominator > 0 else (lower + upper) / 2)
        alpha = float(np.clip(alpha, lower, upper))
        estimate = float(left + alpha * (right - left))
        score = float(np.sum(((start + alpha * slope - observations) / tolerance) ** 2))
        segments.append((float(left + lower * (right - left)), float(left + upper * (right - left)), estimate, score, denominator == 0))
    regions = []
    for lower, upper, estimate, score, flat in segments:
        if regions and lower <= regions[-1][1] + 1e-12 * max(1, abs(lower)):
            regions[-1][1] = max(regions[-1][1], upper)
        else:
            regions.append([lower, upper])
    result["compatible_dose_regions_uM"] = regions
    if not regions:
        outside = ((observations < table.ratios.min(axis=0) - tolerance)
                   | (observations > table.ratios.max(axis=0) + tolerance)).any()
        return {**result, "status": "out_of_response_domain" if outside else "no_continuous_match"}
    if len(regions) > 1:
        return {**result, "status": "ambiguous"}
    if any(segment[4] for segment in segments):
        return {**result, "status": "flat_response"}
    best = min(segments, key=lambda segment: segment[3])
    minimizers = sorted(segment[2] for segment in segments
                        if np.isclose(segment[3], best[3], rtol=1e-9, atol=1e-12))
    distinct = []
    for estimate in minimizers:
        if not distinct or not np.isclose(estimate, distinct[-1], rtol=1e-10, atol=1e-10):
            distinct.append(estimate)
    if len(distinct) > 1:
        return {**result, "status": "ambiguous", "equally_good_estimates_uM": distinct}
    return {
        **result, "status": "single_compatible_region", "dose_estimate_uM": best[2],
        "estimate_is_interpolated": True,
        "region_touches_grid_boundary": bool(regions[0][0] == table.doses_uM[0] or regions[0][1] == table.doses_uM[-1]),
    }
