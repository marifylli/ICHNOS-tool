import json

import numpy as np
import pytest

from ichnos.decoder import DoseTable, build_dose_table
from ichnos.continuous_decoder import (
    decode_continuous_dose, load_interpolation_validation,
    save_interpolation_validation, validate_interpolation,
)


def table(ratios=None):
    return DoseTable([0, 50, 100], [.5, 1],
        ratios if ratios is not None else [[.2, .4], [.5, 1], [.8, 1.6]],
        {"source_kind": "model_generated", "ratio_direction": "red/green", "experimentally_validated": False})


def analytic_validation(model, error=0, allowed=1e-9):
    # Controlled linear/nonmonotonic fixtures isolate the mathematical decoder.
    doses = (model.doses_uM[:-1, None] + np.diff(model.doses_uM)[:, None] * [.25, .5, .75]).ravel()
    predicted = np.column_stack([
        np.interp(doses, model.doses_uM, model.ratios[:, index]) for index in range(2)
    ])
    return dict(schema_version=1, table_artifact_id=model.to_dict()["artifact_id"],
        source_kind="model_generated", experimentally_validated=False,
        tested_doses_uM=doses.tolist(), simulated_ratios=(predicted + error).tolist(),
        max_abs_error_per_segment_time=[[error, error]] * 2,
        allowed_error=[allowed, allowed], passed=bool(error <= allowed))


def decode(model=None, observations=(.425, .85), validation=None, **kwargs):
    model = model or table()
    options = dict(times_hours=[.5, 1], calibrated_ratios=observations,
        ratio_tolerance=.01, interpolation_validation=validation or analytic_validation(model))
    options.update(kwargs)
    return decode_continuous_dose(model, **options)


def test_recovers_between_grid_points_and_returns_compatibility_region():
    result = decode()
    assert result["status"] == "single_compatible_region"
    assert result["dose_estimate_uM"] == pytest.approx(37.5)
    assert len(result["compatible_dose_regions_uM"]) == 1
    lower, upper = result["compatible_dose_regions_uM"][0]
    assert lower < 37.5 < upper
    assert result["discrete_grid_only"] is False
    assert result["tolerance_is_confidence_interval"] is False


def test_nonmonotonic_response_returns_separate_regions():
    model = table([[.2, .4], [.8, 1.6], [.2, .4]])
    result = decode(model, observations=(.5, 1))
    assert result["status"] == "ambiguous"
    assert len(result["compatible_dose_regions_uM"]) == 2
    assert result["dose_estimate_uM"] is None


def test_equal_minima_are_not_arbitrarily_selected_even_in_connected_region():
    model = table([[.2, .4], [.8, 1.6], [.2, .4]])
    result = decode(model, observations=(.5, 1), ratio_tolerance=2)
    assert len(result["compatible_dose_regions_uM"]) == 1
    assert result["status"] == "ambiguous"
    assert result["equally_good_estimates_uM"] == pytest.approx([25, 75])


def test_flat_response_has_no_point_estimate():
    model = table([[.5, 1], [.5, 1], [.5, 1]])
    result = decode(model, observations=(.5, 1))
    assert result["status"] == "flat_response"
    assert result["compatible_dose_regions_uM"] == [[0, 100]]
    assert result["dose_estimate_uM"] is None


def test_no_extrapolation():
    result = decode(observations=(2, 4))
    assert result["status"] == "out_of_response_domain"
    assert result["dose_estimate_uM"] is None


def test_inconsistent_multitime_signal_abstains():
    assert decode(observations=(.2, 1.6))["status"] == "no_continuous_match"


def test_failed_interpolation_check_abstains():
    model = table()
    result = decode(model, validation=analytic_validation(model, error=.1, allowed=.01))
    assert result["status"] == "interpolation_check_failed"
    assert result["dose_estimate_uM"] is None


def test_detection_floor_abstention_is_preserved():
    assert decode(green_measurements=[0, 100], green_floor=10)["status"] == "below_detection_floor"


def test_wrong_validation_table_is_rejected():
    validation = analytic_validation(table())
    validation["table_artifact_id"] = "wrong"
    with pytest.raises(ValueError, match="different table"):
        decode(validation=validation)


def test_incomplete_or_inconsistent_validation_is_rejected():
    validation = analytic_validation(table())
    validation["tested_doses_uM"].pop()
    with pytest.raises(ValueError, match="holdout"):
        decode(validation=validation)
    validation = analytic_validation(table())
    validation["max_abs_error_per_segment_time"] = [[1, 1], [1, 1]]
    with pytest.raises(ValueError, match="inconsistent"):
        decode(validation=validation)


def test_validation_artifact_integrity_and_exclusive_write(tmp_path):
    path = tmp_path / "validation.json"
    validation = analytic_validation(table())
    save_interpolation_validation(validation, path)
    assert load_interpolation_validation(path) == validation
    with pytest.raises(FileExistsError):
        save_interpolation_validation(validation, path)
    data = json.loads(path.read_text())
    data["passed"] = False
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="modified"):
        load_interpolation_validation(path)


def test_independent_offgrid_simulation_is_recovered_with_checked_interpolation():
    model = build_dose_table(
        variant="ox", doses_uM=[25, 30, 35, 40, 45, 50], times_hours=[.5, 1],
        initialization="finite-preincubation", preincubation_hours=2,
    )
    validation = validate_interpolation(model, allowed_error=.001)
    assert validation["passed"]
    independent = build_dose_table(
        variant="ox", doses_uM=[37.5, 42.5], times_hours=[.5, 1],
        initialization="finite-preincubation", preincubation_hours=2,
    )
    # These doses are absent from the calibration grid.
    for dose, observation in zip(independent.doses_uM, independent.ratios):
        result = decode_continuous_dose(
            model, times_hours=[.5, 1], calibrated_ratios=observation,
            ratio_tolerance=1e-5, interpolation_validation=validation,
        )
        assert result["status"] == "single_compatible_region"
        assert abs(result["dose_estimate_uM"] - dose) < .5
        lower, upper = result["compatible_dose_regions_uM"][0]
        assert lower <= dose <= upper
