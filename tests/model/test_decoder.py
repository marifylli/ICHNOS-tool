import json

import numpy as np
import pytest

from ichnos.decoder import DoseTable, build_dose_table, decode_dose, load_dose_table, save_dose_table


def table(ratios=None):
    return DoseTable(
        [0, 25, 75], [.5, 1],
        ratios if ratios is not None else [[.1, .2], [.4, .5], [.7, .8]],
        {"source_kind": "model_generated", "ratio_direction": "red/green",
         "experimentally_validated": False},
    )


def decode(model=None, ratios=(.4, .5), **kwargs):
    options = dict(times_hours=[.5, 1], calibrated_ratios=ratios, ratio_tolerance=.01)
    options.update(kwargs)
    return decode_dose(model or table(), **options)


def test_recovers_unique_grid_dose():
    result = decode()
    assert result["status"] == "unique_grid_match"
    assert result["dose_estimate_uM"] == 25
    assert result["candidate_doses_uM"] == [25]
    assert result["experimentally_validated"] is False
    assert result["time_estimated"] is False


def test_duplicate_response_is_ambiguous_not_arbitrary_pick():
    result = decode(table([[.4, .5], [.4, .5], [.7, .8]]))
    assert result["status"] == "ambiguous"
    assert result["candidate_doses_uM"] == [0, 25]
    assert result["dose_estimate_uM"] is None


def test_large_tolerance_preserves_all_candidates():
    result = decode(ratio_tolerance=1)
    assert result["status"] == "ambiguous"
    assert result["candidate_doses_uM"] == [0, 25, 75]


def test_outside_response_envelope_is_rejected_without_extrapolation():
    result = decode(ratios=(10, 20))
    assert result["status"] == "out_of_response_domain"
    assert result["dose_estimate_uM"] is None


def test_no_joint_match_within_envelope_abstains():
    result = decode(ratios=(.1, .8))
    assert result["status"] == "no_grid_match"
    assert result["candidate_doses_uM"] == []


def test_offgrid_ratio_does_not_get_a_continuous_dose_estimate():
    result = decode(ratios=(.25, .35))
    assert result["status"] == "no_grid_match"
    assert result["dose_estimate_uM"] is None


def test_below_explicit_green_floor_abstains_even_with_exact_ratio():
    result = decode(green_measurements=[0, 100], green_floor=10)
    assert result["status"] == "below_detection_floor"
    assert result["dose_estimate_uM"] is None
    assert result["detection_floor_checked"] is True


def test_no_automatic_experimental_detection_floor():
    assert decode()["detection_floor_checked"] is False


@pytest.mark.parametrize("kwargs", [
    {"times_hours": [.5, 2]}, {"times_hours": [1, .5]},
    {"calibrated_ratios": [.5]}, {"calibrated_ratios": [np.nan, .5]},
    {"calibrated_ratios": [-1, .5]}, {"calibrated_ratios": [True, .5]},
    {"ratio_tolerance": 0}, {"ratio_tolerance": -1},
    {"ratio_tolerance": np.inf}, {"ratio_tolerance": True},
    {"ratio_tolerance": [.01]}, {"green_floor": 10},
    {"green_measurements": [100, 100]},
    {"green_floor": 10, "green_measurements": [100]},
])
def test_rejects_invalid_or_mismatched_observations(kwargs):
    with pytest.raises(ValueError):
        decode(**kwargs)


def test_artifact_roundtrip_and_integrity(tmp_path):
    path = tmp_path / "table.json"
    model = table()
    save_dose_table(model, path)
    restored = load_dose_table(path)
    assert restored.to_dict() == model.to_dict()
    with pytest.raises(FileExistsError):
        save_dose_table(model, path)
    data = json.loads(path.read_text())
    data["ratios"][0][0] = 999
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="modified"):
        load_dose_table(path)


@pytest.mark.parametrize("variant", ["ox", "er"])
def test_model_generated_grid_matches_independent_protocol_simulation(variant):
    import libsbml
    from ichnos import naming, protocol
    from ichnos.fluorescence import model_fluorescence

    model = build_dose_table(
        variant=variant, doses_uM=[0, 25, 75], times_hours=[.5, 1],
        initialization="finite-preincubation", preincubation_hours=2,
    )
    exposure = protocol.StressProtocol(variant=variant, dose=25)
    runner, loaded_model = protocol.load_protocol_model(model.provenance["baseline_sbml"], exposure)
    result, _ = protocol.run_protocol_after_preincubation(
        runner, exposure, times_hours=[.5, 1], preincubation_hours=2,
        id_to_name=naming.id_to_name_map(loaded_model),
    )
    recovered = decode_dose(
        model, times_hours=[.5, 1],
        calibrated_ratios=model_fluorescence(result)["ratio_red_green"],
        ratio_tolerance=1e-9,
    )
    assert recovered["status"] == "unique_grid_match"
    assert recovered["dose_estimate_uM"] == 25
    assert all(run["exposure"]["model"] == "constant" for run in model.provenance["runs"])
    assert model.provenance["f"] > 0
    assert libsbml.readSBMLFromString(model.provenance["baseline_sbml"]).getModel() is not None


@pytest.mark.parametrize("kwargs", [
    {"doses_uM": [25, 0]}, {"doses_uM": [0]},
    {"initialization": "finite-preincubation", "preincubation_hours": None},
    {"initialization": "equilibrium", "preincubation_hours": 2},
])
def test_builder_rejects_invalid_grid_or_initialization(kwargs):
    options = dict(variant="ox", doses_uM=[0, 25], times_hours=[.5, 1], initialization="equilibrium")
    options.update(kwargs)
    with pytest.raises(ValueError):
        build_dose_table(**options)
