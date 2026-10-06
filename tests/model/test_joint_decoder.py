import json
import subprocess
import sys
import numpy as np
import pytest
from ichnos.decoder import DoseTable, build_dose_table, save_dose_table
from ichnos.joint_decoder import decode_joint_dose_time


def table(values=None):
    return DoseTable([0, 25], [.5, 1, 1.5], values or [[.1, .2, .3], [.4, .8, 1.2]],
                     dict(source_kind="model_generated", ratio_direction="red/green", experimentally_validated=False))


def decode(model=None, **kwargs):
    options = dict(relative_times_hours=[0, .5], elapsed_time_grid_hours=[.5, 1],
                   calibrated_ratios=[.8, 1.2], ratio_tolerance=.001)
    options.update(kwargs)
    return decode_joint_dose_time(model or table(), **options)


def test_unique_pair_and_time_origin():
    result = decode()
    assert result["status"] == "unique_grid_pair"
    assert result["dose_estimate_uM"] == 25
    assert result["elapsed_time_estimate_hours"] == 1
    assert result["time_estimated"] is True
    assert result["experimentally_validated"] is False
    assert result["tolerance_is_confidence_interval"] is False


def test_single_measurement_dose_time_confounding():
    result = decode(table([[.1, .4, .9], [.4, .8, 1.2]]), relative_times_hours=[0], calibrated_ratios=[.4])
    assert result["status"] == "ambiguous"
    assert len(result["candidate_pairs"]) == 2
    assert result["dose_estimate_uM"] is None
    assert result["elapsed_time_estimate_hours"] is None


def test_flat_time_response_keeps_both_times():
    result = decode(table([[.1, .1, .1], [.8, .8, .8]]), calibrated_ratios=[.8, .8])
    assert result["status"] == "ambiguous"
    assert len(result["candidate_pairs"]) == 2


@pytest.mark.parametrize("options", [
    dict(relative_times_hours=[.1, .5]), dict(relative_times_hours=[0, 0]),
    dict(elapsed_time_grid_hours=[1, .5]), dict(elapsed_time_grid_hours=[.5, .5]),
    dict(elapsed_time_grid_hours=[.75]), dict(calibrated_ratios=[np.nan, 1.2]),
    dict(ratio_tolerance=0), dict(green_measurements=[1, 2]),
])
def test_invalid_inputs_rejected(options):
    with pytest.raises(ValueError):
        decode(**options)


def test_no_pair_inside_envelope():
    assert decode(calibrated_ratios=[.1, 1.2])["status"] == "no_grid_pair_match"


def test_outside_envelope_and_floor():
    assert decode(calibrated_ratios=[5, 6])["status"] == "out_of_response_domain"
    result = decode(green_measurements=[0, 100], green_floor=1)
    assert result["status"] == "below_detection_floor"
    assert result["candidate_pairs"] == []


def test_independent_model_run_recovers_grid_pair():
    options = dict(variant="ox", initialization="finite-preincubation", preincubation_hours=2)
    grid = build_dose_table(**options, doses_uM=[25, 75], times_hours=[.5, 1, 1.5])
    independent = build_dose_table(**options, doses_uM=[75, 100], times_hours=[1, 1.5])
    result = decode(grid, calibrated_ratios=independent.ratios[0], ratio_tolerance=1e-8)
    assert result["status"] == "unique_grid_pair"
    assert result["dose_estimate_uM"] == 75
    assert result["elapsed_time_estimate_hours"] == 1


def test_cli_records_provenance_and_preserves_existing_output(tmp_path):
    model = tmp_path / "table.json"
    save_dose_table(table(), model)
    observations = tmp_path / "ratios.csv"
    observations.write_text("relative_time_hours,calibrated_ratio_red_green\n0,0.8\n0.5,1.2\n")
    out = tmp_path / "result.json"
    command = [sys.executable, "scripts/run_joint_decoder.py", "--table", str(model),
               "--observations", str(observations), "--elapsed-times-hours", ".5", "1",
               "--ratio-tolerance", ".001", "--data-kind", "synthetic", "--out", str(out)]
    run = subprocess.run(command, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    data = json.loads(out.read_text())
    assert data["status"] == "unique_grid_pair"
    assert data["artifact_id"] == table().to_dict()["artifact_id"]
    assert len(data["observations_csv_sha256"]) == 64
    before = out.read_bytes()
    assert subprocess.run(command, capture_output=True).returncode != 0
    assert out.read_bytes() == before
