import hashlib
import json
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
from types import SimpleNamespace

import pandas as pd
import pytest

from ichnos.calibration import fit_session_calibration, model_reference_id, save_session_calibration
from ichnos.decoder import build_dose_table, save_dose_table
from ichnos.sample_decoder import decode_sample, check_reference_compatibility


ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def baseline(tmp_path_factory):
    directory = tmp_path_factory.mktemp("sample-linkage")
    reference_dir = directory / "reference"
    completed = subprocess.run([
        sys.executable, str(ROOT / "scripts/run_protocol.py"),
        "--variant", "ox", "--dose", "25", "--dose-units", "uM",
        "--times-hours", ".25", ".5", ".75", "1", "1.5", "2",
        "--initialization", "finite-preincubation", "--preincubation-hours", "2",
        "--out-dir", str(reference_dir),
    ], cwd=ROOT, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr
    reference = json.loads((reference_dir / "metadata.json").read_text())
    predicted = pd.read_csv(reference_dir / "fluorescence.csv").ratio_red_green.to_numpy()
    calibration = fit_session_calibration(
        predicted, 1.5 * predicted, session_id="s1",
        reference_id=model_reference_id(reference), source_kind="synthetic",
        source="model-generated known-scale reference",
    )
    save_session_calibration(calibration, directory / "calibration.json")
    pd.DataFrame([{
        "session_id": "s1", "reference_id": calibration.reference_id,
        "calibration_path": "calibration.json",
    }]).to_csv(directory / "calibrations.csv", index=False)
    table = build_dose_table(
        variant="ox", doses_uM=[0, 25, 75], times_hours=[.5, 1],
        initialization="finite-preincubation", preincubation_hours=2,
    )
    save_dose_table(table, directory / "table.json")
    # Use real segmentation/correction/extraction of simplified synthetic images.
    import numpy as np
    from ichnos_image.pipeline import ImageSet, process_image_set
    from ichnos_image.export import export_csv
    yy, xx = np.mgrid[:128, :128]
    mask = np.zeros((128, 128), dtype=bool)
    for cy, cx in ((32, 32), (64, 96), (96, 32)):
        mask |= (yy - cy) ** 2 + (xx - cx) ** 2 < 10 ** 2
    records = []
    for index, (time, ratio) in enumerate(zip(table.times_hours, table.ratios[2])):
        green = np.where(mask, 1000., 0.)
        images = ImageSet(
            green=green, red=1.5 * ratio * green,
            bright_field=np.where(mask, .2, .5),
            session_id="s1", sample_id="sample-a", condition_id="dose75",
            timepoint=index, acquisition_order=index,
            exposure_ms_green=100, exposure_ms_red=200,
            nd_filter_green=0, nd_filter_red=0, objective="40X",
            burner_hours=1, lamp_warmup_minutes=30,
            sampling_time_hours=float(time), measurement_time_hours=float(time + .25),
        )
        cells = process_image_set(images, bleed_green_to_red=0)
        assert len(cells) == 3 and all(cell.qc_pass for cell in cells)
        records.extend(cells)
    export_csv(records, directory / "cells.csv")
    run = runpy.run_path(str(ROOT / "scripts/summarize_cells.py"))["run"]
    run(SimpleNamespace(
        cells=directory / "cells.csv", out_dir=directory / "summary",
        data_kind="synthetic", calibrations=directory / "calibrations.csv",
        min_cells=3, green_floor=0,
    ))
    return directory


@pytest.fixture
def dataset(baseline, tmp_path):
    directory = tmp_path / "data"
    shutil.copytree(baseline, directory)
    return directory


def options(directory):
    return dict(
        samples_csv=directory / "summary/samples.csv",
        summary_metadata=directory / "summary/metadata.json",
        dose_table=directory / "table.json",
        calibration_reference=directory / "reference/metadata.json",
        session_id="s1", sample_id="sample-a", time_field="sampling_time_hours",
        ratio_tolerance=1e-9,
    )


def edit_summary(directory, change):
    path = directory / "summary/samples.csv"
    frame = pd.read_csv(path)
    change(frame)
    frame.to_csv(path, index=False)
    metadata_path = directory / "summary/metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["samples_csv_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    metadata_path.write_text(json.dumps(metadata))


def test_population_to_decoder_recovers_unknown_sample_dose(dataset):
    result = decode_sample(**options(dataset))
    assert result["status"] == "unique_grid_match"
    assert result["dose_estimate_uM"] == 75
    linkage = result["sample_linkage"]
    assert linkage["time_field"] == "sampling_time_hours"
    assert linkage["c_session"] == pytest.approx(1.5)
    assert linkage["calibration_applied_again"] is False
    assert linkage["reference_dose_is_not_assumed_to_be_sample_dose"] is True
    assert linkage["n_cells_used"] == [3, 3]


@pytest.mark.parametrize("change", [
    {"session_id": "missing"}, {"sample_id": "missing"},
    {"condition_id": "missing"}, {"time_field": "timepoint"},
    {"time_field": "measurement_time_hours"},
])
def test_explicit_selection_and_time_contract(dataset, change):
    kwargs = options(dataset)
    kwargs.update(change)
    with pytest.raises(ValueError):
        decode_sample(**kwargs)


def test_rejects_modified_csv_without_matching_metadata(dataset):
    path = dataset / "summary/samples.csv"
    path.write_text(path.read_text() + "\n")
    with pytest.raises(ValueError, match="CSV differs"):
        decode_sample(**options(dataset))


@pytest.mark.parametrize("field,value", [
    ("sampling_time_hours", None), ("n_cells_used", 2),
    ("calibrated_ratio_red_green_median", None), ("c_session", 2),
    ("calibration_reference_id", "wrong-reference"),
    ("calibrated_ratio_red_green_median", 99),
])
def test_rejects_invalid_or_inconsistent_sample_rows(dataset, field, value):
    edit_summary(dataset, lambda frame: frame.__setitem__(field, value))
    with pytest.raises(ValueError):
        decode_sample(**options(dataset))


def test_multiple_conditions_require_explicit_selection(dataset):
    path = dataset / "summary/samples.csv"
    frame = pd.read_csv(path)
    other = frame.copy()
    other["condition_id"] = "other"
    combined = pd.concat([frame, other])
    path.write_text(combined.to_csv(index=False))
    meta_path = dataset / "summary/metadata.json"
    metadata = json.loads(meta_path.read_text())
    metadata["samples_csv_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    meta_path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="several conditions"):
        decode_sample(**options(dataset))
    assert decode_sample(**options(dataset), condition_id="dose75")["dose_estimate_uM"] == 75


def test_reference_metadata_must_match_recorded_calibration_id(dataset):
    path = dataset / "reference/metadata.json"
    metadata = json.loads(path.read_text())
    metadata["fluorescence"]["f"] = 2
    path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="reference_id"):
        decode_sample(**options(dataset))


@pytest.mark.parametrize("change,match", [
    ("model", "prepared model"), ("profile", "profile"),
    ("solver", "solver"), ("exposure", "exposure"),
    ("initialization", "initialization"), ("f", "f differs"),
])
def test_forward_model_compatibility_checks(dataset, change, match):
    from ichnos.decoder import load_dose_table
    table = load_dose_table(dataset / "table.json")
    reference = json.loads((dataset / "reference/metadata.json").read_text())
    if change == "model":
        reference["model_sbml_sha256"] = "different"
    elif change == "profile":
        reference["parameter_profile"]["name"] = "different"
    elif change == "solver":
        reference["solver_used"]["relative_tolerance"] = .1
    elif change == "exposure":
        reference["exposure"]["model"] = "first_order_decay"
    elif change == "initialization":
        reference["initialization"]["zero_stress_duration_hours"] = 3
    else:
        reference["fluorescence"]["f"] = 2
    with pytest.raises(ValueError, match=match):
        check_reference_compatibility(table, reference)


def test_old_summaries_require_regeneration(dataset):
    path = dataset / "summary/metadata.json"
    data = json.loads(path.read_text())
    data.pop("samples_csv_sha256")
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="regenerate"):
        decode_sample(**options(dataset))


def test_sample_decoder_cli_and_overwrite_guard(dataset):
    output = dataset / "decoded.json"
    command = [sys.executable, str(ROOT / "scripts/decode_sample.py")]
    for name, value in options(dataset).items():
        command.extend(["--" + name.replace("_", "-"), str(value)])
    command.extend(["--out", str(output)])
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(output.read_text())["dose_estimate_uM"] == 75
    previous = output.read_bytes()
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert completed.returncode != 0
    assert output.read_bytes() == previous


def test_continuous_mode_uses_the_same_checked_sample_linkage(dataset):
    from ichnos.decoder import load_dose_table
    from ichnos.continuous_decoder import validate_interpolation, save_interpolation_validation
    model = load_dose_table(dataset / "table.json")
    validation = validate_interpolation(model, allowed_error=.1)
    assert validation["passed"]
    path = dataset / "validation.json"
    save_interpolation_validation(validation, path)
    result = decode_sample(**options(dataset), mode="continuous", interpolation_validation=path)
    assert result["status"] == "single_compatible_region"
    assert result["dose_estimate_uM"] == pytest.approx(75)
    assert result["sample_linkage"]["decoder_mode"] == "continuous"
    assert result["sample_linkage"]["calibration_applied_again"] is False


@pytest.mark.parametrize("change", [
    {"mode": "continuous"}, {"mode": "unknown"},
    {"mode": "discrete", "interpolation_validation": "unused.json"},
])
def test_invalid_continuous_mode_contract_is_rejected(dataset, change):
    with pytest.raises(ValueError):
        decode_sample(**options(dataset), **change)
