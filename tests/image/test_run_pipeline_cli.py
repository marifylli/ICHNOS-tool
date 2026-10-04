import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from PIL import Image


@pytest.fixture
def cli():
    script = (
        Path(__file__).resolve().parents[2]
        / "scripts"
        / "run_pipeline.py"
    )
    spec = importlib.util.spec_from_file_location(
        "run_pipeline_under_test", script
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def inputs(tmp_path):
    green = np.zeros((8, 8, 3), dtype=np.uint8)
    green[:] = [10, 20, 30]
    green[0, 0, 2] = 200

    red = np.zeros((8, 8, 3), dtype=np.uint8)
    red[:] = [40, 50, 60]
    red[1, 1, 2] = 201

    paths = {}
    for name, array in (
        ("green", green),
        ("red", red),
        ("control_green", green),
        ("control_red", red),
    ):
        path = tmp_path / f"{name}.png"
        Image.fromarray(array).save(path)
        paths[name] = str(path)

    manifest = tmp_path / "manifest.csv"
    pd.DataFrame([{
        "session_id": "session_a",
        "timepoint": 30,
        "acquisition_order": 1,
        "green_path": paths["green"],
        "red_path": paths["red"],
        "exposure_ms_green": 100,
        "exposure_ms_red": 100,
        "nd_filter_green": 1,
        "nd_filter_red": 1,
        "objective": "40x",
        "burner_hours": 10,
        "lamp_warmup_minutes": 30,
    }]).to_csv(manifest, index=False)

    controls = tmp_path / "controls.csv"
    pd.DataFrame([{
        "session_id": "session_a",
        "control_green_path": paths["control_green"],
        "control_red_path": paths["control_red"],
    }]).to_csv(controls, index=False)

    return manifest, controls, paths, green, red


@pytest.mark.parametrize("method", ["G", "sum"])
def test_main_preserves_rgb_and_saturation_options(
    cli, inputs, tmp_path, monkeypatch, method
):
    manifest, controls, paths, green, red = inputs
    output = tmp_path / "output.csv"
    loaded = []
    captured = {}
    real_load = cli.load_image

    def tracked_load(path):
        loaded.append(str(path))
        return real_load(path)

    def fake_calibration(green_plane, red_plane):
        captured["control_green"] = green_plane.copy()
        captured["control_red"] = red_plane.copy()
        return {
            "bleed_green_to_red": 0.05,
            "residual_check": 0.0,
        }

    def fake_process(image_sets, out_path, **kwargs):
        captured["image_sets"] = image_sets
        captured["options"] = kwargs
        pd.DataFrame([{
            "session_id": "session_a",
            "cell_id": 1,
            "qc_pass": True,
        }]).to_csv(out_path, index=False)
        return out_path

    monkeypatch.setattr(cli, "load_image", tracked_load)
    monkeypatch.setattr(
        cli, "calibrate_crosstalk_from_control", fake_calibration
    )
    monkeypatch.setattr(cli, "process_experiment", fake_process)
    monkeypatch.setattr(sys, "argv", [
        "run_pipeline.py",
        "--manifest", str(manifest),
        "--controls", str(controls),
        "--out", str(output),
        "--green-extraction", method,
        "--red-extraction", "R",
        "--saturation-value", "200",
    ])

    cli.main()

    expected_green = (
        green[..., 1]
        if method == "G"
        else green.astype(float).sum(axis=-1)
    )
    expected_red = red[..., 0]

    assert len(captured["image_sets"]) == 1
    image_set = captured["image_sets"][0]
    np.testing.assert_allclose(image_set.green, expected_green)
    np.testing.assert_allclose(image_set.red, expected_red)
    np.testing.assert_allclose(
        captured["control_green"], expected_green
    )
    np.testing.assert_allclose(
        captured["control_red"], expected_red
    )

    expected_mask = np.zeros((8, 8), dtype=bool)
    expected_mask[0, 0] = True
    expected_mask[1, 1] = True
    np.testing.assert_array_equal(
        image_set.raw_saturation_mask, expected_mask
    )
    assert image_set.saturation_value == 200
    assert captured["options"]["bleed_green_to_red"] == {
        "session_a": 0.05
    }

    # Each sample and control image is read exactly once.
    assert sorted(loaded) == sorted(paths.values())


def test_manifest_rejects_mismatched_image_shapes(cli, inputs):
    manifest, _, paths, _, _ = inputs
    Image.fromarray(
        np.zeros((7, 8, 3), dtype=np.uint8)
    ).save(paths["red"])

    with pytest.raises(ValueError, match="matching shapes"):
        cli._build_image_sets(
            manifest,
            green_extraction="G",
            red_extraction="R",
        )