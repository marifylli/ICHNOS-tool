"""A blank cell in the acquisition log must not destroy the run.

The log is kept by hand at the microscope and has gaps: an exposure nobody
wrote down, a burner-hours field left empty. Twice on the team's 2026-10-05
session the run processed every image and then raised while writing its own
provenance record, because a NaN cannot be written with allow_nan=False.
The gap is real and must stay visible -- as null, not as a substituted
number -- but it is not a reason to throw the work away.
"""
import json
import runpy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from PIL import Image

from ichnos_image import pipeline
from ichnos_image.pipeline import ImageSet

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def build_image_sets():
    return runpy.run_path(str(ROOT / "scripts" / "run_pipeline.py"))["_build_image_sets"]


@pytest.fixture
def images(tmp_path):
    paths = []
    for name, value in (("green.png", 100), ("red.png", 120)):
        path = tmp_path / name
        Image.fromarray(np.full((8, 8), value, dtype=np.uint8)).save(path)
        paths.append(path)
    return paths


def manifest_with(tmp_path, images, **gaps):
    row = {
        "session_id": "s", "timepoint": 1, "acquisition_order": 1,
        "green_path": str(images[0]), "red_path": str(images[1]),
        "exposure_ms_green": 100.0, "exposure_ms_red": 100.0,
        "nd_filter_green": 0.0, "nd_filter_red": 0.0, "objective": "60X",
        "burner_hours": 5.0, "lamp_warmup_minutes": 30.0,
    }
    row.update(gaps)
    path = tmp_path / "manifest.csv"
    pd.DataFrame([row]).to_csv(path, index=False)
    return path


@pytest.mark.parametrize("missing", ["exposure_ms_green", "exposure_ms_red",
                                     "nd_filter_green", "nd_filter_red"])
def test_a_missing_setting_is_recorded_as_null_not_a_crash(
    tmp_path, images, build_image_sets, missing
):
    image_set = build_image_sets(manifest_with(tmp_path, images, **{missing: ""}))[0]
    recorded = json.loads(image_set.acquisition_json)
    assert recorded[missing] is None
    assert recorded["exposure_ms_green"] in (None, 100.0)


def test_nothing_is_substituted_for_the_missing_value(tmp_path, images, build_image_sets):
    """A zero or a default here would be indistinguishable from a reading,
    and every stage that divides by an exposure would use it."""
    image_set = build_image_sets(manifest_with(tmp_path, images, exposure_ms_green=""))[0]
    assert json.loads(image_set.acquisition_json)["exposure_ms_green"] is None
    assert np.isnan(image_set.exposure_ms_green)


def cells_image_set(**overrides):
    rng = np.random.default_rng(0)
    yy, xx = np.mgrid[:120, :120]
    green = 30 + rng.normal(0, 2, (120, 120))
    for cy, cx in [(30, 30), (30, 90), (90, 30), (90, 90)]:
        green[(yy - cy) ** 2 + (xx - cx) ** 2 <= 64] += 60
    fields = dict(
        session_id="s", timepoint=0, acquisition_order=0, green=green,
        red=green * 0.8 + 5, exposure_ms_green=100, exposure_ms_red=100,
        nd_filter_green=0, nd_filter_red=0, objective="40X", burner_hours=1,
        lamp_warmup_minutes=30, sample_id="a",
    )
    fields.update(overrides)
    return ImageSet(**fields)


@pytest.mark.parametrize("gap", [dict(burner_hours=float("nan")),
                                 dict(lamp_warmup_minutes=float("nan")),
                                 dict(exposure_ms_red=float("nan"))])
def test_a_gap_does_not_lose_an_already_processed_run(tmp_path, gap):
    """The expensive failure: every image processed, then the run discarded
    at the final write."""
    out = pipeline.process_experiment(
        [cells_image_set(**gap)], tmp_path / "out.csv", bleed_green_to_red=0.05,
        segmentation_method="sparse", segmentation_kwargs=dict(min_size=30),
    )
    assert out.exists() and sum(1 for _ in out.open()) > 1
    metadata = json.loads(Path(str(out) + ".manifest.json").read_text())
    recorded = metadata["images"][0]
    name, _ = next(iter(gap.items()))
    assert recorded[name] is None


def test_the_gap_is_visible_in_the_manifest_rather_than_filled_in(tmp_path):
    out = pipeline.process_experiment(
        [cells_image_set(burner_hours=float("nan"))], tmp_path / "out.csv",
        bleed_green_to_red=0.05, segmentation_method="sparse",
        segmentation_kwargs=dict(min_size=30),
    )
    text = Path(str(out) + ".manifest.json").read_text()
    assert '"burner_hours": null' in text
    assert "NaN" not in text  # invalid JSON, and what the old code tried to write
