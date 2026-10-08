"""The registration offsets have to survive the trip from the manifest into
the ImageSet, including the shapes a hand-edited CSV actually takes: missing
column, empty cell, negative number, text."""
import runpy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def build_image_sets():
    return runpy.run_path(str(ROOT / "scripts" / "run_pipeline.py"))["_build_image_sets"]


@pytest.fixture
def images(tmp_path):
    paths = []
    for name, value in (("green.png", 100), ("red.png", 120), ("bf.png", 140)):
        path = tmp_path / name
        Image.fromarray(np.full((8, 8), value, dtype=np.uint8)).save(path)
        paths.append(path)
    return paths


def manifest_row(green, red, bright):
    return {
        "session_id": "s", "timepoint": 1, "acquisition_order": 1,
        "green_path": str(green), "red_path": str(red), "bright_field_path": str(bright),
        "exposure_ms_green": 100.0, "exposure_ms_red": 100.0,
        "nd_filter_green": 0.0, "nd_filter_red": 0.0, "objective": "60X",
        "burner_hours": 5.0, "lamp_warmup_minutes": 30.0,
    }


def build(tmp_path, build_image_sets, images, **extra):
    row = manifest_row(*images)
    row.update(extra)
    path = tmp_path / "manifest.csv"
    pd.DataFrame([row]).to_csv(path, index=False)
    return build_image_sets(path)[0]


def test_shifts_are_read(tmp_path, build_image_sets, images):
    image_set = build(tmp_path, build_image_sets, images,
                      brightfield_shift_dy=-14, brightfield_shift_dx=11)
    assert (image_set.brightfield_shift_dy, image_set.brightfield_shift_dx) == (-14.0, 11.0)


@pytest.mark.parametrize("extra", [{}, dict(brightfield_shift_dy="", brightfield_shift_dx="")])
def test_an_absent_or_blank_shift_means_no_shift(tmp_path, build_image_sets, images, extra):
    """Not an error: most manifests predate the column, and an unmeasured
    shift is usually small. Where it is not, the masks land on background
    and the per-cell QC says so."""
    image_set = build(tmp_path, build_image_sets, images, **extra)
    assert (image_set.brightfield_shift_dy, image_set.brightfield_shift_dx) == (0.0, 0.0)


def test_a_nonsense_shift_is_rejected_rather_than_silently_zero(tmp_path, build_image_sets, images):
    """Zero is the default for 'not measured', so a typo that fell through
    as zero would be indistinguishable from a measured alignment."""
    with pytest.raises(ValueError, match="brightfield_shift_dx"):
        build(tmp_path, build_image_sets, images, brightfield_shift_dx="left a bit")
