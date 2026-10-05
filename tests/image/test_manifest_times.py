import runpy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def build_image_sets():
    namespace = runpy.run_path(
        str(ROOT / "scripts" / "run_pipeline.py")
    )
    return namespace["_build_image_sets"]


def _manifest_row(green_path, red_path):
    return {
        "session_id": "session-1",
        "timepoint": 1,
        "acquisition_order": 1,
        "green_path": str(green_path),
        "red_path": str(red_path),
        "exposure_ms_green": 100.0,
        "exposure_ms_red": 100.0,
        "nd_filter_green": 0.0,
        "nd_filter_red": 0.0,
        "objective": "60X",
        "burner_hours": 5.0,
        "lamp_warmup_minutes": 30.0,
    }


@pytest.fixture
def image_paths(tmp_path):
    paths = []
    for name in ("green.png", "red.png"):
        path = tmp_path / name
        Image.fromarray(
            np.full((8, 8), 100, dtype=np.uint8)
        ).save(path)
        paths.append(path)
    return paths


def test_manifest_preserves_fractional_times(
    tmp_path, image_paths, build_image_sets
):
    row = _manifest_row(*image_paths)
    row.update(
        sampling_time_hours=0.5,
        measurement_time_hours=0.75,
    )
    manifest = tmp_path / "manifest.csv"
    pd.DataFrame([row]).to_csv(manifest, index=False)

    image_sets = build_image_sets(manifest)

    assert len(image_sets) == 1
    assert image_sets[0].timepoint == 1
    assert image_sets[0].sampling_time_hours == 0.5
    assert image_sets[0].measurement_time_hours == 0.75


@pytest.mark.parametrize("include_empty_columns", [False, True])
def test_manifest_missing_times_remain_unknown(
    tmp_path, image_paths, build_image_sets, include_empty_columns
):
    row = _manifest_row(*image_paths)
    if include_empty_columns:
        row.update(
            sampling_time_hours="",
            measurement_time_hours="",
        )

    manifest = tmp_path / "manifest.csv"
    pd.DataFrame([row]).to_csv(manifest, index=False)

    image_set = build_image_sets(manifest)[0]

    assert image_set.sampling_time_hours is None
    assert image_set.measurement_time_hours is None


@pytest.mark.parametrize(
    "field",
    ["sampling_time_hours", "measurement_time_hours"],
)
@pytest.mark.parametrize("value", [-1.0, float("inf"), "invalid"])
def test_manifest_invalid_time_rejected_before_image_loading(
    tmp_path, build_image_sets, field, value
):
    row = _manifest_row(
        tmp_path / "missing_green.png",
        tmp_path / "missing_red.png",
    )
    row[field] = value

    manifest = tmp_path / "manifest.csv"
    pd.DataFrame([row]).to_csv(manifest, index=False)

    with pytest.raises(ValueError, match=field):
        build_image_sets(manifest)