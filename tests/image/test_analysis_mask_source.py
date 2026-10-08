"""The session analysis script and the pipeline must agree on which image the
masks are cut from, or a disagreement between the two paths is a difference in
mask source rather than independent verification of anything."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from ichnos_image import pipeline


@pytest.fixture(scope="module")
def script():
    path = (
        Path(__file__).resolve().parents[2]
        / "analysis" / "h2o2-2026-10-05" / "analyze_h2o2_session.py"
    )
    spec = importlib.util.spec_from_file_location("analyze_h2o2_under_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def field():
    """One population bright in green, one bright in red, on a dark level."""
    rng = np.random.default_rng(0)
    shape = (160, 160)
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    green = 40 + rng.normal(0, 1.5, shape)
    red = 40 + rng.normal(0, 1.5, shape)
    for cy, cx in [(40, 40), (40, 80), (40, 120)]:
        green[(yy - cy) ** 2 + (xx - cx) ** 2 <= 81] += 60
    for cy, cx in [(120, 40), (120, 80), (120, 120)]:
        red[(yy - cy) ** 2 + (xx - cx) ** 2 <= 81] += 60
    flat = np.ones(shape)
    return dict(green=green, red=red, dark=10.0, flat_g=flat, flat_r=flat)


def source(script, field, choice):
    return script.segmentation_source(
        field["green"], field["red"], field["dark"], field["flat_g"], field["flat_r"], choice
    )


def test_the_default_is_the_symmetric_one(script):
    """Not 'red', which is what the script did before: a red mask lets mCherry
    pick the cells whose red/green ratio is then reported."""
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--segmentation-source", choices=script.SEGMENTATION_SOURCES,
                        default="sum")
    assert parser.parse_args([]).segmentation_source == "sum"


def test_the_choices_match_the_pipelines(script):
    """Same words, so a run of each can be compared without a translation
    table. The pipeline additionally offers brightfield; this session has none.
    """
    assert set(script.SEGMENTATION_SOURCES) <= set(pipeline.SEGMENTATION_SOURCES)


def test_the_source_name_is_reported_not_inferred(script, field):
    assert source(script, field, "green")[1] == "green"
    assert source(script, field, "red")[1] == "red"
    assert source(script, field, "sum")[1] == "green+red"


def test_an_unknown_source_is_rejected(script, field):
    with pytest.raises(ValueError, match="segmentation source"):
        source(script, field, "gren")


def test_the_dark_level_is_removed_before_summing(script, field):
    """Summing two uncorrected channels would add the dark level twice and
    shift the high-pass threshold, so the correction has to come first."""
    summed, _ = source(script, field, "sum")
    green, _ = source(script, field, "green")
    red, _ = source(script, field, "red")
    assert np.allclose(summed, green + red)
    assert summed.min() < field["green"].min() + field["red"].min() - field["dark"]


def test_each_single_channel_source_finds_only_its_own_population(script, field):
    """The mechanism of the bias, stated as a test: this is why the ratio moves
    with the mask source."""
    top = lambda labels: (labels[:80] > 0).sum()  # noqa: E731 - green population
    bottom = lambda labels: (labels[80:] > 0).sum()  # noqa: E731 - red population
    kwargs = dict(k=3.0, min_area=60, max_area=5000, smooth=1.5)

    from_green, _, _ = script.segment(source(script, field, "green")[0], **kwargs)
    from_red, _, _ = script.segment(source(script, field, "red")[0], **kwargs)
    from_sum, _, _ = script.segment(source(script, field, "sum")[0], **kwargs)

    assert top(from_green) > 10 * max(bottom(from_green), 1)
    assert bottom(from_red) > 10 * max(top(from_red), 1)
    assert min(top(from_sum), bottom(from_sum)) > 0
