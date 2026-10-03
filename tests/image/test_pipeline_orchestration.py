"""End-to-end orchestration: process_experiment() over multiple synthetic
image sets should produce one combined CSV covering every timepoint.
"""
import numpy as np
import pandas as pd

from ichnos_image import CrosstalkCalibrationWarning, ImageSet, process_experiment
from ichnos_image.synthesize import make_synthetic_pair


def _make_bf_and_channels(seed, bleed=0.05, true_ratio=0.3):
    size = 150
    rng = np.random.default_rng(seed)
    bf = np.full((size, size), 0.5)
    base_green = np.full((size, size), 50.0)
    yy, xx = np.mgrid[0:size, 0:size]
    for cy, cx in rng.integers(20, size - 20, size=(10, 2)):
        blob = ((yy - cy) ** 2 + (xx - cx) ** 2) < 8**2
        bf[blob] -= 0.3
        base_green[blob] = 3000.0

    green, red, _ = make_synthetic_pair(
        base_green, bleed_green_to_red=bleed, true_ratio_red_green=true_ratio,
        noise_std=5.0, rng=rng,
    )
    return bf, green, red


def test_process_experiment_covers_every_timepoint(tmp_path):
    bleed = 0.05
    image_sets = []
    for timepoint in range(3):
        bf, green, red = _make_bf_and_channels(seed=timepoint, bleed=bleed)
        image_sets.append(
            ImageSet(
                green=green,
                red=red,
                bright_field=bf,
                session_id="session-1",
                timepoint=timepoint,
                acquisition_order=timepoint,
                exposure_ms_green=800.0,
                exposure_ms_red=2000.0,
                nd_filter_green=0.0,
                nd_filter_red=0.0,
                objective="60x",
                burner_hours=5.0,
                lamp_warmup_minutes=30.0,
            )
        )

    out_csv = process_experiment(image_sets, tmp_path / "experiment.csv", bleed_green_to_red=bleed)

    df = pd.read_csv(out_csv)
    assert set(df["timepoint"]) == {0, 1, 2}
    assert len(df) > 0
    assert (df.groupby("timepoint").size() > 0).all()


def test_process_experiment_per_session_calibration(tmp_path):
    """bleed_green_to_red as a {session_id: value} dict applies a different
    calibrated coefficient per session, as a real multi-session experiment
    would need."""
    bleed_by_session = {"session-A": 0.02, "session-B": 0.15}
    image_sets = []
    for session_id, bleed in bleed_by_session.items():
        bf, green, red = _make_bf_and_channels(seed=0, bleed=bleed)
        image_sets.append(
            ImageSet(
                green=green, red=red, bright_field=bf,
                session_id=session_id, timepoint=0, acquisition_order=0,
                exposure_ms_green=800.0, exposure_ms_red=2000.0,
                nd_filter_green=0.0, nd_filter_red=0.0, objective="60x",
                burner_hours=5.0, lamp_warmup_minutes=30.0,
            )
        )

    out_csv = process_experiment(image_sets, tmp_path / "multi_session.csv", bleed_green_to_red=bleed_by_session)
    df = pd.read_csv(out_csv)
    assert set(df["session_id"]) == set(bleed_by_session)


def test_scalar_bleed_across_multiple_sessions_warns(tmp_path):
    """A single bleed_green_to_red reused across distinct sessions should
    warn -- crosstalk calibration is meant to be per-session (a real
    GFP-only control per imaging session), not shared across sessions.
    Reusing it silently across sessions is exactly the mistake the team's
    protocol note (one control per session) is meant to prevent.
    """
    image_sets = []
    for session_id in ("session-A", "session-B"):
        bf, green, red = _make_bf_and_channels(seed=0, bleed=0.05)
        image_sets.append(
            ImageSet(
                green=green, red=red, bright_field=bf,
                session_id=session_id, timepoint=0, acquisition_order=0,
                exposure_ms_green=800.0, exposure_ms_red=2000.0,
                nd_filter_green=0.0, nd_filter_red=0.0, objective="60x",
                burner_hours=5.0, lamp_warmup_minutes=30.0,
            )
        )

    import pytest

    # Asserting on the specific category keeps its negative counterpart
    # (test_scalar_bleed_within_one_session_does_not_warn) meaningful: if the
    # warning were ever dropped, this test fails rather than that one silently
    # passing.
    with pytest.warns(CrosstalkCalibrationWarning, match="per-session"):
        process_experiment(image_sets, tmp_path / "warned.csv", bleed_green_to_red=0.05)


def test_scalar_bleed_within_one_session_does_not_warn(tmp_path):
    image_sets = []
    for timepoint in range(2):
        bf, green, red = _make_bf_and_channels(seed=timepoint, bleed=0.05)
        image_sets.append(
            ImageSet(
                green=green, red=red, bright_field=bf,
                session_id="session-1", timepoint=timepoint, acquisition_order=timepoint,
                exposure_ms_green=800.0, exposure_ms_red=2000.0,
                nd_filter_green=0.0, nd_filter_red=0.0, objective="60x",
                burner_hours=5.0, lamp_warmup_minutes=30.0,
            )
        )

    import warnings

    # Scoped to CrosstalkCalibrationWarning on purpose. simplefilter("error")
    # would also turn unrelated third-party deprecations (e.g. scikit-image)
    # into failures, which says nothing about this function's behaviour.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        process_experiment(image_sets, tmp_path / "not_warned.csv", bleed_green_to_red=0.05)

    crosstalk_warnings = [w for w in caught if issubclass(w.category, CrosstalkCalibrationWarning)]
    assert not crosstalk_warnings, [str(w.message) for w in crosstalk_warnings]


def test_process_experiment_with_rolling_ball_background(tmp_path):
    """background_method="rolling_ball" should run end-to-end, auto-picking
    a radius from each image set's own objective via
    correct.rolling_ball_radius_for_objective() -- the wiring gap flagged
    after the build report (the option existed in correct.py but nothing in
    pipeline.py ever called it)."""
    bleed = 0.05
    bf, green, red = _make_bf_and_channels(seed=0, bleed=bleed)
    image_sets = [
        ImageSet(
            green=green, red=red, bright_field=bf,
            session_id="session-1", timepoint=0, acquisition_order=0,
            exposure_ms_green=800.0, exposure_ms_red=2000.0,
            nd_filter_green=0.0, nd_filter_red=0.0, objective="100X",
            burner_hours=5.0, lamp_warmup_minutes=30.0,
        )
    ]

    out_csv = process_experiment(
        image_sets, tmp_path / "rolling_ball.csv", bleed_green_to_red=bleed, background_method="rolling_ball"
    )
    df = pd.read_csv(out_csv)
    assert len(df) > 0


def test_process_experiment_rolling_ball_unknown_objective_raises(tmp_path):
    bf, green, red = _make_bf_and_channels(seed=0)
    image_sets = [
        ImageSet(
            green=green, red=red, bright_field=bf,
            session_id="session-1", timepoint=0, acquisition_order=0,
            exposure_ms_green=800.0, exposure_ms_red=2000.0,
            nd_filter_green=0.0, nd_filter_red=0.0, objective="40X",  # not in the reference table
            burner_hours=5.0, lamp_warmup_minutes=30.0,
        )
    ]

    import pytest

    with pytest.raises(KeyError):
        process_experiment(
            image_sets, tmp_path / "should_fail.csv", bleed_green_to_red=0.05, background_method="rolling_ball"
        )
