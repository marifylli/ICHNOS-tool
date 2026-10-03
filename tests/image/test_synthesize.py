"""Synthetic dual-channel data generation, and a sanity check that a correct
crosstalk-unmixing step (given the true bleed coefficient) recovers the known
ground-truth red signal -- validates the additive model synthesize.py and
correct.unmix_crosstalk() share.
"""
import numpy as np

from ichnos_image import correct, synthesize


def _synthetic_base(size=96, n_cells=4, seed=0, intensity_jitter=0.0):
    rng = np.random.default_rng(seed)
    green = np.full((size, size), 50.0)
    yy, xx = np.mgrid[0:size, 0:size]
    for cy, cx in rng.integers(15, size - 15, size=(n_cells, 2)):
        blob = ((yy - cy) ** 2 + (xx - cx) ** 2) < 8**2
        intensity = 4000.0 * (1 + rng.uniform(-intensity_jitter, intensity_jitter))
        green[blob] = intensity
    return green


def test_make_synthetic_pair_recovers_known_bleed():
    real_green = _synthetic_base()
    bleed = 0.08
    true_ratio = 0.4
    rng = np.random.default_rng(1)

    green, red, true_red = synthesize.make_synthetic_pair(
        real_green, bleed_green_to_red=bleed, true_ratio_red_green=true_ratio,
        noise_std=5.0, rng=rng,
    )

    green_bg, _ = correct.subtract_background(green, percentile=1.0)
    red_bg, _ = correct.subtract_background(red, percentile=1.0)
    _, red_unmixed = correct.unmix_crosstalk(green_bg, red_bg, bleed_green_to_red=bleed)

    cell_mask = real_green > 1000
    assert np.abs(red_unmixed[cell_mask].mean() - true_red[cell_mask].mean()) < 50.0


def test_estimate_crosstalk_coefficient_recovers_known_bleed():
    real_green = _synthetic_base(size=160, n_cells=10, seed=2)
    bleed = 0.07
    rng = np.random.default_rng(3)

    # pure crosstalk (true_ratio=0) -- the easy case, should recover bleed tightly
    green, red, _ = synthesize.make_synthetic_pair(
        real_green, bleed_green_to_red=bleed, true_ratio_red_green=0.0,
        noise_std=5.0, rng=rng,
    )
    estimated_bleed, _ = correct.estimate_crosstalk_coefficient(green, red)
    assert abs(estimated_bleed - bleed) < 0.02


def test_ratio_jitter_reduces_crosstalk_estimate_contamination():
    """With a single exact ratio shared by every cell, true_red is exactly
    proportional to green -- the same shape as crosstalk -- so no estimator
    can recover bleed alone from a single image; it can only recover
    bleed + true_ratio (this is what an earlier, non-jittered version of
    synthesize.make_synthetic_pair produced, confirmed via
    scripts/validate_crosstalk.py: error == true_ratio, exactly, for every
    nonzero ratio tested). Per-cell jitter breaks that degeneracy. This test
    locks in that jitter measurably helps; it doesn't demand high absolute
    precision, since a ~30-cell toy image is far sparser than the thousands
    of cells in the real reference images scripts/validate_crosstalk.py runs
    against.
    """
    real_green = _synthetic_base(size=200, n_cells=30, seed=2, intensity_jitter=0.5)
    bleed = 0.07
    true_ratio = 0.3

    green_flat, red_flat, _ = synthesize.make_synthetic_pair(
        real_green, bleed_green_to_red=bleed, true_ratio_red_green=true_ratio,
        ratio_jitter_frac=0.0, noise_std=5.0, rng=np.random.default_rng(3),
    )
    estimated_flat, _ = correct.estimate_crosstalk_coefficient(green_flat, red_flat)
    assert abs(estimated_flat - (bleed + true_ratio)) < 0.03

    green_jit, red_jit, _ = synthesize.make_synthetic_pair(
        real_green, bleed_green_to_red=bleed, true_ratio_red_green=true_ratio,
        ratio_jitter_frac=0.5, noise_std=5.0, rng=np.random.default_rng(3),
    )
    estimated_jit, _ = correct.estimate_crosstalk_coefficient(green_jit, red_jit)
    assert abs(estimated_jit - bleed) < abs(estimated_flat - bleed)


def test_calibrate_from_control_then_apply_to_experimental_image():
    """The intended Stage 4 production workflow: calibrate once from a
    genuine GFP-only control (true_ratio=0), then apply that coefficient to
    a *different* experimental image (true_ratio>0) -- rather than
    re-estimating blindly from the experimental image itself, which
    test_ratio_jitter_reduces_crosstalk_estimate_contamination shows is
    unreliable. This should recover the experimental image's true red signal
    well, since bleed_green_to_red is now known rather than estimated
    alongside an unknown ratio.
    """
    control_base = _synthetic_base(size=150, n_cells=20, seed=5, intensity_jitter=0.4)
    bleed = 0.06

    control_green, control_red, _ = synthesize.make_synthetic_pair(
        control_base, bleed_green_to_red=bleed, true_ratio_red_green=0.0,
        noise_std=5.0, rng=np.random.default_rng(10),
    )
    calibration = correct.calibrate_crosstalk_from_control(control_green, control_red)
    assert abs(calibration["bleed_green_to_red"] - bleed) < 0.01
    assert abs(calibration["residual_check"]) < 10.0

    experimental_base = _synthetic_base(size=150, n_cells=20, seed=6, intensity_jitter=0.4)
    true_ratio = 0.5
    exp_green, exp_red, exp_true_red = synthesize.make_synthetic_pair(
        experimental_base, bleed_green_to_red=bleed, true_ratio_red_green=true_ratio,
        ratio_jitter_frac=0.3, noise_std=5.0, rng=np.random.default_rng(11),
    )
    _, exp_red_unmixed = correct.unmix_crosstalk(
        exp_green, exp_red, bleed_green_to_red=calibration["bleed_green_to_red"]
    )

    cell_mask = experimental_base > 1000
    recovered = exp_red_unmixed[cell_mask].mean()
    truth = exp_true_red[cell_mask].mean()
    assert abs(recovered - truth) / truth < 0.15


def test_build_dataset(tmp_path):
    real_green = _synthetic_base()
    base_path = tmp_path / "base_green.png"
    synthesize.save_png16(real_green, base_path)

    labels_path = synthesize.build_dataset(
        {"toy": base_path},
        tmp_path / "synthetic_out",
        bleed_values=(0.05, 0.1),
        true_ratio_values=(0.0, 0.5),
        seed=0,
    )

    import pandas as pd
    df = pd.read_csv(labels_path)
    assert len(df) == 4  # 1 base x 2 bleed x 2 ratio
    for row in df.itertuples():
        assert (labels_path.parent / row.green_file).exists()
        assert (labels_path.parent / row.red_file).exists()
