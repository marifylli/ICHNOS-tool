"""Stage 1 flat-field estimation: recovers a known illumination profile from
a stack of images with different, randomly-placed cell content.
"""
import numpy as np

from ichnos_image import correct


def _vignetting_profile(size=150):
    yy, xx = np.mgrid[0:size, 0:size]
    cy, cx = size / 2, size / 2
    r2 = (yy - cy) ** 2 + (xx - cx) ** 2
    profile = 1.0 - 0.5 * (r2 / r2.max())  # bright center, dim corners
    return profile


def _cell_field(size=150, n_cells=15, seed=0):
    rng = np.random.default_rng(seed)
    field = np.full((size, size), 100.0)
    yy, xx = np.mgrid[0:size, 0:size]
    for cy, cx in rng.integers(10, size - 10, size=(n_cells, 2)):
        blob = ((yy - cy) ** 2 + (xx - cx) ** 2) < 6**2
        field[blob] = 3000.0
    return field


def test_estimate_flat_field_recovers_known_vignetting():
    true_profile = _vignetting_profile()
    images = [_cell_field(seed=i) * true_profile for i in range(8)]

    estimated = correct.estimate_flat_field(images, smooth_sigma=20.0)

    true_normalized = true_profile / true_profile.mean()
    # compare shape, not exact scale/edge effects from smoothing
    correlation = np.corrcoef(estimated.ravel(), true_normalized.ravel())[0, 1]
    assert correlation > 0.9


def test_flat_field_correct_flattens_illumination():
    true_profile = _vignetting_profile()
    flat_field_only = np.full((150, 150), 500.0) * true_profile  # no cells, pure shading

    images = [_cell_field(seed=i) * true_profile for i in range(8)]
    estimated = correct.estimate_flat_field(images, smooth_sigma=20.0)

    corrected = correct.flat_field_correct(flat_field_only, estimated)
    # a flat-field-only image should become much more uniform after correction
    uncorrected_cv = flat_field_only.std() / flat_field_only.mean()
    corrected_cv = corrected.std() / corrected.mean()
    assert corrected_cv < uncorrected_cv * 0.5


def test_estimate_flat_field_requires_multiple_images():
    import pytest

    with pytest.raises(ValueError):
        correct.estimate_flat_field([_cell_field()])


def _field_with_density(n_cells, size=150, true_bg=200.0, seed=0):
    rng = np.random.default_rng(seed)
    field = rng.normal(true_bg, 8, size=(size, size))
    yy, xx = np.mgrid[0:size, 0:size]
    for cy, cx in rng.integers(10, size - 10, size=(n_cells, 2)):
        blob = ((yy - cy) ** 2 + (xx - cx) ** 2) < 6**2
        field[blob] = 4000.0
    return field


def test_raw_histogram_mode_breaks_at_high_cell_density():
    """Documents *why* subtract_background needs the sanity-percentile
    fallback: the bare mode estimator alone jumps to the cell intensity
    once cells cover most of the field."""
    dense_field = _field_with_density(n_cells=300)  # ~67% cell coverage
    mode_level = correct._histogram_mode_background(dense_field)
    assert mode_level > 1000  # latched onto the cell peak (~4000), not background (~200)


def test_subtract_background_mode_falls_back_at_high_density():
    true_bg = 200.0
    sparse = _field_with_density(n_cells=5, true_bg=true_bg)
    dense = _field_with_density(n_cells=300, true_bg=true_bg)

    _, bg_sparse = correct.subtract_background(sparse, method="mode")
    _, bg_dense = correct.subtract_background(dense, method="mode")

    assert abs(bg_sparse - true_bg) < 15
    assert abs(bg_dense - true_bg) < 15  # fallback kept this near truth, not ~4000


def test_rolling_ball_handles_spatial_gradient_better_than_scalar_mode():
    """method="mode"/"percentile" fit one scalar for the whole image, which
    can't track a spatially-varying background (e.g. an illumination
    gradient) -- method="rolling_ball" (skimage.restoration.rolling_ball,
    the algorithm behind ImageJ/FIJI's "Subtract Background") estimates a
    full per-pixel background map instead, and should leave a much smaller,
    more even residual across the gradient.
    """
    size = 200
    yy, xx = np.mgrid[0:size, 0:size]
    true_bg = 50 + (xx / size) * 250.0  # gradient: 50 (left) -> 300 (right)

    rng = np.random.default_rng(0)
    image = true_bg.copy()
    cell_mask = np.zeros((size, size), dtype=bool)
    for cy, cx in rng.integers(15, size - 15, size=(15, 2)):
        blob = ((yy - cy) ** 2 + (xx - cx) ** 2) < 8**2
        image[blob] += 3000
        cell_mask |= blob
    bg_only = ~cell_mask

    left = np.zeros((size, size), dtype=bool)
    left[:, : size // 4] = True
    right = np.zeros((size, size), dtype=bool)
    right[:, -size // 4 :] = True

    corr_mode, _ = correct.subtract_background(image, method="mode")
    corr_rb, _ = correct.subtract_background(image, method="rolling_ball", rolling_ball_radius=60)

    mode_spread = abs(corr_mode[bg_only & left].mean() - corr_mode[bg_only & right].mean())
    rb_spread = abs(corr_rb[bg_only & left].mean() - corr_rb[bg_only & right].mean())

    assert rb_spread < mode_spread * 0.5


def test_suggest_rolling_ball_radius_scales_with_pixel_size():
    # smaller pixel size (higher magnification) -> a cell spans more pixels -> larger radius needed
    radius_100x = correct.suggest_rolling_ball_radius(pixel_size_um=0.13)
    radius_60x = correct.suggest_rolling_ball_radius(pixel_size_um=0.1076)
    assert radius_60x > radius_100x > 0


def test_rolling_ball_radius_for_objective_uses_config_lookup():
    from ichnos.config import OBJECTIVE_PIXEL_SIZE_UM_REFERENCE

    radius = correct.rolling_ball_radius_for_objective("100X")
    expected = correct.suggest_rolling_ball_radius(OBJECTIVE_PIXEL_SIZE_UM_REFERENCE["100X"])
    assert radius == expected

    # case-insensitive: manifests/microscope software aren't consistently cased
    assert correct.rolling_ball_radius_for_objective("100x") == expected

    import pytest

    with pytest.raises(KeyError):
        correct.rolling_ball_radius_for_objective("unknown-objective")


def test_undersized_rolling_ball_radius_eats_cell_signal():
    """Validates the formula empirically (mirrors
    scripts/calibrate_rolling_ball_radius.py's real-image finding): a radius
    well below the suggested value lets the ball dip into cells and erode
    their signal, while the suggested radius has already captured nearly all
    of the signal a much larger radius would.
    """
    size = 300
    pixel_size_um = 0.13
    cell_diameter_um = 5.0
    rng = np.random.default_rng(0)

    true_bg = 200.0
    image = np.full((size, size), true_bg)
    yy, xx = np.mgrid[0:size, 0:size]
    cell_radius_px = (cell_diameter_um / 2) / pixel_size_um
    cell_mask = np.zeros((size, size), dtype=bool)
    for cy, cx in rng.integers(30, size - 30, size=(20, 2)):
        blob = ((yy - cy) ** 2 + (xx - cx) ** 2) < cell_radius_px**2
        image[blob] += 3000
        cell_mask |= blob

    suggested = correct.suggest_rolling_ball_radius(pixel_size_um, cell_diameter_um=cell_diameter_um)

    corr_small, _ = correct.subtract_background(image, method="rolling_ball", rolling_ball_radius=suggested * 0.25)
    corr_suggested, _ = correct.subtract_background(image, method="rolling_ball", rolling_ball_radius=suggested)
    corr_large, _ = correct.subtract_background(image, method="rolling_ball", rolling_ball_radius=suggested * 3)

    signal_small = corr_small[cell_mask].mean()
    signal_suggested = corr_suggested[cell_mask].mean()
    signal_large = corr_large[cell_mask].mean()

    assert signal_small < signal_suggested * 0.9  # too-small radius measurably erodes signal
    assert signal_suggested > signal_large * 0.95  # suggested radius is already near the plateau


def test_pixel_size_at_sample_um_matches_sc30_manual_formula():
    # SC30 manual: pixel_size_um = sensor_pixel_size_um / (objective_mag * adapter_mag)
    assert correct.pixel_size_at_sample_um(60, 1.0, sensor_pixel_size_um=3.2) == 3.2 / 60
    assert correct.pixel_size_at_sample_um(100, 0.5, sensor_pixel_size_um=3.2) == 3.2 / 50

    # defaults to the SC30's own sensor pixel size
    from ichnos.config import CAMERA_SC30_SENSOR_PIXEL_SIZE_UM

    assert correct.pixel_size_at_sample_um(40) == CAMERA_SC30_SENSOR_PIXEL_SIZE_UM / 40


def test_pixel_size_at_sample_um_scales_with_binning():
    no_binning = correct.pixel_size_at_sample_um(20, 0.5, sensor_pixel_size_um=3.2, binning_factor=1)
    binned_2x = correct.pixel_size_at_sample_um(20, 0.5, sensor_pixel_size_um=3.2, binning_factor=2)
    assert binned_2x == no_binning * 2
