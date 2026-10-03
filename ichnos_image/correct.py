"""Per-channel corrections (Stages 1/3/4): background subtraction, flat-field
illumination, channel co-registration, crosstalk/bleed-through unmixing, and
photobleaching correction.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import shift as ndi_shift
from scipy.stats import theilslopes
from skimage.registration import phase_cross_correlation

from ichnos.config import (
    BACKGROUND_MODE_SANITY_PERCENTILE,
    CAMERA_SC30_SENSOR_PIXEL_SIZE_UM,
    OBJECTIVE_PIXEL_SIZE_UM_REFERENCE,
    YEAST_CELL_DIAMETER_UM,
)


def pixel_size_at_sample_um(
    objective_magnification: float,
    camera_adapter_magnification: float = 1.0,
    sensor_pixel_size_um: float = CAMERA_SC30_SENSOR_PIXEL_SIZE_UM,
    binning_factor: int = 1,
) -> float:
    """Real µm/pixel at the sample, from the camera sensor's own pixel size
    and the optical magnification in front of it -- straight from the SC30
    manual's formula: pixel_size_um = sensor_pixel_size_um / (objective_mag *
    adapter_mag).

    Defaults to the SC30's sensor pixel size (3.2um), the team's real
    camera. camera_adapter_magnification should be one of
    ichnos.config.KNOWN_CAMERA_ADAPTER_MAGNIFICATIONS (0.25-1.0, per the
    CKX41/Camera Adapter System manuals) once the team confirms which C-mount
    adapter they actually use -- printed on the adapter barrel itself.
    binning_factor (1/2/3/4, see ichnos.config.CAMERA_SC30_BINNING_MODES)
    scales the effective sensor pixel size: at Nx binning, N physical pixels
    combine into one, so it's N times larger (and resolution/exposure change
    accordingly -- see that dict).

    This is the number to feed into suggest_rolling_ball_radius() or a new
    OBJECTIVE_PIXEL_SIZE_UM_REFERENCE entry once the team confirms objective
    + adapter + binning mode -- OBJECTIVE_PIXEL_SIZE_UM_REFERENCE currently
    holds unrelated placeholder values from the public YRC dataset's
    (different) camera, not this one.
    """
    effective_sensor_pixel_size_um = sensor_pixel_size_um * binning_factor
    return effective_sensor_pixel_size_um / (objective_magnification * camera_adapter_magnification)


def suggest_rolling_ball_radius(
    pixel_size_um: float, cell_diameter_um: float = YEAST_CELL_DIAMETER_UM, safety_factor: float = 2.5
) -> float:
    """Suggest a rolling_ball radius (px) for subtract_background(method="rolling_ball"),
    from the objective's pixel size rather than a fixed default.

    The ball's radius has to stay well above the cell's radius in *pixels*,
    not micrometers, or it dips into cells and eats real signal. Pixel size
    depends on the objective + camera sensor, so the right radius is
    objective-specific -- compute it once per objective you use and lock it
    in, don't reuse one global default across a 40x and a 100x session.

    safety_factor=2.5 (cell radius x this) is a reasonable starting point,
    not a theoretically exact multiple -- there's no principled "correct"
    value, only an empirical one. Verify it on a few real images per
    objective (see scripts/calibrate_rolling_ball_radius.py, which does
    exactly that on the real reference images this project has so far) once
    real data for that objective exists, rather than trusting the formula
    alone.
    """
    cell_radius_px = (cell_diameter_um / 2) / pixel_size_um
    return cell_radius_px * safety_factor


def rolling_ball_radius_for_objective(objective: str, **kwargs) -> float:
    """suggest_rolling_ball_radius(), looked up by objective name against
    ichnos.config.OBJECTIVE_PIXEL_SIZE_UM_REFERENCE (case-insensitive, e.g.
    "60x" and "60X" both match -- objective naming isn't consistently
    cased across manifests/microscope software).

    Those reference pixel sizes come from the public YRC dataset (a
    different microscope), not the team's Olympus -- this exists so the
    mechanism (objective -> radius) is ready to wire in real Olympus values
    the moment they're known, not as a source of real Olympus numbers today.
    """
    normalized = {key.upper(): value for key, value in OBJECTIVE_PIXEL_SIZE_UM_REFERENCE.items()}
    pixel_size_um = normalized.get(objective.upper())
    if pixel_size_um is None:
        raise KeyError(
            f"no pixel size on record for objective {objective!r}; add it to "
            "ichnos.config.OBJECTIVE_PIXEL_SIZE_UM_REFERENCE, or call "
            "suggest_rolling_ball_radius(pixel_size_um=...) directly"
        )
    return suggest_rolling_ball_radius(pixel_size_um, **kwargs)


def subtract_background(
    channel: np.ndarray,
    method: str = "mode",
    percentile: float = 5.0,
    n_bins: int = 256,
    sanity_percentile: float = BACKGROUND_MODE_SANITY_PERCENTILE,
    rolling_ball_radius: float = 100.0,
) -> tuple[np.ndarray, float]:
    """Estimate + remove background/dark noise for one channel (Stage 3),
    per-image and dynamic (not a fixed value carried over from a control),
    since autofluorescence shifts with metabolic state.

    method="mode" (default): the background level is a single scalar, the
    peak of a smoothed pixel-intensity histogram -- more principled than an
    arbitrary percentile cutoff when the true background isn't simply "the
    darkest N% of pixels" (e.g. it's skewed, or its width varies with
    metabolic state). It has its own failure mode though: in a dense field
    where cell pixels start to outnumber background pixels, the *cells*
    become the histogram's tallest peak and the estimate jumps to a
    cell-intensity level, badly wrong -- confirmed empirically at ~60%+ cell
    coverage on synthetic test fields (tests/test_correct.py). To guard
    against that without giving up mode's accuracy in normal (low/moderate-
    density, realistic-for-yeast) fields, the mode estimate is only trusted
    if it falls at or below the image's `sanity_percentile`-th percentile;
    otherwise this falls back to the plain percentile estimate
    (method="percentile" forces that path always).

    method="rolling_ball": skimage.restoration.rolling_ball, the standard
    algorithm behind ImageJ/FIJI's "Subtract Background" -- estimates a
    *spatially-varying* background (not one scalar), rolling a ball of
    `rolling_ball_radius` under the intensity landscape. Handles uneven
    illumination within a single image directly, unlike the other two
    methods here or estimate_flat_field() (which needs a stack of several
    same-session images). Assumes dark pixels are background, which holds
    for fluorescence channels (GFP/mCherry: background is dark, cells are
    bright). radius must be well above the largest cell's radius in pixels,
    or the ball dips into cells and eats real signal -- default 100 is
    comfortably above yeast cell sizes (~10-40px radius) at the pixel sizes
    in the team's reference images, but re-check if the real objective/
    binning gives much larger cells in pixels.

    Returns (corrected_channel, estimated_background_level) -- for
    method="rolling_ball" the second value is the *mean* of the per-pixel
    background map (a summary for logging/QC), even though the actual
    subtraction used the full spatially-varying map.
    """
    percentile_level = float(np.percentile(channel, percentile))

    if method == "percentile":
        background_level = percentile_level
        corrected = np.clip(channel.astype(float) - background_level, 0, None)
    elif method == "mode":
        mode_level = _histogram_mode_background(channel, n_bins=n_bins)
        sanity_bound = float(np.percentile(channel, sanity_percentile))
        background_level = mode_level if mode_level <= sanity_bound else percentile_level
        corrected = np.clip(channel.astype(float) - background_level, 0, None)
    elif method == "rolling_ball":
        from skimage.restoration import rolling_ball

        background_map = rolling_ball(channel.astype(float), radius=rolling_ball_radius)
        background_level = float(background_map.mean())
        corrected = np.clip(channel.astype(float) - background_map, 0, None)
    else:
        raise ValueError(f"unknown background method: {method!r}")

    return corrected, background_level


def _histogram_mode_background(channel: np.ndarray, n_bins: int = 256) -> float:
    from scipy.ndimage import gaussian_filter1d

    values = channel.ravel().astype(float)
    counts, edges = np.histogram(values, bins=n_bins)
    smoothed_counts = gaussian_filter1d(counts.astype(float), sigma=2)
    peak_bin = int(np.argmax(smoothed_counts))
    return float((edges[peak_bin] + edges[peak_bin + 1]) / 2)


def estimate_flat_field(images: list[np.ndarray], smooth_sigma: float = 50.0) -> np.ndarray:
    """Retrospectively estimate an illumination/vignetting profile from a
    stack of same-session, same-channel images (Stage 1).

    Not a full CIDRE (Smith et al. 2015) or BaSiC (Peng et al. 2017) fit --
    those decompose the stack into low-rank illumination + sparse
    foreground via iterative optimization. This is the simple version of the
    same idea: cells sit at different positions in each image, so a robust
    per-pixel statistic (median) across the stack averages most cell content
    out, leaving the illumination pattern common to all of them; a large
    Gaussian blur then removes whatever high-frequency cell signal survived
    the median, keeping only the smooth, low-frequency shading. Good enough
    as a stand-in, but needs more images than CIDRE/BaSiC to be reliable
    (those use more sophisticated fits precisely to work with fewer) -- with
    very few images or dense/uneven cell coverage, prefer a dedicated
    prospective flat-field slide image instead.

    Returns a mean-normalized flat field; pass straight into
    flat_field_correct().
    """
    if len(images) < 2:
        raise ValueError("need at least 2 images to estimate a shared illumination profile")
    stack = np.stack([im.astype(float) for im in images], axis=0)
    median_image = np.median(stack, axis=0)
    from skimage.filters import gaussian

    smoothed = gaussian(median_image, sigma=smooth_sigma)
    return smoothed / np.mean(smoothed)


def flat_field_correct(channel: np.ndarray, flat_field: np.ndarray) -> np.ndarray:
    """Divide out an illumination/vignetting profile (CIDRE/BaSiC-style flat field, Stage 1).

    flat_field is a per-pixel illumination reference (mean-normalized here),
    from estimate_flat_field() on a same-session image stack, or a
    prospective flat-field slide image if too few images are available.
    """
    flat = flat_field.astype(float)
    flat = flat / np.mean(flat)
    return channel.astype(float) / np.clip(flat, 1e-6, None)


def estimate_registration_shift(
    reference_channel: np.ndarray, moving_channel: np.ndarray, upsample_factor: int = 10
) -> tuple[float, float]:
    """Sub-pixel (row, col) shift of moving_channel relative to reference_channel
    via phase cross-correlation (Stage 1 co-registration / Stage 7 registration QC).
    """
    shift_rc, *_ = phase_cross_correlation(reference_channel, moving_channel, upsample_factor=upsample_factor)
    return float(shift_rc[0]), float(shift_rc[1])


def apply_shift(channel: np.ndarray, shift_rc: tuple[float, float]) -> np.ndarray:
    """Apply a (row, col) sub-pixel shift, e.g. from estimate_registration_shift."""
    return ndi_shift(channel.astype(float), shift=shift_rc, order=1, mode="nearest")


def unmix_crosstalk(green: np.ndarray, red: np.ndarray, bleed_green_to_red: float) -> tuple[np.ndarray, np.ndarray]:
    """Linear unmixing for GFP -> mCherry channel bleed-through (Stage 4).

    bleed_green_to_red is the fraction of green signal leaking into the red
    detection channel. Get it from calibrate_crosstalk_from_control() run on
    a genuine single-fluorophore GFP-only control image -- see that
    function's docstring for why that is the only reliable source for a
    tandem-timer reporter (not a per-experimental-image estimate). Extend to
    a full 2x2 mixing matrix if red -> green leak also proves non-negligible.
    """
    red_unmixed = np.clip(red.astype(float) - bleed_green_to_red * green.astype(float), 0, None)
    return green.astype(float), red_unmixed


def estimate_crosstalk_coefficient(
    green: np.ndarray,
    red: np.ndarray,
    n_bins: int = 50,
    quantile: float = 0.05,
    min_pixels_per_bin: int = 50,
) -> tuple[float, float]:
    """Estimate the green -> red bleed-through coefficient from a single dual-
    channel image, without a dedicated GFP-only control (Stage 4).

    Model: observed_red = true_red + bleed*green + noise, with true_red >= 0
    always. The key assumption is that true_red is not present at every
    pixel for every green intensity -- i.e. within each green-intensity bin,
    at least a few pixels/cells are crosstalk-only (true_red ~= 0). Under
    that assumption the *low quantile* of red within each bin traces the
    crosstalk floor, and Theil-Sen regression of that floor against green
    gives a robust (bleed, intercept) estimate.

    Validated on synthesize.build_dataset's 75-sample synthetic set (real
    GFP images, per-cell-independent true ratios): essentially exact
    (MAE ~0.0001) when true_ratio_red_green == 0 -- i.e. on a genuine
    single-fluorophore (GFP-only) control image, which is the intended use.
    On images with real red signal present, the estimate is biased upward
    and the bias grows with the amount of true signal (empirically, roughly
    80-90% of true_ratio_red_green leaks into the bleed estimate even with
    substantial per-cell ratio variation) -- the "enough independent cells
    at every green level" assumption is optimistic for typical field-of-view
    cell counts and density. Re-run scripts/validate_crosstalk.py after any
    change here; do not call this on a real dual-fluorophore experimental
    image and trust the result -- calibrate from an actual GFP-only control
    via calibrate_crosstalk_from_control() instead. For the team's tandem-
    timer reporter specifically, no blind per-image method (NMF, SSASU,
    or otherwise) can substitute for a control: true red and green signal
    come from the same fused molecule at the same location, so true_red is,
    by construction, exactly proportional to green within a cell -- the
    same degenerate case demonstrated above, not something a fancier
    algorithm can see past without outside information. A control-free
    blind-separation approach only has a chance when the two channels track
    genuinely different, spatially-independent proteins (as in the public
    reference dataset's GFP-target + mCherry-NUP49 pairs), which is not the
    team's real experimental setup.

    Returns (bleed_green_to_red, intercept); pass the first straight into
    unmix_crosstalk().
    """
    green_flat = green.ravel().astype(float)
    red_flat = red.ravel().astype(float)

    # Equal-width bins over the green intensity range, not equal-population
    # (quantile) bins: real images are often bimodal (background vs. cells),
    # and quantile bins would pack dozens of near-duplicate points into
    # whichever cluster holds most pixels (usually background), swamping the
    # few points from the other cluster in the robust (Theil-Sen) fit below.
    bin_edges = np.unique(np.linspace(green_flat.min(), green_flat.max(), n_bins + 1))
    if len(bin_edges) < 3:
        raise ValueError("not enough distinct green intensities to bin")

    bin_idx = np.digitize(green_flat, bin_edges[1:-1])
    xs, ys = [], []
    for b in range(len(bin_edges) - 1):
        sel = bin_idx == b
        if sel.sum() < min_pixels_per_bin:
            continue
        xs.append(green_flat[sel].mean())
        ys.append(float(np.quantile(red_flat[sel], quantile)))
    if len(xs) < 2:
        raise ValueError("not enough populated intensity bins to fit a crosstalk line")

    slope, intercept, _lo, _hi = theilslopes(ys, xs)
    return max(float(slope), 0.0), float(intercept)


def calibrate_crosstalk_from_control(
    control_green: np.ndarray, control_red: np.ndarray, cell_percentile: float = 90.0, **kwargs
) -> dict:
    """Stage 4 production entrypoint: calibrate the green -> red bleed
    coefficient once per imaging session from a genuine single-fluorophore
    (GFP-only, no mCherry) control image, for use on every experimental
    image from that session via unmix_crosstalk().

    This is the only reliable calibration route for the team's tandem-timer
    reporter -- see estimate_crosstalk_coefficient()'s docstring for why no
    blind per-image method can substitute for it: true red signal in a real
    experimental image is, by construction, exactly proportional to green
    within each cell (same fused molecule, same location), which is
    mathematically indistinguishable from crosstalk without a control.

    Returns a dict with:
    - bleed_green_to_red, intercept: pass bleed_green_to_red into
      unmix_crosstalk() for every experimental image from this session.
    - residual_check: mean unmixed-red level in this control's own brightest
      (cell) pixels after applying the estimated coefficient back to it --
      should be close to 0 (the control has no real red signal by
      definition). A large residual is a warning sign that this image was
      not actually a clean GFP-only control (e.g. mCherry leaked in from
      somewhere, or it was accidentally an experimental image).
    """
    bleed, intercept = estimate_crosstalk_coefficient(control_green, control_red, **kwargs)
    _, red_unmixed = unmix_crosstalk(control_green, control_red, bleed_green_to_red=bleed)

    cell_mask = control_green > np.percentile(control_green, cell_percentile)
    residual_check = float(red_unmixed[cell_mask].mean() if cell_mask.any() else red_unmixed.mean())

    return {
        "bleed_green_to_red": bleed,
        "intercept": intercept,
        "residual_check": residual_check,
    }


def correct_photobleaching(intensity_series: np.ndarray, timepoints: np.ndarray) -> np.ndarray:
    """Fit + remove a per-session exponential decay from a channel's time series (Stage 4/6).

    Falls back to returning the series unchanged if the fit doesn't converge
    (e.g. too few timepoints), rather than raising mid-pipeline.
    """
    from scipy.optimize import curve_fit

    def decay(t, a, k, c):
        return a * np.exp(-k * t) + c

    try:
        params, _ = curve_fit(decay, timepoints, intensity_series, p0=(intensity_series[0], 1e-3, 0), maxfev=5000)
    except RuntimeError:
        return intensity_series
    baseline = decay(timepoints, *params)
    return intensity_series / np.clip(baseline / baseline[0], 1e-6, None)
