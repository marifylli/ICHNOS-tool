"""bright-field/DIC -> per-cell label masks (Stage 2), plus a couple of
Stage-2/7 helpers (border-touching cells, focus score) that naturally operate
on the same phase/DIC image used for segmentation.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from skimage import filters, morphology, segmentation, measure
from skimage.feature import peak_local_max
from inspect import signature

_REMOVE_SMALL_OBJECTS_HAS_MAX_SIZE = (
    "max_size"
    in signature(morphology.remove_small_objects).parameters
)


def _remove_objects_below_size(
    binary: np.ndarray, min_size: int
) -> np.ndarray:
    """Keep connected objects containing at least min_size pixels."""
    if min_size <= 0:
        return binary.copy()

    if _REMOVE_SMALL_OBJECTS_HAS_MAX_SIZE:
        return morphology.remove_small_objects(
            binary, max_size=min_size - 1
        )

    # Compatibility with scikit-image versions before 0.26.
    return morphology.remove_small_objects(
        binary, min_size=min_size
    )

def foreground_fraction(mask: np.ndarray) -> float:
    """Share of the frame a segmentation claims as cells (Stage 2 QC).

    Accepts a boolean mask or an integer label mask. The point of measuring
    it is that the two ways segmentation fails on a real fluorescence frame
    look nothing alike in the output but are both obvious here: a threshold
    that lands inside the background claims tens of percent of the frame and
    the watershed then shatters it into thousands of cell-sized fragments,
    while a threshold above the signal claims almost nothing and silently
    returns a handful of the brightest cells. Object count alone separates
    neither case from a real field -- the fragments are cell-sized, and a
    sparse field genuinely has few cells.
    """
    return float((np.asarray(mask) > 0).mean())


def segment_cells(
    bf_image: np.ndarray, method: str = "otsu", min_size: int = 30, resize_factor: float = 1.0, **cellpose_kwargs
) -> np.ndarray:
    """Return an integer label mask (0 = background) for a bright-field/DIC image.

    method="sparse" is for fluorescence frames, where cells occupy a few
    percent of the pixels and the histogram has one mode (background), not
    two. Both other methods assume otherwise; see _segment_sparse.

    method="cellpose" needs the optional cellpose+torch dependencies (deep
    learning segmentation, preferred for stress-distorted morphology per the
    team's protocol) -- ~110s/image on CPU for a real 535x512 DIC image, so
    batch runs may need either GPU or a resize_factor < 1.0 (see below --
    unlike "otsu", cellpose tolerates this well). See
    scripts/estimate_cellpose_throughput.py to turn that per-image number
    into a "how long will my real batch take" estimate.

    method="transmitted" is for bright-field frames, where cells are dark
    rims with bright halos rather than bright objects; see
    _segment_transmitted. "sparse" on such a frame returns the halos.

    method="otsu" (default) is a classical Otsu + watershed fallback with no
    heavy dependencies, good enough for tests/prototyping.

    resize_factor < 1.0 downsamples before segmenting and scales the label
    mask back up (nearest-neighbor, so label ids stay intact). Measured
    effect differs sharply by method -- always check
    scripts/evaluate_resize_accuracy_tradeoff.py / a real-image cell-count
    sanity check for your own case before relying on it, don't assume either
    number below transfers:
    - "cellpose": degrades gracefully on a real DIC image -- resize_factor=0.5
      kept 75/79 cells (95%) at 2.6x the speed (42s vs 110s); 0.25 kept
      65/79 (82%) at 5.7x the speed (19s). A real, usable throughput lever.
    - "otsu": breaks down hard well before that -- F1 against synthetic
      ground truth collapsed from 0.83 (resize_factor=1.0) to 0.33 at 0.75
      and to 0.0 by 0.5, because its edge-detection kernels (segment.py's
      internal sigma values) are fixed in *pixels* and stop resolving
      nearby cell boundaries once the image shrinks -- they don't scale
      with resize_factor. Don't downsample "otsu" runs; it isn't the
      resize/upsample bookkeeping that fails (label positions stay
      correct), it's this method's fixed-scale kernels.
    """
    if resize_factor != 1.0:
        return _segment_resized(bf_image, method=method, min_size=min_size, resize_factor=resize_factor, **cellpose_kwargs)
    if method == "cellpose":
        return _segment_cellpose(bf_image, **cellpose_kwargs)
    if method == "otsu":
        return _segment_otsu(bf_image, min_size=min_size)
    if method == "sparse":
        return _segment_sparse(bf_image, min_size=min_size, **cellpose_kwargs)
    if method == "transmitted":
        return _segment_transmitted(bf_image, min_size=min_size, **cellpose_kwargs)
    raise ValueError(f"unknown segmentation method: {method!r}")


def _segment_transmitted(
    image: np.ndarray,
    *,
    min_size: int = 200,
    texture_window_px: int = 9,
    noise_sigmas: float = 2.5,
    closing_px: int = 2,
    erosion_px: int = 8,
) -> np.ndarray:
    """Segment a bright-field frame, where cells are neither bright nor dark.

    A yeast cell in transmitted light is a dark rim with a bright halo on a
    bright, unevenly lit background. Its interior is often the same grey as
    the background, so thresholding intensity in either direction finds rims
    and halos rather than cells, and "sparse" -- which looks for bright
    blobs -- returns the halos. What separates cell from background here is
    not brightness but structure: the background is smooth and the cells are
    not, at the scale of a few pixels. This thresholds local standard
    deviation, which measures exactly that and is indifferent to the slow
    illumination gradient across the frame.

    The masks come out wider than the cells, because the halo is textured
    too. On synthetic cells of known size the raw mask was three times their
    area; erosion brings it back to roughly cell-sized, which matters when
    the mask is used to register against another image -- it sharpened the
    recovered displacement from (+20, -20) to the true (+17, -23). Erosion
    is backed off rather than applied blindly, since a cell smaller than the
    radius asked for would simply disappear.
    """
    if texture_window_px < 2 or noise_sigmas <= 0 or closing_px < 0 or erosion_px < 0:
        raise ValueError(
            "texture_window_px must be >= 2, noise_sigmas > 0, "
            "closing_px and erosion_px >= 0"
        )
    data = np.asarray(image, dtype=float)
    if not np.isfinite(data).all():
        raise ValueError("image must be finite")

    texture = _local_std(ndi.median_filter(data, 3), texture_window_px)
    middle = float(np.median(texture))
    noise = 1.4826 * float(np.median(np.abs(texture - middle)))
    if noise == 0:
        return np.zeros(data.shape, dtype=np.int32)

    binary = texture > middle + noise_sigmas * noise
    if closing_px:
        binary = ndi.binary_closing(binary, morphology.disk(closing_px))
    binary = ndi.binary_fill_holes(binary)
    binary = _remove_objects_below_size(binary, min_size)
    binary = _erode_keeping_objects(binary, erosion_px)
    return _split_touching(binary)


def _local_std(image: np.ndarray, window: int) -> np.ndarray:
    mean = ndi.uniform_filter(image, window)
    return np.sqrt(np.maximum(ndi.uniform_filter(image * image, window) - mean * mean, 0.0))


def _erode_keeping_objects(binary: np.ndarray, radius: int, keep: float = 0.5) -> np.ndarray:
    """Erode as far as possible while at least `keep` of the objects survive."""
    if radius <= 0:
        return binary
    before = ndi.label(binary)[1]
    if before == 0:
        return binary
    for current in range(radius, 0, -1):
        candidate = ndi.binary_erosion(binary, morphology.disk(current))
        if ndi.label(candidate)[1] >= keep * before:
            return candidate
    return binary


def _segment_resized(bf_image: np.ndarray, *, method: str, min_size: int, resize_factor: float, **cellpose_kwargs) -> np.ndarray:
    from skimage.transform import rescale, resize

    small = rescale(bf_image, resize_factor, anti_aliasing=True, preserve_range=True)
    small_min_size = max(1, int(round(min_size * resize_factor**2)))
    small_labels = segment_cells(small, method=method, min_size=small_min_size, resize_factor=1.0, **cellpose_kwargs)

    labels = resize(small_labels, bf_image.shape, order=0, preserve_range=True, anti_aliasing=False)
    return labels.astype(np.int32)


def load_cellpose_model(gpu: bool = False):
    """Load a CellposeModel once, to reuse across many segment_cells(method="cellpose", ...)
    calls (e.g. in pipeline.process_experiment over many images) instead of
    reloading the ~1.1GB pretrained model from disk every single call.
    """
    try:
        from cellpose import models
    except ImportError as exc:
        raise ImportError(
            "cellpose (and torch) are required for method='cellpose'. "
            "Install the 'cellpose' extra, or use method='otsu'."
        ) from exc
    return models.CellposeModel(gpu=gpu)


def _segment_cellpose(bf_image: np.ndarray, cellpose_model=None, diameter: float | None = None, gpu: bool = False) -> np.ndarray:
    """Segment with Cellpose (cellpose>=4: models.CellposeModel, not the
    removed models.Cellpose(model_type=...) class from cellpose<4 -- verified
    against a real install: cellpose==4.2.1.1 downloads a ~1.1GB pretrained
    model (cpsam_v2) on first use and takes ~2 min/image on CPU, but produces
    plausible masks (79 cells on the real 179997/TRX2 DIC image, close to the
    68 the "otsu" method finds on the same image).
    """
    model = cellpose_model if cellpose_model is not None else load_cellpose_model(gpu=gpu)
    masks, *_ = model.eval(bf_image, diameter=diameter, channels=None)
    return masks.astype(np.int32)


def _segment_otsu(bf_image: np.ndarray, min_size: int = 30) -> np.ndarray:
    image = bf_image.astype(float)
    image = (image - image.min()) / (np.ptp(image) + 1e-9)

    # Threshold on local edge/texture "activity" rather than raw intensity.
    # DIC/phase cells stand out by internal structure and boundary contrast,
    # not by being reliably brighter or darker than one global level -- a
    # plain global-intensity Otsu threshold (even after flattening a slow
    # illumination gradient) misclassified large low-contrast background
    # regions as foreground on real DIC test images, producing a spurious
    # web of "cells" across empty background. Edge activity stays low there
    # regardless of absolute brightness or gentle shading.
    smoothed = filters.gaussian(image, sigma=1.0)
    activity = filters.gaussian(filters.sobel(smoothed), sigma=3.0)
    threshold = filters.threshold_otsu(activity)
    binary = activity > threshold

    # edges alone are a ring, not a filled cell -- close small gaps in the
    # ring then fill it completely (not just small holes) to get solid masks
    binary = morphology.closing(
    binary, morphology.disk(3), mode="ignore"
    )
    binary = ndi.binary_fill_holes(binary)
    binary = _remove_objects_below_size(binary, min_size)

    return _split_touching(binary)


def _segment_sparse(
    image: np.ndarray,
    *,
    min_size: int = 30,
    noise_sigmas: float = 3.0,
    background_sigma_px: float = 30.0,
    smoothing_sigma_px: float = 1.5,
    boundary_erosion_px: int = 0,
) -> np.ndarray:
    """Segment a fluorescence frame where cells are a small minority of pixels.

    Why not "otsu" here: that method thresholds edge activity with Otsu's
    rule, which splits a histogram into two classes of comparable weight.
    On the team's 2026-10-05 frames the cells are 1-4% of the pixels, so
    there is no second class to find and the split falls wherever the
    background happens to spread. Measured on those frames, the resulting
    foreground ran from 0.1% (18 objects where ~90 cells are visible) to
    87% (shattered into 1301 fragments) across images from one session --
    the swing tracked how dim the frame was, not how well it was segmented.

    This method makes the sparse assumption explicit instead. It removes the
    background with a high-pass (anything varying more slowly than
    background_sigma_px is illumination, not a cell), estimates the noise of
    what is left from its MAD, and keeps what rises noise_sigmas above it.
    The threshold is therefore set by each frame's own noise rather than by
    its histogram shape, which is what makes it hold as frames get dimmer.

    noise_sigmas is a choice, not a calibration: 3.0 is the usual
    detection-threshold convention, but the value that suits a given
    microscope should be set against fields whose cells have been counted by
    eye. Raise it to drop faint cells, lower it to admit more noise.

    Masks come out larger than the cells, because the smoothing that makes
    the threshold robust also spreads each cell outward: against synthetic
    discs the masks recover every cell but run about 1.5x their area, so a
    per-cell mean is pulled toward the local background. That dilution is
    similar across a session and so largely cancels in a ratio or in a
    comparison against a time-matched control, but it is not negligible in
    an absolute intensity. boundary_erosion_px shrinks the masks back:
    measured on the same discs, 2 px took precision from 0.65 to 0.98 while
    keeping recall at 0.98. It is 0 by default because the same 2 px cost a
    third of the objects on the team's real dim frames (77 -> 41) -- real
    cells are not discs, and the faint ones erode away entirely. Tighten it
    only against fields counted by eye.
    """
    if image.ndim != 2:
        raise ValueError("image must be 2D")
    if noise_sigmas <= 0 or background_sigma_px <= 0:
        raise ValueError("noise_sigmas and background_sigma_px must be positive")
    if boundary_erosion_px < 0:
        raise ValueError("boundary_erosion_px must not be negative")

    # The median filter takes out single-pixel spikes, which would otherwise
    # inflate the MAD and so raise the threshold for the whole frame.
    smoothed = filters.gaussian(ndi.median_filter(image.astype(float), size=3), sigma=smoothing_sigma_px)
    high_pass = smoothed - filters.gaussian(smoothed, sigma=background_sigma_px)

    # MAD, not the standard deviation: the cells are in this image too, and
    # a standard deviation would count them as noise and threshold them away.
    noise = 1.4826 * float(np.median(np.abs(high_pass - np.median(high_pass))))
    if noise <= 0:
        return np.zeros(image.shape, dtype=np.int32)

    binary = ndi.binary_fill_holes(high_pass > noise_sigmas * noise)
    if boundary_erosion_px:
        binary = morphology.erosion(binary, morphology.disk(boundary_erosion_px))
    # After erosion, so that min_size is applied to the masks actually
    # returned rather than to the inflated ones.
    binary = _remove_objects_below_size(binary, min_size)
    labels = _split_touching(binary)

    # Watershed can split retained components into smaller labels.
    ids, areas = np.unique(labels, return_counts=True)
    small_ids = ids[(ids != 0) & (areas < min_size)]
    if small_ids.size:
        labels[np.isin(labels, small_ids)] = 0

    return labels


def _split_touching(binary: np.ndarray) -> np.ndarray:
    """Watershed on the distance transform, so touching cells stay separate."""
    distance = ndi.distance_transform_edt(binary)
    peaks = peak_local_max(distance, min_distance=5, labels=binary.astype(int))
    markers = np.zeros_like(binary, dtype=np.int32)
    for i, (r, c) in enumerate(peaks, start=1):
        markers[r, c] = i
    markers = ndi.label(markers)[0]
    return segmentation.watershed(-distance, markers, mask=binary).astype(np.int32)


def border_touching_labels(label_mask: np.ndarray) -> set[int]:
    """Cell IDs whose mask touches the image border (candidates for edge_flag, Stage 7)."""
    border = np.concatenate([
        label_mask[0, :], label_mask[-1, :], label_mask[:, 0], label_mask[:, -1],
    ])
    return {int(label) for label in np.unique(border) if label != 0}


def focus_score(gray_image: np.ndarray) -> float:
    """Variance-of-Laplacian focus metric for the phase/DIC channel (Stage 7 focus QC).

    Normalized to the image's own dynamic range (min/max) before computing
    the Laplacian, so the score -- and any threshold calibrated against it
    -- is comparable across bit depths and exposure levels, not just within
    one camera's raw pixel-value scale. This matters concretely here: the
    team's real Olympus camera (SC30, per its manual) is 8-bit/channel
    (0-255), while the public reference images this project's threshold was
    originally calibrated on, and most of its tests run against, are 16-bit
    (0-65535) -- raw (non-normalized) variance-of-Laplacian would differ by
    roughly (65535/255)^2 =~ 66000x between them for similar-looking images,
    making a threshold calibrated on one meaningless on the other. Still
    re-run scripts/calibrate_focus_threshold.py once real Olympus images
    exist -- normalization fixes the *scale* mismatch, not necessarily the
    absolute pass/fail cutoff for genuinely different optics/noise.

    Calibrate the pass/fail threshold on known good/bad fields for the
    current microscope, then compare against it in export.build_records.
    """
    from scipy.ndimage import laplace
    image = gray_image.astype(float)
    value_range = image.max() - image.min()
    normalized = (image - image.min()) / value_range if value_range > 0 else image - image.min()
    return float(np.var(laplace(normalized)))
