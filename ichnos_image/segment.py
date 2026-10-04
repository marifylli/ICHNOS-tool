"""bright-field/DIC -> per-cell label masks (Stage 2), plus a couple of
Stage-2/7 helpers (border-touching cells, focus score) that naturally operate
on the same phase/DIC image used for segmentation.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from skimage import filters, morphology, segmentation, measure
from skimage.feature import peak_local_max


def segment_cells(
    bf_image: np.ndarray, method: str = "otsu", min_size: int = 30, resize_factor: float = 1.0, **cellpose_kwargs
) -> np.ndarray:
    """Return an integer label mask (0 = background) for a bright-field/DIC image.

    method="cellpose" needs the optional cellpose+torch dependencies (deep
    learning segmentation, preferred for stress-distorted morphology per the
    team's protocol) -- ~110s/image on CPU for a real 535x512 DIC image, so
    batch runs may need either GPU or a resize_factor < 1.0 (see below --
    unlike "otsu", cellpose tolerates this well). See
    scripts/estimate_cellpose_throughput.py to turn that per-image number
    into a "how long will my real batch take" estimate.

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
    raise ValueError(f"unknown segmentation method: {method!r}")


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
    binary = morphology.remove_small_objects(binary, min_size=min_size)

    distance = ndi.distance_transform_edt(binary)
    peaks = peak_local_max(distance, min_distance=5, labels=binary.astype(int))
    markers = np.zeros_like(binary, dtype=np.int32)
    for i, (r, c) in enumerate(peaks, start=1):
        markers[r, c] = i
    markers = ndi.label(markers)[0]

    labels = segmentation.watershed(-distance, markers, mask=binary)
    return labels.astype(np.int32)


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
