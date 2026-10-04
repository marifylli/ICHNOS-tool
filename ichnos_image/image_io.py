"""Reading image files off disk, and turning what a colour camera stores into
the 2D fluorescence plane the rest of the pipeline works on.

Extracted out of synthesize.py, which is a synthetic-data helper: the runtime
pipeline should not import its loader from a demo module.

The part that needs care is the second job. The camera on this setup is a
*colour* camera (QImaging MicroPublisher 3.3 RTV) saving 24-bit RGB. A GFP
acquisition and an mCherry acquisition are two separate exposures taken
through two different filter cubes, and each one lands on disk as a
three-component RGB image, not as a single-channel measurement. Deciding
which component -- or which calibrated combination of components -- stands
for "the green signal" is a measurement decision that changes every number
downstream, so this module refuses to make it silently.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
from PIL import Image

from . import instrument


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------


def load_image(path: str | Path) -> np.ndarray:
    """Load an image file as float64, shape unchanged.

    An RGB file comes back with its trailing axis of size 3 intact. Nothing
    is collapsed, rescaled or reordered here -- see
    extract_fluorescence_plane() for going from RGB to a 2D plane.
    """
    return np.array(Image.open(path)).astype(np.float64)


def load_png16(path: str | Path) -> np.ndarray:
    """Historical name for load_image(), kept so existing callers and the
    synthetic fixtures keep working. New code should call load_image().

    The name is misleading twice over: it loads any format PIL can read, and
    it returns float64 rather than uint16. That float64 return is exactly why
    saturation must come from the instrument profile and never be inferred
    from the array's dtype -- see saturation_value_for_file().
    """
    return load_image(path)


def save_png16(array: np.ndarray, path: str | Path) -> None:
    """Write a 16-bit PNG. Used for synthetic fixtures, not for real data."""
    Image.fromarray(np.clip(array, 0, 65535).astype(np.uint16)).save(path)


def saturation_value_for_file(array: np.ndarray) -> float:
    """The value at which this setup's stored files clip.

    Deliberately ignores the array's dtype. load_image() returns float64 for
    everything, and the previous pipeline inferred saturation from dtype,
    which gave 65535 for an 8-bit file loaded into a wider container -- so
    genuinely clipped pixels at 255 went unflagged. The stored depth is a
    property of the acquisition preset, not of the array in memory.
    """
    del array  # signature kept for call-site clarity; the answer is the rig's
    return instrument.SATURATION_VALUE


# --------------------------------------------------------------------------
# RGB -> one fluorescence plane
# --------------------------------------------------------------------------

#: How to reduce a stored RGB acquisition to one 2D plane.
#:
#: "R", "G", "B"   take that single component. For a colour camera behind a
#:                 single-band emission filter this is usually the honest
#:                 choice: the signal lands predominantly on the Bayer
#:                 component matching the emission band, and the other two
#:                 carry mostly crosstalk and noise.
#: "sum"           add the three components. Keeps all the photons, but the
#:                 per-component white balance (R5/G3/B2 on this preset) is
#:                 then baked in unevenly, so it is only meaningful if that
#:                 balance is fixed and documented.
#: "luminance"     the usual 0.299/0.587/0.114 perceptual weighting. Designed
#:                 for how displays look to people, NOT for photometry.
EXTRACTION_METHODS = ("R", "G", "B", "sum", "luminance")

_LUMINANCE_WEIGHTS = (0.299, 0.587, 0.114)


class UncalibratedExtractionError(RuntimeError):
    """Raised when a 2D plane is requested from an RGB image without saying
    how to reduce it."""


def extract_fluorescence_plane(
    image: np.ndarray,
    *,
    method: Optional[str] = None,
    cube: Optional[str] = None,
) -> np.ndarray:
    """Reduce a stored acquisition to the 2D fluorescence plane.

    A 2D input passes straight through -- synthetic fixtures and any
    monochrome source are already a plane.

    For a 3-component input, `method` is required. There is no default, on
    purpose: which component represents the signal is an open measurement
    decision for this setup (instrument.UNRESOLVED["rgb_extraction"]), and a
    default would make an unvalidated choice look like a settled one. A
    generic RGB-to-grayscale conversion in particular must never be the
    silent fallback for quantitative comparison.

    `cube` ("B" for green, "G" for red) is recorded for provenance and
    cross-checked against instrument.FILTER_CUBES. It does not select the
    method -- the mapping from cube to component is exactly what has to be
    measured, not assumed.
    """
    arr = np.asarray(image, dtype=np.float64)

    if arr.ndim == 2:
        return arr
    if arr.ndim != 3 or arr.shape[-1] not in (3, 4):
        raise ValueError(
            f"expected a 2D plane or an RGB(A) image, got shape {arr.shape}"
        )

    if cube is not None and cube not in instrument.FILTER_CUBES:
        raise ValueError(
            f"{cube!r} is not a cube position on this microscope; "
            f"available: {sorted(instrument.FILTER_CUBES)}"
        )

    if method is None:
        raise UncalibratedExtractionError(
            "this is a 3-component RGB acquisition and no extraction method was "
            "given. Pass method= one of "
            f"{EXTRACTION_METHODS}. There is no default: which component "
            "carries the signal for each filter cube is an unresolved "
            "measurement question for this setup, and picking one silently "
            "would make an unvalidated choice look measured."
        )
    if method not in EXTRACTION_METHODS:
        raise ValueError(f"unknown method {method!r}; expected one of {EXTRACTION_METHODS}")

    rgb = arr[..., :3]
    if method in ("R", "G", "B"):
        return rgb[..., "RGB".index(method)]
    if method == "sum":
        return rgb.sum(axis=-1)
    return rgb @ np.asarray(_LUMINANCE_WEIGHTS)


def plane_for_channel(
    image: np.ndarray, channel: str, *, method: Optional[str] = None
) -> np.ndarray:
    """extract_fluorescence_plane() for a named biological channel.

    channel is "green" (GFP, cube B) or "red" (mCherry, cube G). This only
    records which cube the acquisition used; it still will not choose an
    extraction method for you.
    """
    if channel not in instrument.CUBE_FOR_CHANNEL:
        raise ValueError(
            f"unknown channel {channel!r}; expected one of "
            f"{sorted(instrument.CUBE_FOR_CHANNEL)}"
        )
    return extract_fluorescence_plane(
        image, method=method, cube=instrument.CUBE_FOR_CHANNEL[channel]
    )


def saturation_mask_for_image(
    image: np.ndarray,
    saturation_value: float,
) -> np.ndarray:
    """Check clipping in the stored image, before RGB extraction."""
    if not np.isfinite(saturation_value) or saturation_value <= 0:
        raise ValueError("saturation_value must be finite and positive")

    arr = np.asarray(image)

    if arr.ndim == 2:
        return arr >= saturation_value

    if arr.ndim == 3 and arr.shape[-1] in (3, 4):
        return np.any(arr[..., :3] >= saturation_value, axis=-1)

    raise ValueError(f"unsupported image shape: {arr.shape}")