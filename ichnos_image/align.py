"""Put a bright-field mask on the fluorescence frame it belongs with.

Two questions, one measurement. Does the bright-field frame show the same
field of view as the fluorescence pair, and if so, by how much are they
displaced? Both matter before bright-field masks can be used at all: a mask
from a different field reports background as cells, and a mask off by tens
of pixels does the same thing more subtly.

The measurement is: lay the mask over the fluorescence frame and take the
mean fluorescence under it. Do that for every possible displacement -- which
one correlation gives at once -- and the best displacement is the
registration offset, while how far it stands above the rest of that surface
is the evidence that the field matches at all.

Why not correlate the two images directly. A bright-field frame and a
fluorescence frame of one field do not look alike: cells are dark rims on a
bright background in one and bright blobs on a dark background in the other,
and most cells that are visible in the first are invisible in the second. An
earlier attempt compared their local texture and judged the result against a
reference measured between two fluorescence frames, which is a comparison of
a different kind; it scored true matches at 0.2-0.35 against a 0.64
reference and read 24 of them as misses. Asking whether the masks land on
signal avoids the comparison entirely, and works when only a handful of the
cells fluoresce -- the case this exists for.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi
from scipy.signal import fftconvolve


@dataclass(frozen=True)
class MaskAlignment:
    """Where a mask best sits on an image, and how clearly."""

    #: Displacement to apply to the mask, in pixels (row, column).
    dy: float
    dx: float
    #: Height of the best displacement above the rest of the surface, in
    #: robust standard deviations. Compare it against the same number
    #: computed for a frame the mask certainly does not belong to; there is
    #: no threshold that transfers between sessions.
    z: float
    #: Mean image value inside the placed mask over mean outside. 1.0 means
    #: the mask has found nothing.
    contrast: float
    foreground_fraction: float

    @property
    def shift_px(self) -> float:
        return float(np.hypot(self.dy, self.dx))


def shift_mask(mask: np.ndarray, dy: int, dx: int) -> np.ndarray:
    """Translate without wrapping.

    A rolled mask brings cells that leave one edge back in at the opposite
    one, where they would be scored against unrelated pixels -- which
    inflates exactly the agreement being measured.
    """
    dy, dx = int(dy), int(dx)
    out = np.zeros_like(mask)
    rows, cols = mask.shape
    if abs(dy) >= rows or abs(dx) >= cols:
        return out
    out[max(dy, 0):rows + min(dy, 0), max(dx, 0):cols + min(dx, 0)] = mask[
        max(-dy, 0):rows - max(dy, 0), max(-dx, 0):cols - max(dx, 0)
    ]
    return out


def align_mask(
    mask: np.ndarray,
    image: np.ndarray,
    *,
    max_shift: int = 200,
    background_sigma_px: float = 30.0,
    guard_px: int = 25,
) -> MaskAlignment:
    """Best displacement of `mask` onto `image`, and the evidence for it."""
    mask = np.asarray(mask) > 0
    image = np.asarray(image, dtype=float)
    if mask.shape != image.shape:
        raise ValueError("mask and image must have the same shape")
    if max_shift < 1 or guard_px < 0:
        raise ValueError("max_shift must be >= 1 and guard_px >= 0")
    area = float(mask.sum())
    if area == 0:
        return MaskAlignment(np.nan, np.nan, np.nan, np.nan, 0.0)

    # Remove the illumination background first, so that "mean under the
    # mask" means signal rather than whichever part of the vignette the mask
    # happens to cover.
    flat = image - ndi.gaussian_filter(image, background_sigma_px)
    flipped = mask[::-1, ::-1].astype(float)
    total = fftconvolve(flat, flipped, mode="same")
    # How much of the mask is still inside the frame at each displacement.
    # Without dividing by this, a large displacement scores low merely
    # because part of the mask has left the image, which drags the answer
    # toward zero: on a synthetic pair with a known (+17, -23) shift the
    # undivided version reported (+11, -16).
    overlap = fftconvolve(np.ones_like(flat), flipped, mode="same")
    with np.errstate(invalid="ignore", divide="ignore"):
        surface = np.where(overlap > 0.5 * area, total / np.maximum(overlap, 1.0), np.nan)

    centre = np.asarray(surface.shape) // 2
    lo = np.maximum(centre - max_shift, 0)
    hi = np.minimum(centre + max_shift + 1, surface.shape)
    window = surface[lo[0]:hi[0], lo[1]:hi[1]]
    if not np.isfinite(window).any():
        return MaskAlignment(np.nan, np.nan, np.nan, np.nan, float(mask.mean()))

    peak = np.unravel_index(int(np.nanargmax(window)), window.shape)
    # The null is the rest of this very surface, with the peak's own
    # neighbourhood removed: scores next to a true peak are themselves high,
    # and leaving them in would inflate the spread and hide the peak.
    rest = window.copy()
    rest[max(peak[0] - guard_px, 0):peak[0] + guard_px + 1,
         max(peak[1] - guard_px, 0):peak[1] + guard_px + 1] = np.nan
    middle = float(np.nanmedian(rest))
    spread = 1.4826 * float(np.nanmedian(np.abs(rest - middle)))
    z = float((window[peak] - middle) / spread) if spread > 0 else np.nan

    dy = float(peak[0] + lo[0] - centre[0])
    dx = float(peak[1] + lo[1] - centre[1])
    placed = shift_mask(mask, int(dy), int(dx))
    outside = image[~placed]
    contrast = (float(image[placed].mean() / outside.mean())
                if placed.any() and outside.size and outside.mean() else np.nan)
    return MaskAlignment(dy, dx, z, contrast, float(mask.mean()))
