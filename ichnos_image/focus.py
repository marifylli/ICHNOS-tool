"""Per-cell focus QC (Stage 7) and cross-channel focus agreement.

Why this exists, beyond segment.focus_score(): that one scores a whole
field off the phase/DIC channel. Two things it cannot catch, both present
in the team's first real H2O2 session (2026-10-05, where the wet lab
reported difficulty focusing):

1. Defocus varies *within* a field — cells sit at different heights under
   the coverslip, so a field that scores well overall still contains cells
   whose light has spread out of the mask into their own background ring.
   extract.extract_per_cell() then under-measures the cell and
   over-measures its local background, and the two errors push the same
   way: every response is compressed toward no-change. A field-level score
   cannot flag those cells, and the compression is invisible in the output.

2. Green and red are separate exposures, and focus can differ between them
   (refocusing between channels, chromatic focal shift). The red/green
   ratio is then a focus artefact rather than biology. Nothing in the
   pipeline currently compares the two channels' sharpness.

The metric is an edge-contrast ratio: the gradient magnitude on the cell's
boundary band divided by the gradient magnitude in the surrounding
background ring, both after a small Gaussian smoothing. A sharp cell has a
steep boundary and scores well above 1; as it defocuses the boundary
spreads and its gradient falls toward the background's own texture, so the
score falls toward 1.0 (the floor: the boundary is indistinguishable from
background).

Two properties this shape buys, both needed for this data:

  - It is a ratio of two gradients of the same image, so it does not move
    when exposure or gain changes. That session varied exposure between
    1.2 s and 2.4 s, so any absolute metric would be unusable across it.
  - The denominator is measured per cell, so a dim cell is judged against
    its own noise rather than against a bright cell's. On an 8-bit sensor
    with cells only a few grey levels above background, a metric that does
    not do this reads noise as sharpness.

The score is relative: rank and threshold cells within a dataset, and
calibrate the cutoff per microscope on known good and bad fields (see
tools/calibrate_focus_threshold.py). It does not measure how far from
focus a cell is, and it cannot repair a blurred cell — a flagged cell is
dropped, not corrected.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import binary_dilation, binary_erosion, gaussian_filter
from skimage.measure import regionprops
from skimage.morphology import disk

#: Floor of the metric: at 1.0 the cell boundary is no steeper than the
#: background's own texture.
NO_EDGE_SCORE = 1.0


@dataclass(frozen=True)
class CellFocus:
    """Per-cell focus QC for one channel."""

    cell_id: int
    #: Edge gradient / background gradient. 1.0 = no detectable edge.
    #: NaN when `scored` is False.
    focus_score: float
    #: False when the cell has too little boundary or background ring to
    #: measure, so the score is NaN rather than a silent 0.0 or 1.0.
    scored: bool


@dataclass(frozen=True)
class ChannelFocusAgreement:
    """How similar two channels' focus is over the same cells."""

    n_cells: int
    median_green: float
    median_red: float
    #: (green - 1) / (red - 1) folded to <= 1, i.e. comparing the two
    #: channels' edge contrast above the no-edge floor. 1.0 means equally
    #: sharp; lower means one channel is blurrier than the other.
    agreement: float
    agrees: bool


def _patch(region, pad: int, shape) -> tuple[slice, slice]:
    r0, c0, r1, c1 = region.bbox
    return (
        slice(max(r0 - pad, 0), min(r1 + pad, shape[0])),
        slice(max(c0 - pad, 0), min(c1 + pad, shape[1])),
    )


def score_cells(
    label_mask: np.ndarray,
    channel: np.ndarray,
    *,
    smoothing_sigma_px: float = 1.5,
    boundary_band_px: int = 2,
    ring_width_px: int = 6,
    min_ring_px: int = 20,
    min_boundary_px: int = 5,
) -> list[CellFocus]:
    """Per-cell focus score for one channel.

    Pass the channel the intensities are measured from, so the score
    describes the same pixels. `smoothing_sigma_px` suppresses pixel noise
    before the gradient is taken; raise it for noisier sensors, but keep it
    well below the cell radius or genuine edges get smoothed away too.

    Cells touching the image border are scored like any other — combine
    with the existing edge_flag downstream rather than dropping them here.
    """
    if label_mask.shape != channel.shape:
        raise ValueError("label_mask and channel must have the same shape")
    if smoothing_sigma_px < 0:
        raise ValueError("smoothing_sigma_px must be non-negative")

    smoothed = gaussian_filter(channel.astype(float), smoothing_sigma_px)
    gy, gx = np.gradient(smoothed)
    gradient = np.hypot(gy, gx)
    occupied = label_mask > 0
    pad = ring_width_px + 2
    band, ring_se = disk(boundary_band_px), disk(ring_width_px)

    scores: list[CellFocus] = []
    for region in regionprops(label_mask):
        cell_id = int(region.label)
        sl = _patch(region, pad, label_mask.shape)
        cell = label_mask[sl] == cell_id
        boundary = binary_dilation(cell, band) & ~binary_erosion(cell, band)
        ring = binary_dilation(cell, ring_se) & ~occupied[sl]

        if boundary.sum() < min_boundary_px or ring.sum() < min_ring_px:
            scores.append(CellFocus(cell_id, float("nan"), False))
            continue

        patch = gradient[sl]
        background = float(np.median(patch[ring]))
        if background <= 0:
            scores.append(CellFocus(cell_id, float("nan"), False))
            continue

        # 75th percentile, not the mean: a cell is rarely in focus on one
        # side and out on the other, but it is routinely overlapped by a
        # neighbour or clipped by the patch edge along part of its
        # boundary, and those stretches drag a mean down.
        edge = float(np.percentile(patch[boundary], 75))
        scores.append(CellFocus(cell_id, max(edge / background, NO_EDGE_SCORE), True))

    return scores


def channel_agreement(
    green_scores: list[CellFocus],
    red_scores: list[CellFocus],
    *,
    min_agreement: float = 0.6,
    min_cells: int = 5,
) -> ChannelFocusAgreement:
    """Compare two channels' per-cell focus over the cells scored in both.

    Low agreement invalidates the red/green ratio for that image set: the
    channels were not acquired at the same focal plane, so their
    intensities are not comparable however well each one is measured
    individually. The response is to drop that image set's ratios, not to
    correct them. `agrees` is False when fewer than `min_cells` cells were
    scored in both channels, since the comparison is not informative then.

    Scores are compared above the no-edge floor rather than as raw values:
    a channel at 1.1 and one at 1.3 differ threefold in edge contrast, but
    their raw ratio of 0.85 would read as near-agreement.
    """
    by_id = {s.cell_id: s.focus_score for s in green_scores if s.scored}
    pairs = [
        (by_id[s.cell_id], s.focus_score)
        for s in red_scores
        if s.scored and s.cell_id in by_id
    ]
    if len(pairs) < min_cells:
        return ChannelFocusAgreement(
            len(pairs), float("nan"), float("nan"), float("nan"), False
        )

    green = float(np.median([g for g, _ in pairs]))
    red = float(np.median([r for _, r in pairs]))
    above_floor = (green - NO_EDGE_SCORE, red - NO_EDGE_SCORE)
    high = max(above_floor)
    agreement = max(min(above_floor), 0.0) / high if high > 0 else 0.0
    return ChannelFocusAgreement(
        len(pairs), green, red, agreement, agreement >= min_agreement
    )
