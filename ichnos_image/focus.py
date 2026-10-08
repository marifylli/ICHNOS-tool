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

A caution before the metric: it reads *measured* sharpness, which a noisy
channel loses as surely as a defocused one. Verified on this session's own
images — dropping the red channel's contrast to the green channel's while
leaving the focal plane untouched moved its median score from 1.86 to
1.59, most of the way to the green channel's 1.36. So a lower score in the
dimmer channel is not evidence that it was focused differently, and
`channel_agreement` refuses to draw that conclusion at all when a channel's
cells are too close to the noise to carry a measurable edge (see
`comparable` there). On this session the green channel sits at a
contrast-to-noise ratio near 0.85 — its cells are within the noise — so the
session's apparent focus difference is reported as "cannot tell" rather
than as a finding.

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
    #: (cell mean - background median) / background noise, from the same
    #: ring. How much of focus_score the cell's contrast alone could
    #: explain: two channels' scores are only comparable at similar values.
    contrast_to_noise: float
    #: False when the cell has too little boundary or background ring to
    #: measure, so the score is NaN rather than a silent 0.0 or 1.0.
    scored: bool


@dataclass(frozen=True)
class ChannelFocusAgreement:
    """How similar two channels' focus is over the same cells."""

    n_cells: int
    median_green: float
    median_red: float
    median_cnr_green: float
    median_cnr_red: float

    #: (green - 1) / (red - 1) folded to <= 1, i.e. comparing the two
    #: channels' edge contrast above the no-edge floor. 1.0 means equally
    #: sharp; lower means one channel is blurrier than the other.
    agreement: float
    #: False when either channel's cells sit too close to the noise for its
    #: focus score to mean anything. Then `agrees` is False as "cannot
    #: tell", not as "the channels were focused differently".
    comparable: bool
    #: True only when the channels are comparable AND their focus agrees.
    #: False with comparable=False means "cannot tell", not "disagrees".
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

    image = channel.astype(float)
    smoothed = gaussian_filter(image, smoothing_sigma_px)
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
            scores.append(CellFocus(cell_id, float("nan"), float("nan"), False))
            continue

        patch = gradient[sl]
        background = float(np.median(patch[ring]))
        if background <= 0:
            scores.append(CellFocus(cell_id, float("nan"), float("nan"), False))
            continue

        ring_values = image[sl][ring]
        ring_level = float(np.median(ring_values))
        # MAD, not standard deviation: a stray bright speck in the ring
        # would inflate an SD and make a clean cell look noisy.
        noise = 1.4826 * float(np.median(np.abs(ring_values - ring_level)))
        cnr = (float(image[sl][cell].mean()) - ring_level) / noise if noise > 0 else float("nan")

        # 75th percentile, not the mean: a cell is rarely in focus on one
        # side and out on the other, but it is routinely overlapped by a
        # neighbour or clipped by the patch edge along part of its
        # boundary, and those stretches drag a mean down.
        edge = float(np.percentile(patch[boundary], 75))
        scores.append(
            CellFocus(cell_id, max(edge / background, NO_EDGE_SCORE), cnr, True)
        )

    return scores


def channel_agreement(
    green_scores: list[CellFocus],
    red_scores: list[CellFocus],
    *,
    min_agreement: float = 0.6,
    min_cells: int = 5,
    min_cnr: float = 1.0,
) -> ChannelFocusAgreement:
    """Compare two channels' per-cell focus over the cells scored in both.

    Low agreement flags inconsistent measured sharpness and makes the ratio
    suspect; it does not prove different focal planes. Biological distribution,
    contrast and acquisition differences can also change edge scores. An
    enforced policy rejects these ratios instead of inventing a correction. `agrees` is False when fewer than `min_cells` cells were
    scored in both channels, since the comparison is not informative then.

    Scores are compared above the no-edge floor rather than as raw values:
    a channel at 1.1 and one at 1.3 differ threefold in edge contrast, but
    their raw ratio of 0.85 would read as near-agreement.

    `comparable` guards the whole comparison. A channel whose cells barely
    rise above their own noise scores low on focus whatever its focal
    plane, because there is no edge left to measure through the noise; its
    score cannot be read as sharpness at all, let alone compared with
    another channel's. Defocus and underexposure both destroy contrast, and
    one image cannot tell them apart, so the function does not try: below
    `min_cnr` in either channel it reports "cannot tell" and the fix is in
    the acquisition (more exposure on the dim channel), not in the
    analysis.
    """
    by_id = {s.cell_id: s for s in green_scores if s.scored}
    pairs = [
        (by_id[s.cell_id], s)
        for s in red_scores
        if s.scored and s.cell_id in by_id
    ]
    nan = float("nan")
    if len(pairs) < min_cells:
        return ChannelFocusAgreement(
            len(pairs), nan, nan, nan, nan, nan, False, False
        )

    green = float(np.median([g.focus_score for g, _ in pairs]))
    red = float(np.median([r.focus_score for _, r in pairs]))
    def finite_median(values):
        values = np.asarray(values)
        return float(np.median(values[np.isfinite(values)])) if np.isfinite(values).any() else nan
    cnr_green = finite_median([g.contrast_to_noise for g, _ in pairs])
    cnr_red = finite_median([r.contrast_to_noise for _, r in pairs])

    comparable = bool(np.isfinite([cnr_green, cnr_red]).all() and min(cnr_green, cnr_red) >= min_cnr)

    above_floor = (green - NO_EDGE_SCORE, red - NO_EDGE_SCORE)
    high = max(above_floor)
    agreement = max(min(above_floor), 0.0) / high if high > 0 else 0.0
    return ChannelFocusAgreement(
        len(pairs), green, red, cnr_green, cnr_red, agreement,
        comparable, comparable and agreement >= min_agreement,
    )


@dataclass(frozen=True)
class FocusPolicy:
    """Explicit QC thresholds; report mode measures without changing eligibility.

    Enforce mode rejects unknown/low-CNR measurements as indeterminate, never
    labels low contrast as proven defocus. Thresholds require local validation.
    """
    mode: str = 'report'
    min_score: float | None = None
    min_cnr: float = 1.0
    min_agreement: float = 0.6

    def __post_init__(self):
        if self.mode not in {'report', 'enforce'}:
            raise ValueError('focus mode must be report or enforce')
        if self.mode == 'enforce' and self.min_score is None:
            raise ValueError('enforced focus QC requires an explicit min_score')
        for name in ('min_cnr', 'min_agreement', 'min_score'):
            value = getattr(self, name)
            if value is not None and (isinstance(value, bool) or not np.isfinite(value) or value <= 0):
                raise ValueError(f'{name} must be finite and positive')
        if self.min_agreement > 1:
            raise ValueError('min_agreement must be <= 1')


def evaluate_cell(green: CellFocus, red: CellFocus, policy: FocusPolicy):
    """Return diagnostics and eligibility separately, preserving unknown states."""
    if not green.scored or not red.scored:
        status, agreement = 'focus_unscorable', None
    elif not np.isfinite([green.contrast_to_noise, red.contrast_to_noise]).all():
        status, agreement = 'focus_cnr_unknown', None
    elif min(green.contrast_to_noise, red.contrast_to_noise) < policy.min_cnr:
        status, agreement = 'focus_low_cnr', None
    else:
        comparison = channel_agreement([green], [red], min_cells=1,
            min_cnr=policy.min_cnr, min_agreement=policy.min_agreement)
        agreement = comparison.agreement
        if policy.min_score is not None and min(green.focus_score, red.focus_score) < policy.min_score:
            status = 'focus_low_score'
        elif not comparison.agrees:
            status = 'focus_channel_disagreement'
        else:
            status = 'focus_pass'
    return status, agreement, policy.mode != 'enforce' or status == 'focus_pass'
