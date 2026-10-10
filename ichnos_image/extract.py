"""masks + corrected channels -> per-cell integrated features (Stage 5), plus
the raw red/green ratio (Stage 6 — geometric ratio only; FRET/maturation-
kinetics correction and the per-session calibration factor c_session belong
in the core ichnos ratio/decoder package, not here).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from skimage.measure import regionprops
from skimage.morphology import dilation, erosion, disk



@dataclass
class CellFeatures:
    cell_id: int
    area_px: int
    raw_mean_green: float
    raw_mean_red: float
    corrected_mean_green: float
    corrected_mean_red: float
    integrated_green: float
    integrated_red: float
    sat_flag: bool
    local_background_green: float
    local_background_red: float
    sat_flag_legacy: bool | None = None
    focus_score_green: float | None = None
    focus_score_red: float | None = None
    contrast_to_noise_green: float | None = None
    contrast_to_noise_red: float | None = None
    focus_agreement: float | None = None
    focus_status: str = 'not_evaluated'
    focus_qc_pass: bool = True
    focus_qc_mode: str = 'report'


def extract_per_cell(
    label_mask: np.ndarray,
    raw_green: np.ndarray,
    raw_red: np.ndarray,
    corrected_green: np.ndarray,
    corrected_red: np.ndarray,
    saturation_value: float | None = None,
    erosion_px: int = 1,
    annulus_width_px: int = 3,
    subtract_local_background: bool = True,
    raw_saturation_mask: np.ndarray | None = None,
    legacy_saturation_mask: np.ndarray | None = None,
    cell_focus: dict | None = None,
) -> list[CellFeatures]:
    """Per-cell mean/integrated intensity (Stage 5).

    Measurement region: the whole-cell mask eroded by `erosion_px` (settles
    the earlier open question -- shrinks the region away from the cell
    boundary, where partial-volume mixing with neighbors or background is
    worst), or the full mask if erosion would wipe out a very small cell.

    subtract_local_background=True (default) additionally estimates a
    per-cell local background from a ring (`annulus_width_px` wide) of
    non-cell pixels immediately surrounding the cell, and subtracts its
    median from that cell's (already Stage-3-corrected) intensities. This
    complements, not duplicates, Stage 3's single global background level:
    it captures *local* residual variation across the field (e.g. diffuse
    light near a dense cluster of bright cells) that one global estimate
    can't. Always run Stage 3's subtract_background() first and pass its
    output as corrected_green/corrected_red -- this only corrects the
    remaining local offset, not the main background level.

    saturation_value defaults to the raw channel dtype's max (e.g. 65535 for
    16-bit, matching the real downloaded reference images); a cell with any
    raw pixel at/above it is flagged sat_flag=True and should be dropped
    downstream — a saturated green channel biases the red/green ratio
    toward a falsely large "time since stress". Only auto-inferred for
    integer dtypes, where it reflects a real sensor bit depth: for a float
    channel (e.g. already background-subtracted, or synthetic test data)
    there's no way to recover the sensor's true saturation point from the
    array's own values, and defaulting to the array's own max would trivially
    flag whichever cell happens to be brightest, saturated or not (this was
    an actual bug here previously -- caught by tests/test_pipeline.py once
    it started asserting on qc_pass, not just running without crashing).
    Pass saturation_value explicitly for float channels.

    Don't trust dtype-based auto-inference blindly for the team's real
    camera either: the Olympus SC30 (per its manual) is 8-bit/channel --
    real saturation is 255 (ichnos.config.CAMERA_SC30_SATURATION_VALUE), NOT
    65535. If a real SC30 file ever ends up loaded into a wider container
    (e.g. an 8-bit source saved/read back as uint16), dtype-based inference
    would silently pick 65535 and never flag genuine saturation. Pass
    saturation_value=CAMERA_SC30_SATURATION_VALUE explicitly for real
    Olympus/SC30 data rather than relying on this default.
    """
    if raw_saturation_mask is not None:
        if raw_saturation_mask.shape != label_mask.shape:
            raise ValueError(
                "raw_saturation_mask must match label_mask shape"
            )
        if raw_saturation_mask.dtype != np.bool_:
            raise ValueError("raw_saturation_mask must be boolean")

    if raw_saturation_mask is None and saturation_value is None:
        if not np.issubdtype(raw_green.dtype, np.integer):
            raise ValueError(
                "saturation_value must be given explicitly for a non-integer raw_green "
                "dtype; it can't be inferred from the array's own max (see docstring)."
            )
        saturation_value = float(np.iinfo(raw_green.dtype).max)

    features = []
    # Every morphological step below reaches at most erosion_px + annulus_width_px
    # pixels from the cell, so work in the cell's bounding box padded by that much.
    # Running erosion/dilation on the full frame for each cell made extraction
    # take ~70 s per 2048x1536 field; the cropped version gives identical values.
    pad = erosion_px + annulus_width_px + 1
    height, width = label_mask.shape
    for region in regionprops(label_mask):
        cell_id = int(region.label)
        r0, c0, r1, c1 = region.bbox
        box = (slice(max(r0 - pad, 0), min(r1 + pad, height)), slice(max(c0 - pad, 0), min(c1 + pad, width)))
        labels_box = label_mask[box]
        mask = labels_box == cell_id

        measure_mask = (
            erosion(mask, disk(erosion_px), mode="ignore")
            if erosion_px > 0
            else mask
        )
        if not measure_mask.any():
            measure_mask = mask  # too small to erode -- fall back to the full cell

        raw_g, raw_r = raw_green[box][measure_mask], raw_red[box][measure_mask]
        corr_g, corr_r = corrected_green[box][measure_mask], corrected_red[box][measure_mask]

        local_bg_g = local_bg_r = 0.0
        if subtract_local_background:
            outer = dilation(
                mask,
                disk(erosion_px + annulus_width_px),
                mode="ignore",
            )
            inner = dilation(
                mask,
                disk(erosion_px),
                mode="ignore",
            )
            annulus = outer & ~inner & (labels_box == 0)  # exclude other cells' territory too
            if annulus.any():
                local_bg_g = float(np.median(corrected_green[box][annulus]))
                local_bg_r = float(np.median(corrected_red[box][annulus]))
                corr_g = np.clip(corr_g - local_bg_g, 0, None)
                corr_r = np.clip(corr_r - local_bg_r, 0, None)



        if raw_saturation_mask is not None:
            sat = bool(raw_saturation_mask[box][mask].any())
        else:
            sat = bool(
                (raw_green[box][mask] >= saturation_value).any()
                or (raw_red[box][mask] >= saturation_value).any()
            )
        features.append(
            CellFeatures(
                cell_id=cell_id,
                **((cell_focus or {}).get(cell_id, {})),
                area_px=int(region.area),
                raw_mean_green=float(raw_g.mean()),
                raw_mean_red=float(raw_r.mean()),
                corrected_mean_green=float(corr_g.mean()),
                corrected_mean_red=float(corr_r.mean()),
                integrated_green=float(corr_g.sum()),
                integrated_red=float(corr_r.sum()),
                sat_flag=sat,
                local_background_green=local_bg_g,
                local_background_red=local_bg_r,
                sat_flag_legacy=(bool(legacy_saturation_mask[box][mask].any())
                                 if legacy_saturation_mask is not None else sat),
            )
        )
    return features


def compute_ratio(corrected_mean_green: float, corrected_mean_red: float) -> float | None:
    """Raw corrected red/green ratio; None when green is non-positive (undefined ratio)."""
    if corrected_mean_green <= 0:
        return None
    return corrected_mean_red / corrected_mean_green
