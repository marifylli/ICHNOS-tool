"""extract_per_cell works in each cell's padded bounding box; it must give
exactly what the full-frame computation gave, including for cells touching
the image border and cells whose annuli overlap a neighbour."""
import numpy as np
from skimage.measure import regionprops
from skimage.morphology import dilation, disk, erosion

from ichnos_image import extract


def _full_frame_reference(label_mask, raw_g_img, raw_r_img, corr_g_img, corr_r_img, saturation_value,
                          erosion_px=1, annulus_width_px=3):
    """The previous implementation: every step on the whole frame."""
    out = {}
    for region in regionprops(label_mask):
        mask = label_mask == region.label
        measure = erosion(mask, disk(erosion_px), mode="ignore")
        if not measure.any():
            measure = mask
        corr_g, corr_r = corr_g_img[measure], corr_r_img[measure]
        outer = dilation(mask, disk(erosion_px + annulus_width_px), mode="ignore")
        inner = dilation(mask, disk(erosion_px), mode="ignore")
        annulus = outer & ~inner & (label_mask == 0)
        bg_g = bg_r = 0.0
        if annulus.any():
            bg_g, bg_r = float(np.median(corr_g_img[annulus])), float(np.median(corr_r_img[annulus]))
            corr_g, corr_r = np.clip(corr_g - bg_g, 0, None), np.clip(corr_r - bg_r, 0, None)
        sat = bool((raw_g_img[mask] >= saturation_value).any() or (raw_r_img[mask] >= saturation_value).any())
        out[region.label] = (float(raw_g_img[measure].mean()), float(raw_r_img[measure].mean()),
                             float(corr_g.mean()), float(corr_r.mean()), float(corr_g.sum()),
                             float(corr_r.sum()), sat, bg_g, bg_r)
    return out


def test_cropped_extraction_equals_full_frame():
    rng = np.random.default_rng(20261010)
    labels = np.zeros((60, 80), dtype=np.int32)
    labels[0:6, 0:7] = 1          # top-left corner
    labels[20:28, 30:38] = 2      # interior
    labels[20:28, 39:45] = 3      # neighbour one pixel away from 2
    labels[54:60, 70:80] = 4      # bottom-right corner
    labels[10, 50] = 5            # single pixel, too small to erode
    raw_g = rng.integers(0, 256, labels.shape).astype(float)
    raw_r = rng.integers(0, 256, labels.shape).astype(float)
    corr_g, corr_r = raw_g - 20.0, raw_r - 30.0

    got = {f.cell_id: f for f in extract.extract_per_cell(labels, raw_g, raw_r, corr_g, corr_r, saturation_value=255)}
    expected = _full_frame_reference(labels, raw_g, raw_r, corr_g, corr_r, 255)

    assert set(got) == set(expected)
    for cell_id, values in expected.items():
        f = got[cell_id]
        assert (f.raw_mean_green, f.raw_mean_red, f.corrected_mean_green, f.corrected_mean_red,
                f.integrated_green, f.integrated_red, f.sat_flag, f.local_background_green,
                f.local_background_red) == values
