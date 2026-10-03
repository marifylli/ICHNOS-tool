"""Stage 5 feature extraction: local-background-annulus correction should
recover a cell's true signal even where Stage 3's single global background
level misses a local residual offset (e.g. diffuse glow near a bright spot).
"""
import numpy as np
from scipy import ndimage as ndi

from ichnos_image import extract


def _two_cell_scene(size=150):
    yy, xx = np.mgrid[0:size, 0:size]

    global_bg = 50.0
    glow_offset = 80.0  # local residual a single global background scalar can't catch
    cell_signal = 2000.0

    raw_green = np.full((size, size), global_bg)
    glow_region = ((yy - 40) ** 2 + (xx - 40) ** 2) < 30**2
    raw_green[glow_region] += glow_offset

    cell_a = ((yy - 40) ** 2 + (xx - 40) ** 2) < 8**2  # inside the glow
    cell_b = ((yy - 110) ** 2 + (xx - 110) ** 2) < 8**2  # clean region, no glow
    raw_green[cell_a] += cell_signal
    raw_green[cell_b] += cell_signal

    raw_red = np.zeros((size, size))
    label_mask, _ = ndi.label(cell_a | cell_b)

    # simulate Stage 3's global-only correction: removes the flat 50, leaves
    # the local +80 glow untouched -- that's what local-background-annulus
    # correction in extract_per_cell is meant to clean up
    corrected_green = np.clip(raw_green - global_bg, 0, None)
    corrected_red = raw_red.copy()

    return label_mask, raw_green, raw_red, corrected_green, corrected_red, cell_signal


def test_local_background_annulus_corrects_residual_local_offset():
    label_mask, raw_green, raw_red, corrected_green, corrected_red, cell_signal = _two_cell_scene()

    with_local = extract.extract_per_cell(
        label_mask, raw_green, raw_red, corrected_green, corrected_red,
        subtract_local_background=True, saturation_value=65535.0,
    )
    without_local = extract.extract_per_cell(
        label_mask, raw_green, raw_red, corrected_green, corrected_red,
        subtract_local_background=False, saturation_value=65535.0,
    )

    by_id_with = {f.cell_id: f for f in with_local}
    by_id_without = {f.cell_id: f for f in without_local}
    glow_id, clean_id = sorted(by_id_with)  # label order: glow cell (near origin) first

    # without local-background correction, the glow cell's mean is inflated
    # by the ~80 residual offset on top of the true cell_signal
    assert by_id_without[glow_id].corrected_mean_green > cell_signal + 40
    # with it, the glow cell's mean lands close to the true cell_signal
    assert abs(by_id_with[glow_id].corrected_mean_green - cell_signal) < 40
    assert by_id_with[glow_id].local_background_green > 40

    # the clean cell (no local glow nearby) should be about the same either way
    assert abs(by_id_with[clean_id].corrected_mean_green - by_id_without[clean_id].corrected_mean_green) < 15
    assert by_id_with[clean_id].local_background_green < 15


def test_saturation_value_required_for_float_channels():
    """Regression test: saturation_value used to silently default to the
    array's own max for float channels, which trivially flagged whichever
    cell happened to be brightest as "saturated" -- caught via
    tests/test_pipeline.py once it started checking qc_pass. Now it's a
    clear error instead of a silent wrong answer.
    """
    import pytest

    label_mask, raw_green, raw_red, corrected_green, corrected_red, _ = _two_cell_scene()
    assert raw_green.dtype == np.float64
    with pytest.raises(ValueError):
        extract.extract_per_cell(label_mask, raw_green, raw_red, corrected_green, corrected_red)
