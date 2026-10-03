import numpy as np
import pytest

from ichnos_image import image_io, extract


@pytest.mark.parametrize("method", ["G", "sum", "luminance"])
@pytest.mark.parametrize(
    "level, expected",
    [(254, False), (255, True)],
)
def test_rgb_saturation_is_checked_before_extraction(
    method, level, expected
):
    rgb = np.zeros((32, 32, 3), dtype=float)
    rgb[8:24, 8:24] = 100

    # A boundary pixel in the blue component.
    rgb[8, 8, 2] = level

    labels = np.zeros((32, 32), dtype=int)
    labels[8:24, 8:24] = 1

    plane = image_io.plane_for_channel(
        rgb, "green", method=method
    )
    mask = image_io.saturation_mask_for_image(rgb, 255)

    features = extract.extract_per_cell(
        labels,
        plane,
        plane,
        plane,
        plane,
        saturation_value=255,
        raw_saturation_mask=mask,
    )

    assert features[0].sat_flag is expected