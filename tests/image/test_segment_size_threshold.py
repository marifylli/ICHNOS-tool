import numpy as np
import pytest

from ichnos_image.segment import _remove_objects_below_size


def test_objects_at_minimum_size_are_kept():
    mask = np.zeros((5, 33), dtype=bool)
    mask[0, :29] = True
    mask[2, :30] = True
    mask[4, :31] = True

    result = _remove_objects_below_size(mask, min_size=30)

    expected = mask.copy()
    expected[0] = False

    np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize("min_size", [0, 1])
def test_zero_or_one_minimum_preserves_single_pixel(min_size):
    mask = np.zeros((3, 3), dtype=bool)
    mask[1, 1] = True

    result = _remove_objects_below_size(mask, min_size=min_size)

    np.testing.assert_array_equal(result, mask)