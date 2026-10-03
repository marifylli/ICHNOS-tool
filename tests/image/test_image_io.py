"""Loading, and the RGB-to-plane decision.

The pipeline's inputs from this setup are 24-bit RGB files, one per filter
cube. Collapsing three components into one measurement is a choice that
changes every intensity downstream, so these tests check that the choice has
to be made explicitly and that saturation is not inferred from the array.
"""
import numpy as np
import pytest

from ichnos_image import image_io, instrument


def _rgb(r, g, b, size=8):
    img = np.zeros((size, size, 3), dtype=np.float64)
    img[..., 0] = r
    img[..., 1] = g
    img[..., 2] = b
    return img


def test_two_dimensional_input_passes_through():
    plane = np.arange(16, dtype=np.float64).reshape(4, 4)
    assert np.array_equal(image_io.extract_fluorescence_plane(plane), plane)


def test_rgb_without_a_method_raises_rather_than_guessing():
    """No default. A silent luminance conversion would turn an open
    measurement question into an answer nobody checked.
    """
    with pytest.raises(image_io.UncalibratedExtractionError) as exc:
        image_io.extract_fluorescence_plane(_rgb(10, 200, 5))
    assert "no extraction method" in str(exc.value)


@pytest.mark.parametrize("method, expected", [("R", 10.0), ("G", 200.0), ("B", 5.0)])
def test_single_component_extraction(method, expected):
    plane = image_io.extract_fluorescence_plane(_rgb(10, 200, 5), method=method)
    assert plane.shape == (8, 8)
    assert np.allclose(plane, expected)


def test_sum_and_luminance_differ_from_the_green_component():
    """Recorded because it is the point: these three reductions of the same
    file are three different measurements, and the right one for this rig is
    not yet decided.
    """
    img = _rgb(10, 200, 5)
    green = image_io.extract_fluorescence_plane(img, method="G")
    total = image_io.extract_fluorescence_plane(img, method="sum")
    luma = image_io.extract_fluorescence_plane(img, method="luminance")

    assert np.allclose(total, 215.0)
    assert np.allclose(luma, 0.299 * 10 + 0.587 * 200 + 0.114 * 5)
    assert not np.allclose(green, total)
    assert not np.allclose(green, luma)


def test_unknown_method_and_unknown_cube_raise():
    with pytest.raises(ValueError):
        image_io.extract_fluorescence_plane(_rgb(1, 2, 3), method="grayscale")
    with pytest.raises(ValueError):
        image_io.extract_fluorescence_plane(_rgb(1, 2, 3), method="G", cube="Z")


def test_plane_for_channel_uses_the_right_cube():
    img = _rgb(10, 200, 5)
    # The cube is recorded, not used to pick the method -- that mapping is
    # exactly what still has to be measured.
    assert np.allclose(image_io.plane_for_channel(img, "green", method="G"), 200.0)
    with pytest.raises(image_io.UncalibratedExtractionError):
        image_io.plane_for_channel(img, "red")
    with pytest.raises(ValueError):
        image_io.plane_for_channel(img, "blue", method="B")


def test_saturation_comes_from_the_rig_not_the_dtype():
    """load_image() returns float64 for everything. Inferring saturation from
    the container gave 65535 for an 8-bit file, so clipped pixels at 255 went
    unflagged.
    """
    eight_bit_content = np.full((4, 4), 255.0)  # float64 container, 8-bit content
    assert image_io.saturation_value_for_file(eight_bit_content) == instrument.SATURATION_VALUE
    assert image_io.saturation_value_for_file(eight_bit_content) == 255.0


def test_round_trip_through_disk(tmp_path):
    path = tmp_path / "frame.png"
    original = (np.arange(64, dtype=np.float64).reshape(8, 8) * 500)
    image_io.save_png16(original, path)
    assert np.allclose(image_io.load_image(path), original)
    # historical alias still works for existing callers
    assert np.allclose(image_io.load_png16(path), original)
