"""The instrument profile has to describe the microscope the lab actually
owns. These tests pin the values the September 2026 equipment survey
recorded, so a future edit that reintroduces a plausible-but-wrong default
fails here rather than silently rescaling every measurement downstream.
"""
import pytest

from ichnos_image import instrument


def test_camera_is_the_qimaging_not_the_sc30():
    assert "QImaging" in instrument.CAMERA.model
    assert instrument.CAMERA.sensor_pixel_size_um == 3.45  # not the SC30's 3.2
    assert instrument.CAMERA.full_resolution_px == (2048, 1536)  # SC30 is 2048x1532
    assert instrument.CAMERA.is_colour
    assert instrument.CAMERA.actively_cooled


@pytest.mark.parametrize(
    "objective, expected_um_per_px",
    [("4X", 1.725), ("10X", 0.690), ("40X", 0.1725), ("60X", 0.115)],
)
def test_pixel_sizes_match_the_survey(objective, expected_um_per_px):
    assert instrument.pixel_size_um(objective) == pytest.approx(expected_um_per_px)


def test_a_yeast_cell_spans_a_usable_number_of_pixels_at_40x():
    """Sanity check on the geometry, not just the arithmetic: a 5 µm cell
    should be tens of pixels across at the objective proposed for
    fluorescence, and only a few at 4X.
    """
    assert 25 < 5.0 / instrument.pixel_size_um("40X") < 35
    assert 5.0 / instrument.pixel_size_um("4X") < 5


def test_objective_not_on_this_microscope_raises():
    # 100X was in the old reference table, borrowed from an unrelated public
    # dataset. This turret has 4X/10X/40X/60X.
    with pytest.raises(KeyError):
        instrument.pixel_size_um("100X")


def test_binning_scales_the_sample_pixel():
    assert instrument.pixel_size_um("40X", binning=2) == pytest.approx(
        2 * instrument.pixel_size_um("40X")
    )


def test_saturation_is_the_stored_depth_not_the_adc_depth():
    """255 (8 bit per stored RGB component) and 1023 (10 bit ADC) are
    different numbers describing different things. QC tests against the
    stored file.
    """
    assert instrument.SATURATION_VALUE == 255.0
    assert instrument.ADC_MAX_VALUE == 1023
    assert instrument.SATURATION_VALUE != instrument.ADC_MAX_VALUE
    # Clipping in the stored file does not prove the sensor saturated.
    assert instrument.SATURATION_MEANS_SENSOR_SATURATED is False


def test_cube_assignment():
    assert instrument.CUBE_FOR_CHANNEL["green"] == "B"
    assert instrument.CUBE_FOR_CHANNEL["red"] == "G"
    # U is on the slider but unused in this experiment.
    assert instrument.FILTER_CUBES["U"].used_for is None


def test_unknowns_are_none_not_guessed():
    """A missing measurement stays missing. Filling any of these with a
    plausible default makes an unmeasured quantity indistinguishable from a
    measured one.
    """
    assert instrument.LAMP_INTENSITY_PERCENT is None
    assert instrument.LAMP_BURNER_HOURS is None
    assert instrument.FLUORESCENCE_PRESET.post_snap_macro is None
    assert instrument.FLUORESCENCE_PRESET.post_snap_macro_confirmed is False


def test_fluorescence_objective_is_the_confirmed_one():
    """40X was confirmed for the October 2026 sessions; it must be a turret objective."""
    assert instrument.FLUORESCENCE_OBJECTIVE == "40X"
    assert instrument.FLUORESCENCE_OBJECTIVE in instrument.OBJECTIVES


def test_profile_does_not_claim_quantitative_calibration():
    """gamma=1 does not make the response validated-linear while a post-snap
    LUT is unconfirmed, and no dark frames, PTC or flat fields exist yet.
    """
    assert instrument.RESPONSE_IS_VALIDATED_LINEAR is False
    assert instrument.BACKGROUND_IS_CLIPPED_AT_ZERO is True
    assert instrument.is_quantitatively_calibrated() is False
    for key in ("photon_transfer_curve", "flat_fields", "rgb_extraction"):
        assert key in instrument.UNRESOLVED


def test_exposure_is_read_from_the_file_not_recomputed():
    """"Adjust Exp for Binning" already produced the 4x capture exposure that
    is written into the TIFF. Applying the factor again would quadruple every
    exposure.
    """
    assert instrument.EXPOSURE_READ_FROM_FILE is True
