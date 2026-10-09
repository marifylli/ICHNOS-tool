"""The imaging hardware the team actually uses, as data rather than as
constants scattered through the pipeline.

Everything here comes from the equipment survey (September 2026), which
recorded the setup from photographs of the hardware, of the acquisition
software, and from two TIFF files the lab sent. Where that survey marked
something as unknown, the field here is None -- a missing value is recorded
as missing and never filled with a plausible default, because a wrong
default is indistinguishable from a measurement downstream.

This module replaces the Olympus SC30 constants that the image pipeline
previously carried. The SC30 is not this lab's camera; the real camera is a
QImaging MicroPublisher 3.3 RTV. The old values are in this repository's git
history, and in DryLabTool, if that camera ever needs supporting.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


# --------------------------------------------------------------------------
# Objectives
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Objective:
    name: str
    magnification: float
    numerical_aperture: float
    phase_ring: Optional[str]  # phase-contrast annulus this objective needs
    correction_collar_mm: Optional[tuple[float, float]]


# The four objectives physically on the turret, by their colour band.
# The 60X NA is 0.70; an earlier survey draft recorded 0.20 in error.
OBJECTIVES: dict[str, Objective] = {
    "4X": Objective("UPlanFLN 4X", 4.0, 0.13, "PhL", None),
    "10X": Objective("CPlanFLN 10X", 10.0, 0.30, "PhC", None),
    "40X": Objective("LUCPlanFLN 40X", 40.0, 0.60, "Ph2", (0.0, 2.0)),
    "60X": Objective("LUCPlanFLN 60X", 60.0, 0.70, None, (0.1, 1.3)),
}

# Objective used for the fluorescence sessions of 5-9 October 2026 (H2O2 and
# DTT, transfected and non-transfected): 40X, confirmed by the team on
# 2026-10-10. The correction-collar position for the vessel is still open.
FLUORESCENCE_OBJECTIVE: Optional[str] = "40X"


# --------------------------------------------------------------------------
# Filter cubes
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class FilterCube:
    position: str
    used_for: Optional[str]


# Slider positions U, B, G. GFP is imaged through B, mCherry through G.
# U (DAPI-type dyes) is not used in this experiment.
FILTER_CUBES: dict[str, FilterCube] = {
    "U": FilterCube("U", None),
    "B": FilterCube("B", "green"),
    "G": FilterCube("G", "red"),
}

CUBE_FOR_CHANNEL = {"green": "B", "red": "G"}

# A D460/50M emission filter (435-485 nm, visibly worn) was found loose in a
# holder. It is not a standard CKX41 filter. The working assumption is that
# it belongs to the U cube, in which case it does not sit in the GFP or
# mCherry light path -- but this is unconfirmed, and until it is, any claim
# that the emission path is fully characterised is unsupported.
UNPLACED_FILTER = "D460/50M"


# --------------------------------------------------------------------------
# Camera
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Camera:
    model: str
    sensor: str
    sensor_pixel_size_um: float
    full_resolution_px: tuple[int, int]
    adc_bit_depth: int
    actively_cooled: bool
    is_colour: bool


CAMERA = Camera(
    model="QImaging MicroPublisher 3.3 RTV",
    sensor="Sony ICX252 CCD",
    sensor_pixel_size_um=3.45,
    full_resolution_px=(2048, 1536),
    adc_bit_depth=10,
    actively_cooled=True,  # Peltier, roughly 10 C below ambient
    is_colour=True,
)

# Olympus U-TV0.5XC-3, confirmed from the label on the adapter barrel.
CAMERA_ADAPTER_MAGNIFICATION = 0.5


# --------------------------------------------------------------------------
# Digitisation vs. storage -- these are two different numbers
# --------------------------------------------------------------------------

# The camera digitises at 10 bit (0-1023). The acquisition preset saves
# 24-bit colour, i.e. 8 bit per RGB component (0-255), so two bits are
# discarded on the way to disk. The value a loaded file can reach is
# therefore 255, not 1023 and not 65535.
#
# These must not be used interchangeably. SATURATION_VALUE is the one QC
# should test against, because it describes the stored file. ADC_MAX_VALUE
# describes the sensor and is kept for noise/dynamic-range reasoning.
SATURATION_VALUE = 255.0
SAVED_BIT_DEPTH_PER_CHANNEL = 8
ADC_MAX_VALUE = 2 ** CAMERA.adc_bit_depth - 1  # 1023

# A pixel at SATURATION_VALUE means the *stored* signal clipped. It does not
# prove the sensor saturated: gain, white balance and any post-capture LUT
# can clip earlier in the chain.
SATURATION_MEANS_SENSOR_SATURATED = False


# --------------------------------------------------------------------------
# Acquisition preset 6, "QI Fluoro (Gain)"
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class AcquisitionPreset:
    name: str
    gain: float
    gamma: float
    offset: float
    white_balance_rgb: tuple[float, float, float]
    binning: int
    capture_resolution_px: tuple[int, int]
    saturation_warning_enabled: bool
    post_snap_macro: Optional[str]
    post_snap_macro_confirmed: bool


FLUORESCENCE_PRESET = AcquisitionPreset(
    name="QI Fluoro (Gain)",  # preset 6
    gain=3.0,
    gamma=1.0,
    offset=0.0,
    white_balance_rgb=(5.0, 3.0, 2.0),
    binning=1,
    capture_resolution_px=(2048, 1536),
    saturation_warning_enabled=False,
    # CAM_APPLY_LUT is confirmed on preset 3 (phase contrast). Whether it
    # also runs on preset 6 was never photographed. Recorded as unknown.
    post_snap_macro=None,
    post_snap_macro_confirmed=False,
)

# gamma=1 is necessary but not sufficient for a linear response. In the two
# sample TIFFs, one intensity value in five never occurs in the fluorescence
# image, which is a signature of rescaling before saving -- consistent with
# a post-snap LUT. Until the preset 6 macro tab is photographed, files from
# this setup are NOT validated raw linear measurements.
RESPONSE_IS_VALIDATED_LINEAR = False

# Offset 0 clips the dark background at zero: 13% of red and 43% of blue
# pixels in the sample files sit at exactly 0. Clipped background biases any
# background estimate upward and breaks the assumption that read noise is
# symmetric about the background level.
BACKGROUND_IS_CLIPPED_AT_ZERO = True

# "Adjust Exp for Binning" is ticked, and preview runs at 2x2 while capture
# runs at 1x1, so the capture exposure comes out as 4x the preview exposure
# (100/400 ms and 300/1200 ms were both observed in the same preset). The
# capture exposure is written into each TIFF.
#
# The pipeline reads the exposure from the file. It must NOT apply this
# factor again -- the value in the TIFF is already the capture exposure.
EXPOSURE_IS_FOUR_TIMES_PREVIEW = True
EXPOSURE_READ_FROM_FILE = True

# Image-Pro Plus writes a default channel name of "Alexa Fluor 560" into
# every file regardless of the cube in use. It is meaningless here and must
# never be used to identify a channel; the image manifest says which cube
# each acquisition used.
UNRELIABLE_TIFF_CHANNEL_NAME = "Alexa Fluor 560"


# --------------------------------------------------------------------------
# Illumination
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class FluorescenceLamp:
    model: str
    bulb: str
    watts: float
    fibre_coupled: bool
    intensity_settings_percent: tuple[int, ...]


FLUORESCENCE_LAMP = FluorescenceLamp(
    model="Olympus U-HGLGPS",
    bulb="SHI-130 OL",
    watts=130.0,
    fibre_coupled=True,
    intensity_settings_percent=(0, 3, 6, 12, 25, 50, 100),
)

# The intensity dial changes fluorescence brightness directly, so it has to
# be recorded per session. Which position the lab actually shoots at is not
# known; the dial was at 0 when photographed, which is the off position.
LAMP_INTENSITY_PERCENT: Optional[int] = None

# The hour counter display was unlit in the photograph and could not be read.
LAMP_BURNER_HOURS: Optional[float] = None

TRANSMITTED_LIGHT_LAMP = "Olympus U-LS30-3 (6V30WHAL)"

MICROSCOPE = "Olympus CKX41"
PHASE_SLIDER = "Olympus IX2-SL"  # positions PhL and Ph2; needs centring
ACQUISITION_SOFTWARE = "Image-Pro Plus 7.0.1.658"
CAMERA_DRIVER = "QImaging 7.0.3.23"


# --------------------------------------------------------------------------
# Derived geometry
# --------------------------------------------------------------------------


def pixel_size_um(objective: str, binning: int = 1) -> float:
    """Sample-plane µm per pixel for one of this microscope's objectives.

    sensor_pixel_size * binning / (objective_magnification * adapter). The
    adapter is the confirmed 0.5X U-TV0.5XC-3, so the sample pixel is twice
    as large as the objective magnification alone would suggest.
    """
    key = objective.upper()
    if key not in OBJECTIVES:
        raise KeyError(
            f"{objective!r} is not an objective on this microscope; "
            f"available: {sorted(OBJECTIVES)}"
        )
    effective_sensor_pixel_um = CAMERA.sensor_pixel_size_um * binning
    return effective_sensor_pixel_um / (
        OBJECTIVES[key].magnification * CAMERA_ADAPTER_MAGNIFICATION
    )


OBJECTIVE_PIXEL_SIZE_UM: dict[str, float] = {
    name: pixel_size_um(name) for name in OBJECTIVES
}


# --------------------------------------------------------------------------
# What is still missing before this profile supports quantitative claims
# --------------------------------------------------------------------------

# Hardware is identified. Measurement calibration is a separate question and
# is not complete. Each entry below blocks a specific quantitative claim, not
# the pipeline running.
UNRESOLVED: dict[str, str] = {
    "fluorescence_objective": "which objective fluorescence is taken with; sets the spatial scale",
    "correction_collar": "collar position for the chosen vessel",
    "lamp_intensity": "dial position used for acquisition; scales brightness directly",
    "lamp_burner_hours": "hour counter could not be read",
    "exposure_per_channel": "agreed fixed exposure for green and red",
    "gain_mcherry": "gain used for the red channel",
    "preset6_macro": "whether CAM_APPLY_LUT runs after a fluorescence snap",
    "save_format": "the Save As format and capture-depth options actually used",
    "d460_filter": "where D460/50M sits, and that it is out of the GFP/mCherry path",
    "rgb_extraction": "which 2D component or calibrated combination is used per cube",
    "dark_frames": "dark frames per exposure",
    "photon_transfer_curve": "gain, read noise and detection floor",
    "flat_fields": "per-session illumination reference",
    "crosstalk_controls": "GFP-only and mCherry-only strains",
}


def is_quantitatively_calibrated() -> bool:
    """False while anything in UNRESOLVED is open.

    Absolute intensity comparisons across sessions, and any uncertainty
    reported as a real experimental accuracy rather than a provisional
    figure, depend on this being True.
    """
    return not UNRESOLVED


# --------------------------------------------------------------------------
# Calibrated QC thresholds
# --------------------------------------------------------------------------
# Moved here from ichnos/config.py in Step 3, when that module became the
# model/variant configuration. These are properties of the image analysis as
# calibrated against a specific acquisition setup, so they belong on the
# image side and must not collide with model config.
#
# All of them were calibrated on 16-bit public reference images from a
# different microscope. They are the best available starting points, not
# measurements of this rig, and each needs recalibrating once real images
# from the QImaging camera exist.

# Stage 7 QC: segment.focus_score() (variance of Laplacian, normalized to the
# image's own dynamic range -- see that function's docstring for why: raw,
# non-normalized variance would be ~66000x different between an 8-bit and a
# 16-bit camera for similar-looking images, which matters concretely here
# since the team's real Olympus camera (SC30, see CAMERA_SC30_* below) is
# 8-bit/channel while the public reference images this was calibrated on are
# 16-bit) below this is flagged. Calibrated by
# ichnos_image/scripts/calibrate_focus_threshold.py: applies synthetic blur
# to the 3 real DIC reference images (179997/165478/182391) and finds where
# segmentation mask IoU vs. the sharp baseline first drops below 0.8; this is
# the mean focus_score at that break point across the 3 images.
# Small sample (3 images, one objective) -- recalibrate with real wet-lab
# images once available. Re-run the calibration script (not just rescale this
# number) after any further change to focus_score()'s definition.
FOCUS_SCORE_THRESHOLD = 5.580033091722408e-06

# Stage 7 QC: correct.estimate_registration_shift()'s magnitude (px) above
# this is flagged. Calibrated by
# ichnos_image/scripts/calibrate_registration_threshold.py: applies a known
# synthetic green/red shift to real GFP images (179997/165478) and finds
# where the resulting per-cell ratio's error vs. true ratio first exceeds
# 10%; this is the more conservative (smaller) of the two images' break
# points. The 10% error tolerance is a reasonable default, not a requirement
# from the team -- tighten it if a stricter accuracy target is set.
REGISTRATION_SHIFT_THRESHOLD_PX = 3.5

# Stage 7 QC: images taken less than this many minutes after lamp ignition
# are flagged (mercury/xenon burners drift in intensity while warming up).
# Not calibrated here -- this is the team's own stated protocol figure,
# applied as a QC rule rather than re-derived.
LAMP_WARMUP_THRESHOLD_MINUTES = 15.0

# Stage 3 background estimation: only trust the histogram-mode background
# estimate if it falls at or below this percentile of the image; above that,
# fall back to the plain low-percentile estimate. Calibrated empirically in
# ichnos_image/tests/test_correct.py: the bare mode estimator latches onto
# the cell-intensity peak instead of the background peak once cells cover
# roughly 60%+ of the field, and this sanity bound catches that before it
# happens (background stays below the 20th percentile in every density
# tested up to that point).
BACKGROUND_MODE_SANITY_PERCENTILE = 20.0

# Typical yeast (S. cerevisiae) cell diameter, µm -- the default in
# correct.suggest_rolling_ball_radius(). Per the team's own figure.
YEAST_CELL_DIAMETER_UM = 5.0
