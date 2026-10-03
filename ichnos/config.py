"""Calibrated pipeline constants, with their provenance, in one place instead
of scattered as literals inside function signatures.

Everything here was measured, not guessed -- see each constant's comment for
the script that derived it and its caveats (sample size, what it assumes).
Re-run the cited script and update the constant (and this comment) if the
underlying assumptions change -- e.g. once real wet-lab images replace the
public reference dataset these were calibrated on.
"""
from __future__ import annotations

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

# Stage 1/3: pixel size (µm/pixel) by objective, for
# correct.suggest_rolling_ball_radius(). These are REFERENCE values measured
# from the real YRC public-dataset images in data/raw/metadata.csv (a
# different microscope, not the team's Olympus) -- placeholders to unblock
# building the objective-aware mechanism now. Replace/extend with the real
# Olympus objective(s) + sensor pixel size once known, and re-run
# ichnos_image/scripts/calibrate_rolling_ball_radius.py against real Olympus
# images to confirm the suggested radius actually holds -- see that script's
# findings for why an empirical check matters, not just the formula.
OBJECTIVE_PIXEL_SIZE_UM_REFERENCE = {
    "60X": 0.10760,  # YRC Tkach-screen dual-channel entries (e.g. 210071, 201943)
    "100X": 0.13000,  # YRC single-GFP entries (e.g. 179997, 165478, 182391)
}

# Typical yeast (S. cerevisiae) cell diameter, µm -- used as the default in
# correct.suggest_rolling_ball_radius(). Per the team's own figure.
YEAST_CELL_DIAMETER_UM = 5.0

# Olympus SC30 camera (the team's real camera, per its manual -- extracted
# 2026-09-21). Sensor pixel size at the chip, NOT at the sample -- combine
# with an objective + camera adapter magnification via
# correct.pixel_size_at_sample_um() to get the real µm/pixel for a given
# setup. Distinct from OBJECTIVE_PIXEL_SIZE_UM_REFERENCE above, which is
# borrowed from the unrelated YRC public dataset's own (different) camera.
CAMERA_SC30_SENSOR_PIXEL_SIZE_UM = 3.2
CAMERA_SC30_CHIP_SIZE = '1/2"'
CAMERA_SC30_EFFECTIVE_AREA_MM = (6.55, 4.92)

# 24-bit total RGB = 8-bit/channel -- i.e. saturation is 255, NOT 1023
# (10-bit) or 65535 (16-bit, the public reference images' depth). If a real
# SC30 file ever gets loaded into a wider container (e.g. an 8-bit source
# saved/read as uint16), extract.extract_per_cell()'s dtype-based
# saturation_value auto-detection would silently infer 65535 and be wrong --
# pass saturation_value=CAMERA_SC30_SATURATION_VALUE explicitly for real
# Olympus/SC30 data rather than trusting auto-detection.
CAMERA_SC30_BIT_DEPTH = 8
CAMERA_SC30_SATURATION_VALUE = 255.0

# No thermoelectric cooling on the SC30 (manual states only a 0-45C
# operating range, no active cooling) -- more dark current/thermal noise
# than a cooled camera, especially at longer exposures (e.g. the 2.00s
# mCherry exposure seen on the public reference images). Doesn't change any
# code here directly (Stage 3's background subtraction is already per-image
# and dynamic, not a fixed assumed value), but raises the priority of
# measuring this camera's real Photon Transfer Curve (read noise, gain, QE
# are all still unknown) before trusting absolute intensity comparisons
# across sessions/exposures.
CAMERA_SC30_ACTIVE_COOLING = False

# Confirmed independently from the actual SC30 installation manual
# (2026-09-21, not just the colleague's summary of it) -- matches exactly:
# CMOS, 1/2" chip, 3.2x3.2um pixels, 6.55x4.92mm effective area, 24-bit RGB
# (8-bit/channel), C-mount, USB 2.0, 0-45C operating range (no cooling
# mentioned). New from this manual:
CAMERA_SC30_MAX_RESOLUTION_PX = (2048, 1532)

# Binning combines neighboring sensor pixels into blocks -- trades
# resolution for sensitivity/speed. This changes the EFFECTIVE sensor pixel
# size fed into pixel_size_at_sample_um(): at 2x binning a "pixel" is
# physically 2x as large (6.4um, not 3.2um), etc. {binning_factor:
# (resolution_px, max_fps, exposure_range_s)}, per the manual's table.
CAMERA_SC30_BINNING_MODES = {
    1: {"resolution_px": (2048, 1532), "max_fps": 10, "exposure_s": (0.00061, 1.74)},
    2: {"resolution_px": (1024, 768), "max_fps": 28, "exposure_s": (0.00044, 0.97)},
    3: {"resolution_px": (680, 512), "max_fps": 37, "exposure_s": (0.00049, 0.99)},
    4: {"resolution_px": (508, 384), "max_fps": 49, "exposure_s": (0.00049, 0.94)},
}

# Real C-mount camera adapters (sit between objective and the C-mount SC30)
# listed in the CKX41/CKX31 and Camera Adapter System manuals as usable on
# this microscope. This is the list of *plausible* options -- confirm which
# one the team's physical adapter actually is (usually printed on the
# adapter barrel itself) before trusting a specific magnification.
KNOWN_CAMERA_ADAPTER_MAGNIFICATIONS = {
    "U-TV0.25XC": 0.25,
    "U-TV0.35XC-2": 0.35,
    "U-TV0.5XC-3": 0.5,
    "U-TV0.63XC": 0.63,
    "U-TV1XC": 1.0,
}

# CKX41/CKX31 culture microscope -- the team's real microscope, confirmed
# from its manual (an inverted, phase-contrast-only scope for culture
# dishes -- no DIC). UIS2 objectives actually available on this model, with
# NA / working distance (mm) / resolution (um) as stated. "Ph"-suffixed
# objectives are for the IX2-SL/IX2-SLP phase-contrast sliders.
# LUCPlanFLN (long-working-distance, correction-collar) objectives exist
# specifically because standard culture vessels are much thicker than a
# coverslip: PlanCN60X/100XO's 0.2/0.13mm working distance is impractical
# through a normal dish bottom -- the LUC objectives (6.6-7.8mm at 20X) are
# the ones actually usable for imaging through a culture dish.
CKX41_OBJECTIVES = {
    "PlanCN4X": {"na": 0.10, "wd_mm": 18.5, "resolution_um": 3.36},
    "PlanCN10X": {"na": 0.25, "wd_mm": 10.5, "resolution_um": 1.30},
    "PlanCN20X": {"na": 0.40, "wd_mm": 1.2, "resolution_um": 0.84},
    "PlanCN40X": {"na": 0.65, "wd_mm": 0.6, "resolution_um": 0.54},
    "PlanCN60X": {"na": 0.80, "wd_mm": 0.2, "resolution_um": 0.42},
    "PlanCN100XO": {"na": 1.25, "wd_mm": 0.13, "resolution_um": 0.27},  # oil immersion
    "LUCPlanFLN20X": {"na": 0.45, "wd_mm": (6.6, 7.8), "resolution_um": 0.75},  # correction collar
    "LUCPlanFLN40X": {"na": 0.60, "wd_mm": (2.7, 4.0), "resolution_um": 0.56},  # correction collar
    "LUCPlanFLN60X": {"na": 0.70, "wd_mm": (1.5, 2.2), "resolution_um": 0.48},  # correction collar
}

# Reflected fluorescence system for CKX41 (mercury HBO50W burner) -- the
# team's real B/G excitation filter cubes, confirmed from its manual.
# B-excitation matches GFP-family fluorophores (EGFP, S65T, RSGFP);
# G-excitation matches mCherry/RFP-family -- this is the real filter set
# standing in for what the public YRC reference images' 488/540 (GFP) and
# 561/600 (mCherry) approximate.
CKX41_FLUORESCENCE_FILTERS = {
    "B": {"excitation_nm": (460, 490), "dichroic_nm": 500, "barrier_nm": 520,
          "applications": "FITC/GFP-family (EGFP, S65T, RSGFP), acridine orange, auramine"},
    "G": {"excitation_nm": (480, 550), "dichroic_nm": 570, "barrier_nm": 590,
          "applications": "Rhodamine/TRITC, RFP/mCherry-family, propidium iodide"},
}
