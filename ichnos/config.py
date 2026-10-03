"""Calibrated QC constants for the image pipeline, with their provenance, in
one place instead of scattered as literals inside function signatures.

Hardware description -- camera, objectives, filter cubes, lamp, acquisition
preset -- is NOT here. It lives in ichnos_image/instrument.py, which
describes one specific microscope. This module holds thresholds that are
properties of the analysis, not of the instrument.

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

# Typical yeast (S. cerevisiae) cell diameter, µm -- used as the default in
# correct.suggest_rolling_ball_radius(). Per the team's own figure.
YEAST_CELL_DIAMETER_UM = 5.0
