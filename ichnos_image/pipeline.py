"""End-to-end orchestration over a batch of image sets (Stages 2-5+8) --
the "glue" that turns segment/correct/extract/export from library functions
into something that processes a real experiment and writes one combined CSV.

Deliberately does NOT do Stage 1/4 *calibration* itself (which control image
is the session's crosstalk reference, or estimating a flat field, are human
decisions) -- calibrate those once per session with correct.calibrate_crosstalk_from_control()
and correct.estimate_flat_field(), and pass the results in.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import correct, export, extract, segment
from .schema import CellRecord


class CrosstalkCalibrationWarning(UserWarning):
    """Raised when one scalar crosstalk coefficient is applied across several
    imaging sessions. A subclass of UserWarning, so existing UserWarning
    filters still catch it; having its own category lets a test assert on
    *this* warning rather than on "any warning at all", which would otherwise
    break on unrelated third-party deprecations.
    """



@dataclass
class ImageSet:
    """One bright-field/DIC (optional) + GFP + mCherry image set to process."""

    green: np.ndarray
    red: np.ndarray
    session_id: str
    timepoint: int
    acquisition_order: int
    exposure_ms_green: float
    exposure_ms_red: float
    nd_filter_green: float
    nd_filter_red: float
    objective: str
    burner_hours: float
    lamp_warmup_minutes: float
    bright_field: np.ndarray | None = None  # None -> segment off the green channel instead
    saturation_value: float = 65535.0
    raw_saturation_mask: np.ndarray | None = None


def process_image_set(
    image_set: ImageSet,
    *,
    bleed_green_to_red: float,
    flat_field_green: np.ndarray | None = None,
    flat_field_red: np.ndarray | None = None,
    segmentation_method: str = "otsu",
    segmentation_kwargs: dict | None = None,
    background_method: str = "mode",
    rolling_ball_radius: float | None = None,
) -> list[CellRecord]:
    """Run Stages 2-5+8 on one image set, given a pre-calibrated crosstalk
    coefficient (from correct.calibrate_crosstalk_from_control(), once per
    session) and optional pre-estimated flat fields
    (correct.estimate_flat_field(), once per session/channel).

    segmentation_kwargs is forwarded to segment.segment_cells() -- e.g. pass
    {"cellpose_model": segment.load_cellpose_model()} once, from the caller,
    when segmentation_method="cellpose" over many image sets, instead of
    reloading the ~1.1GB pretrained model on every single call; or
    {"resize_factor": 0.5} for a speed/accuracy tradeoff (see
    segment.segment_cells()'s docstring for measured numbers -- safe for
    "cellpose", not for "otsu").

    background_method="rolling_ball" switches Stage 3 to
    skimage.restoration.rolling_ball instead of the scalar "mode" default --
    handles a spatially-varying background within one image, which "mode"
    can't. rolling_ball_radius defaults to
    correct.rolling_ball_radius_for_objective(image_set.objective) when not
    given explicitly; that raises KeyError if the objective isn't in
    ichnos.config.OBJECTIVE_PIXEL_SIZE_UM_REFERENCE (currently only "60X"
    and "100X", from the public reference dataset -- add the team's real
    Olympus objective(s) there, or pass rolling_ball_radius explicitly, once
    known).
    """
    segmentation_source = image_set.bright_field if image_set.bright_field is not None else image_set.green
    labels = segment.segment_cells(segmentation_source, method=segmentation_method, **(segmentation_kwargs or {}))

    green, red = image_set.green, image_set.red
    if flat_field_green is not None:
        green = correct.flat_field_correct(green, flat_field_green)
    if flat_field_red is not None:
        red = correct.flat_field_correct(red, flat_field_red)

    shift_rc = correct.estimate_registration_shift(green, red)
    registration_shift_px = float(np.hypot(*shift_rc))
    red = correct.apply_shift(red, shift_rc)

    background_kwargs = {}
    if background_method == "rolling_ball":
        radius = (
            rolling_ball_radius
            if rolling_ball_radius is not None
            else correct.rolling_ball_radius_for_objective(image_set.objective)
        )
        background_kwargs["rolling_ball_radius"] = radius

    green_bg, _ = correct.subtract_background(green, method=background_method, **background_kwargs)
    red_bg, _ = correct.subtract_background(red, method=background_method, **background_kwargs)
    green_corr, red_corr = correct.unmix_crosstalk(green_bg, red_bg, bleed_green_to_red=bleed_green_to_red)

    features = extract.extract_per_cell(
        labels, image_set.green, image_set.red, green_corr, red_corr,
        saturation_value=image_set.saturation_value,
        raw_saturation_mask=image_set.raw_saturation_mask,
    )

    edge_ids = segment.border_touching_labels(labels)
    focus = segment.focus_score(segmentation_source)

    return export.build_records(
        features,
        session_id=image_set.session_id,
        timepoint=image_set.timepoint,
        edge_flagged_ids=edge_ids,
        focus_score=focus,
        registration_shift_px=registration_shift_px,
        exposure_ms_green=image_set.exposure_ms_green,
        exposure_ms_red=image_set.exposure_ms_red,
        nd_filter_green=image_set.nd_filter_green,
        nd_filter_red=image_set.nd_filter_red,
        objective=image_set.objective,
        burner_hours=image_set.burner_hours,
        lamp_warmup_minutes=image_set.lamp_warmup_minutes,
        acquisition_order=image_set.acquisition_order,
    )


def process_experiment(
    image_sets: list[ImageSet],
    out_csv: str | Path,
    *,
    bleed_green_to_red: float | dict[str, float],
    flat_field_green: np.ndarray | dict[str, np.ndarray] | None = None,
    flat_field_red: np.ndarray | dict[str, np.ndarray] | None = None,
    segmentation_method: str = "otsu",
    segmentation_kwargs: dict | None = None,
    background_method: str = "mode",
    rolling_ball_radius: float | None = None,
) -> Path:
    """Process every image set and write one combined CSV (Stage 8, final
    dataset -- every cell, every timepoint, every session in one file).

    bleed_green_to_red / flat_field_* accept either one value used for every
    image set, or a {session_id: value} dict for a multi-session experiment
    where each session was calibrated separately (the normal case -- see
    module docstring on why calibration itself happens outside this
    function).

    Warns (doesn't raise -- this function doesn't know whether that single
    value came from a real calibration or is deliberate, e.g. a single-
    session experiment or an illustrative demo) when bleed_green_to_red is
    one scalar reused across multiple distinct sessions: per-session
    crosstalk calibration exists because lamp drift, exposure changes, etc.
    make one session's coefficient unreliable for another (see
    correct.calibrate_crosstalk_from_control()'s docstring, and the team's
    own protocol requirement that every imaging session needs its own real
    GFP-only control). scripts/run_pipeline.py enforces this harder, at the
    CLI/manifest level, by refusing to run a session with neither a control
    nor an explicit fallback.
    """
    out_csv = Path(out_csv)
    if out_csv.exists():
        out_csv.unlink()

    sessions = {s.session_id for s in image_sets}
    if not isinstance(bleed_green_to_red, dict) and len(sessions) > 1:
        import warnings

        warnings.warn(
            f"bleed_green_to_red={bleed_green_to_red} applied to all {len(sessions)} sessions "
            f"({sorted(sessions)}) -- crosstalk calibration is meant to be per-session (a real "
            "GFP-only control per imaging session), not shared across sessions. Pass a "
            "{session_id: bleed} dict from calibrate_crosstalk_from_control() per session unless "
            "this is deliberate (e.g. an illustrative demo).",
            CrosstalkCalibrationWarning,
            stacklevel=2,
        )

    for image_set in image_sets:
        bleed = (
            bleed_green_to_red[image_set.session_id]
            if isinstance(bleed_green_to_red, dict)
            else bleed_green_to_red
        )
        flat_g = (
            flat_field_green.get(image_set.session_id)
            if isinstance(flat_field_green, dict)
            else flat_field_green
        )
        flat_r = (
            flat_field_red.get(image_set.session_id) if isinstance(flat_field_red, dict) else flat_field_red
        )

        records = process_image_set(
            image_set,
            bleed_green_to_red=bleed,
            flat_field_green=flat_g,
            flat_field_red=flat_r,
            segmentation_method=segmentation_method,
            segmentation_kwargs=segmentation_kwargs,
            background_method=background_method,
            rolling_ball_radius=rolling_ball_radius,
        )
        export.export_csv(records, out_csv, append=True)

    return out_csv
