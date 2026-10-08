"""End-to-end orchestration over a batch of image sets (Stages 2-5+8) --
the "glue" that turns segment/correct/extract/export from library functions
into something that processes a real experiment and writes one combined CSV.

Deliberately does NOT do Stage 1/4 *calibration* itself (which control image
is the session's crosstalk reference, or estimating a flat field, are human
decisions) -- calibrate those once per session with correct.calibrate_crosstalk_from_control()
and correct.estimate_flat_field(), and pass the results in.
"""
from __future__ import annotations

from dataclasses import dataclass, fields, asdict, is_dataclass
from pathlib import Path

import numpy as np
import hashlib
import inspect

from ichnos.artifacts import publish_csv_bundle, code_provenance, file_fingerprints

from . import correct, export, extract, segment
from .focus import FocusPolicy, score_cells, evaluate_cell
from .schema import CellRecord

from ichnos.schema import validate_elapsed_hours, validate_optional_identifier


class SegmentationQCError(RuntimeError):
    """Segmentation claimed an implausible share of the frame (Stage 2 QC).

    Raised rather than flagged per cell, because there are no cells to flag:
    when the threshold lands inside the background every "cell" downstream is
    noise, and when it lands above the signal the few that survive are a
    brightness-biased sample of the field. Either way the per-cell rows are
    not a worse measurement of the right thing, they are a measurement of
    something else, and QC that let them through would be read as a pass.

    process_image_set raises it; _process_experiment catches it per frame,
    records the refusal and carries on. That split is deliberate: the one
    frame is unusable, the other frames are not, and an aborted batch
    discards good work that has already been paid for -- on the team's
    2026-10-05 session an abort on the last frame threw away 24 processed
    frames and half an hour.
    """


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
    green_saturation_mask: np.ndarray | None = None
    red_saturation_mask: np.ndarray | None = None

    sampling_time_hours: float | None = None
    measurement_time_hours: float | None = None
    sample_id: str | None = None
    condition_id: str | None = None
    specimen_id: str | None = None
    biological_replicate_id: str | None = None
    #: False for a frame that should contain no cells at all -- a
    #: medium-only control, a blank. The segmentation QC then inverts: an
    #: empty result is the expected one, and finding cells is the failure,
    #: because it means either contamination or a threshold low enough to
    #: be reading noise as cells on every other frame in the session.
    expect_cells: bool = True
    acquisition_json: str | None = None
    source_paths: tuple[str, ...] = ()

    def __post_init__(self):
        for name in ("sample_id", "condition_id", "specimen_id", "biological_replicate_id"):
            validate_optional_identifier(getattr(self, name), field_name=name)
        for name in ("sampling_time_hours", "measurement_time_hours"):
            setattr(
                self,
                name,
                validate_elapsed_hours(
                    getattr(self, name),
                    field_name=name,
                ),
            )


SEGMENTATION_SOURCES = ("sum", "green", "red", "brightfield")


def _segmentation_source(image_set, choice: str):
    """Pick the image the cell masks are cut from, and name it for the record.

    "sum" is the default because the headline measurement is a ratio of the
    two channels, and choosing cells by their brightness in either one tilts
    that ratio: on the team's 2026-10-05 frames the median red/green went
    from 0.9 with masks from green to 3.9 with masks from red, on the same
    field. Adding the channels treats them alike, so the selection no longer
    pushes the ratio one way.

    It is not a neutral choice, only a symmetric one -- a cell visible in
    neither channel is still missed. "brightfield" is the genuinely
    independent option and is preferable when a bright-field frame of the
    same field of view exists, which in that session it did not.
    """
    if choice not in SEGMENTATION_SOURCES:
        raise ValueError(f"segmentation_source must be one of {SEGMENTATION_SOURCES}, got {choice!r}")
    if choice == "brightfield":
        if image_set.bright_field is None:
            raise ValueError(
                "segmentation_source='brightfield' needs ImageSet.bright_field, which is None. "
                "Supply the bright-field frame of the same field of view, or choose 'sum'."
            )
        return image_set.bright_field, "brightfield"
    if choice == "green":
        return image_set.green, "green"
    if choice == "red":
        return image_set.red, "red"
    return image_set.green + image_set.red, "green+red"


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
    focus_policy: FocusPolicy | None = None,
    segmentation_source: str = "sum",
    min_foreground_fraction: float = 0.001,
    max_foreground_fraction: float = 0.2,
    blank_max_foreground_fraction: float = 0.01,
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
    if (image_set.green.ndim != 2 or image_set.red.shape != image_set.green.shape
            or not np.isfinite(image_set.green).all() or not np.isfinite(image_set.red).all()):
        raise ValueError("green/red must be finite matching 2D planes")
    if not np.isfinite(image_set.saturation_value) or image_set.saturation_value <= 0:
        raise ValueError("saturation_value must be finite and positive")
    source_image, mask_source = _segmentation_source(image_set, segmentation_source)
    labels = segment.segment_cells(source_image, method=segmentation_method, **(segmentation_kwargs or {}))

    # Checked here, before anything runs per cell: this is the cheapest
    # moment to find out, and every later stage scales with the cell count
    # that a failure inflates. The bounds are deliberately wide -- a real
    # yeast field sits at a few percent -- so that passing says only "not
    # obviously broken", never "segmented well".
    fraction = segment.foreground_fraction(labels)
    where = (f"session {image_set.session_id!r}, sample {image_set.sample_id!r}, "
             f"method {segmentation_method!r}")
    if not image_set.expect_cells:
        if fraction > blank_max_foreground_fraction:
            raise SegmentationQCError(
                f"a frame declared cell-free claimed {fraction:.1%} of itself as cells "
                f"({labels.max()} objects), above {blank_max_foreground_fraction:.1%}; {where}. "
                "Either the frame is not blank, or the threshold is low enough to be "
                "reading noise as cells -- in which case every other frame in this "
                "session is affected too."
            )
    elif not min_foreground_fraction <= fraction <= max_foreground_fraction:
        raise SegmentationQCError(
            f"segmentation claimed {fraction:.1%} of the frame as cells "
            f"({labels.max()} objects), outside the plausible "
            f"{min_foreground_fraction:.1%}-{max_foreground_fraction:.1%}; {where}. "
            "Too high means the threshold fell inside the background and the objects "
            "are noise; too low means it fell above the signal and only the brightest "
            "cells survived. For fluorescence frames try segmentation_method='sparse'; "
            "for a frame that genuinely holds no cells set ImageSet.expect_cells=False."
        )

    green, red = image_set.green, image_set.red
    if flat_field_green is not None:
        green = correct.flat_field_correct(green, flat_field_green)
    if flat_field_red is not None:
        red = correct.flat_field_correct(red, flat_field_red)

    shift_rc = correct.estimate_registration_shift(green, red)
    registration_shift_px = float(np.hypot(*shift_rc))
    red = correct.apply_shift(red, shift_rc)

    # The green/segmentation frame is fixed; red pixels and their clipping
    # support must move together. Legacy combined masks cannot be separated:
    # conservatively treat them as potentially belonging to either channel.
    masks = []
    for name, supplied, channel in (
        ("green", image_set.green_saturation_mask, image_set.green),
        ("red", image_set.red_saturation_mask, image_set.red),
    ):
        mask = supplied if supplied is not None else channel >= image_set.saturation_value
        if mask.dtype != np.bool_ or mask.shape != image_set.green.shape:
            raise ValueError(f"{name} saturation mask must be boolean and match image shape")
        if supplied is None and image_set.raw_saturation_mask is not None:
            legacy = image_set.raw_saturation_mask
            if legacy.dtype != np.bool_ or legacy.shape != mask.shape:
                raise ValueError("raw_saturation_mask must be boolean and match image shape")
            mask = mask | legacy
        masks.append(mask)
    legacy_saturation = masks[0] | masks[1]
    registered_saturation = masks[0] | correct.register_saturation_mask(masks[1], shift_rc)

    background_kwargs = {}
    if background_method == "rolling_ball":
        radius = (
            rolling_ball_radius
            if rolling_ball_radius is not None
            else correct.rolling_ball_radius_for_objective(image_set.objective)
        )
        background_kwargs["rolling_ball_radius"] = radius

    policy = focus_policy or FocusPolicy()
    green_focus, red_focus = score_cells(labels, green), score_cells(labels, red)
    cell_focus = {}
    def finite_or_none(value):
        return float(value) if value is not None and np.isfinite(value) else None
    for g, r in zip(green_focus, red_focus):
        status, agreement, passed = evaluate_cell(g, r, policy)
        cell_focus[g.cell_id] = dict(focus_score_green=finite_or_none(g.focus_score),
            focus_score_red=finite_or_none(r.focus_score),
            contrast_to_noise_green=finite_or_none(g.contrast_to_noise),
            contrast_to_noise_red=finite_or_none(r.contrast_to_noise),
            focus_agreement=finite_or_none(agreement), focus_status=status,
            focus_qc_pass=passed, focus_qc_mode=policy.mode)

    green_bg, _ = correct.subtract_background(green, method=background_method, **background_kwargs)
    red_bg, _ = correct.subtract_background(red, method=background_method, **background_kwargs)
    green_corr, red_corr = correct.unmix_crosstalk(green_bg, red_bg, bleed_green_to_red=bleed_green_to_red)

    features = extract.extract_per_cell(
        labels, image_set.green, image_set.red, green_corr, red_corr,
        saturation_value=image_set.saturation_value,
        raw_saturation_mask=registered_saturation,
        legacy_saturation_mask=legacy_saturation,
        cell_focus=cell_focus,
    )

    edge_ids = segment.border_touching_labels(labels)
    focus = segment.focus_score(source_image)

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
        mask_source=mask_source,
        sampling_time_hours=image_set.sampling_time_hours,
        measurement_time_hours=image_set.measurement_time_hours,
        sample_id=image_set.sample_id,
        condition_id=image_set.condition_id,
        specimen_id=image_set.specimen_id, biological_replicate_id=image_set.biological_replicate_id,
        acquisition_json=image_set.acquisition_json,
    )


def _process_experiment(
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
    focus_policy: FocusPolicy | None = None,
    segmentation_source: str = "sum",
    min_foreground_fraction: float = 0.001,
    max_foreground_fraction: float = 0.2,
    blank_max_foreground_fraction: float = 0.01,
    refusals: list | None = None,
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
        raise FileExistsError(f"output path already exists: {out_csv}")

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

    refused = refusals if refusals is not None else []
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

        try:
            records = process_image_set(
                image_set,
                bleed_green_to_red=bleed,
                flat_field_green=flat_g,
                flat_field_red=flat_r,
                segmentation_method=segmentation_method,
                segmentation_kwargs=segmentation_kwargs,
                segmentation_source=segmentation_source,
                background_method=background_method,
                rolling_ball_radius=rolling_ball_radius,
                focus_policy=focus_policy,
                min_foreground_fraction=min_foreground_fraction,
                max_foreground_fraction=max_foreground_fraction,
                blank_max_foreground_fraction=blank_max_foreground_fraction,
            )
        except SegmentationQCError as error:
            # One unusable frame is not a reason to discard the frames that
            # already processed. The refusal is carried out to the manifest
            # instead, so a run that silently has fewer frames than its
            # manifest cannot be mistaken for a complete one.
            refused.append(dict(session_id=image_set.session_id, sample_id=image_set.sample_id,
                                acquisition_order=image_set.acquisition_order, reason=str(error)))
            continue
        export.export_csv(records, out_csv, append=True)

    if len(refused) == len(image_sets):
        raise SegmentationQCError(
            f"every one of the {len(image_sets)} frames was refused, so there is no "
            f"output to publish. First: {refused[0]['reason']}"
        )
    return out_csv


def process_experiment(image_sets, out_csv, *, input_paths=(), **options):
    """Process into a temporary CSV and publish with a reproducibility manifest.

    Existing CSVs/manifests are never overwritten. Array hashes cover API inputs;
    file hashes additionally cover CLI manifests, controls and source images.
    """
    out_csv = Path(out_csv)
    if out_csv.exists() or out_csv.with_suffix(out_csv.suffix + '.manifest.json').exists():
        raise FileExistsError(f"output CSV or manifest already exists: {out_csv}")
    image_sets = list(image_sets)
    if not image_sets:
        raise ValueError("image_sets must not be empty")
    sources = [*input_paths, *(p for item in image_sets for p in item.source_paths)]
    def describe(value):
        if is_dataclass(value):
            return describe(asdict(value))
        if isinstance(value, np.ndarray):
            return dict(shape=list(value.shape), dtype=str(value.dtype),
                        sha256=hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest())
        if isinstance(value, dict):
            return {str(k): describe(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [describe(v) for v in value]
        if value is None or isinstance(value, (str, bool, int, float)):
            return value
        return dict(type=type(value).__qualname__, configuration_not_serializable=True)
    def write(csv):
        fingerprints = file_fingerprints(sources)
        bound = inspect.signature(_process_experiment).bind(image_sets, csv, **options)
        bound.apply_defaults()
        effective = {k:v for k,v in bound.arguments.items() if k not in {'image_sets','out_csv'}}
        effective['focus_policy'] = effective['focus_policy'] or FocusPolicy()
        from .instrument import FOCUS_SCORE_THRESHOLD, REGISTRATION_SHIFT_THRESHOLD_PX, LAMP_WARMUP_THRESHOLD_MINUTES
        metadata = dict(schema_version=2, kind="image_processing_run", code=code_provenance(),
            qc_thresholds=dict(field_focus=FOCUS_SCORE_THRESHOLD,registration_shift_px=REGISTRATION_SHIFT_THRESHOLD_PX,
                lamp_warmup_minutes=LAMP_WARMUP_THRESHOLD_MINUTES),
            inputs=fingerprints, images=[{f.name: describe(getattr(item, f.name))
                                         for f in fields(item)} for item in image_sets],
            options=describe({k: v for k, v in effective.items() if k != "refusals"}),
            saturation_registration="linear-support, red-to-green frame",
            photobleaching_corrected=False)
        refusals = []
        _process_experiment(image_sets, csv, refusals=refusals, **options)
        metadata["refused_image_sets"] = refusals
        if file_fingerprints(sources) != fingerprints:
            raise ValueError("input files changed during processing")
        return metadata
    return publish_csv_bundle(out_csv, write)
