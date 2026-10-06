"""General-purpose CLI: run the full ichnos_image pipeline (segment ->
correct -> extract -> export) over an arbitrary experiment described by a
manifest CSV, and a per-session crosstalk-control manifest, producing one
combined CSV.

Image files are read with PIL (ichnos_image.image_io.load_image despite
its name -- Image.open handles PNG, TIFF, and most other common formats;
raw microscope formats like .czi/.nd2 would need a dedicated reader, e.g.
Bio-Formats via python-bioformats, not included here).

--manifest CSV columns (one row per image set):
    session_id, timepoint, acquisition_order, green_path, red_path,
    bright_field_path (optional, blank if none), exposure_ms_green,
    exposure_ms_red, nd_filter_green, nd_filter_red, objective, burner_hours,
    lamp_warmup_minutes
Optional timing columns:
sampling_time_hours, measurement_time_hours
Optional sample identity columns:
sample_id, condition_id (sample_id is required by the separate summary CLI)

Both are actual elapsed hours from stress onset. Missing values remain
unknown. timepoint is a point identifier, not an elapsed time in hours.

--controls CSV columns (one row per session_id needing crosstalk
calibration; a session missing from this file must instead get
--bleed-default):
    session_id, control_green_path, control_red_path

Usage:
    python scripts/run_pipeline.py --manifest images.csv --controls controls.csv --out experiment.csv
    python scripts/run_pipeline.py --manifest images.csv --bleed-default 0.05 --out experiment.csv
"""
import argparse
from pathlib import Path
from skimage.color import rgb2gray

import pandas as pd
from ichnos_image.instrument import SATURATION_VALUE
from ichnos.schema import validate_elapsed_hours, validate_optional_identifier

from ichnos_image import ImageSet, process_experiment
from ichnos_image.correct import calibrate_crosstalk_from_control
from ichnos_image.image_io import (
    EXTRACTION_METHODS,
    load_image,
    plane_for_channel,
    saturation_mask_for_image,
)

def _load_bright_field(path: str | Path):
    """Convert bright-field images to a 2D segmentation plane."""
    image = load_image(path)

    if image.ndim == 2:
        return image

    if image.ndim == 3 and image.shape[-1] in (3, 4):
        return rgb2gray(image[..., :3])

    raise ValueError(f"unsupported bright-field shape: {image.shape}")

def _optional_elapsed_hours(row, name):
    value = getattr(row, name, None)
    if value is None or pd.isna(value):
        return None

    return validate_elapsed_hours(value, field_name=name)


def _optional_identifier(row, name):
    value = getattr(row, name, None)
    if value is None or pd.isna(value):
        return None
    return validate_optional_identifier(str(value), field_name=name)

def _build_image_sets(
    manifest_path: Path,
    *,
    green_extraction: str | None = None,
    red_extraction: str | None = None,
    saturation_value: float = SATURATION_VALUE,
) -> list[ImageSet]:
    manifest = pd.read_csv(
        manifest_path,
        dtype={"session_id": str, "sample_id": str, "condition_id": str},
    )
    image_sets = []

    for row in manifest.itertuples():
        sample_id = _optional_identifier(row, "sample_id")
        condition_id = _optional_identifier(row, "condition_id")
        sampling_time_hours = _optional_elapsed_hours(
        row, "sampling_time_hours"
        )
        measurement_time_hours = _optional_elapsed_hours(
            row, "measurement_time_hours"
        )
        bright_field = (
            _load_bright_field(row.bright_field_path)
            if getattr(row, "bright_field_path", "")
            and pd.notna(row.bright_field_path)
            else None
        )
#
        raw_green = load_image(row.green_path)
        raw_red = load_image(row.red_path)

        green_saturation = saturation_mask_for_image(
            raw_green, saturation_value
        )
        red_saturation = saturation_mask_for_image(
            raw_red, saturation_value
        )

        if green_saturation.shape != red_saturation.shape:
            raise ValueError("green and red images must have matching shapes")

        raw_saturation = green_saturation | red_saturation

        image_sets.append(
            ImageSet(
                green=plane_for_channel(
                    raw_green,
                    "green",
                    method=green_extraction,
                ),
                red=plane_for_channel(
                    raw_red,
                    "red",
                    method=red_extraction,
                ),
                saturation_value=saturation_value,
                raw_saturation_mask=raw_saturation,
                bright_field=bright_field,
                session_id=str(row.session_id),
                timepoint=int(row.timepoint),
                acquisition_order=int(row.acquisition_order),
                exposure_ms_green=float(row.exposure_ms_green),
                exposure_ms_red=float(row.exposure_ms_red),
                nd_filter_green=float(row.nd_filter_green),
                nd_filter_red=float(row.nd_filter_red),
                objective=str(row.objective),
                burner_hours=float(row.burner_hours),
                lamp_warmup_minutes=float(row.lamp_warmup_minutes),
                sampling_time_hours=sampling_time_hours,
                measurement_time_hours=measurement_time_hours,
                sample_id=sample_id,
                condition_id=condition_id,
            )
        )

    return image_sets


def _calibrate_bleed_per_session(
    controls_path: Path | None,
    sessions: set[str],
    bleed_default: float | None,
    *,
    green_extraction: str | None = None,
    red_extraction: str | None = None,
) -> dict:
    bleed_by_session = {}
    if controls_path is not None:
        controls = pd.read_csv(controls_path)
        for row in controls.itertuples():
            control_green = plane_for_channel(
                load_image(row.control_green_path),
                "green",
                method=green_extraction,
            )
            control_red = plane_for_channel(
                load_image(row.control_red_path),
                "red",
                method=red_extraction,
            )
            calibration = calibrate_crosstalk_from_control(control_green, control_red)
            bleed_by_session[str(row.session_id)] = calibration["bleed_green_to_red"]
            print(
                f"calibrated session={row.session_id}: "
                f"bleed_green_to_red={calibration['bleed_green_to_red']:.4f} "
                f"residual_check={calibration['residual_check']:.2f}"
            )

    missing = sessions - set(bleed_by_session)
    if missing:
        if bleed_default is None:
            raise SystemExit(
                f"no control and no --bleed-default for session(s): {sorted(missing)}. "
                "Add them to --controls or pass --bleed-default."
            )
        for session_id in missing:
            bleed_by_session[session_id] = bleed_default

    return bleed_by_session


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", required=True, type=Path, help="CSV describing image sets to process")
    parser.add_argument("--controls", type=Path, default=None, help="CSV of per-session crosstalk controls")
    parser.add_argument(
        "--bleed-default", type=float, default=None,
        help="fallback bleed_green_to_red for sessions not in --controls (illustrative use only)",
    )
    parser.add_argument("--out", required=True, type=Path, help="output combined CSV path")
    parser.add_argument("--segmentation-method", default="otsu", choices=["otsu", "cellpose"])
    parser.add_argument(
        "--cellpose-gpu", action="store_true",
        help="use GPU for cellpose (--segmentation-method cellpose only)",
    )
    parser.add_argument(
        "--resize-factor", type=float, default=1.0,
        help="downsample before segmenting for speed (see segment.segment_cells() docstring: "
        "safe for --segmentation-method cellpose, NOT for otsu -- otsu's accuracy collapses well "
        "before 0.5)",
    )
    parser.add_argument(
        "--background-method", default="mode", choices=["mode", "percentile", "rolling_ball"],
        help="Stage 3 background subtraction method (default: mode)",
    )
    parser.add_argument(
        "--rolling-ball-radius", type=float, default=None,
        help="override the auto-computed rolling_ball radius (px); by default it's looked up per "
        "image set from its objective via ichnos.config.OBJECTIVE_PIXEL_SIZE_UM_REFERENCE "
        "(--background-method rolling_ball only)",
    )
    parser.add_argument(
        "--green-extraction",
        choices=EXTRACTION_METHODS,
        default=None,
        help="Explicit RGB extraction for green samples and controls",
    )
    parser.add_argument(
        "--red-extraction",
        choices=EXTRACTION_METHODS,
        default=None,
        help="Explicit RGB extraction for red samples and controls",
    )
    parser.add_argument(
        "--saturation-value",
        type=float,
        default=SATURATION_VALUE,
        help=(
            "Clipping threshold per stored component; "
            "default comes from the instrument profile"
        ),
    )
    args = parser.parse_args()

    image_sets = _build_image_sets(
        args.manifest,
        green_extraction=args.green_extraction,
        red_extraction=args.red_extraction,
        saturation_value=args.saturation_value,
    )

    sessions = {s.session_id for s in image_sets}

    bleed_by_session = _calibrate_bleed_per_session(
        args.controls,
        sessions,
        args.bleed_default,
        green_extraction=args.green_extraction,
        red_extraction=args.red_extraction,
    )

    segmentation_kwargs = {}
    if args.resize_factor != 1.0:
        segmentation_kwargs["resize_factor"] = args.resize_factor
    if args.segmentation_method == "cellpose":
        from ichnos_image.segment import load_cellpose_model

        print("loading cellpose model once for the whole batch...")
        segmentation_kwargs["cellpose_model"] = load_cellpose_model(gpu=args.cellpose_gpu)

    out_path = process_experiment(
        image_sets, args.out, bleed_green_to_red=bleed_by_session,
        segmentation_method=args.segmentation_method, segmentation_kwargs=segmentation_kwargs or None,
        background_method=args.background_method, rolling_ball_radius=args.rolling_ball_radius,
    )

    df = pd.read_csv(out_path)
    print(f"\n{len(image_sets)} image set(s) -> {len(df)} cell records -> {out_path}")
    print(df.groupby("session_id").agg(n_cells=("cell_id", "count"), qc_pass_frac=("qc_pass", "mean")))


if __name__ == "__main__":
    main()
