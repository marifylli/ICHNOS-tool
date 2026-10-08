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
expect_cells: "false" for a frame that should hold no cells at all, such as
a medium-only control. Segmentation QC then inverts for that frame: empty is
the expected result, and finding cells is the failure. Defaults to true.

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
import json
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

from ichnos_image import manifest as manifest_io
from ichnos_image.manifest import _optional_identifier, _optional_elapsed_hours, _load_bright_field


def _build_image_sets(*args, **kwargs):
    return manifest_io._build_image_sets(*args, image_loader=load_image, **kwargs)


def _calibrate_bleed_per_session(*args, **kwargs):
    return manifest_io._calibrate_bleed_per_session(*args, image_loader=load_image,
        calibrator=calibrate_crosstalk_from_control, **kwargs)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", required=True, type=Path, help="CSV describing image sets to process")
    parser.add_argument("--controls", type=Path, default=None, help="CSV of per-session crosstalk controls")
    parser.add_argument(
        "--bleed-default", type=float, default=None,
        help="fallback bleed_green_to_red for sessions not in --controls (illustrative use only)",
    )
    parser.add_argument("--out", required=True, type=Path, help="output combined CSV path")
    parser.add_argument("--segmentation-method", default="otsu",
                        choices=["otsu", "sparse", "transmitted", "cellpose"],
                        help="'sparse' for fluorescence frames where cells are a few percent of "
                             "the pixels; 'transmitted' for bright-field, where cells are dark "
                             "rims with bright halos rather than bright objects")
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
    parser.add_argument('--focus-mode', choices=['report', 'enforce'], default='report')
    parser.add_argument('--focus-min-score', type=float)
    parser.add_argument('--focus-min-cnr', type=float, default=1.0)
    parser.add_argument('--focus-min-agreement', type=float, default=0.6)
    parser.add_argument("--segmentation-source", default="sum",
                        choices=["sum", "green", "red", "brightfield"],
                        help="which image the cell masks are cut from; 'sum' keeps the red/green "
                             "ratio from being tilted by the channel that chose the cells. "
                             "'brightfield' is the only genuinely independent choice: it needs a "
                             "bright_field_path per row, is normally paired with "
                             "--segmentation-method transmitted, and uses the optional "
                             "brightfield_shift_dy/dx columns to put the masks on the "
                             "fluorescence frames (scripts/align_brightfield.py measures both "
                             "those offsets and whether the frame shows the same field at all)")
    parser.add_argument('--min-foreground-fraction', type=float, default=0.001,
                        help="refuse a frame whose segmentation claims less of it than this")
    parser.add_argument('--max-foreground-fraction', type=float, default=0.2,
                        help="refuse a frame whose segmentation claims more of it than this")
    args = parser.parse_args()
    from ichnos_image.focus import FocusPolicy
    focus_policy = FocusPolicy(args.focus_mode, args.focus_min_score,
                               args.focus_min_cnr, args.focus_min_agreement)

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
        image_sets, args.out, bleed_green_to_red=bleed_by_session, focus_policy=focus_policy,
        input_paths=[args.manifest] + ([args.controls] + [args.controls.parent / p for p in pd.read_csv(args.controls)[["control_green_path", "control_red_path"]].to_numpy().ravel()] if args.controls else []),
        segmentation_method=args.segmentation_method, segmentation_kwargs=segmentation_kwargs or None,
        segmentation_source=args.segmentation_source,
        background_method=args.background_method, rolling_ball_radius=args.rolling_ball_radius,
        min_foreground_fraction=args.min_foreground_fraction,
        max_foreground_fraction=args.max_foreground_fraction,
    )

    # Absent when process_experiment is stubbed out (tests), so treat a
    # missing manifest as "nothing refused" rather than failing the run
    # after it has already produced its CSV.
    run_manifest = Path(str(out_path) + ".manifest.json")
    refused = (
        json.loads(run_manifest.read_text()).get("refused_image_sets", [])
        if run_manifest.exists() else []
    )
    if refused:
        print(f"\n{len(refused)} image set(s) refused by segmentation QC and left out of the CSV:")
        for item in refused:
            print(f"  - {item['sample_id']}: {item['reason']}")

    df = pd.read_csv(out_path)
    print(f"\n{len(image_sets) - len(refused)}/{len(image_sets)} image set(s) -> {len(df)} cell records -> {out_path}")
    print(df.groupby("session_id").agg(n_cells=("cell_id", "count"), qc_pass_frac=("qc_pass", "mean")))


if __name__ == "__main__":
    main()
