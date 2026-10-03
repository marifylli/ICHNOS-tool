"""General-purpose CLI: run the full ichnos_image pipeline (segment ->
correct -> extract -> export) over an arbitrary experiment described by a
manifest CSV, and a per-session crosstalk-control manifest, producing one
combined CSV.

Image files are read with PIL (ichnos_image.synthesize.load_png16 despite
its name -- Image.open handles PNG, TIFF, and most other common formats;
raw microscope formats like .czi/.nd2 would need a dedicated reader, e.g.
Bio-Formats via python-bioformats, not included here).

--manifest CSV columns (one row per image set):
    session_id, timepoint, acquisition_order, green_path, red_path,
    bright_field_path (optional, blank if none), exposure_ms_green,
    exposure_ms_red, nd_filter_green, nd_filter_red, objective, burner_hours,
    lamp_warmup_minutes

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

import pandas as pd

from ichnos_image import ImageSet, process_experiment
from ichnos_image.correct import calibrate_crosstalk_from_control
from ichnos_image.synthesize import load_png16


def _build_image_sets(manifest_path: Path) -> list[ImageSet]:
    manifest = pd.read_csv(manifest_path)
    image_sets = []
    for row in manifest.itertuples():
        bright_field = (
            load_png16(row.bright_field_path)
            if getattr(row, "bright_field_path", "") and pd.notna(row.bright_field_path)
            else None
        )
        image_sets.append(
            ImageSet(
                green=load_png16(row.green_path),
                red=load_png16(row.red_path),
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
            )
        )
    return image_sets


def _calibrate_bleed_per_session(controls_path: Path | None, sessions: set[str], bleed_default: float | None) -> dict:
    bleed_by_session = {}
    if controls_path is not None:
        controls = pd.read_csv(controls_path)
        for row in controls.itertuples():
            control_green = load_png16(row.control_green_path)
            control_red = load_png16(row.control_red_path)
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
    args = parser.parse_args()

    image_sets = _build_image_sets(args.manifest)
    sessions = {s.session_id for s in image_sets}
    bleed_by_session = _calibrate_bleed_per_session(args.controls, sessions, args.bleed_default)

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
