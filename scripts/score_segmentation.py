"""Score a segmentation method against the hand-marked reference crops.

    python3 scripts/score_segmentation.py --reference data/segmentation_reference \\
        --method brightfield_rim --out outputs/segscore/rim_default.csv \\
        --overlay-dir outputs/segscore/rim_default

    # a different setting
    python3 scripts/score_segmentation.py ... --param rim_sigmas=2.5 --param opening_um=0.5

    # Cellpose (needs `pip install cellpose`; downloads its model on first use)
    python3 scripts/score_segmentation.py ... --method cellpose --cellpose-gpu

    # the masks the full run of 2026-10-10 used, as the baseline
    python3 scripts/score_segmentation.py ... --method sparse --source sum \\
        --manifest data/acquisition/manifest/images.csv

Reads crops.csv (which frame each crop was cut from, and where) and, per
crop, <crop stem>.points.csv and optionally <crop stem>.unsure.csv. Crops
without a points file are skipped and listed. The whole frame is segmented,
as the pipeline does, and scored inside the crop's window; see
ichnos_image/segmentation_score.py for what counts as what.

Writes one row per crop plus a pooled row (sums of the counts, recall and
precision recomputed from them). With --overlay-dir, also one PNG per crop:
masks with one point in green, merged in orange, false in red, unsure in
grey; missed points as yellow crosses.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from skimage import segmentation as sk_segmentation

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ichnos_image import segment
from ichnos_image.image_io import load_image
from ichnos_image.manifest import _load_bright_field
from ichnos_image.segmentation_score import read_points, score_labels

COLOURS = {"separated": (40, 200, 60), "merged": (255, 150, 0), "false": (230, 30, 30), "unsure": (150, 150, 150)}


def parse_param(text: str) -> tuple[str, object]:
    key, sep, value = text.partition("=")
    if not sep or not key:
        raise argparse.ArgumentTypeError(f"--param needs key=value, got {text!r}")
    try:
        return key, float(value)
    except ValueError:
        return key, value


def source_image(crop: pd.Series, source: str, manifest: pd.DataFrame | None) -> np.ndarray:
    if source == "brightfield":
        return _load_bright_field(crop.bright_field_path)
    if manifest is None:
        raise SystemExit("--source sum needs --manifest, to find the fluorescence pair of each crop")
    rows = manifest[(manifest.session_id == crop.session_id) & (manifest.sample_id == crop.sample_id)
                    & (manifest.bright_field_path == crop.bright_field_path)]
    if len(rows) != 1:
        raise SystemExit(f"{crop.crop}: {len(rows)} manifest rows match; expected 1")
    green, red = load_image(rows.iloc[0].green_path), load_image(rows.iloc[0].red_path)
    if green.ndim != 2 or red.ndim != 2:
        raise SystemExit(f"{crop.crop}: fluorescence frames are not single-plane; extract a plane first")
    return green + red


def overlay(crop_png: Path, labels: np.ndarray, window, detail: dict, out: Path) -> None:
    y0, x0, height, width = window
    base = np.asarray(Image.open(crop_png).convert("L"), dtype=np.uint8)
    rgb = np.dstack([base] * 3)
    crop_labels = labels[y0:y0 + height, x0:x0 + width]
    edges = sk_segmentation.find_boundaries(crop_labels, mode="inner")
    for mask_id, status in detail["masks"].items():
        rgb[edges & (crop_labels == mask_id)] = COLOURS[status]
    for row, col in detail["missed"]:
        r, c = int(row), int(col)
        for d in range(-5, 6):
            for rr, cc in ((r + d, c), (r, c + d)):
                if 0 <= rr < height and 0 <= cc < width:
                    rgb[rr, cc] = (255, 230, 0)
    Image.fromarray(rgb).save(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", type=Path, required=True, help="folder with crops.csv and the points files")
    ap.add_argument("--method", required=True, choices=["brightfield_rim", "transmitted", "sparse", "otsu", "cellpose"])
    ap.add_argument("--cellpose-gpu", action="store_true", help="--method cellpose: use the GPU (CUDA, or Apple MPS)")
    ap.add_argument("--source", default="brightfield", choices=["brightfield", "sum"])
    ap.add_argument("--manifest", type=Path, help="images.csv; needed for --source sum")
    ap.add_argument("--param", action="append", type=parse_param, default=[],
                    help="setting passed to the method, key=value (repeatable)")
    ap.add_argument("--tolerance-um", type=float, default=0.5,
                    help="how far from a point a mask may be and still count as under it")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--overlay-dir", type=Path)
    args = ap.parse_args(argv)

    if args.out.exists():
        raise SystemExit(f"output path already exists: {args.out}")
    crops = pd.read_csv(args.reference / "crops.csv", dtype=str)
    manifest = pd.read_csv(args.manifest, dtype=str) if args.manifest else None
    params = dict(args.param)
    if args.method == "cellpose":
        # Load once: the model is ~1.1 GB and would otherwise be read per crop.
        params["cellpose_model"] = segment.load_cellpose_model(gpu=args.cellpose_gpu)

    rows, skipped = [], []
    for _, crop in crops.iterrows():
        stem = Path(crop.crop).stem
        points_path = args.reference / f"{stem}.points.csv"
        if not points_path.exists():
            skipped.append(stem)
            continue
        unsure_path = args.reference / f"{stem}.unsure.csv"
        pixel_size = float(crop.pixel_size_um)
        kwargs = dict(params)
        if args.method == "brightfield_rim":
            kwargs.setdefault("pixel_size_um", pixel_size)
        labels = segment.segment_cells(source_image(crop, args.source, manifest), method=args.method, **kwargs)
        window = tuple(int(crop[k]) for k in ("y0", "x0", "height", "width"))
        score, detail = score_labels(
            labels, read_points(points_path),
            read_points(unsure_path) if unsure_path.exists() else None,
            tolerance_px=args.tolerance_um / pixel_size, window=window,
        )
        rows.append({"crop": stem, "session_id": crop.session_id, **score.as_dict()})
        print(f"{stem:20s} marked {score.marked:4d}  recall {score.recall:.2f}  precision {score.precision:.2f}  "
              f"merged {score.merged_masks:3d}  false {score.false_masks:3d}")
        if args.overlay_dir:
            args.overlay_dir.mkdir(parents=True, exist_ok=True)
            overlay(args.reference / crop.crop, labels, window, detail, args.overlay_dir / f"{stem}.png")

    if not rows:
        raise SystemExit(f"no points files in {args.reference}; skipped: {skipped}")
    table = pd.DataFrame(rows)
    total = table[["marked", "unsure", "masks_scored", "detected", "separated", "merged_masks", "false_masks"]].sum()
    # A mask holding exactly one point and a point alone in its mask are the
    # same pairs, so "separated" counts both.
    recall = total.separated / total.marked
    precision = total.separated / total.masks_scored if total.masks_scored else float("nan")
    pooled = {"crop": "ALL", "session_id": "", **total.to_dict(), "recall": recall, "precision": precision,
              "f1": 2 * recall * precision / (recall + precision) if recall + precision else 0.0}
    table = pd.concat([table, pd.DataFrame([pooled])], ignore_index=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", newline="") as handle:
        table.to_csv(handle, index=False)
    Path(str(args.out) + ".json").write_text(json.dumps(
        {"method": args.method, "source": args.source,
         "params": {k: v for k, v in params.items() if k != "cellpose_model"},
         "tolerance_um": args.tolerance_um, "skipped_crops": skipped}, indent=2, ensure_ascii=False))
    print(f"ALL{'':17s} marked {int(total.marked):4d}  recall {recall:.2f}  precision {precision:.2f}"
          + (f"   (skipped, no points: {', '.join(skipped)})" if skipped else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
