"""Build the image-pipeline manifest by matching TIFF timestamps to the acquisition log.

    python3 scripts/build_image_manifest.py \\
        --log data/acquisition/acquisition_log.csv \\
        --root A=~/Downloads/wetransfer_6-10_2026-10-06_1906 \\
        --root "B=~/Downloads/igem h2o2 non transfected" \\
        --root C=~/Downloads/wetransfer_igem-dtt-vol2_2026-10-07_1803 \\
        --root "D=~/Downloads/igem DTT non transfected" \\
        --objective 40X \\
        --out-dir data/acquisition/manifest

Writes three files into a new --out-dir:
  images.csv        manifest for scripts/run_pipeline.py (absolute image paths)
  match_report.csv  every image and every log row, matched or not, with residuals
  offsets.csv       the camera-clock offset estimated for each day
Images stay where they are; nothing is renamed or moved.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ichnos_image.acquisition_match import (
    MANIFEST_COLUMNS,
    estimate_offsets,
    manifest_rows,
    match,
    scan_root,
    write_csv,
)
from ichnos_image.instrument import OBJECTIVES


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--log", type=Path, required=True)
    ap.add_argument("--root", action="append", required=True, metavar="CODE=PATH",
                    help="image folder for a folder code of the log (A, B, C, D); repeat")
    ap.add_argument("--objective", required=True, choices=sorted(OBJECTIVES),
                    help="objective used for fluorescence; not recorded in the log, so it must be stated")
    ap.add_argument("--tolerance-min", type=float, default=2.5,
                    help="max time difference for an image matched to a row of another folder")
    ap.add_argument("--same-folder-tolerance-min", type=float, default=5.0,
                    help="max time difference for an image matched to a row of its own folder")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)

    if args.out_dir.exists():
        raise SystemExit(f"output path already exists: {args.out_dir}")
    with args.log.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    images = []
    for item in args.root:
        code, sep, path = item.partition("=")
        if not sep:
            raise SystemExit(f"--root needs CODE=PATH, got {item!r}")
        found = scan_root(code, Path(path).expanduser())
        print(f"{code}: {len(found)} red/green pairs in {path}")
        images.extend(found)

    offsets = estimate_offsets(images, rows)
    report = match(images, rows, offsets, tolerance_min=args.tolerance_min,
                   same_folder_tolerance_min=args.same_folder_tolerance_min)
    manifest = manifest_rows(images, rows, offsets, objective=args.objective)

    args.out_dir.mkdir(parents=True)
    write_csv(manifest, args.out_dir / "images.csv", MANIFEST_COLUMNS)
    write_csv(report, args.out_dir / "match_report.csv",
              ["status", "root", "image_folder", "log_folder", "log_dose_uM", "log_tp", "log_rep",
               "log_index", "residual_min", "camera_red_time", "red_path", "green_path"])
    write_csv([{"date": d, "offset_min": f"{o:.2f}"} for d, o in sorted(offsets.items())],
              args.out_dir / "offsets.csv", ["date", "offset_min"])

    for day, offset in sorted(offsets.items()):
        print(f"clock offset {day}: {offset:+.1f} min")
    counts = Counter(r["status"] for r in report)
    print(dict(counts))
    for r in report:
        if r["status"] == "matched_other_folder":
            print(f"  saved in {r['image_folder']} but taken at the time of "
                  f"{r['log_dose_uM']} uM {r['log_tp']} {r['log_rep']} (log folder {r['log_folder']})")
    print(f"{len(manifest)} manifest rows -> {args.out_dir / 'images.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
