"""Step 2b: add exposure-normalised intensities (grey levels per second) to a cells CSV.

    python3 scripts/normalize_exposure.py \\
        --cells outputs/full_20261010/cells.csv \\
        --out outputs/full_20261010/cells_per_s.csv

Adds corrected_mean_green_per_s, corrected_mean_red_per_s,
integrated_green_per_s and integrated_red_per_s. Nothing else changes; the
red/green ratio is already exposure-free because both channels share one
exposure per sample (checked on every row).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ichnos_image.exposure import add_per_second_columns


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    if args.out.exists():
        raise SystemExit(f"output path already exists: {args.out}")
    cells = pd.read_csv(args.cells, low_memory=False)
    out = add_per_second_columns(cells)
    out.to_csv(args.out, index=False)
    summary = out.groupby("session_id").agg(
        exposure_ms=("exposure_ms_green", lambda s: " ".join(str(int(v)) for v in sorted(set(s)))),
        green=("corrected_mean_green", "median"), green_per_s=("corrected_mean_green_per_s", "median"))
    print(summary.round(2).to_string())
    print(f"{len(out)} cells -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
