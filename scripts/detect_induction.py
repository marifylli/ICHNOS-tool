"""Does the reporter respond to the stress? Each dose against the same day's 0 µM.

    python3 scripts/detect_induction.py --cells outputs/full_20261012_cellpose/cells_all_per_s.csv \\
        --out-dir outputs/induction_20261013

Writes into --out-dir (which must not exist yet):
  effects.csv   per day, dose and readout (green, red, positive share): mean excess over the
                same day's 0 µM, 95 % interval, minimum detectable effect, bootstrap p and
                Benjamini-Hochberg q over all rows
  trend.csv     per day and readout: slope of the effect against log10(dose), 95 % interval, p, q
  run.json      settings and input fingerprint

See ichnos_image/induction.py for why every comparison stays within one day.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ichnos_image import induction
from ichnos_image import timecourse as tc


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--positive-quantile", type=float, default=0.9,
                    help="a cell is reporter-positive above this quantile of non-transfected cells, in both channels")
    args = ap.parse_args(argv)

    if args.out_dir.exists():
        raise SystemExit(f"output directory already exists: {args.out_dir}")
    cells = pd.read_csv(args.cells, low_memory=False)
    needed = [f"corrected_mean_{ch}_per_s" for ch in tc.CHANNELS]
    if missing := [c for c in needed if c not in cells.columns]:
        raise SystemExit(f"{args.cells} lacks {missing}; run scripts/normalize_exposure.py first")

    cells = tc.annotate(cells)
    cells = tc.subtract_baseline(cells, tc.autofluorescence_baseline(cells))
    effects, trend = induction.induction_table(cells, n_boot=args.n_boot, seed=args.seed,
                                               quantile=args.positive_quantile)
    args.out_dir.mkdir(parents=True)
    effects.to_csv(args.out_dir / "effects.csv", index=False)
    trend.to_csv(args.out_dir / "trend.csv", index=False)
    digest = hashlib.sha256(args.cells.read_bytes()).hexdigest()
    (args.out_dir / "run.json").write_text(json.dumps({
        "cells": str(args.cells), "cells_sha256": digest, "n_boot": args.n_boot, "seed": args.seed,
        "positive_quantile": args.positive_quantile,
        "reference": "same day's 0 µM curve, linearly interpolated to each time point"}, indent=2))

    pd.set_option("display.width", 200)
    shown = effects[effects["detected_bh"]][["day", "stressor", "dose_uM", "readout", "effect",
                                             "ci_low", "ci_high", "mde", "q_bh"]]
    print(f"{len(effects)} comparisons; {int(effects['detected_bh'].sum())} with q < 0.05:")
    print(shown.round(3).to_string(index=False) if len(shown) else "  none")
    print("\ndose trend (slope per 10-fold dose):")
    print(trend[["day", "stressor", "readout", "slope_per_log10_dose", "ci_low", "ci_high", "q_bh"]]
          .round(3).to_string(index=False))
    print(f"\nwrote {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
