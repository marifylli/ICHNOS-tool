"""The model's predicted induction against the measured one (fold over the same day's 0 µM).

    python3 scripts/compare_model.py --cells outputs/full_20261012_cellpose/cells_all_per_s.csv \\
        --out-dir outputs/model_compare_20261013 --clearance none 1 3

Writes into --out-dir (which must not exist yet):
  data_fold.csv    reporter-positive cells: fold over 0 µM per day, dose and channel, 95 % interval,
                   fold MDE, q (from the within-day test)
  model_fold.csv   model fold at the same doses and imaging times, per clearance scenario
  compare.csv      the two side by side; model_inside_ci, model_over_data
  compare.png      fold against dose: days as points with intervals, model scenarios as lines
  run.json         settings and input fingerprint

The model is the default (uncalibrated) profile; clearance rates are assumed,
not fitted. See ichnos_image/model_compare.py for the scale both sides share.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ichnos_image import model_compare as mc
from ichnos_image import timecourse as tc

DAY_STYLE = [("#2a78d6", "o"), ("#eb6834", "s"), ("#1baf7a", "^")]
SCENARIO_STYLE = ["-", "--", ":", "-."]


def clearance(text: str):
    if text.lower() == "none":
        return None
    value = float(text)
    if not value > 0:
        raise argparse.ArgumentTypeError("clearance must be > 0 per hour, or 'none'")
    return value


def plot(table: pd.DataFrame, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker  # noqa: F401  (minor tick labels off)

    stressors = sorted(table["stressor"].unique())
    fig, axes = plt.subplots(len(stressors), 2, figsize=(10, 3.6 * len(stressors)), squeeze=False)
    for row, stressor in enumerate(stressors):
        part = table[table["stressor"] == stressor]
        for col, readout in enumerate(mc.READOUTS):
            ax = axes[row, col]
            sub = part[part["readout"] == readout]
            data = sub.drop_duplicates(["session_id", "dose_uM"])
            for k, (day, d) in enumerate(data.groupby("day")):
                colour, marker = DAY_STYLE[k % len(DAY_STYLE)]
                d = d.sort_values("dose_uM")
                ax.errorbar(d["dose_uM"] * (1 + 0.04 * (k - 1)), d["fold"],
                            yerr=[d["fold"] - d["fold_ci_low"], d["fold_ci_high"] - d["fold"]],
                            fmt=marker, color=colour, capsize=3, label=f"data {day}")
            for j, (scenario, m) in enumerate(sub.groupby(sub["clearance_per_h"].fillna(0))):
                m = m.groupby("dose_uM")["model_fold"].mean().sort_index()
                name = "model, constant stress" if scenario == 0 else f"model, clearance {scenario:g}/h"
                ax.plot(m.index, m.values, SCENARIO_STYLE[j % 4], color="0.25", label=name)
            ax.axhline(1, color="0.7", lw=0.8)
            ax.set_xscale("log")
            doses = sorted(sub["dose_uM"].unique())
            ax.set_xticks(doses, [f"{d:g}" for d in doses])
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            ax.set_xlabel(f"{stressor.upper()} dose (µM)")
            ax.set_ylabel(f"{readout} fold over same-day 0 µM")
            ax.set_title(f"{stressor.upper()}, {readout}")
            ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--clearance", type=clearance, nargs="+", default=[None],
                    help="first-order stress clearance per hour; 'none' = constant stress (repeatable list)")
    ap.add_argument("--profile", default="default")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--positive-quantile", type=float, default=0.9)
    args = ap.parse_args(argv)

    if args.out_dir.exists():
        raise SystemExit(f"output directory already exists: {args.out_dir}")
    cells = pd.read_csv(args.cells, low_memory=False)
    cells = tc.annotate(cells)
    cells = tc.subtract_baseline(cells, tc.autofluorescence_baseline(cells))

    data = mc.data_fold(cells, n_boot=args.n_boot, seed=args.seed, quantile=args.positive_quantile)
    times = mc.imaging_times(cells)
    stressor_of = cells.drop_duplicates("session_id").set_index("session_id")["stressor"].to_dict()
    model = mc.model_fold(times, stressor_of, clearances=args.clearance, profile=args.profile)
    table = mc.compare(data, model)

    args.out_dir.mkdir(parents=True)
    data.to_csv(args.out_dir / "data_fold.csv", index=False)
    model.to_csv(args.out_dir / "model_fold.csv", index=False)
    table.to_csv(args.out_dir / "compare.csv", index=False)
    plot(table, args.out_dir / "compare.png")
    (args.out_dir / "run.json").write_text(json.dumps({
        "cells": str(args.cells), "cells_sha256": hashlib.sha256(args.cells.read_bytes()).hexdigest(),
        "profile": args.profile, "clearance_per_h": [c for c in args.clearance],
        "n_boot": args.n_boot, "seed": args.seed, "positive_quantile": args.positive_quantile,
        "model_calibrated": False, "clearance_is_assumed": True}, indent=2))

    pd.set_option("display.width", 200)
    shown = table.assign(clearance=table["clearance_per_h"].fillna(0))
    wide = shown.pivot_table(index=["stressor", "day", "readout", "dose_uM", "fold", "fold_ci_low", "fold_ci_high"],
                             columns="clearance", values="model_fold").reset_index()
    wide.columns = [c if not isinstance(c, float) else ("model const" if c == 0 else f"model k={c:g}/h")
                    for c in wide.columns]
    print(wide.round(2).to_string(index=False))
    print(f"\nmodel inside the measured 95 % interval: {int(table['model_inside_ci'].sum())} of {len(table)} rows")
    print(f"wrote {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
