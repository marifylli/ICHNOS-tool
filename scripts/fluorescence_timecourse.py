"""Reporter time courses per stressor, dose and day, and when the fluorescence peaks.

    python3 scripts/fluorescence_timecourse.py --cells outputs/full_20261012_cellpose/cells_per_s.csv \\
        --out-dir outputs/timecourse_20261013

Input: the per-cell CSV after exposure normalisation (``scripts/normalize_exposure.py``),
i.e. with ``corrected_mean_green_per_s`` and ``corrected_mean_red_per_s``.

Writes into --out-dir (which must not exist yet):
  baseline.csv              autofluorescence per stressor and channel, from non-transfected cells
  baseline_diagnostics.csv  non-transfected medians per dose and time, to judge the constant baseline
  cells_af.csv              the input with green_af / red_af (baseline subtracted)
  fields.csv                one row per field: median green_af, red_af over usable cells
  timecourse.csv            one row per sample and time point (fields pooled)
  peaks.csv                 per transfected curve, green and red: peak time, 5-95 % interval,
                            whether it lies between samples, and the excess over the same
                            day's 0 µM curve at that time
  timecourse_<stressor>.png curves per dose, one line per day, peak marked
  run.json                  settings and input fingerprint

See ichnos_image/timecourse.py for the choices (constant baseline per stressor,
field medians, parabola refinement, resampling interval). Time is
``measurement_time_hours``: hours from stress addition to imaging.
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

from ichnos_image import timecourse as tc

# Validated categorical slots 1-3 (light surface); days keep their slot in every figure.
DAY_COLOURS = ("#2a78d6", "#eb6834", "#1baf7a")
DAY_MARKERS = ("o", "s", "^")
CHANNEL_LABELS = {"green": "green (sfGFP), grey/s above baseline",
                  "red": "red (mCherry), grey/s above baseline",
                  "ratio": "red / green"}


def fingerprint(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def plot(stressor: str, curves: pd.DataFrame, fields: pd.DataFrame, peaks: pd.DataFrame, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    doses = sorted(curves["dose_uM"].unique())
    days = sorted(curves["day"].unique())
    channels = ["green", "red", "ratio"]
    fig, axes = plt.subplots(len(channels), len(doses), figsize=(3.0 * len(doses), 2.6 * len(channels)),
                             sharex=True, squeeze=False)
    for col, dose in enumerate(doses):
        for row, channel in enumerate(channels):
            ax = axes[row][col]
            for k, day in enumerate(days):
                colour, marker = DAY_COLOURS[k % len(DAY_COLOURS)], DAY_MARKERS[k % len(DAY_MARKERS)]
                curve = curves[(curves["dose_uM"] == dose) & (curves["day"] == day)].sort_values("time_h")
                if curve.empty:
                    continue
                if channel != "ratio":
                    pts = fields[(fields["dose_uM"] == dose) & (fields["day"] == day)]
                    ax.scatter(pts["time_h"], pts[f"{channel}_af"], s=10, color=colour, alpha=0.35,
                               marker=marker, linewidths=0)
                ax.plot(curve["time_h"], curve[channel], color=colour, marker=marker, markersize=5,
                        linewidth=2, label=f"{day[6:8]}/{day[4:6]}")
                peak = peaks[(peaks["dose_uM"] == dose) & (peaks["day"] == day) & (peaks["channel"] == channel)]
                if not peak.empty and np.isfinite(peak["t_peak"].iloc[0]):
                    p = peak.iloc[0]
                    ax.axvline(p["t_peak"], color=colour, linewidth=1, linestyle=":" if not p["bracketed"] else "-")
            ax.axhline(0, color="#b0afaa", linewidth=0.8)
            ax.grid(color="#e6e5e1", linewidth=0.6)
            ax.tick_params(labelsize=8, colors="#52514e")
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            if row == 0:
                ax.set_title(f"{stressor.upper()} {dose:g} µM", fontsize=10)
            if col == 0:
                ax.set_ylabel(CHANNEL_LABELS[channel], fontsize=8, color="#52514e")
            if row == len(channels) - 1:
                ax.set_xlabel("hours since stress added (at imaging)", fontsize=8, color="#52514e")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, title="day", loc="upper right", fontsize=8, title_fontsize=8, frameon=False)
    fig.text(0.01, 0.005, "Line: median of field medians per time point. Dots: single fields. "
             "Vertical line: peak (solid = between samples; dotted = at first/last sample, only a bound).",
             fontsize=7, color="#52514e")
    fig.tight_layout(rect=(0, 0.02, 0.95, 1))
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--ratio-floor", type=float, default=0.0,
                    help="leave the ratio empty where either baseline-subtracted channel is at or below this")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)

    if args.out_dir.exists():
        raise SystemExit(f"output directory already exists: {args.out_dir}")
    cells = pd.read_csv(args.cells, low_memory=False)
    needed = [f"corrected_mean_{ch}_per_s" for ch in tc.CHANNELS]
    if missing := [c for c in needed if c not in cells.columns]:
        raise SystemExit(f"{args.cells} lacks {missing}; run scripts/normalize_exposure.py first")

    cells = tc.annotate(cells)
    baseline = tc.autofluorescence_baseline(cells)
    diagnostics = tc.baseline_diagnostics(cells)
    cells = tc.subtract_baseline(cells, baseline)
    fields = tc.field_summaries(cells)
    curves = tc.timecourse(fields, ratio_floor=args.ratio_floor)

    peak_rows = []
    transfected = curves[curves["strain"] == "tr"]
    for keys, curve in transfected.groupby(tc.CURVE_KEYS, dropna=False):
        table = tc.bootstrap_peaks(cells, curve, fields, n_boot=args.n_boot, seed=args.seed)
        for name, value in reversed(list(zip(tc.CURVE_KEYS, keys))):
            table.insert(0, name, value)
        peak_rows.append(table)
    peaks = pd.concat(peak_rows, ignore_index=True)
    peaks["excess_over_0uM"] = tc.excess_over_control(peaks, transfected)
    peaks = peaks[tc.CURVE_KEYS + [
        "channel", "t_peak", "t_peak_p05", "t_peak_p95", "t_max_sampled", "position", "bracketed",
        "share_at_sampled_max", "rise_from_first", "excess_over_0uM", "n_timepoints", "t_first", "t_last"]]

    args.out_dir.mkdir(parents=True)
    baseline.to_csv(args.out_dir / "baseline.csv", index=False)
    diagnostics.to_csv(args.out_dir / "baseline_diagnostics.csv", index=False)
    cells.to_csv(args.out_dir / "cells_af.csv", index=False)
    fields.to_csv(args.out_dir / "fields.csv", index=False)
    curves.to_csv(args.out_dir / "timecourse.csv", index=False)
    peaks.to_csv(args.out_dir / "peaks.csv", index=False)
    if not args.no_figures:
        for stressor in sorted(transfected["stressor"].unique()):
            plot(stressor, transfected[transfected["stressor"] == stressor],
                 fields[(fields["stressor"] == stressor) & (fields["strain"] == "tr")],
                 peaks[peaks["stressor"] == stressor], args.out_dir / f"timecourse_{stressor}.png")
    (args.out_dir / "run.json").write_text(json.dumps({
        "cells": str(args.cells), "cells_sha256": fingerprint(args.cells), "n_boot": args.n_boot,
        "seed": args.seed, "ratio_floor": args.ratio_floor, "time_axis": "measurement_time_hours",
        "baseline": "median of usable non-transfected cells per stressor and channel",
    }, indent=2))

    print(baseline[["stressor", "n_cells", "n_fields", "days", "baseline_green", "baseline_red"]]
          .round(2).to_string(index=False))
    shown = peaks[peaks["channel"] == "green"][["stressor", "day", "dose_uM", "t_peak", "t_peak_p05",
                                                "t_peak_p95", "position", "excess_over_0uM", "t_first", "t_last"]]
    print("\ngreen peak (h since stress added):")
    print(shown.round(2).to_string(index=False))
    print(f"\nwrote {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
