"""Control-normalized H2O2 session data against model predictions.

Reads conditions.csv from analyze_h2o2_session.py; does not re-segment.

Each stressed point is divided by the 0-dose control interpolated to the
same elapsed time. The session ratio scale cancels, so no calibration is
needed. The model is normalized the same way: ratio(d, t) / ratio(0, t),
and likewise for green.

Not a fit: model curves use the stated profile and protocol unchanged.
Optional clearance rates are overlaid as alternative exposure scenarios.

Rows: ratio, green/exposure and red/exposure, each relative to control.
Red is recomputed from cells.csv (usable cells, red_signal / exposure).

Caveats written to report_compare.json:
  - control uncertainty is not propagated into the relative error bars;
  - points outside the control time range use the nearest control value;
  - dose labels are assumed to be micromolar unless stated otherwise.

Usage
  python compare_to_model.py --analysis-dir h2o2_analysis
  python compare_to_model.py --analysis-dir h2o2_analysis --clearance-rates 0.5 2
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ichnos.snapshot_decoder import build_observable_table


def control_series(data, column):
    control = data[(data.dose_label == "0") & data[column].notna()]
    grouped = control.groupby("time_label").agg(t=("elapsed_h_at_red", "mean"), v=(column, "mean"))
    grouped = grouped.sort_values("t")
    return grouped.t.to_numpy(), grouped.v.to_numpy()


def relative(data, column, low, high):
    t_c, v_c = control_series(data, column)
    if len(t_c) < 2:
        raise SystemExit(f"need at least two control points with {column}")
    out = data[(data.dose_label != "0") & data[column].notna()].copy()
    ref = np.interp(out.elapsed_h_at_red, t_c, v_c)
    out[f"{column}_rel"] = out[column] / ref
    out[f"{column}_rel_ci_low"] = out[low] / ref
    out[f"{column}_rel_ci_high"] = out[high] / ref
    out[f"{column}_control_extrapolated"] = (out.elapsed_h_at_red < t_c.min()) | (out.elapsed_h_at_red > t_c.max())
    return out


def red_summary(cells_path, rng, n_boot=2000):
    """Median per-cell red/exposure per image, with a cell bootstrap CI."""
    cells = pd.read_csv(cells_path)
    cells = cells[cells.usable & cells.exposure_acq_ms.notna()]
    rows = []
    for image, group in cells.groupby("image"):
        values = (group.red_signal / (group.exposure_acq_ms / 1000)).to_numpy()
        low = high = np.nan
        if len(values) >= 5:
            meds = np.median(rng.choice(values, (n_boot, len(values))), axis=1)
            low, high = np.percentile(meds, [2.5, 97.5])
        rows.append(dict(image=image, red_per_s_median=np.median(values),
                         red_per_s_ci_low=low, red_per_s_ci_high=high))
    return pd.DataFrame(rows)


def model_relative(doses, times, *, variant, profile, initialization, clearance):
    table = build_observable_table(
        variant=variant, doses_uM=[0.0] + doses, times_hours=times,
        initialization=initialization, profile_name=profile,
        clearance_rate_per_hour=clearance,
    )
    ratio, green = table.ratio_table.ratios, table.green
    ratio_rel, green_rel = ratio[1:] / ratio[0], green[1:] / green[0]
    # Model ratio = f * Reporter_red / Observed_Green, so red_rel = ratio_rel * green_rel.
    return ratio_rel, green_rel, ratio_rel * green_rel


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--analysis-dir", type=Path, required=True)
    ap.add_argument("--out-name", default="model_comparison")
    ap.add_argument("--variant", default="ox")
    ap.add_argument("--profile", default="default")
    ap.add_argument("--initialization", default="equilibrium")
    ap.add_argument("--dose-units-assumption", default="uM",
                    help="recorded only; the notes do not state units (open question 2)")
    ap.add_argument("--clearance-rates", type=float, nargs="*", default=[],
                    help="extra model scenarios, first-order clearance in 1/h")
    ap.add_argument("--seed", type=int, default=20261006)
    ap.add_argument("--exclude", nargs="*", default=["300:3h:01a"],
                    help="dose:time:field entries left out of trends (default: cluster field, notes [22])")
    args = ap.parse_args()
    out = args.analysis_dir / args.out_name
    out.mkdir(exist_ok=False)

    data = pd.read_csv(args.analysis_dir / "conditions.csv", dtype={"dose_label": str, "time_label": str})
    red = red_summary(args.analysis_dir / "cells.csv", np.random.default_rng(args.seed))
    data = data.merge(red, on="image", how="left", validate="one_to_one")
    data = data[data.dose_label != "medium"]
    key = data.dose_label + ":" + data.time_label + ":" + data.field
    excluded = data[key.isin(args.exclude)]
    data = data[~key.isin(args.exclude)].dropna(subset=["ratio_median"])

    ratio_rel = relative(data, "ratio_median", "ratio_ci_low", "ratio_ci_high")
    green_rel = relative(data, "green_per_s_median", "green_per_s_ci_low", "green_per_s_ci_high")
    red_rel = relative(data, "red_per_s_median", "red_per_s_ci_low", "red_per_s_ci_high")
    ratio_rel.to_csv(out / "ratio_relative.csv", index=False)
    green_rel.to_csv(out / "green_relative.csv", index=False)
    red_rel.to_csv(out / "red_relative.csv", index=False)

    doses = sorted({float(d) for d in data.dose_label if d != "0"})
    times = list(np.round(np.arange(0.25, data.elapsed_h_at_red.max() + 0.3, 0.05), 3))
    scenarios = {"constant exposure": None, **{f"clearance {k:g}/h": k for k in args.clearance_rates}}
    models = {}
    for name, rate in scenarios.items():
        print(f"simulating model: {name}")
        models[name] = model_relative(doses, times, variant=args.variant, profile=args.profile,
                                      initialization=args.initialization, clearance=rate)
    rows = []
    for name, (r, g, rd) in models.items():
        for i, dose in enumerate(doses):
            for j, t in enumerate(times):
                rows.append(dict(scenario=name, dose=dose, time_h=t, ratio_rel=r[i, j],
                                 green_rel=g[i, j], red_rel=rd[i, j]))
    pd.DataFrame(rows).to_csv(out / "model_relative.csv", index=False)

    fig, axes = plt.subplots(3, len(doses), figsize=(3.4 * len(doses), 9), sharex=True, sharey="row")
    styles = ["-", "--", ":", "-."]
    for i, dose in enumerate(doses):
        label = f"{dose:g}"
        for row, (frame, col, title) in enumerate([
            (ratio_rel, "ratio_median", "red/green vs control"),
            (green_rel, "green_per_s_median", "green/exposure vs control"),
            (red_rel, "red_per_s_median", "red/exposure vs control"),
        ]):
            ax = axes[row, i]
            for k, (name, curves) in enumerate(models.items()):
                ax.plot(times, curves[row][i], "k" + styles[k % 4], lw=1.2,
                        label=f"model, {name}" if i == 0 else None)
            pts = frame[frame.dose_label == label].sort_values("elapsed_h_at_red")
            extra = pts[f"{col}_control_extrapolated"]
            y = pts[f"{col}_rel"]
            yerr = [y - pts[f"{col}_rel_ci_low"], pts[f"{col}_rel_ci_high"] - y]
            ax.errorbar(pts.elapsed_h_at_red, y, yerr=yerr, fmt="o", color="C3", capsize=3,
                        label="data (one field)" if i == 0 else None)
            ax.scatter(pts.elapsed_h_at_red[extra], y[extra], facecolor="white", edgecolor="C3", zorder=3,
                       label="control extrapolated" if i == 0 else None)
            ax.axhline(1, color="grey", lw=0.8)
            if row == 0:
                ax.set_title(f"{label} ({args.dose_units_assumption} assumed)")
            if i == 0:
                ax.set_ylabel(title)
            if row == 2:
                ax.set_xlabel("hours from H2O2 addition to image")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle(f"Pilot, one field per condition; model {args.variant}/{args.profile}, not fitted; "
                 "error bars = cell bootstrap of the target only", fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "relative_vs_model.png", dpi=150)

    early = ratio_rel[ratio_rel.elapsed_h_at_red <= 2.6]
    report = dict(
        scope="pilot descriptive comparison; not a fit, not validation",
        model=dict(variant=args.variant, profile=args.profile, initialization=args.initialization,
                   scenarios={k: v for k, v in scenarios.items()}),
        dose_units_assumption=args.dose_units_assumption, excluded=excluded[["dose_label", "time_label", "field"]]
            .to_dict("records"),
        early_window_hours=2.6,
        early_points_below_control=int((early.ratio_median_rel < 1).sum()), early_points=len(early),
        caveats=["single control field per time; all comparisons share it",
                 "control uncertainty not propagated",
                 "elapsed time = addition to image; exposure history unresolved (Q1)",
                 "red background did not scale with recorded exposure; red/exposure is indicative only",
                 "model red_rel derived as ratio_rel * green_rel"],
    )
    (out / "report_compare.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
