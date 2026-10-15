"""Decode the dose of single fields, validated leave-one-day-out.

    python3 scripts/evaluate_decoder_lodo.py --cells data/reference/cells_all_per_s.csv.gz \\
        --out-dir outputs/decoder_lodo_20261013

Writes into --out-dir (which must not exist yet):
  fields.csv        features of every transfected field (see ichnos_image/lodo_decoder.py)
  predictions.csv   posterior over the doses for every field and pooled sample of every held-out day
  folds.csv         which day was held out, how many training days, whether between-day spread was estimable
  metrics.csv       accuracy vs chance, mean p(true dose), information gain (bits), entropy vs prior,
                    80 % credible-set coverage and size; per held-out day and the mean over days
  posterior.png     mean posterior (single fields) per true dose: a flat row means the decoder learnt nothing
  run.json          settings and input fingerprint

Decoders: absolute (one image, nothing else), paired (minus the same day's
0 µM), model (paired log folds, class means from the ICHNOS model, default
profile, uncalibrated; --clearance sets the assumed H2O2/DTT clearance).
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

from ichnos_image import lodo_decoder as ld
from ichnos_image import model_compare as mc
from ichnos_image import timecourse as tc

# sequential blue ramp, steps 100-700 of the reference palette (light -> dark)
RAMP = ["#ffffff", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#1c5cab", "#104281", "#0d366b"]
GRID = np.round(np.arange(0.5, 8.01, 0.25), 2)


def model_log_fold_table(fields: pd.DataFrame, clearance: float | None, profile: str) -> pd.DataFrame:
    """Model log fold over 0 µM for every stressor and dose of the data, on a time grid."""
    rows = []
    for stressor, part in fields.groupby("stressor"):
        variant = mc.STRESSOR_VARIANT.get(stressor)
        if variant is None:
            continue
        control = mc._run_model(variant, 0.0, list(GRID), clearance, profile)
        for dose in sorted(part["dose_uM"].unique()):
            run = control if dose == 0 else mc._run_model(variant, float(dose), list(GRID), clearance, profile)
            for i, t in enumerate(GRID):
                rows.append({"stressor": stressor, "dose_uM": float(dose), "time_h": float(t),
                             **{f"log_fold_{ch}": float(np.log(run[ch][i] / control[ch][i])) for ch in mc.READOUTS}})
    return pd.DataFrame(rows)


def interpolator(table: pd.DataFrame):
    def model_log_fold(stressor: str, dose: float, time: float) -> np.ndarray:
        curve = table[(table["stressor"] == stressor) & (table["dose_uM"] == dose)].sort_values("time_h")
        return np.array([np.interp(time, curve["time_h"], curve[f"log_fold_{ch}"]) for ch in mc.READOUTS])
    return model_log_fold


def plot(predictions: pd.DataFrame, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("seq_blue", RAMP)
    fields = predictions[predictions["unit"] == "field"]
    panels = list(fields.groupby(["stressor", "decoder"], sort=True))
    stressors = sorted(fields["stressor"].unique())
    decoders = [d for d in ld.FEATURES if d in set(fields["decoder"])]
    fig, axes = plt.subplots(len(stressors), len(decoders), figsize=(3.6 * len(decoders), 3.4 * len(stressors)),
                             squeeze=False)
    for ax in axes.ravel():
        ax.set_visible(False)
    for (stressor, decoder), part in panels:
        ax = axes[stressors.index(stressor), decoders.index(decoder)]
        ax.set_visible(True)
        cols = [c for c in part.columns if c.startswith("p_") and c.endswith("uM") and part[c].notna().any()]
        doses = [float(c[2:-2]) for c in cols]
        order = np.argsort(doses)
        cols, doses = [cols[i] for i in order], [doses[i] for i in order]
        matrix = np.array([part[part["true_dose_uM"] == d][cols].mean().to_numpy() for d in doses])
        ax.imshow(matrix, cmap=cmap, vmin=0, vmax=1)
        for r in range(len(doses)):
            for c in range(len(doses)):
                v = matrix[r, c]
                ax.text(c, r, f"{v:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if v > 0.55 else "#1f1f1f")
            ax.add_patch(plt.Rectangle((r - 0.5, r - 0.5), 1, 1, fill=False, ec="#1f1f1f", lw=1.2))
        labels = [f"{d:g}" for d in doses]
        ax.set_xticks(range(len(doses)), labels, fontsize=8)
        ax.set_yticks(range(len(doses)), labels, fontsize=8)
        ax.set_xlabel("decoded dose (µM)", fontsize=8)
        ax.set_ylabel("true dose (µM)", fontsize=8)
        acc = (part["map_dose_uM"] == part["true_dose_uM"]).mean()
        ax.set_title(f"{stressor.upper()} · {decoder}\nall held-out fields: accuracy {acc:.2f} "
                     f"(chance {1 / len(doses):.2f})", fontsize=9)
    fig.suptitle("Mean posterior per true dose, single fields, leave-one-day-out "
                 "(flat rows = no information; prior = 1/K)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def clearance(text: str):
    return None if text.lower() == "none" else float(text)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--decoders", nargs="+", default=list(ld.FEATURES), choices=list(ld.FEATURES))
    ap.add_argument("--clearance", type=clearance, default=None, help="model decoder: clearance per hour, or none")
    ap.add_argument("--profile", default="default")
    ap.add_argument("--min-cells", type=int, default=20)
    ap.add_argument("--positive-quantile", type=float, default=0.9)
    args = ap.parse_args(argv)

    if args.out_dir.exists():
        raise SystemExit(f"output directory already exists: {args.out_dir}")
    cells = tc.annotate(pd.read_csv(args.cells, low_memory=False))
    cells = tc.subtract_baseline(cells, tc.autofluorescence_baseline(cells))
    fields = ld.add_paired(ld.field_features(cells, quantile=args.positive_quantile, min_cells=args.min_cells))

    model_fn = None
    if "model" in args.decoders:
        model_fn = interpolator(model_log_fold_table(fields, args.clearance, args.profile))
    predictions, folds = [], []
    for decoder in args.decoders:
        p, f = ld.leave_one_day_out(fields, decoder, model_fn)
        predictions.append(p)
        folds.append(f)
    predictions = pd.concat(predictions, ignore_index=True)
    folds = pd.concat(folds, ignore_index=True)
    metrics = ld.summarize(predictions)

    args.out_dir.mkdir(parents=True)
    fields.to_csv(args.out_dir / "fields.csv", index=False)
    predictions.to_csv(args.out_dir / "predictions.csv", index=False)
    folds.to_csv(args.out_dir / "folds.csv", index=False)
    metrics.to_csv(args.out_dir / "metrics.csv", index=False)
    plot(predictions, args.out_dir / "posterior.png")
    (args.out_dir / "run.json").write_text(json.dumps({
        "cells": str(args.cells), "cells_sha256": hashlib.sha256(args.cells.read_bytes()).hexdigest(),
        "decoders": args.decoders, "clearance_per_h": args.clearance, "profile": args.profile,
        "min_cells": args.min_cells, "positive_quantile": args.positive_quantile,
        "prior": "uniform over the doses of the experiment", "credible_mass": ld.CREDIBLE,
        "validation": "leave-one-day-out per stressor; non-transfected sessions only in the baseline"}, indent=2))

    pd.set_option("display.width", 220)
    shown = metrics[metrics["held_out"] == "mean over days"]
    print("mean over held-out days:")
    print(shown[["stressor", "decoder", "unit", "n", "accuracy", "chance", "mean_p_true", "info_gain_bits",
                 "entropy_vs_prior", "credible80_coverage", "credible80_size", "outside_training"]].round(2).to_string(index=False))
    single = folds[folds["n_training_days"] == 1]
    if len(single):
        print("\nnote: one usable training day, so the between-day spread is set to 0 (overconfident): "
              + ", ".join(f"{r.decoder}/{r.held_out}" for r in single.itertuples()))
    empty = folds[folds["n_test_fields"] == 0]
    if len(empty):
        print("note: no scorable fields (too few reporter-positive cells): "
              + ", ".join(f"{r.decoder}/{r.held_out}" for r in empty.itertuples()))
    print(f"\nwrote {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
