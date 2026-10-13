"""Decode the dose of a new culture from its images, with the decoder's own track record attached.

A user images a culture exposed to an unknown dose, ideally next to an
unstressed (0 µM) culture of the same day, runs the image pipeline, and asks
which dose it was. This module answers with a posterior over the doses of a
reference data set (``lodo_decoder``), trained on every transfected day of
the stressor in that set, and with a verdict in words:

  dose estimate         the posterior is informative and the decoder beat
                        chance on held-out reference days
  not validated         the posterior points somewhere, but on held-out
                        reference days this decoder did no better than
                        chance; the dose is not supported
  no dose information   the posterior is close to the prior (entropy at
                        least 90 % of it): the image cannot tell the doses apart
  unlike reference      the field is farther from every reference dose than
                        the reference data would ever put it (chi2 99th
                        percentile): different reporter, setup or stressor
  needs a 0 µM image    paired decoder without a same-day 0 µM culture
  too few cells         fewer usable cells than the reference fields had

The reference contributes everything fitted: the autofluorescence baseline
and the reporter-positive threshold (its non-transfected cells) and the
dose classes (its transfected days). New images must therefore be taken
with the same reporter, objective and imaging settings (docs/data_contract.md).

The new cell table needs: session_id, condition_id, sample_id,
measurement_time_hours, qc_pass, edge_flag, corrected_mean_green_per_s,
corrected_mean_red_per_s (scripts/run_pipeline.py, then
scripts/normalize_exposure.py). A condition_id of the form
``<stressor>_<dose>uM_<label>`` gives that field a known dose (``h2o2_0uM_ctrl``
marks the same-day control); any other name means the dose is unknown.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
from scipy.stats import chi2

from . import induction
from . import lodo_decoder as ld
from . import timecourse as tc

UNINFORMATIVE_ENTROPY = 0.9
REQUIRED = ("session_id", "condition_id", "sample_id", "measurement_time_hours", "qc_pass", "edge_flag",
            "corrected_mean_green_per_s", "corrected_mean_red_per_s")
OFFERED = ("absolute", "paired")  # the model decoder failed validation (docs/limitations.md)


def prepare_new_cells(cells: pd.DataFrame, stressor: str, baseline: pd.DataFrame) -> pd.DataFrame:
    """Label the user's cells and subtract the reference autofluorescence baseline."""
    missing = [c for c in REQUIRED if c not in cells.columns]
    if missing:
        raise ValueError(f"cell table lacks {missing}; run run_pipeline.py and normalize_exposure.py first")
    out = cells.copy()
    out["session_id"] = out["session_id"].astype(str)
    out["condition_id"] = out["condition_id"].astype(str)
    out["sample_id"] = out["sample_id"].astype(str)
    out["stressor"], out["strain"], out["day"] = stressor, "tr", out["session_id"]
    pattern = re.compile(rf"^{re.escape(stressor)}_(\d+(?:\.\d+)?)uM(?:_|$)")
    out["dose_uM"] = [float(m.group(1)) if (m := pattern.match(c)) else np.nan for c in out["condition_id"]]
    return tc.subtract_baseline(out, baseline)


def train(reference: pd.DataFrame, stressor: str, decoder: str = "paired", quantile: float = 0.9,
          min_cells: int = 20, exclude_sessions: tuple = ()) -> dict:
    """Fit on every transfected reference day of ``stressor``; validate leave-one-day-out on the same days.

    ``reference`` is a raw exposure-normalised cell table (transfected and non-transfected sessions).
    """
    if decoder not in OFFERED:
        raise ValueError(f"decoder must be one of {OFFERED}; the model decoder is not offered "
                         "because it was confidently wrong on held-out days")
    cells = tc.annotate(reference)
    cells = cells[~cells["session_id"].isin(exclude_sessions)]
    cells = cells[cells["stressor"] == stressor]
    if cells.empty or not (cells["strain"] == "nt").any():
        raise ValueError(f"reference has no {stressor} sessions with non-transfected cells")
    baseline = tc.autofluorescence_baseline(cells)
    cells = tc.subtract_baseline(cells, baseline)
    limits = induction.positive_thresholds(cells, quantile)
    fields = ld.add_paired(ld.field_features(cells, quantile=quantile, min_cells=min_cells, limits=limits))
    if fields.empty:
        raise ValueError(f"no transfected {stressor} fields with ≥ {min_cells} usable cells in the reference")
    features = ld.FEATURES[decoder]
    doses = sorted(fields["dose_uM"].unique())
    classes = ld.fit_classes(fields, features, doses)
    predictions, _ = ld.leave_one_day_out(fields, decoder)
    validation = ld.summarize(predictions)  # empty with a single reference day: nothing to hold out
    mean = (validation[(validation["held_out"] == "mean over days") & (validation["unit"] == "field")]
            if len(validation) else validation)
    return {"stressor": stressor, "decoder": decoder, "features": features, "doses": doses, "classes": classes,
            "baseline": baseline, "limits": limits, "quantile": quantile, "min_cells": min_cells,
            "reference_days": sorted(fields["session_id"].unique()), "reference_fields": fields,
            "validation": validation,
            "validated_gain_bits": float(mean["info_gain_bits"].iloc[0]) if len(mean) else np.nan,
            "validated_accuracy": float(mean["accuracy"].iloc[0]) if len(mean) else np.nan}


def _verdict(row: dict, trained: dict) -> str:
    if row["n_cells"] < trained["min_cells"]:
        return "too few cells"
    if not row["features_available"]:
        return "needs a 0 µM image" if trained["decoder"] == "paired" else "too few cells"
    if row["outside_reference"]:
        return "unlike reference"
    if row["entropy_vs_prior"] >= UNINFORMATIVE_ENTROPY:
        return "no dose information"
    gain = trained["validated_gain_bits"]
    if not np.isfinite(gain) or gain <= 0:
        return "not validated"
    return "dose estimate"


def decode(cells: pd.DataFrame, trained: dict) -> pd.DataFrame:
    """One row per new field and per new sample (its fields pooled), with posterior and verdict."""
    prepared = prepare_new_cells(cells, trained["stressor"], trained["baseline"])
    fields = ld.field_features(prepared, quantile=trained["quantile"], min_cells=1, limits=trained["limits"])
    if fields.empty:
        raise ValueError("no usable cells in the new table")
    fields = ld.add_paired(fields)
    features, doses, classes = trained["features"], trained["doses"], trained["classes"]
    limit = float(chi2.ppf(ld.FIT_QUANTILE, len(features)))
    k = len(doses)

    units = [("field", r["session_id"], r["condition_id"], r["sample_id"], r["dose_uM"], r["time_h"],
              int(r["n_cells"]), 1, r[features].to_numpy(float)) for _, r in fields.iterrows()]
    for (session, condition), g in fields.groupby(["session_id", "condition_id"]):
        if len(g) > 1:
            units.append(("sample", session, condition, ",".join(g["sample_id"]), g["dose_uM"].iloc[0],
                          float(g["time_h"].mean()), int(g["n_cells"].min()), len(g),
                          g[features].mean().to_numpy(float)))
    rows = []
    for unit, session, condition, sample, dose, time_h, n_cells, n_fields, x in units:
        row = {"unit": unit, "session_id": session, "condition_id": condition, "sample_id": sample,
               "known_dose_uM": dose, "time_h": time_h, "n_cells": n_cells, "n_fields": n_fields,
               "features_available": bool(np.isfinite(x).all())}
        if row["features_available"]:
            p, d2 = ld.posterior(x, classes, n_fields)
            entropy = float(-(p[p > 0] * np.log2(p[p > 0])).sum())
            cs = ld.credible_set(p)
            row.update({"map_dose_uM": doses[int(np.argmax(p))],
                        "credible80_uM": " ".join(f"{doses[i]:g}" for i in sorted(cs)),
                        "entropy_vs_prior": entropy / np.log2(k), "min_mahalanobis2": d2,
                        "outside_reference": bool(d2 > limit),
                        **{f"p_{d:g}uM": float(v) for d, v in zip(doses, p)}})
        else:
            row.update({"map_dose_uM": np.nan, "credible80_uM": "", "entropy_vs_prior": np.nan,
                        "min_mahalanobis2": np.nan, "outside_reference": False})
        row["verdict"] = _verdict(row, trained)
        rows.append(row)
    return pd.DataFrame(rows)
