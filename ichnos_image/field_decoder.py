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

Time since the stress (the second half of the question) is decoded the same
way into time windows (default 0–2, 2–4, 4–6, 6–8 h), from the tandem timer:
log(red/green) of the reporter-positive cells, which by design does not
depend on how much reporter a cell carries, plus their green level. Trained on
the stressed reference fields of all doses, validated leave-one-day-out, with
its own verdicts (time estimate / not validated / no time information /
unlike reference / too few reporter cells / unstressed control / time not
available). A held-out field whose time window no other reference day covers
cannot be scored and is counted as unscorable, not as correct or wrong.
The paired dose decoder needs a time only to line the 0 µM culture up with
the sample; with no time given, the day's 0 µM images are averaged.

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
TIME_FEATURES = ["log_timer_ratio", "log_green_pos"]
DEFAULT_TIME_EDGES = (0.0, 2.0, 4.0, 6.0, 8.0)


def add_timer_features(fields: pd.DataFrame) -> pd.DataFrame:
    """log(red/green) and log(green) of each field's median reporter-positive cell (NaN if not both > 0)."""
    out = fields.copy()
    green, red = out["green_pos"].astype(float), out["red_pos"].astype(float)
    ok = (green > 0) & (red > 0)
    out["log_timer_ratio"] = np.log(red.where(ok) / green.where(ok))
    out["log_green_pos"] = np.log(green.where(ok))
    return out


def time_labels(edges) -> list[str]:
    return [f"{a:g}–{b:g} h" for a, b in zip(edges[:-1], edges[1:])]


def time_bin(hours, edges) -> np.ndarray:
    """Index of the window holding each time (the last window includes its upper edge); -1 outside or unknown."""
    t = np.asarray(hours, float)
    index = np.searchsorted(np.asarray(edges, float), t, side="right") - 1
    index[np.isclose(t, edges[-1])] = len(edges) - 2
    index[~np.isfinite(t) | (index < 0) | (index > len(edges) - 2)] = -1
    return index


def _time_track_record(stressed: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Leave-one-day-out for the time decoder, scoring only fields whose window the other days cover."""
    rows, per_day = [], []
    for held in sorted(stressed["session_id"].unique()):
        train, test = stressed[stressed["session_id"] != held], stressed[stressed["session_id"] == held]
        bins = sorted(train["dose_uM"].unique())
        try:
            classes = ld.fit_classes(train, TIME_FEATURES, bins) if len(bins) >= 2 else None
        except ValueError:
            classes = None
        scored = test[test["dose_uM"].isin(bins)] if classes is not None else test.iloc[0:0]
        day_rows = []
        for _, r in scored.iterrows():
            p, _ = ld.posterior(r[TIME_FEATURES].to_numpy(float), classes)
            truth, k = bins.index(r["dose_uM"]), len(bins)
            entropy = float(-(p[p > 0] * np.log2(p[p > 0])).sum())
            day_rows.append({"correct": int(np.argmax(p)) == truth, "chance": 1 / k,
                             "gain": float(np.log2(max(p[truth], 1e-300)) + np.log2(k)),
                             "entropy_rel": entropy / np.log2(k), "covered": truth in ld.credible_set(p)})
        rows += day_rows
        d = pd.DataFrame(day_rows)
        per_day.append({"held_out": held, "n_scored": len(d), "n_unscorable": len(test) - len(d),
                        "windows_in_training": len(bins),
                        **({"accuracy": d["correct"].mean(), "chance": d["chance"].mean(),
                            "info_gain_bits": d["gain"].mean(), "entropy_vs_prior": d["entropy_rel"].mean(),
                            "credible80_coverage": d["covered"].mean()} if len(d) else {})})
    table = pd.DataFrame(per_day)
    pooled = pd.DataFrame(rows)
    summary = ({"n_scored": len(pooled), "n_unscorable": int(table["n_unscorable"].sum()) if len(table) else 0,
                "accuracy": pooled["correct"].mean(), "chance": pooled["chance"].mean(),
                "info_gain_bits": pooled["gain"].mean(), "entropy_vs_prior": pooled["entropy_rel"].mean(),
                "credible80_coverage": pooled["covered"].mean()} if len(pooled) else
               {"n_scored": 0, "n_unscorable": int(table["n_unscorable"].sum()) if len(table) else 0})
    if len(table):
        table = pd.concat([table, pd.DataFrame([{"held_out": "all held-out fields", **summary}])], ignore_index=True)
    return table, summary


def train_time(fields: pd.DataFrame, edges=DEFAULT_TIME_EDGES) -> dict:
    """Time-window classes from the stressed reference fields; unavailable (with the reason) when they cannot be fitted."""
    edges = tuple(float(e) for e in edges)
    if len(edges) < 3 or any(b <= a for a, b in zip(edges[:-1], edges[1:])) or edges[0] < 0:
        raise ValueError("time windows need at least 3 increasing edges ≥ 0, e.g. 0,2,4,6,8")
    labels = time_labels(edges)
    stressed = add_timer_features(fields[fields["dose_uM"] > 0])
    stressed = stressed.assign(window=time_bin(stressed["time_h"], edges)).dropna(subset=TIME_FEATURES)
    stressed = stressed[stressed["window"] >= 0].assign(dose_uM=lambda f: f["window"].astype(float))
    base = {"edges": edges, "labels": labels}
    bins = sorted(stressed["dose_uM"].unique())
    if len(bins) < 2:
        return {**base, "available": False,
                "reason": "the reference's stressed fields with reporter-positive cells cover fewer than 2 time windows"}
    try:
        classes = ld.fit_classes(stressed, TIME_FEATURES, bins)
    except ValueError as exc:
        return {**base, "available": False, "reason": str(exc)}
    validation, record = _time_track_record(stressed)
    return {**base, "available": True, "classes": classes, "windows": [int(b) for b in bins],
            "n_fields": len(stressed), "validation": validation, "record": record}


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
          min_cells: int = 20, exclude_sessions: tuple = (), time_edges=DEFAULT_TIME_EDGES) -> dict:
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
            "validated_accuracy": float(mean["accuracy"].iloc[0]) if len(mean) else np.nan,
            "time": train_time(fields, time_edges)}


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


def _time_verdict(row: dict, trained: dict) -> str:
    time = trained["time"]
    if row["known_dose_uM"] == 0:
        return "unstressed control"
    if not time["available"]:
        return "time not available"
    if row["n_cells"] < trained["min_cells"]:
        return "too few cells"
    if not row["time_features_available"]:
        return "too few reporter cells"
    if row["time_outside_reference"]:
        return "unlike reference"
    if row["time_entropy_vs_prior"] >= UNINFORMATIVE_ENTROPY:
        return "no time information"
    gain = time["record"].get("info_gain_bits", np.nan)
    if not np.isfinite(gain) or gain <= 0:
        return "not validated"
    return "time estimate"


def _decode_time(x: np.ndarray, n_fields: int, time: dict) -> dict:
    if not time["available"]:
        return {"time_features_available": bool(np.isfinite(x).all()), "time_map": "", "time_credible80": "",
                "time_entropy_vs_prior": np.nan, "time_outside_reference": False}
    labels = [time["labels"][w] for w in time["windows"]]
    if not np.isfinite(x).all():
        return {"time_features_available": False, "time_map": "", "time_credible80": "",
                "time_entropy_vs_prior": np.nan, "time_outside_reference": False}
    p, d2 = ld.posterior(x, time["classes"], n_fields)
    entropy = float(-(p[p > 0] * np.log2(p[p > 0])).sum())
    return {"time_features_available": True, "time_map": labels[int(np.argmax(p))],
            "time_credible80": ", ".join(labels[i] for i in sorted(ld.credible_set(p))),
            "time_entropy_vs_prior": entropy / np.log2(len(labels)),
            "time_outside_reference": bool(d2 > chi2.ppf(ld.FIT_QUANTILE, len(TIME_FEATURES))),
            **{f"pt_{label.replace('–', '-').replace(' ', '')}": float(v) for label, v in zip(labels, p)}}


def decode(cells: pd.DataFrame, trained: dict) -> pd.DataFrame:
    """One row per new field and per new sample (its fields pooled), with posterior and verdict."""
    prepared = prepare_new_cells(cells, trained["stressor"], trained["baseline"])
    fields = ld.field_features(prepared, quantile=trained["quantile"], min_cells=1, limits=trained["limits"])
    if fields.empty:
        raise ValueError("no usable cells in the new table")
    fields = add_timer_features(ld.add_paired(fields))
    features, doses, classes = trained["features"], trained["doses"], trained["classes"]
    limit = float(chi2.ppf(ld.FIT_QUANTILE, len(features)))
    k = len(doses)

    units = [("field", r["session_id"], r["condition_id"], r["sample_id"], r["dose_uM"], r["time_h"],
              int(r["n_cells"]), 1, r[features].to_numpy(float), r[TIME_FEATURES].to_numpy(float))
             for _, r in fields.iterrows()]
    for (session, condition), g in fields.groupby(["session_id", "condition_id"]):
        if len(g) > 1:
            units.append(("sample", session, condition, ",".join(g["sample_id"]), g["dose_uM"].iloc[0],
                          float(g["time_h"].mean()), int(g["n_cells"].min()), len(g),
                          g[features].mean().to_numpy(float), g[TIME_FEATURES].mean().to_numpy(float)))
    rows = []
    for unit, session, condition, sample, dose, time_h, n_cells, n_fields, x, x_time in units:
        row = {"unit": unit, "session_id": session, "condition_id": condition, "sample_id": sample,
               "known_dose_uM": dose, "known_time_h": time_h, "n_cells": n_cells, "n_fields": n_fields,
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
        row.update(_decode_time(x_time, n_fields, trained["time"]))
        row["time_verdict"] = _time_verdict(row, trained)
        if row["time_verdict"] == "unstressed control":  # no stress, so no time since it: show no answer
            row.update({k: ("" if isinstance(v, str) else np.nan) for k, v in row.items()
                        if k.startswith(("time_map", "time_credible", "pt_"))})
        rows.append(row)
    return pd.DataFrame(rows)
