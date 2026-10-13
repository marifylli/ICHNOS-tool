"""Decode the dose of one field, validated leave-one-day-out.

The question a user asks is: from one image of my culture, which dose was
it exposed to? This module answers it with a posterior over the doses of
the experiment and scores that answer honestly: each transfected day of a
stressor is held out in turn, everything fitted from data is fitted on the
other days of that stressor, and the held-out day is decoded one field at a
time (and, as a second unit, the fields of one sample pooled).

Features of a field (usable transfected cells, after the non-transfected
autofluorescence baseline; see timecourse.py):
  green, red      median baseline-subtracted intensity (grey/s)
  positive        share of reporter-positive cells (green AND red above the
                  non-transfected 90th percentile of that stressor)
  green_pos,      median of the reporter-positive cells only, used by the
  red_pos         model decoder

Three decoders, all with a uniform prior over the doses:
  absolute  the field's (green, red, positive) as they are. What a user gets
            with one image and nothing else.
  paired    the same minus the same day's 0 µM, interpolated to the field's
            time (fields of the scored sample left out of the reference).
            What a user gets who images an unstressed culture alongside.
  model     log fold of (green_pos, red_pos) over the same day's 0 µM, with
            class means taken from the ICHNOS model instead of the data;
            only the spread is fitted on the training days.

Likelihood: Gaussian, mean per dose, covariance = between-day + within-day
/ n_fields. Class means are the average of the training days' means, the
between-day covariance their spread around it; with a single training day
(DTT) it cannot be estimated and is set to zero, which makes that fold's
posteriors overconfident. Fields of one sample are not independent: pooling
averages their features instead of multiplying their likelihoods.

Scores, per held-out day: accuracy of the most probable dose against chance
(1/K), mean posterior on the true dose, information gain in bits over the
prior (log2 p_true + log2 K: positive is better than guessing, negative is
confidently wrong), entropy relative to the prior's (1 = no information), and
the coverage and size of the 80 % credible set, and the share of units
flagged as unlike every dose of the training days (squared Mahalanobis
distance to the nearest dose above the 99th percentile of chi2).
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import chi2

from . import induction
from . import timecourse as tc

FEATURES = {
    "absolute": ["green", "red", "positive"],
    "paired": ["d_green", "d_red", "d_positive"],
    "model": ["log_fold_green_pos", "log_fold_red_pos"],
}
CREDIBLE = 0.8
# A field farther than this quantile of chi2(k) from every dose looks like none of
# the training data; the posterior is still reported, but flagged.
FIT_QUANTILE = 0.99


def field_features(cells: pd.DataFrame, quantile: float = 0.9, min_cells: int = 20,
                   min_positive: int = 5, limits: pd.DataFrame | None = None) -> pd.DataFrame:
    """One row per transfected field with at least ``min_cells`` usable cells.

    ``cells`` must be annotated and carry ``green_af``/``red_af``.
    green_pos/red_pos are NaN when the field has fewer than ``min_positive`` positive cells.
    ``limits``: reporter-positive thresholds per stressor (``induction.positive_thresholds``);
    taken from the non-transfected cells of ``cells`` when not given.
    """
    if limits is None:
        limits = induction.positive_thresholds(cells, quantile)
    use = cells[tc.usable(cells) & (cells["strain"] == "tr") & cells["stressor"].isin(limits.index)]
    rows = []
    for (session, condition, sample), g in use.groupby(["session_id", "condition_id", "sample_id"]):
        if len(g) < min_cells:
            continue
        limit = limits.loc[g["stressor"].iloc[0]]
        pos = ((g["corrected_mean_green_per_s"] > limit["green"])
               & (g["corrected_mean_red_per_s"] > limit["red"])).to_numpy()
        enough = pos.sum() >= min_positive
        rows.append({
            "session_id": session, "day": g["day"].iloc[0], "stressor": g["stressor"].iloc[0],
            "dose_uM": float(g["dose_uM"].iloc[0]), "condition_id": condition, "sample_id": sample,
            "time_h": float(g["measurement_time_hours"].median()), "n_cells": len(g), "n_positive": int(pos.sum()),
            "green": float(g["green_af"].median()), "red": float(g["red_af"].median()),
            "positive": float(pos.mean()),
            "green_pos": float(g.loc[pos, "green_af"].median()) if enough else np.nan,
            "red_pos": float(g.loc[pos, "red_af"].median()) if enough else np.nan,
        })
    return pd.DataFrame(rows)


def _reference(day: pd.DataFrame, column: str, time: float, exclude_condition: str) -> float:
    """The day's 0 µM value of ``column`` at ``time``: per-condition means, linear in time, flat outside."""
    control = day[(day["dose_uM"] == 0) & (day["condition_id"] != exclude_condition)].dropna(subset=[column])
    if control.empty:
        return np.nan
    points = control.groupby("condition_id").agg(t=("time_h", "mean"), v=(column, "mean")).sort_values("t")
    return float(np.interp(time, points["t"], points["v"]))


def add_paired(fields: pd.DataFrame) -> pd.DataFrame:
    """Add the paired features and the log folds of positive cells over the same day's 0 µM."""
    out = fields.copy()
    for column in ("d_green", "d_red", "d_positive", "log_fold_green_pos", "log_fold_red_pos"):
        out[column] = np.nan
    for session, day in out.groupby("session_id"):
        for i, row in day.iterrows():
            for name in ("green", "red", "positive"):
                out.at[i, f"d_{name}"] = row[name] - _reference(day, name, row["time_h"], row["condition_id"])
            for name in ("green_pos", "red_pos"):
                ref = _reference(day, name, row["time_h"], row["condition_id"])
                if row[name] > 0 and ref > 0:
                    out.at[i, f"log_fold_{name}"] = float(np.log(row[name] / ref))
    return out


def fit_classes(train: pd.DataFrame, features: list[str], doses: list[float]) -> dict:
    """Per-dose means (average of training-day means) and the between- and within-day covariances."""
    train = train.dropna(subset=features)
    k = len(features)
    day_means = train.groupby(["session_id", "dose_uM"])[features].mean()
    missing = [d for d in doses if d not in day_means.index.get_level_values("dose_uM")]
    if missing:
        raise ValueError(f"no training fields at dose(s) {missing}")
    means = np.array([day_means.xs(d, level="dose_uM").mean().to_numpy() for d in doses])

    within = train[features].to_numpy() - day_means.loc[list(zip(train["session_id"], train["dose_uM"]))].to_numpy()
    dof_within = len(train) - len(day_means)
    cov_within = within.T @ within / dof_within if dof_within > 0 else np.zeros((k, k))

    between = day_means.to_numpy() - means[[doses.index(d) for d in day_means.index.get_level_values("dose_uM")]]
    dof_between = len(day_means) - len(doses)
    cov_between = between.T @ between / dof_between if dof_between > 0 else np.zeros((k, k))
    return {"doses": list(doses), "features": features, "means": means, "within": cov_within,
            "between": cov_between, "between_estimated": dof_between > 0,
            "n_training_days": int(train["session_id"].nunique())}


def mahalanobis2(x: np.ndarray, means: np.ndarray, covariance: np.ndarray) -> tuple[np.ndarray, float]:
    """Squared Mahalanobis distance of x to each row of ``means``, and the log-determinant term.

    A tiny ridge keeps a singular covariance usable.
    """
    cov = covariance + 1e-9 * np.eye(len(x)) * max(1.0, float(np.trace(covariance)))
    chol = np.linalg.cholesky(cov)
    z = np.linalg.solve(chol, (x[None, :] - means).T)
    return np.sum(z**2, axis=0), float(np.log(np.diag(chol)).sum())


def posterior(x: np.ndarray, classes: dict, n_fields: int = 1, means: np.ndarray | None = None
              ) -> tuple[np.ndarray, float]:
    """Posterior over the doses (uniform prior), and the smallest squared distance to any dose."""
    covariance = classes["between"] + classes["within"] / n_fields
    d2, logdet = mahalanobis2(np.asarray(x, float), classes["means"] if means is None else means, covariance)
    ll = -0.5 * d2 - logdet
    return np.exp(ll - logsumexp(ll)), float(d2.min())


def credible_set(p: np.ndarray, mass: float = CREDIBLE) -> np.ndarray:
    """Indices of the most probable doses until their mass reaches ``mass``."""
    order = np.argsort(p)[::-1]
    stop = int(np.searchsorted(np.cumsum(p[order]), mass - 1e-12)) + 1
    return order[:stop]


def _units(test: pd.DataFrame, features: list[str]) -> list[dict]:
    """Single fields, and the fields of one sample at one time pooled (feature means)."""
    units = []
    for _, row in test.iterrows():
        units.append({"unit": "field", "condition_id": row["condition_id"], "sample_id": row["sample_id"],
                      "dose_uM": row["dose_uM"], "time_h": row["time_h"], "n_fields": 1,
                      "x": row[features].to_numpy(float)})
    for condition, g in test.groupby("condition_id"):
        units.append({"unit": "sample", "condition_id": condition, "sample_id": ",".join(g["sample_id"]),
                      "dose_uM": g["dose_uM"].iloc[0], "time_h": float(g["time_h"].mean()), "n_fields": len(g),
                      "x": g[features].mean().to_numpy(float)})
    return units


ModelLogFold = Callable[[str, float, float], np.ndarray]


def leave_one_day_out(fields: pd.DataFrame, decoder: str, model_log_fold: ModelLogFold | None = None
                      ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Posterior for every field (and pooled sample) of every held-out day.

    Returns (predictions, folds). ``model_log_fold(stressor, dose, time)``
    gives the model's (log fold green, log fold red) and is required for
    decoder="model".
    """
    features = FEATURES[decoder]
    if decoder == "model" and model_log_fold is None:
        raise ValueError("decoder 'model' needs model_log_fold")
    rows, folds = [], []
    for stressor, part in fields.groupby("stressor"):
        doses = sorted(part["dose_uM"].unique())
        days = sorted(part["session_id"].unique())
        if len(days) < 2:
            continue
        for held_out in days:
            train = part[part["session_id"] != held_out]
            test = part[part["session_id"] == held_out].dropna(subset=features)
            try:
                classes = fit_classes(train, features, doses)
            except ValueError as exc:
                folds.append({"stressor": stressor, "held_out": held_out, "decoder": decoder, "skipped": str(exc)})
                continue
            folds.append({"stressor": stressor, "held_out": held_out, "decoder": decoder, "skipped": "",
                          "n_training_days": classes["n_training_days"],
                          "between_day_estimated": classes["between_estimated"],
                          "n_test_fields": len(test)})
            fit_limit = float(chi2.ppf(FIT_QUANTILE, len(features)))
            for unit in _units(test, features):
                means = None
                if decoder == "model":
                    means = np.array([model_log_fold(stressor, d, unit["time_h"]) for d in doses])
                p, d2 = posterior(unit["x"], classes, unit["n_fields"], means)
                truth = doses.index(unit["dose_uM"])
                cs = credible_set(p)
                rows.append({
                    "stressor": stressor, "held_out": held_out, "decoder": decoder, "unit": unit["unit"],
                    "condition_id": unit["condition_id"], "sample_id": unit["sample_id"],
                    "time_h": unit["time_h"], "n_fields": unit["n_fields"], "true_dose_uM": unit["dose_uM"],
                    "map_dose_uM": doses[int(np.argmax(p))], "p_true": float(p[truth]),
                    "entropy_bits": float(-(p[p > 0] * np.log2(p[p > 0])).sum()), "n_doses": len(doses),
                    "credible_set_uM": " ".join(f"{doses[i]:g}" for i in sorted(cs)),
                    "credible_covers_truth": bool(truth in cs), "credible_size": len(cs),
                    "min_mahalanobis2": d2, "outside_training": bool(d2 > fit_limit),
                    **{f"p_{d:g}uM": float(v) for d, v in zip(doses, p)},
                })
    return pd.DataFrame(rows), pd.DataFrame(folds)


def summarize(predictions: pd.DataFrame) -> pd.DataFrame:
    """Scores per stressor, decoder, unit and held-out day, plus the mean over held-out days."""
    if predictions.empty:
        return pd.DataFrame()
    p = predictions.assign(
        correct=predictions["map_dose_uM"] == predictions["true_dose_uM"],
        gain_bits=np.log2(predictions["p_true"].clip(lower=1e-300)) + np.log2(predictions["n_doses"]),
        entropy_rel=predictions["entropy_bits"] / np.log2(predictions["n_doses"]))
    keys = ["stressor", "decoder", "unit", "held_out"]
    table = p.groupby(keys).agg(
        n=("correct", "size"), accuracy=("correct", "mean"), chance=("n_doses", lambda k: 1 / k.iloc[0]),
        mean_p_true=("p_true", "mean"), info_gain_bits=("gain_bits", "mean"),
        entropy_vs_prior=("entropy_rel", "mean"), credible80_coverage=("credible_covers_truth", "mean"),
        credible80_size=("credible_size", "mean"), outside_training=("outside_training", "mean"),
        n_doses=("n_doses", "first")).reset_index()
    mean = (table.groupby(["stressor", "decoder", "unit"])
            .agg({c: "mean" for c in table.columns if c not in keys})
            .reset_index().assign(held_out="mean over days"))
    return pd.concat([table, mean], ignore_index=True).sort_values(keys[:3] + ["held_out"]).reset_index(drop=True)
