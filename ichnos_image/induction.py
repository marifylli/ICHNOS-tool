"""Does the reporter respond to the stress, and how large a response could these data have shown?

Every comparison is made within one day, against that day's 0 µM culture,
because each day used a different transformation and the non-transfected
cultures were imaged on other days. Within a day both share the transformation,
the slide preparation and the imaging session, so neither confound enters.

Three readouts per field, all on usable cells (``timecourse.usable``):
  green, red  median baseline-subtracted intensity (grey/s), as in timecourse.py
  positive    share of cells whose green AND red both exceed the 90th percentile
              of the non-transfected cells of that stressor: the reporter-
              expressing subpopulation, which the population median dilutes
              when only some cells carry the construct

For a stressed sample, at each of its time points, the excess is its value
minus the 0 µM curve of the same day interpolated to that time (left out
outside the control's sampled range). The mean excess over time points is the
effect. Its uncertainty comes from resampling cells within fields and fields
within time points of both the stressed and the control curve.

Reported per day, dose and readout: the effect, its 95 % interval, and the
minimum detectable effect, MDE = (1.96 + 0.84) x SE: the true effect that a
two-sided 5 % test would detect with 80 % probability given this noise. "No
effect detected, and an effect larger than the MDE would have been" is the
statement these data support when the interval covers zero.

What the resampling does not cover is day-to-day biology; days are reported
separately for that reason.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import timecourse as tc

READOUTS = ("green", "red", "positive")
Z_MDE = 1.96 + 0.84


def positive_thresholds(cells: pd.DataFrame, quantile: float = 0.9) -> pd.DataFrame:
    """Per stressor, the non-transfected quantile of each per-second channel."""
    nt = cells[(cells["strain"] == "nt") & tc.usable(cells)]
    rows = [{"stressor": s, **{ch: float(g[f"corrected_mean_{ch}_per_s"].quantile(quantile)) for ch in tc.CHANNELS}}
            for s, g in nt.groupby("stressor")]
    return pd.DataFrame(rows).set_index("stressor")


def _field_arrays(cells: pd.DataFrame, thresholds: pd.DataFrame) -> dict:
    use = cells[tc.usable(cells) & (cells["strain"] == "tr")]
    out = {}
    for key, group in use.groupby(["session_id", "condition_id", "sample_id"]):
        limit = thresholds.loc[group["stressor"].iloc[0]]
        positive = ((group["corrected_mean_green_per_s"] > limit["green"])
                    & (group["corrected_mean_red_per_s"] > limit["red"])).to_numpy(float)
        out[key] = {"green": group["green_af"].to_numpy(float), "red": group["red_af"].to_numpy(float),
                    "positive": positive}
    return out


def _point_value(arrays: dict, keys: list, rng=None) -> dict:
    """Median of field medians (green, red) and mean of field shares (positive)."""
    if rng is not None:
        keys = [keys[j] for j in rng.integers(0, len(keys), len(keys))]
    per_field = {r: [] for r in READOUTS}
    for k in keys:
        a = arrays[k]
        idx = rng.integers(0, len(a["green"]), len(a["green"])) if rng is not None else slice(None)
        per_field["green"].append(np.median(a["green"][idx]))
        per_field["red"].append(np.median(a["red"][idx]))
        per_field["positive"].append(np.mean(a["positive"][idx]))
    return {"green": float(np.median(per_field["green"])), "red": float(np.median(per_field["red"])),
            "positive": float(np.mean(per_field["positive"]))}


def _points(cells: pd.DataFrame, arrays: dict) -> pd.DataFrame:
    """One row per transfected sample and time point, with the field keys behind it."""
    use = cells[tc.usable(cells) & (cells["strain"] == "tr")]
    table = (use.groupby(["session_id", "day", "stressor", "dose_uM", "condition_id"])
             .agg(time_h=("measurement_time_hours", "median"),
                  samples=("sample_id", lambda s: sorted(set(s))))
             .reset_index())
    table["keys"] = [[(r.session_id, r.condition_id, s) for s in r.samples if (r.session_id, r.condition_id, s) in arrays]
                     for r in table.itertuples()]
    return table[table["keys"].map(len) > 0].reset_index(drop=True)


def _excess(points: pd.DataFrame, values: dict) -> dict:
    """Mean over time points of (value - same-day 0 µM interpolated), per stressed sample."""
    out = {}
    for session, day_points in points.groupby("session_id"):
        control = day_points[day_points["dose_uM"] == 0].sort_values("time_h")
        if len(control) < 2:
            continue
        t_ctrl = control["time_h"].to_numpy(float)
        for dose, sample in day_points[day_points["dose_uM"] > 0].groupby("dose_uM"):
            inside = sample[(sample["time_h"] >= t_ctrl.min()) & (sample["time_h"] <= t_ctrl.max())]
            if inside.empty:
                continue
            for readout in READOUTS:
                v_ctrl = np.array([values[i][readout] for i in control.index])
                diffs = [values[i][readout] - np.interp(t, t_ctrl, v_ctrl)
                         for i, t in zip(inside.index, inside["time_h"])]
                out[(session, dose, readout)] = (float(np.mean(diffs)), len(diffs),
                                                 float(np.mean(v_ctrl)))
    return out


def induction_table(cells: pd.DataFrame, n_boot: int = 500, seed: int = 0,
                    quantile: float = 0.9) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Effect of each dose against the same day's 0 µM, with interval and MDE; and the dose trend.

    Returns (effects, trend). ``effects`` has one row per day, dose and
    readout with a bootstrap p-value and its Benjamini-Hochberg q over all rows.

    ``cells`` must already carry ``green_af``/``red_af`` (``timecourse.subtract_baseline``).
    """
    thresholds = positive_thresholds(cells, quantile)
    arrays = _field_arrays(cells, thresholds)
    points = _points(cells, arrays)
    observed = _excess(points, {i: _point_value(arrays, k) for i, k in points["keys"].items()})
    rng = np.random.default_rng(seed)
    replicates = []
    for _ in range(n_boot):
        values = {i: _point_value(arrays, k, rng) for i, k in points["keys"].items()}
        replicates.append({key: effect for key, (effect, _, _) in _excess(points, values).items()})
    draws = {key: [r[key] for r in replicates if key in r] for key in observed}

    meta = points.drop_duplicates("session_id").set_index("session_id")[["day", "stressor"]]
    rows = []
    for (session, dose, readout), (effect, n_points, control_level) in observed.items():
        reps = np.asarray(draws[(session, dose, readout)], float)
        se = float(np.std(reps, ddof=1)) if len(reps) > 1 else np.nan
        lo, hi = (np.percentile(reps, [2.5, 97.5]) if len(reps) else (np.nan, np.nan))
        rows.append({
            "session_id": session, "day": meta.loc[session, "day"], "stressor": meta.loc[session, "stressor"],
            "dose_uM": dose, "readout": readout, "n_timepoints": n_points,
            "control_level": control_level, "effect": effect, "ci_low": float(lo), "ci_high": float(hi),
            "se": se, "mde": Z_MDE * se,
            "ci_excludes_0": bool(lo > 0 or hi < 0),
        })
    table = pd.DataFrame(rows).sort_values(["stressor", "day", "readout", "dose_uM"]).reset_index(drop=True)
    table["p_two_sided"] = [_bootstrap_p(np.asarray(draws[(r.session_id, r.dose_uM, r.readout)], float))
                            for r in table.itertuples()]
    table["q_bh"] = _benjamini_hochberg(table["p_two_sided"].to_numpy())
    table["detected_bh"] = table["q_bh"] < 0.05
    trend = _trend(observed, replicates, meta)
    return table, trend


def _bootstrap_p(reps: np.ndarray) -> float:
    """Two-sided p from the share of replicates on the far side of zero (floored at 1/(n+1))."""
    if len(reps) == 0:
        return np.nan
    tail = min(np.mean(reps <= 0), np.mean(reps >= 0))
    return float(min(1.0, max(2 * tail, 1.0 / (len(reps) + 1))))


def _benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """False-discovery-rate adjusted p-values: many doses x days x readouts are tested at once."""
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return out


def _trend(observed: dict, replicates: list, meta: pd.DataFrame) -> pd.DataFrame:
    """Per day and readout: slope of the effect against log10(dose), with a 95 % interval.

    A response to the stress should grow with dose; effects of one sign at
    every dose but without order can come from the control curve alone.
    """
    groups = {}
    for (session, dose, readout) in observed:
        groups.setdefault((session, readout), []).append(dose)
    rows = []
    for (session, readout), doses in groups.items():
        doses = sorted(doses)
        if len(doses) < 3:
            continue
        x = np.log10(np.asarray(doses, float))
        y = np.array([observed[(session, d, readout)][0] for d in doses])
        slope = float(np.polyfit(x, y, 1)[0])
        reps = []
        for r in replicates:
            if all((session, d, readout) in r for d in doses):
                reps.append(np.polyfit(x, [r[(session, d, readout)] for d in doses], 1)[0])
        reps = np.asarray(reps, float)
        lo, hi = np.percentile(reps, [2.5, 97.5]) if len(reps) else (np.nan, np.nan)
        rows.append({"session_id": session, "day": meta.loc[session, "day"], "stressor": meta.loc[session, "stressor"],
                     "readout": readout, "n_doses": len(doses), "slope_per_log10_dose": slope,
                     "ci_low": float(lo), "ci_high": float(hi), "p_two_sided": _bootstrap_p(reps)})
    columns = ["session_id", "day", "stressor", "readout", "n_doses", "slope_per_log10_dose",
               "ci_low", "ci_high", "p_two_sided"]
    # fewer than 3 doses on every day: no trend can be fitted, and the table stays empty
    trend = pd.DataFrame(rows, columns=columns).sort_values(["stressor", "day", "readout"]).reset_index(drop=True)
    trend["q_bh"] = _benjamini_hochberg(trend["p_two_sided"].to_numpy())
    return trend
