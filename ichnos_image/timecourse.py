"""Population time courses of the reporter, the autofluorescence baseline, and when the signal peaks.

Inputs are the per-cell rows of the image pipeline after exposure normalisation
(``scripts/normalize_exposure.py``), so intensities are grey levels per second.

Steps and the choices behind them:

1. Autofluorescence baseline. Non-transfected cells carry no reporter, so their
   corrected intensity is autofluorescence (plus residual background and, for
   red, bleed-through from green autofluorescence). The baseline is the median
   over every QC-passing non-transfected cell of the same stressor, one number
   per channel. A single number per stressor is a choice: the non-transfected
   cultures were imaged on other days (8/10 H2O2, 9/10 DTT) than the transfected
   ones, so a baseline matched by day is impossible, and one matched by dose and
   time would carry that day's noise into every transfected point. How much the
   baseline moves with dose and time is reported (``baseline_diagnostics``) so
   the choice can be judged. Subtracting a constant does not move a peak in
   time; it does change fold changes and the red/green ratio.

2. Field summary. One value per field and channel: the median over its
   QC-passing cells that do not touch the frame edge.

3. Time course. Fields of one sample at one time (R1-R3, nearby fields of the
   same slide) are pooled by the median of their field medians. They are not
   biological replicates; days are.

4. Peak. The time point with the highest value. If it is the first or last
   sampled time the peak is not bracketed and only a bound is reported.
   Otherwise a parabola through it and its two neighbours gives a time between
   samples. Uncertainty comes from resampling cells within fields and fields
   within time points; it covers the measurement of one culture, not
   day-to-day biology, which is judged by comparing days.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

CHANNELS = ("green", "red")
_DOSE = re.compile(r"^(?P<stressor>[a-z0-9]+)_(?P<dose>\d+(?:\.\d+)?)uM(?:_|$)")
_SESSION = re.compile(r"^(?P<day>\d{8})_(?P<stressor>[a-z0-9]+)_(?P<strain>tr|nt)$")


def annotate(cells: pd.DataFrame) -> pd.DataFrame:
    """Add ``day``, ``stressor``, ``strain`` (tr/nt) and ``dose_uM`` parsed from the identifiers.

    Rows whose session is not a stressor session (e.g. the medium controls,
    ``..._none_tr``) get ``stressor = "none"`` and no dose. Anything else that
    does not parse raises, rather than being dropped silently.
    """
    out = cells.copy()
    sessions = out["session_id"].astype(str).str.extract(_SESSION)
    if sessions.isna().any(axis=None):
        bad = sorted(out.loc[sessions.isna().any(axis=1), "session_id"].astype(str).unique())
        raise ValueError(f"session_id not of the form <yyyymmdd>_<stressor>_<tr|nt>: {bad}")
    out["day"], out["stressor"], out["strain"] = sessions["day"], sessions["stressor"], sessions["strain"]
    doses = out["condition_id"].astype(str).str.extract(_DOSE)
    stress = out["stressor"] != "none"
    unparsed = stress & (doses["dose"].isna() | (doses["stressor"] != out["stressor"]))
    if unparsed.any():
        bad = sorted(out.loc[unparsed, "condition_id"].astype(str).unique())
        raise ValueError(f"condition_id does not start with <stressor>_<dose>uM matching its session: {bad}")
    out["dose_uM"] = pd.to_numeric(doses["dose"]).where(stress)
    return out


def usable(cells: pd.DataFrame) -> pd.Series:
    """QC-passing cells that do not touch the frame edge."""
    qc = cells["qc_pass"].astype(str).str.lower().isin({"true", "1"})
    edge = cells["edge_flag"].astype(str).str.lower().isin({"true", "1"})
    return qc & ~edge


def _per_s(channel: str) -> str:
    return f"corrected_mean_{channel}_per_s"


def autofluorescence_baseline(cells: pd.DataFrame) -> pd.DataFrame:
    """One baseline per stressor and channel, from non-transfected cells."""
    nt = cells[(cells["strain"] == "nt") & usable(cells)]
    if nt.empty:
        raise ValueError("no usable non-transfected cells to take the autofluorescence baseline from")
    rows = []
    for stressor, group in nt.groupby("stressor"):
        row = {"stressor": stressor, "n_cells": len(group), "n_fields": group["sample_id"].nunique(),
               "days": ",".join(sorted(group["day"].unique()))}
        for ch in CHANNELS:
            values = group[_per_s(ch)].astype(float)
            row[f"baseline_{ch}"] = float(values.median())
            row[f"baseline_{ch}_q25"] = float(values.quantile(0.25))
            row[f"baseline_{ch}_q75"] = float(values.quantile(0.75))
        rows.append(row)
    return pd.DataFrame(rows)


def baseline_diagnostics(cells: pd.DataFrame) -> pd.DataFrame:
    """Non-transfected medians per dose and time, to judge the constant-baseline choice."""
    nt = cells[(cells["strain"] == "nt") & usable(cells)]
    keys = ["stressor", "dose_uM", "condition_id"]
    agg = {f"median_{ch}": (_per_s(ch), "median") for ch in CHANNELS}
    table = nt.groupby(keys).agg(n_cells=("cell_id", "size"), time_h=("measurement_time_hours", "median"), **agg)
    return table.reset_index().sort_values(["stressor", "dose_uM", "time_h"])


def subtract_baseline(cells: pd.DataFrame, baseline: pd.DataFrame) -> pd.DataFrame:
    """Add ``<channel>_af``: per-second intensity minus the stressor's baseline."""
    out = cells.copy()
    table = baseline.set_index("stressor")
    missing = sorted(set(out.loc[out["stressor"] != "none", "stressor"]) - set(table.index))
    if missing:
        raise ValueError(f"no baseline for stressor(s) {missing}")
    for ch in CHANNELS:
        offset = out["stressor"].map(table[f"baseline_{ch}"])
        out[f"{ch}_af"] = out[_per_s(ch)].astype(float) - offset
    return out


FIELD_KEYS = ["session_id", "day", "stressor", "strain", "dose_uM", "condition_id", "sample_id"]


def field_summaries(cells: pd.DataFrame) -> pd.DataFrame:
    """Median per field of the baseline-subtracted channels, over usable cells."""
    use = cells[usable(cells) & (cells["stressor"] != "none")]
    agg = {f"{ch}_af": (f"{ch}_af", "median") for ch in CHANNELS}
    table = use.groupby(FIELD_KEYS, dropna=False).agg(
        n_cells=("cell_id", "size"), time_h=("measurement_time_hours", "median"), **agg)
    return table.reset_index()


CURVE_KEYS = ["session_id", "day", "stressor", "strain", "dose_uM"]


def _ratio(red: np.ndarray, green: np.ndarray, floor: float) -> np.ndarray:
    """red/green where both are above ``floor``; a ratio of values at or below the
    baseline measures noise, not the timer."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where((green > floor) & (red > floor), red / green, np.nan)


def timecourse(fields: pd.DataFrame, ratio_floor: float = 0.0) -> pd.DataFrame:
    """One point per sample and time: median of its fields' medians.

    ``ratio`` is red/green of those pooled medians, left empty where either
    baseline-subtracted channel is at or below ``ratio_floor``.
    """
    keys = CURVE_KEYS + ["condition_id"]
    table = fields.groupby(keys, dropna=False).agg(
        time_h=("time_h", "median"), n_fields=("sample_id", "size"), n_cells=("n_cells", "sum"),
        green=("green_af", "median"), red=("red_af", "median")).reset_index()
    table["ratio"] = _ratio(table["red"].to_numpy(), table["green"].to_numpy(), ratio_floor)
    return table.sort_values(CURVE_KEYS + ["time_h"]).reset_index(drop=True)


def refine_peak(times: np.ndarray, values: np.ndarray) -> dict:
    """Peak of one curve: the sampled maximum, and a parabola-refined time if bracketed."""
    times, values = np.asarray(times, float), np.asarray(values, float)
    ok = np.isfinite(times) & np.isfinite(values)
    times, values = times[ok], values[ok]
    order = np.argsort(times)
    times, values = times[order], values[order]
    if len(times) < 2:
        return {"t_max_sampled": np.nan, "t_peak": np.nan, "bracketed": False, "position": "too_few_points"}
    i = int(np.argmax(values))
    result = {"t_max_sampled": float(times[i]), "t_peak": float(times[i]), "bracketed": False}
    if i == 0:
        result["position"] = "first_sample"      # falling from the start: peak at or before t_first
    elif i == len(times) - 1:
        result["position"] = "last_sample"       # still rising: peak at or after t_last
    else:
        result["position"] = "bracketed"
        result["bracketed"] = True
        t, v = times[i - 1:i + 2], values[i - 1:i + 2]
        a, b, _ = np.polyfit(t, v, 2)            # exact through three points
        if a < 0:
            vertex = -b / (2 * a)
            result["t_peak"] = float(np.clip(vertex, t[0], t[2]))
    return result


def _cell_arrays(cells: pd.DataFrame) -> dict:
    """Usable cells per field. A field is (session, condition, sample): sample_id
    alone repeats across time points of one series."""
    use = cells[usable(cells) & (cells["stressor"] != "none")]
    return {key: {ch: group[f"{ch}_af"].to_numpy(float) for ch in CHANNELS}
            for key, group in use.groupby(["session_id", "condition_id", "sample_id"])}


def bootstrap_peaks(cells: pd.DataFrame, curve: pd.DataFrame, fields: pd.DataFrame,
                    n_boot: int = 500, seed: int = 0) -> pd.DataFrame:
    """Peak time per channel for one curve, with a resampling interval.

    Each replicate resamples cells within every field and fields within every
    time point, recomputes the curve and its peak. Reported: the peak of the
    observed curve, the 5th-95th percentile of the replicate peaks, and how
    often the replicate maximum fell on the observed sampled maximum.
    """
    rng = np.random.default_rng(seed)
    arrays = _cell_arrays(cells)
    points = []
    for _, point in curve.iterrows():
        here = fields[(fields["session_id"] == point["session_id"])
                      & (fields["condition_id"] == point["condition_id"])]
        keys = [(point["session_id"], point["condition_id"], s) for s in here["sample_id"]]
        points.append((point["time_h"], [k for k in keys if k in arrays]))

    # The ratio is the timer readout and rises with time by design, so a
    # "peak" of it is not asked for; peaks are reported for the two channels.
    observed = {"green": curve["green"].to_numpy(), "red": curve["red"].to_numpy()}
    times = curve["time_h"].to_numpy(float)
    draws = {name: [] for name in observed}
    for _ in range(n_boot):
        green, red = [], []
        for _, samples in points:
            if not samples:
                green.append(np.nan)
                red.append(np.nan)
                continue
            chosen = [samples[j] for j in rng.integers(0, len(samples), len(samples))]
            g, r = [], []
            for s in chosen:
                n = len(arrays[s]["green"])
                idx = rng.integers(0, n, n)
                g.append(np.median(arrays[s]["green"][idx]))
                r.append(np.median(arrays[s]["red"][idx]))
            green.append(np.median(g))
            red.append(np.median(r))
        green, red = np.array(green), np.array(red)
        replicate = {"green": green, "red": red}
        for name, values in replicate.items():
            peak = refine_peak(times, values)
            draws[name].append((peak["t_peak"], peak["t_max_sampled"]))

    rows = []
    for name, values in observed.items():
        peak = refine_peak(times, values)
        reps = np.array(draws[name], float)
        finite = np.isfinite(reps[:, 0])
        rows.append({
            "channel": name, **peak,
            "t_peak_p05": float(np.percentile(reps[finite, 0], 5)) if finite.any() else np.nan,
            "t_peak_p95": float(np.percentile(reps[finite, 0], 95)) if finite.any() else np.nan,
            "share_at_sampled_max": float(np.mean(reps[finite, 1] == peak["t_max_sampled"])) if finite.any() else np.nan,
            "rise_from_first": float(np.nanmax(values) - values[np.argmin(times)]),
            "n_timepoints": int(np.isfinite(values).sum()),
            "t_first": float(np.nanmin(times)), "t_last": float(np.nanmax(times)),
        })
    return pd.DataFrame(rows)


def excess_over_control(peaks: pd.DataFrame, curves: pd.DataFrame) -> pd.Series:
    """Value at the sampled peak minus the same day's 0 µM curve at that time.

    A peak that the unstressed culture shows too (growth, time of day, the
    slide) is not a stress response. The control is interpolated linearly
    between its own samples and left empty outside them.
    """
    out = []
    for _, peak in peaks.iterrows():
        control = curves[(curves["session_id"] == peak["session_id"]) & (curves["dose_uM"] == 0)]
        own = curves[(curves["session_id"] == peak["session_id"]) & (curves["dose_uM"] == peak["dose_uM"])]
        at = own.loc[own["time_h"] == peak["t_max_sampled"], peak["channel"]]
        if peak["dose_uM"] == 0 or control.empty or at.empty:
            out.append(np.nan)
            continue
        control = control.sort_values("time_h")
        t, v = control["time_h"].to_numpy(float), control[peak["channel"]].to_numpy(float)
        inside = t.min() <= peak["t_max_sampled"] <= t.max()
        out.append(float(at.iloc[0] - np.interp(peak["t_max_sampled"], t, v)) if inside else np.nan)
    return pd.Series(out, index=peaks.index)
