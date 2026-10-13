"""Read what the scripts wrote, for display. Nothing here computes a result."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

# verdict -> (colour, what the user should do); order = how strongly a dose is supported
VERDICTS = {
    "dose estimate": ("#1baf7a", "Report the most probable dose with its 80 % credible set."),
    "not validated": ("#eb6834", "The posterior points somewhere, but this decoder did no better than chance "
                                 "on held-out reference days. Do not report a dose."),
    "no dose information": ("#8a8a8a", "The image cannot tell the doses apart. Do not report a dose."),
    "unlike reference": ("#c23b3b", "Unlike any reference image: check reporter, imaging settings and "
                                    "stressor, or extend the reference."),
    "needs a 0 µM image": ("#2a78d6", "Image an unstressed culture on the same day, or use the absolute decoder."),
    "too few cells": ("#2a78d6", "Image more cells or a denser field."),
    # time since the stress
    "time estimate": ("#1baf7a", "Report the most probable time window with its 80 % credible set."),
    "no time information": ("#8a8a8a", "The timer cannot tell the time windows apart here. Do not report a time."),
    "too few reporter cells": ("#2a78d6", "Fewer than 5 reporter-positive cells: the timer cannot be read. "
                                          "Image more cells."),
    "unstressed control": ("#b0b0b0", "A 0 µM culture has no time since stress."),
    "time not available": ("#8a8a8a", "The reference cannot train a time decoder (see the reason above)."),
}


def read_json(path: Path) -> dict:
    try:
        return json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def sessions(cells: Path) -> pd.DataFrame:
    """Session ids of a cell table with stressor and strain parsed where the name allows."""
    ids = pd.read_csv(cells, usecols=["session_id"], dtype=str)["session_id"].dropna().unique()
    parts = pd.Series(sorted(ids)).str.extract(r"^(?P<day>\d{8})_(?P<stressor>[a-z0-9]+)_(?P<strain>tr|nt)$")
    return pd.concat([pd.Series(sorted(ids), name="session_id"), parts], axis=1)


def decode_result(out_dir: Path) -> tuple[dict, pd.DataFrame]:
    out_dir = Path(out_dir)
    decoded = out_dir / "decoded.csv"
    return read_json(out_dir / "summary.json"), (pd.read_csv(decoded) if decoded.exists() else pd.DataFrame())


def track_record_text(summary: dict) -> str:
    acc, gain, chance = (summary.get(k) for k in ("validated_accuracy", "validated_info_gain_bits", "chance"))
    days = len(summary.get("reference_days", []))
    if acc is None or pd.isna(acc):
        return f"Track record unavailable: {days} reference day(s); at least 2 are needed to hold one out."
    verdict = "better than chance" if gain is not None and gain > 0 else "no better than chance"
    return (f"On held-out reference days ({days} day(s)) this decoder scored accuracy {acc:.2f} against "
            f"chance {chance:.2f} and {gain:+.2f} bits per field: {verdict}.")


def time_record_text(summary: dict) -> str:
    if not summary.get("time_available", False):
        reason = summary.get("time_reason") or "not run"
        return f"Time since the stress cannot be decoded with this reference: {reason}."
    acc, gain, chance = (summary.get(k) for k in ("time_validated_accuracy", "time_validated_info_gain_bits",
                                                   "time_chance"))
    windows = ", ".join(summary.get("time_windows_in_reference", []))
    if acc is None or pd.isna(acc):
        return (f"Time windows in the reference: {windows}. Track record unavailable: no held-out reference "
                f"field falls in a window the other days cover.")
    verdict = "better than chance" if gain is not None and gain > 0 else "no better than chance"
    return (f"Time windows in the reference: {windows}. On held-out reference days the time decoder scored "
            f"accuracy {acc:.2f} against chance {chance:.2f} and {gain:+.2f} bits per field "
            f"({summary.get('time_scored_fields')} scored, {summary.get('time_unscorable_fields')} unscorable): {verdict}.")


def analysis_figures(out_dir: Path) -> list[tuple[str, Path]]:
    out_dir = Path(out_dir)
    found = [(f"Time course — {p.stem.split('_', 1)[1].upper()}", p)
             for p in sorted((out_dir / "timecourse").glob("timecourse_*.png"))]
    for title, rel in (("Model against data", "model_compare/compare.png"),
                       ("Single-field decoder, leave-one-day-out", "decoder_lodo/posterior.png")):
        if (out_dir / rel).exists():
            found.append((title, out_dir / rel))
    return found


def analysis_tables(out_dir: Path) -> dict[str, pd.DataFrame]:
    out_dir = Path(out_dir)
    tables = {}
    effects = out_dir / "induction" / "effects.csv"
    if effects.exists():
        e = pd.read_csv(effects)
        tables["Detected effects (q < 0.05)"] = e[e["detected_bh"]][
            ["day", "stressor", "dose_uM", "readout", "effect", "ci_low", "ci_high", "mde", "q_bh"]]
        tables["All comparisons"] = e
    trend = out_dir / "induction" / "trend.csv"
    if trend.exists():
        tables["Dose trend"] = pd.read_csv(trend)
    metrics = out_dir / "decoder_lodo" / "metrics.csv"
    if metrics.exists():
        m = pd.read_csv(metrics)
        tables["Decoder, mean over held-out days"] = m[m["held_out"] == "mean over days"]
    peaks = out_dir / "timecourse" / "peaks.csv"
    if peaks.exists():
        tables["Peak times"] = pd.read_csv(peaks)
    return tables
