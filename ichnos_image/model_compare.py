"""The model's predicted induction against the measured one, as fold over the same day's 0 µM.

Both sides are put on the one scale they share without calibration: how much
a stressed culture exceeds its own 0 µM control, averaged over the times the
stressed culture was imaged,

    fold = 1 + mean_t(stressed(t) - control(t)) / mean_t(control(t))

with control(t) the same day's 0 µM curve interpolated to t (the definition
``induction.induction_table`` uses for its effect). Model units (nM) and
camera units (grey/s) cancel in the ratio, so no fluorescence calibration
is needed, but neither is one tested.

Data side: reporter-positive cells only (green and red above the
non-transfected 90th percentile, as in induction.py), because the model
describes a cell that carries the reporter and the population median is
diluted by cells that do not. Selecting on brightness compresses the fold
towards 1 (cells are kept only above a floor), so a measured fold is a lower
bound on the fold of the expressing cells.

Model side: the default parameter profile, i.e. uncalibrated, at the doses
and the imaging times of each day. Stress is constant unless a first-order
clearance rate is given; H2O2 is degraded by cells within hours, so
clearance rates are assumptions to compare, not fitted values.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import induction
from . import timecourse as tc

STRESSOR_VARIANT = {"h2o2": "ox", "dtt": "er"}
READOUTS = ("green", "red")


def reporter_positive(cells: pd.DataFrame, quantile: float = 0.9) -> pd.DataFrame:
    """Non-transfected cells as they are; transfected cells only if positive in both channels."""
    limits = induction.positive_thresholds(cells, quantile)
    tr = cells["strain"] == "tr"
    known = cells["stressor"].isin(limits.index)
    green = cells["stressor"].map(limits["green"]) if len(limits) else np.nan
    red = cells["stressor"].map(limits["red"]) if len(limits) else np.nan
    positive = (known & (cells["corrected_mean_green_per_s"] > green)
                & (cells["corrected_mean_red_per_s"] > red))
    return cells[~tr | positive].copy()


def imaging_times(cells: pd.DataFrame) -> dict:
    """Per transfected session and dose > 0: the imaging times inside the 0 µM time range.

    These are the times ``induction`` averages over, so the model is sampled at them too.
    """
    use = cells[tc.usable(cells) & (cells["strain"] == "tr") & (cells["stressor"] != "none")]
    points = (use.groupby(["session_id", "dose_uM", "condition_id"])["measurement_time_hours"]
              .median().reset_index())
    out = {}
    for session, day in points.groupby("session_id"):
        control = day[day["dose_uM"] == 0]["measurement_time_hours"]
        if len(control) < 2:
            continue
        lo, hi = control.min(), control.max()
        for dose, sample in day[day["dose_uM"] > 0].groupby("dose_uM"):
            t = sample["measurement_time_hours"]
            t = sorted(set(np.round(t[(t >= lo) & (t <= hi)], 3)))
            if t:
                out[(session, float(dose))] = t
    return out


def data_fold(cells: pd.DataFrame, n_boot: int = 1000, seed: int = 0, quantile: float = 0.9) -> pd.DataFrame:
    """Fold over the same day's 0 µM of reporter-positive cells, per day, dose and channel, with 95 % interval.

    ``cells`` must be annotated and carry ``green_af``/``red_af``.
    """
    effects, _ = induction.induction_table(reporter_positive(cells, quantile), n_boot=n_boot, seed=seed,
                                           quantile=quantile)
    effects = effects[effects["readout"].isin(READOUTS)].copy()
    level = effects["control_level"]
    if (level <= 0).any():
        bad = effects.loc[level <= 0, ["day", "stressor", "readout"]].drop_duplicates().to_dict("records")
        raise ValueError(f"0 µM level of reporter-positive cells is not above the baseline for {bad}")
    effects["fold"] = 1 + effects["effect"] / level
    effects["fold_ci_low"] = 1 + effects["ci_low"] / level
    effects["fold_ci_high"] = 1 + effects["ci_high"] / level
    effects["fold_mde"] = effects["mde"] / level  # smallest detectable change in fold, e.g. 0.3 = 1.3x
    return effects[["session_id", "day", "stressor", "dose_uM", "readout", "n_timepoints", "control_level",
                    "fold", "fold_ci_low", "fold_ci_high", "fold_mde", "q_bh"]].reset_index(drop=True)


def _run_model(variant: str, dose: float, times: list, clearance: float | None, profile: str) -> dict:
    import contextlib
    import io

    import libsbml

    from ichnos import build, naming, params, protocol
    from ichnos.fluorescence import model_fluorescence

    with contextlib.redirect_stdout(io.StringIO()):  # the builder reports every merge step
        sbml = build.build_variant_sbml_string(variant, save_sbml=False)
    document = libsbml.readSBMLFromString(sbml)
    params.apply_profile(document.getModel(), params.load_profile(variant, profile))
    stress = protocol.StressProtocol(variant=variant, dose=float(dose), dose_units="uM",
                                     clears=clearance is not None, clearance_rate_per_hour=clearance)
    runner, model = protocol.load_protocol_model(libsbml.writeSBMLToString(document), stress)
    result, _ = protocol.run_protocol_at_times(runner, stress, times_hours=np.asarray(times, float),
                                               id_to_name=naming.id_to_name_map(model))
    fluorescence = model_fluorescence(result)
    return {ch: np.asarray(fluorescence[ch], float) for ch in READOUTS}


def model_fold(times: dict, stressor_of: dict, clearances=(None,), profile: str = "default",
               run=_run_model) -> pd.DataFrame:
    """Model fold over 0 µM at each (session, dose) and its imaging times, per clearance scenario.

    times: {(session_id, dose_uM): [hours]} (``imaging_times``); stressor_of: {session_id: stressor}.
    ``run`` is the simulator (variant, dose, times, clearance, profile) -> {"green": array, "red": array}.
    """
    rows = []
    for (session, dose), t in sorted(times.items()):
        variant = STRESSOR_VARIANT.get(stressor_of[session])
        if variant is None:
            continue
        for clearance in clearances:
            stressed = run(variant, dose, t, clearance, profile)
            control = run(variant, 0.0, t, clearance, profile)
            for ch in READOUTS:
                fold = 1 + float(np.mean(stressed[ch] - control[ch])) / float(np.mean(control[ch]))
                rows.append({"session_id": session, "stressor": stressor_of[session], "variant": variant,
                             "dose_uM": dose, "readout": ch,
                             "clearance_per_h": np.nan if clearance is None else float(clearance),
                             "times_h": " ".join(f"{x:g}" for x in t), "model_fold": fold})
    return pd.DataFrame(rows)


def compare(data: pd.DataFrame, model: pd.DataFrame) -> pd.DataFrame:
    """One row per day, dose, channel and clearance scenario; is the model fold inside the measured interval?"""
    table = data.merge(model.drop(columns="stressor"), on=["session_id", "dose_uM", "readout"], how="inner")
    table["model_inside_ci"] = (table["model_fold"] >= table["fold_ci_low"]) & (table["model_fold"] <= table["fold_ci_high"])
    table["model_over_data"] = table["model_fold"] / table["fold"]
    return table.sort_values(["stressor", "day", "readout", "clearance_per_h", "dose_uM"],
                             na_position="first").reset_index(drop=True)
