"""Exposure normalisation of per-cell intensities, and the camera-linearity check it relies on.

Exposure changed between samples (Exp 150–600, i.e. 600–2400 ms acquisition),
so corrected intensities in grey levels are not comparable across samples.
Dividing by the exposure time gives grey levels per second, which is valid only
if the camera responds linearly to exposure once its dark offset is removed.

The red/green ratio needs no normalisation: the log records the same exposure
for both channels of a sample, and ``add_per_second_columns`` checks that.

Linearity is checked here with what the October 2026 sessions contain: the
cell-free background of each field against its exposure. That background is
mostly medium autofluorescence plus the camera's dark offset, so within a
session it should rise along a straight line ``background = offset + rate * t``.
This is a necessary condition only. The direct test is one field imaged at
several exposures, which these sessions do not have.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .correct import _histogram_mode_background
from .image_io import load_image

INTENSITY_COLUMNS = ("corrected_mean_green", "corrected_mean_red", "integrated_green", "integrated_red")


def add_per_second_columns(cells: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with ``<column>_per_s`` for every corrected intensity column.

    Raises if a row has a missing or non-positive exposure, or different
    exposures for the two channels (the ratio would then carry the exposure).
    """
    for name in ("exposure_ms_green", "exposure_ms_red"):
        values = pd.to_numeric(cells[name], errors="coerce")
        if values.isna().any() or (values <= 0).any():
            raise ValueError(f"{name} must be a positive number on every row")
    if not np.allclose(cells.exposure_ms_green, cells.exposure_ms_red):
        raise ValueError("green and red exposures differ on some rows; the ratio is not exposure-free")
    out = cells.copy()
    seconds = out.exposure_ms_green.astype(float) / 1000.0
    for name in INTENSITY_COLUMNS:
        out[f"{name}_per_s"] = out[name].astype(float) / seconds
    return out


def field_background(path: str | Path) -> dict:
    """Background statistics of one fluorescence frame, in raw grey levels."""
    image = load_image(path)
    if image.ndim == 3:
        raise ValueError(f"expected a single-channel frame, got shape {image.shape}: {path}")
    return {
        "mode": _histogram_mode_background(image),
        "p01": float(np.percentile(image, 1)),
        "median": float(np.median(image)),
        "saturated_fraction": float((image >= 255).mean()),
    }


def fit_background_vs_exposure(fields: pd.DataFrame, value: str, group: str = "session_id") -> pd.DataFrame:
    """Least-squares line value = offset + rate * exposure_s, per group and pooled.

    Groups with fewer than two distinct exposures cannot be fitted and are
    reported with ``fit = "single_exposure"`` instead of being dropped.
    """
    rows = []
    groups = [(name, part) for name, part in fields.groupby(group)] + [("all", fields)]
    for name, part in groups:
        t = part.exposure_ms.astype(float).to_numpy() / 1000.0
        y = part[value].astype(float).to_numpy()
        record = {group: name, "value": value, "n_fields": len(part),
                  "exposures_ms": " ".join(str(int(e)) for e in sorted(set(part.exposure_ms)))}
        if len(set(t)) < 2:
            record.update(fit="single_exposure", mean=float(y.mean()))
        else:
            rate, offset = np.polyfit(t, y, 1)
            predicted = offset + rate * t
            ss_res = float(((y - predicted) ** 2).sum())
            ss_tot = float(((y - y.mean()) ** 2).sum())
            record.update(fit="linear", offset=float(offset), rate_per_s=float(rate),
                          r2=1 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
                          rmse=float(np.sqrt(ss_res / len(y))))
        rows.append(record)
    return pd.DataFrame(rows)
