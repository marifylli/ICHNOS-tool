"""Check an image manifest against docs/data_contract.md before hours of processing start.

Returns problems as rows (level, row, column, message):
  error    the pipeline would refuse or misread this; fix before running
  warning  the pipeline runs, but an analysis step will be limited or refuse later
Row numbers are 1-based data rows, as a spreadsheet shows them below the header.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from ichnos_image import instrument

REQUIRED = ("session_id", "timepoint", "acquisition_order", "green_path", "red_path", "bright_field_path",
            "exposure_ms_green", "exposure_ms_red", "nd_filter_green", "nd_filter_red", "objective",
            "burner_hours", "lamp_warmup_minutes")
FOR_ANALYSIS = ("sample_id", "condition_id", "measurement_time_hours")
SESSION = re.compile(r"^(?P<day>\d{8})_(?P<stressor>[a-z0-9]+)_(?P<strain>tr|nt)$")
MAX_LISTED = 20


def _issue(level: str, row, column: str, message: str) -> dict:
    return {"level": level, "row": row, "column": column, "message": message}


def check_manifest(path: Path) -> tuple[pd.DataFrame | None, list[dict], dict]:
    """(manifest, problems, summary). The manifest is None if the file cannot be read."""
    path = Path(path)
    try:
        table = pd.read_csv(path, dtype=str, keep_default_na=False)
    except (OSError, pd.errors.ParserError, UnicodeDecodeError, pd.errors.EmptyDataError) as exc:
        return None, [_issue("error", "", "", f"cannot read {path}: {exc}")], {}
    issues = []
    if table.empty:
        issues.append(_issue("error", "", "", "the manifest has no rows"))
    missing = [c for c in REQUIRED if c not in table.columns]
    for column in missing:
        issues.append(_issue("error", "", column, "required column missing (it may be empty, but must exist)"))
    for column in FOR_ANALYSIS:
        if column not in table.columns:
            issues.append(_issue("warning", "", column, "missing: the pipeline runs, the analysis cannot"))

    def rows(mask) -> list[int]:
        return [int(i) + 1 for i in table.index[mask]][:MAX_LISTED]

    # numbers
    for column in ("exposure_ms_green", "exposure_ms_red"):
        if column in table:
            values = pd.to_numeric(table[column], errors="coerce")
            for r in rows(~(values > 0)):
                issues.append(_issue("error", r, column, f"exposure must be a positive number, got {table.at[r - 1, column]!r}"))
    if "measurement_time_hours" in table:
        raw = table["measurement_time_hours"].str.strip()
        values = pd.to_numeric(raw, errors="coerce")
        for r in rows((raw != "") & ~(values >= 0)):
            issues.append(_issue("error", r, "measurement_time_hours",
                                 f"hours since the stress was added must be ≥ 0 or blank, got {table.at[r - 1, 'measurement_time_hours']!r}"))
        for r in rows(raw == ""):
            issues.append(_issue("warning", r, "measurement_time_hours", "blank: this field has no time axis"))
    if "objective" in table:
        known = set(instrument.OBJECTIVES)
        for r in rows(~table["objective"].str.strip().str.upper().isin(known)):
            issues.append(_issue("error", r, "objective",
                                 f"unknown objective {table.at[r - 1, 'objective']!r}; known: {', '.join(sorted(known))}"))

    # files, relative to the manifest
    base = path.parent
    for column, needed in (("green_path", True), ("red_path", True), ("bright_field_path", False)):
        if column not in table:
            continue
        listed = 0
        for i, value in table[column].items():
            value = value.strip()
            if not value:
                if needed:
                    issues.append(_issue("error", int(i) + 1, column, "empty"))
                continue
            if not (base / value).is_file() and listed < MAX_LISTED:
                issues.append(_issue("error", int(i) + 1, column, f"file not found: {base / value}"))
                listed += 1

    # identifiers
    if {"session_id", "condition_id", "sample_id"} <= set(table.columns):
        keys = table[["session_id", "condition_id", "sample_id"]]
        for r in rows(keys.duplicated(keep=False)):
            issues.append(_issue("error", r, "sample_id", "same session, condition and sample as another row"))
    sessions = {}
    if "session_id" in table:
        for session in sorted(table["session_id"].unique()):
            match = SESSION.match(session)
            if not match:
                issues.append(_issue("warning", "", "session_id",
                                     f"{session!r} is not <yyyymmdd>_<stressor>_<tr|nt>: fine for a new culture to "
                                     f"decode, not for a reference experiment"))
                continue
            sessions[session] = match.groupdict()
    if sessions and "condition_id" in table:
        for session, info in sessions.items():
            conditions = table.loc[table["session_id"] == session, "condition_id"]
            pattern = re.compile(rf"^{re.escape(info['stressor'])}_(\d+(?:\.\d+)?)uM(?:_|$)")
            doses = [float(m.group(1)) for c in conditions if (m := pattern.match(c))]
            unknown = int(sum(1 for c in conditions if not pattern.match(c)))
            info.update(fields=int(len(conditions)), doses=sorted(set(doses)), unknown_dose_fields=unknown)
            if unknown:
                issues.append(_issue("warning", "", "condition_id",
                                     f"{session}: {unknown} field(s) without '<stressor>_<dose>uM' — dose unknown"))
            if info["strain"] == "tr" and 0.0 not in doses:
                issues.append(_issue("warning", "", "condition_id",
                                     f"{session}: no 0 µM culture on this day; within-day comparisons need one"))
        stressors = {i["stressor"] for i in sessions.values()}
        for stressor in sorted(stressors):
            if not any(i["stressor"] == stressor and i["strain"] == "nt" for i in sessions.values()):
                issues.append(_issue("warning", "", "session_id",
                                     f"no non-transfected {stressor} session: the autofluorescence baseline "
                                     f"must then come from a reference experiment"))
    summary = {"rows": int(len(table)), "sessions": sessions,
               "errors": sum(i["level"] == "error" for i in issues),
               "warnings": sum(i["level"] == "warning" for i in issues)}
    return table, issues, summary
