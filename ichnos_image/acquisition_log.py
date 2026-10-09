"""Parse the wet-lab microscopy log into one row per imaged field.

The wet lab keeps the acquisition log as a Word document ("Ημερολόγιο
μετρήσεων μικροσκοπίας — H2O2 & DTT"). It is exported to PDF and converted
with ``pdftotext -layout``; this module reads that text.

Sections are found by their headings, not by line number, so a new
paragraph in the document does not shift the parse. Every value that is not
read verbatim from a table row is marked in ``flags``:

- ``acq_derived_4x``: the log gives acq exposure only for Exp 300-600;
  for other codes acq = 4 x Exp, the relation that holds exactly for every
  listed pair.
- ``exp_from_section_note``: the table has no Exp column; the value comes
  from the section text (for example "Exp = 200 σε όλες τις λήψεις").
- ``starred_in_log:...``: the wet lab marked that time with an asterisk.

Nothing is filled in silently. A section whose stated addition time cannot
be found in its own text is reported in ``ParseResult.warnings``.

Times are wall-clock times as written in the log. The camera computer's
clock runs about 50 minutes ahead of them; that offset is applied when
images are matched to rows, not here.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

# Exposure table stated in the log: Exp code (middle three digits of pvw) -> acq ms.
ACQ_MS_BY_EXP = {300: 1200, 350: 1400, 375: 1500, 400: 1600, 500: 2000, 600: 2400}

COLUMNS = [
    "section", "date", "stressor", "strain", "tp_label", "dose_uM", "dose_units", "rep",
    "added", "plating", "bf", "red", "green", "exp", "exposure_acq_ms",
    "sampling_time_hours", "measurement_time_hours", "folder", "flags",
]

_T = r"(\d{1,2}:\d{2})(\*?)"


@dataclass(frozen=True)
class Section:
    """One table of the log and the context its rows inherit."""

    id: str
    heading: str  # regex matched against a whole, whitespace-normalised line
    kind: str
    date: str = ""
    stressor: str = ""
    strain: str = ""
    added: str = ""
    plating: str = ""
    tp_label: str = ""
    folder: str = ""
    exp_default: int | None = None


# Order matters: each section runs until the next heading in this list.
SECTIONS = [
    Section("medium", r"Control σε θρεπτικό", "medium"),
    Section("h2o2_series", r"Οξειδωτικό στρες \(H2O2\) — χρονοσειρά.*", "series_dated",
            stressor="H2O2", strain="transfected"),
    # Heading paragraphs that state the day's addition time; they hold no table rows.
    Section("h2o2_0710_note", r"Οξειδωτικό στρες \(H2O2\) — τριπλέτες 7/10", "note", added="11:00"),
    Section("h2o2_0710_0min", r"Μέτρηση «0 min».*", "dose_only", "2026-10-07", "H2O2", "transfected",
            added="11:00", tp_label="0min", folder="C2:igem h2o2 Vol 2/{dose}/1h"),
    Section("h2o2_0710_s1245", r"Στρώση 12:45", "rep", "2026-10-07", "H2O2", "transfected",
            added="11:00", plating="12:45", tp_label="strosi_12:45", folder="C2:igem h2o2 Vol 2/{dose}/2h"),
    Section("h2o2_0710_s1429", r"Στρώση 14:29", "rep", "2026-10-07", "H2O2", "transfected",
            added="11:00", plating="14:29", tp_label="strosi_14:29", folder="C2:igem h2o2 Vol 2/{dose}/3h"),
    Section("h2o2_0710_s1704", r"Στρώση 17:04", "rep", "2026-10-07", "H2O2", "transfected",
            added="11:00", plating="17:04", tp_label="strosi_17:04", folder="C2:igem h2o2 Vol 2/{dose}/6h"),
    Section("nt_h2o2_0810", r"Μη μετασχηματισμένα κύτταρα \+ H2O2.*", "nt_series", "2026-10-08", "H2O2",
            "non_transfected", added="12:44", plating="12:44", exp_default=200),
    Section("dtt_0610", r"ER στρες \(DTT\) — χρονοσειρά.*", "dtt_series", "2026-10-06", "DTT", "transfected",
            added="12:49"),
    Section("dtt_0710_note", r"ER στρες \(DTT\) — τριπλέτες 7/10", "note", added="12:50"),
    Section("dtt_0710_s1335", r"Στρώση 13:35.*", "rep", "2026-10-07", "DTT", "transfected",
            added="12:50", plating="13:35", tp_label="45min", folder="C1:IGEM DTT VOL2/{dose}/45'"),
    Section("dtt_0710_s1520", r"Στρώση 15:20.*", "rep", "2026-10-07", "DTT", "transfected",
            added="12:50", plating="15:20", tp_label="2h(actual_2.5h)", folder="C1:IGEM DTT VOL2/{dose}/2h"),
    Section("dtt_0710_s1650", r"Στρώση 16:50", "rep", "2026-10-07", "DTT", "transfected",
            added="12:50", plating="16:50", tp_label="4h", folder="C1:IGEM DTT VOL2/{dose}/4h", exp_default=200),
    # "Mη" in the document starts with a Latin M, so the heading is matched from the second word.
    Section("nt_dtt_0910_0h", r".η μετασχηματισμένα \(DTT\).*", "rep", "2026-10-09", "DTT", "non_transfected",
            added="12:20", tp_label="0h", folder="D:igem DTT non transfected/{dose}/0h", exp_default=150),
    Section("nt_dtt_0910_45", r"45λεπτά στρώση DTT", "rep", "2026-10-09", "DTT", "non_transfected",
            added="12:20", tp_label="45min", folder="D:igem DTT non transfected/{dose}/45'", exp_default=150),
    Section("nt_dtt_0910_2h", r"2 ΩΡΕΣ", "rep", "2026-10-09", "DTT", "non_transfected",
            added="12:20", tp_label="2h", folder="D:igem DTT non transfected/{dose}/2h", exp_default=150),
    Section("nt_dtt_0910_4h", r"4 ΩΡΕΣ DTT", "rep", "2026-10-09", "DTT", "non_transfected",
            added="12:20", tp_label="4h", folder="D:igem DTT non transfected/{dose}/4h", exp_default=150),
]

# Addition times per day for the dated H2O2 time series (5/10 and 6/10 share one table).
SERIES_ADDED = {"5": "13:41", "6": "11:55"}
SERIES_FOLDER = {"0h": "A1:igem h2o2/{dose}/0", "30 min": "A1:igem h2o2/{dose}/30", "1h": "A1:igem h2o2/{dose}/1h",
                 "2h": "A1:igem h2o2/{dose}/2h", "3h": "A1:igem h2o2/{dose}/3h",
                 "6h": "A2:igem h2o2 palia/{dose}/6h", "7h": "A2:igem h2o2 palia/{dose}/7h"}
DTT06_FOLDER = {"0 min": "0", "45 min": "45", "1h": "2h", "4h": "4h", "5h": "5h", "6h": "6h"}
NT_H2O2_FOLDER = {"0h": "0h", "30 min": "30'", "1h": "1Η"}

# Images saved in the wrong folder. The log is authoritative (confirmed by the team on
# 2026-10-10): at 2 h the order was 100 -> 50 -> 200 uM. Camera timestamps put the
# 100 uM fields in the 200/2h folder (files misnamed "50-2Η") and the 200 uM fields in
# 100/2h. Raw folders are left untouched; the mapping is applied here.
# (section, dose_uM) -> (folder holding that dose's images, flag)
FOLDER_OVERRIDES = {
    ("nt_dtt_0910_2h", "100"): ("D:igem DTT non transfected/200/2h",
                                "images_in_200/2h_folder(files_named_50-2Η);log_authoritative"),
    ("nt_dtt_0910_2h", "200"): ("D:igem DTT non transfected/100/2h",
                                "images_in_100/2h_folder;log_authoritative"),
}


@dataclass
class ParseResult:
    rows: list[dict] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def _norm(line: str) -> str:
    return " ".join(line.split())


def _hours(start: str, end: str) -> str:
    if not start or not end:
        return ""
    delta = datetime.strptime(end, "%H:%M") - datetime.strptime(start, "%H:%M")
    return f"{delta.total_seconds() / 3600:.3f}"


def _acq_ms(exp: int | None, flags: list[str]) -> str:
    if exp is None:
        return ""
    if exp in ACQ_MS_BY_EXP:
        return str(ACQ_MS_BY_EXP[exp])
    flags.append("acq_derived_4x")
    return str(4 * exp)


def split_sections(lines: list[str]) -> tuple[dict[str, list[str]], list[str]]:
    """Return section id -> its lines, and warnings for headings not found."""
    starts: list[tuple[int, Section]] = []
    warnings = []
    cursor = 0
    for section in SECTIONS:
        pattern = re.compile(rf"^{section.heading}$")
        for index in range(cursor, len(lines)):
            if pattern.match(_norm(lines[index])):
                starts.append((index, section))
                cursor = index + 1
                break
        else:
            warnings.append(f"heading not found: {section.id}")
    chunks = {}
    for position, (index, section) in enumerate(starts):
        end = starts[position + 1][0] if position + 1 < len(starts) else len(lines)
        chunks[section.id] = lines[index:end]
    return chunks, warnings


def _check_added(section: Section, chunk: list[str], warnings: list[str]) -> None:
    text = " ".join(_norm(line) for line in chunk)
    stated = set(re.findall(r"(?:μπ[ηή]κε|προστέθηκε)\D{0,25}?(\d{1,2}:\d{2})", text))
    expected = set(SERIES_ADDED.values()) if section.kind == "series_dated" else {section.added}
    expected.discard("")
    missing = expected - stated
    if missing and section.kind not in ("rep", "dose_only"):  # these inherit the time stated in a note above them
        warnings.append(f"{section.id}: addition time {sorted(missing)} not stated in section text (found {sorted(stated)})")


def _parse_line(section: Section, text: str, next_text: str) -> tuple[dict, list[str]] | None:
    tag = ""
    if "(cluster)" in text or next_text == "(cluster)":
        tag = "cluster"
    elif "(normal)" in text or next_text == "(normal)":
        tag = "normal"
    text = _norm(text.replace("(cluster)", "").replace("(normal)", ""))
    # pdftotext writes "0h", pypdf "0 h" for the same cell.
    text = re.sub(r"^(\d+) h\b", r"\1h", text)
    flags: list[str] = []
    row: dict = {}
    kind = section.kind
    if kind == "medium":
        m = re.fullmatch(rf"(\d+)/10 {_T} (\d+) {_T} {_T}", text)
        if not m:
            return None
        day = m[1]
        row = {"date": f"2026-10-{int(day):02d}", "stressor": "none",
                   "strain": "transfected" if day == "5" else "non_transfected", "tp_label": "medium",
                   "bf": m[2], "exp": m[4], "red": m[5], "green": m[7],
                   "folder": "A1:igem h2o2/θρεπτικό" if day == "5" else "B:igem h2o2 non transfected/ΘΡΕΠΤΙΚΌ"}
        stars = (m[3], m[6], m[8])
    elif kind == "series_dated":
        m = re.fullmatch(rf"(0h|30 min|1h|2h|3h|6h|7h) (\d+)/10 (\d+) {_T} (\d+) {_T} {_T}", text)
        if not m:
            return None
        day = m[2]
        row = {"date": f"2026-10-{int(day):02d}", "tp_label": m[1], "dose_uM": m[3], "rep": tag, "bf": m[4], "exp": m[6],
                   "red": m[7], "green": m[9], "added": SERIES_ADDED.get(day, ""),
                   "folder": SERIES_FOLDER[m[1]].format(dose=m[3])}
        stars = (m[5], m[8], m[10])
    elif kind == "dose_only":
        m = re.fullmatch(rf"(\d+) {_T} (\d+) {_T} {_T}", text)
        if not m:
            return None
        row = {"dose_uM": m[1], "rep": "single", "bf": m[2], "exp": m[4], "red": m[5], "green": m[7]}
        stars = (m[3], m[6], m[8])
    elif kind == "rep":
        m = re.fullmatch(rf"(\d+) (R\d) {_T}(?: (\d+))? {_T} {_T}", text)
        if not m:
            return None
        exp = m[5]
        if not exp:
            exp = str(section.exp_default or "")
            flags.append("exp_from_section_note")
        row = {"dose_uM": m[1], "rep": m[2], "bf": m[3], "exp": exp, "red": m[6], "green": m[8]}
        stars = (m[4], m[7], m[9])
    elif kind == "nt_series":
        m = re.fullmatch(rf"(0h|30 min|1h) (\d+) (R\d) {_T} {_T} {_T}", text)
        if not m:
            return None
        row = {"tp_label": m[1], "dose_uM": m[2], "rep": m[3], "bf": m[4], "exp": str(section.exp_default), "red": m[6],
                   "green": m[8], "folder": f"B:igem h2o2 non transfected/{m[2]}/{NT_H2O2_FOLDER[m[1]]}"}
        flags.append("exp_from_section_note")
        stars = (m[5], m[7], m[9])
    elif kind == "dtt_series":
        m = re.fullmatch(rf"(0 min|45 min|1h|4h|5h|6h) {_T} (\d+) {_T} (\d+) {_T} {_T}", text)
        if not m:
            return None
        row = {"tp_label": m[1], "plating": m[2], "dose_uM": m[4], "rep": tag, "bf": m[5], "exp": m[7], "red": m[8], "green": m[10],
                   "folder": f"A3:igem DTT/{m[4]}/{DTT06_FOLDER[m[1]]}"}
        stars = (m[6], m[9], m[11])
        if m[3]:
            flags.append("plating_time_starred")
        if m[1] == "6h" and m[2] == "19:49":
            # Plated "19:49" but imaged from 19:08; 18:49 is 6.0 h after the 12:49 addition.
            row["plating"] = "18:49"
            flags.append("plating_19:49_read_as_18:49")
    else:
        raise ValueError(f"unknown section kind {kind!r}")
    if any(stars):
        flags.append("starred_in_log:" + ",".join(n for n, s in zip(("bf", "red", "green"), stars) if s))
    return row, flags


def parse_log_text(text: str) -> ParseResult:
    """Parse ``pdftotext -layout`` output of the wet-lab log."""
    # pdftotext separates pages with form feeds; str.splitlines() would treat
    # them as extra line breaks, so split on newlines only.
    lines = text.split("\n")
    chunks, warnings = split_sections(lines)
    result = ParseResult(warnings=warnings)
    for section in SECTIONS:
        chunk = chunks.get(section.id)
        if chunk is None:
            continue
        _check_added(section, chunk, result.warnings)
        if section.kind == "note":
            continue
        for index, raw in enumerate(chunk):
            next_text = _norm(chunk[index + 1]) if index + 1 < len(chunk) else ""
            parsed = _parse_line(section, _norm(raw), next_text)
            if parsed is None:
                continue
            row, flags = parsed
            for key in ("date", "stressor", "strain", "added", "plating", "tp_label", "folder"):
                if not row.get(key):
                    row[key] = getattr(section, key)
            row["folder"] = row["folder"].format(dose=row.get("dose_uM", ""))
            exp = int(row["exp"]) if row.get("exp") else None
            row["exposure_acq_ms"] = _acq_ms(exp, flags)
            row["sampling_time_hours"] = _hours(row["added"], row["plating"])
            row["measurement_time_hours"] = _hours(row["added"], row["red"])
            row["section"] = section.id
            row["flags"] = flags
            result.rows.append(row)
    _flag_known_issues(result.rows)
    for row in result.rows:
        row["dose_units"] = "uM" if row.get("dose_uM") else ""
        row["flags"] = ";".join(row["flags"])
        for column in COLUMNS:
            row.setdefault(column, "")
    return result


def _flag_known_issues(rows: list[dict]) -> None:
    """Record problems of the log itself on the affected rows."""
    seen = set()
    for row in rows:
        key = (row["section"], row["tp_label"], row.get("dose_uM", ""), row.get("rep", ""))
        if row["section"] != "medium" and key in seen:
            row["flags"].append("duplicate_row_in_log")
        seen.add(key)
        if row["section"] == "h2o2_0710_0min":
            row["flags"].append("log_says_~30min_but_11:00_gives_~55min")
        if row["section"] == "h2o2_0710_s1704" and row.get("dose_uM") in ("75", "300", "600"):
            row["flags"].append("waited_~1h_longer_on_slide")
        if row["section"] == "nt_h2o2_0810":
            row["flags"].append("log_covers_0-1h_only")
        override = FOLDER_OVERRIDES.get((row["section"], row.get("dose_uM")))
        if override:
            row["folder"], flag = override
            row["flags"].append(flag)


def write_rows(rows: list[dict], path: Path) -> None:
    """Write rows to a new CSV; refuses to overwrite, like the rest of the tool."""
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"output path already exists: {path}")
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
