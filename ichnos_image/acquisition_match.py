"""Match microscope images to rows of the acquisition log by acquisition time.

Folder and file names on the acquisition computer are typed by hand and are
sometimes wrong (a 100 uM field saved in the 200 uM folder, a 75 uM file in
the 300 uM folder). The TIFF ``DateTime`` tag is written by the camera
software, so it is used as the link between an image and its log row.

The camera computer's clock is offset from the wall clock the log is kept
in (about +50 min in October 2026). The offset is estimated per day as the
median difference between an image's red timestamp and the red time of the
log row for the same folder, then every image is assigned to a log row of
the same day by minimum total time difference (red and green). Pairs that
differ by more than ``tolerance_min`` after the offset are not matched.

Nothing is decided silently: every image and every log row ends up in the
match report, matched or not, with the residual and the reason.
"""
from __future__ import annotations

import csv
import re
import unicodedata
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from statistics import median

import numpy as np
import tifffile
from scipy.optimize import linear_sum_assignment

ROOT_OF_FOLDER_CODE = {"A1": "A", "A2": "A", "A3": "A", "B": "B", "C1": "C", "C2": "C", "D": "D"}
_GREEN = "[Alexa 488]"
_RED = "[Alexa 594]"
# A log time HH:MM covers the whole minute; compare against its midpoint.
_HALF_MINUTE = timedelta(seconds=30)


@dataclass
class ImageField:
    """One imaged field: the red/green pair and the bright-field frames near it."""

    root: str
    red_path: Path
    green_path: Path
    folder: str  # condition folder including the root's name, e.g. "wetransfer_6-10/igem DTT/200/4h"
    red_time: datetime | None
    green_time: datetime | None
    bf_candidates: list[tuple[Path, datetime | None]] = field(default_factory=list)


def tiff_datetime(path: Path) -> datetime | None:
    """Read the TIFF DateTime tag; Image-Pro writes '10/06/2026 06:11:10.138 PM'."""
    with tifffile.TiffFile(path) as tif:
        tag = tif.pages[0].tags.get("DateTime")
        value = None if tag is None else str(tag.value).strip()
    if not value:
        return None
    for fmt in ("%m/%d/%Y %I:%M:%S.%f %p", "%m/%d/%Y %I:%M:%S %p", "%Y:%m:%d %H:%M:%S"):
        try:
            return datetime.strptime(value, fmt)
        except ValueError:
            continue
    return None


def scan_root(code: str, root: Path) -> list[ImageField]:
    """Find every red/green pair under ``root``.

    The layout is ``<condition>/<field folder>/01a/<name> [01a][Alexa 594].tif``
    with bright-field ``Image000*.tif`` frames in ``<condition>`` or in
    ``<field folder>``.
    """
    root = Path(root).expanduser()
    fields_found = []
    # Not rglob(f"*{_RED}.tif"): the brackets in "[Alexa 594]" are a glob character class.
    reds = sorted(p for p in root.rglob("*.tif") if p.name.endswith(f"{_RED}.tif"))
    for red in reds:
        green = red.with_name(red.name.replace(_RED, _GREEN))
        if not green.exists():
            continue
        field_dir = red.parent.parent
        condition_dir = field_dir.parent
        candidates = sorted({*condition_dir.glob("Image000*.tif"), *field_dir.glob("Image000*.tif")})
        fields_found.append(ImageField(
            root=code, red_path=red, green_path=green,
            folder=str(condition_dir.relative_to(root.parent)),
            red_time=tiff_datetime(red), green_time=tiff_datetime(green),
            bf_candidates=[(p, tiff_datetime(p)) for p in candidates],
        ))
    return fields_found


def _log_time(row: dict, column: str) -> datetime | None:
    value = row.get(column, "")
    if not value:
        return None
    return datetime.strptime(f"{row['date']} {value}", "%Y-%m-%d %H:%M") + _HALF_MINUTE


def _row_root_and_folder(row: dict) -> tuple[str, str]:
    code, _, folder = row["folder"].partition(":")
    return ROOT_OF_FOLDER_CODE.get(code, code), folder


def _norm_path(text: str) -> str:
    # macOS may store Greek names decomposed (NFD); the log is composed (NFC).
    return re.sub(r"\s+", " ", unicodedata.normalize("NFC", text)).strip().casefold()


def _same_folder(image_folder: str, row_folder: str) -> bool:
    """True when the image's condition folder is the one the log row names.

    ``image_folder`` includes the root folder's own name, so a log folder
    written either relative to the root or including it matches.
    """
    image, row = _norm_path(image_folder), _norm_path(row_folder)
    return image == row or image.endswith("/" + row)


def estimate_offsets(images: list[ImageField], rows: list[dict]) -> dict[str, float]:
    """Per-day camera-clock offset in minutes (camera minus log), from same-folder pairs."""
    by_day: dict[str, list[float]] = {}
    for image in images:
        if image.red_time is None:
            continue
        day = image.red_time.strftime("%Y-%m-%d")
        diffs = []
        for row in rows:
            root, folder = _row_root_and_folder(row)
            if row["date"] != day or root != image.root or not _same_folder(image.folder, folder):
                continue
            log_red = _log_time(row, "red")
            if log_red is not None:
                diffs.append((image.red_time - log_red).total_seconds() / 60)
        if diffs:
            # Several rows share a folder (R1-R3); the nearest one is this image's.
            by_day.setdefault(day, []).append(min(diffs, key=abs))
    return {day: median(values) for day, values in by_day.items()}


def match(images: list[ImageField], rows: list[dict], offsets: dict[str, float],
          tolerance_min: float = 2.5) -> list[dict]:
    """Assign images to log rows day by day; return one report record per image and per unmatched row."""
    report = []
    matched_rows: set[int] = set()
    days = sorted({i.red_time.strftime("%Y-%m-%d") for i in images if i.red_time} | {r["date"] for r in rows})
    for day in days:
        day_images = [i for i in images if i.red_time and i.red_time.strftime("%Y-%m-%d") == day]
        day_rows = [(k, r) for k, r in enumerate(rows) if r["date"] == day and r.get("red")]
        if not day_images or not day_rows or day not in offsets:
            continue
        shift = timedelta(minutes=offsets[day])
        big = 1e6
        cost = np.full((len(day_images), len(day_rows)), big)
        for a, image in enumerate(day_images):
            for b, (_, row) in enumerate(day_rows):
                if _row_root_and_folder(row)[0] != image.root:
                    continue
                d_red = abs((image.red_time - shift - _log_time(row, "red")).total_seconds()) / 60
                log_green = _log_time(row, "green")
                d_green = d_red if image.green_time is None or log_green is None else \
                    abs((image.green_time - shift - log_green).total_seconds()) / 60
                if max(d_red, d_green) <= tolerance_min:
                    cost[a, b] = d_red + d_green
        for a, b in zip(*linear_sum_assignment(cost)):
            if cost[a, b] >= big:
                continue
            index, row = day_rows[b]
            image = day_images[a]
            matched_rows.add(index)
            image.__dict__["_match"] = (index, row, cost[a, b] / 2)
    for image in images:
        record = {"root": image.root, "image_folder": image.folder, "red_path": str(image.red_path),
                      "green_path": str(image.green_path),
                      "camera_red_time": "" if image.red_time is None else image.red_time.isoformat(sep=" ")}
        found = image.__dict__.get("_match")
        if image.red_time is None:
            record.update(status="no_timestamp")
        elif found is None:
            record.update(status="no_log_row_within_tolerance")
        else:
            index, row, residual = found
            _, row_folder = _row_root_and_folder(row)
            record.update(status="matched", log_index=index, residual_min=f"{residual:.2f}",
                          folder_agrees=_same_folder(image.folder, row_folder))
        report.append(record)
    for index, row in enumerate(rows):
        if index not in matched_rows and row.get("red"):
            report.append({"status": "log_row_without_image", "log_index": index, "root": _row_root_and_folder(row)[0],
                               "image_folder": _row_root_and_folder(row)[1]})
    return report


def pick_bright_field(image: ImageField, row: dict, offset_min: float, tolerance_min: float = 3.0) -> Path | None:
    """The bright-field frame whose corrected time is nearest the log's bright-field time."""
    log_bf = _log_time(row, "bf")
    if log_bf is None:
        return None
    shift = timedelta(minutes=offset_min)
    scored = [(abs((t - shift - log_bf).total_seconds()) / 60, p) for p, t in image.bf_candidates if t is not None]
    if not scored:
        return None
    best, path = min(scored)
    return path if best <= tolerance_min else None


MANIFEST_COLUMNS = [
    "session_id", "sample_id", "condition_id", "specimen_id", "biological_replicate_id", "timepoint",
    "acquisition_order", "green_path", "red_path", "bright_field_path", "exposure_ms_green", "exposure_ms_red",
    "nd_filter_green", "nd_filter_red", "objective", "burner_hours", "lamp_warmup_minutes", "gain_setting",
    "sampling_time_hours", "measurement_time_hours", "expect_cells", "log_flags",
]


def manifest_rows(images: list[ImageField], rows: list[dict], offsets: dict[str, float], *,
                  objective: str, gain_setting: str = "gain_3_gamma_1",
                  nd_filter: str = "", burner_hours: str = "", lamp_warmup_minutes: str = "") -> list[dict]:
    """Pipeline manifest rows for matched images, one per field.

    Each field gets its own sample_id because R1-R3 were imaged minutes apart
    and carry different measurement times. R1-R3 share specimen_id (same
    slide) and every row of one day and stressor shares biological_replicate_id
    (one culture), so fields are never counted as independent replicates.
    """
    out = []
    for image in images:
        found = image.__dict__.get("_match")
        if found is None:
            continue
        _, row, _ = found
        day = row["date"].replace("-", "")
        strain = "nt" if row["strain"] == "non_transfected" else "tr"
        stressor = row["stressor"].lower()
        dose = f"{row['dose_uM']}uM" if row["dose_uM"] else "medium"
        rep = row["rep"] or "F1"
        session = f"{day}_{stressor}_{strain}"
        bf = pick_bright_field(image, row, offsets[row["date"]])
        out.append({
            "session_id": session,
            "sample_id": f"{row['section']}_{dose}_{rep}",
            "condition_id": f"{stressor}_{dose}_{row['tp_label']}",
            "specimen_id": f"{row['section']}_{dose}",
            "biological_replicate_id": f"{day}_{stressor}_{strain}_culture",
            "timepoint": "", "acquisition_order": "",
            "green_path": str(image.green_path), "red_path": str(image.red_path),
            "bright_field_path": "" if bf is None else str(bf),
            "exposure_ms_green": row["exposure_acq_ms"], "exposure_ms_red": row["exposure_acq_ms"],
            "nd_filter_green": nd_filter, "nd_filter_red": nd_filter, "objective": objective,
            "burner_hours": burner_hours, "lamp_warmup_minutes": lamp_warmup_minutes, "gain_setting": gain_setting,
            "sampling_time_hours": row["sampling_time_hours"], "measurement_time_hours": row["measurement_time_hours"],
            "expect_cells": "false" if row["stressor"] == "none" else "true",
            "log_flags": row["flags"], "_red_time": image.red_time, "_tp": row["tp_label"],
        })
    # timepoint: index of the time label within a session, ordered by first acquisition;
    # acquisition_order: order of the field within its session.
    for session in {r["session_id"] for r in out}:
        members = sorted((r for r in out if r["session_id"] == session), key=lambda r: r["_red_time"])
        labels = list(dict.fromkeys(r["_tp"] for r in members))
        for order, r in enumerate(members):
            r["acquisition_order"] = order
            r["timepoint"] = labels.index(r["_tp"])
    for r in out:
        del r["_red_time"], r["_tp"]
    return sorted(out, key=lambda r: (r["session_id"], r["acquisition_order"]))


def write_csv(records: list[dict], path: Path, columns: list[str] | None = None) -> None:
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"output path already exists: {path}")
    columns = columns or sorted({k for r in records for k in r})
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
