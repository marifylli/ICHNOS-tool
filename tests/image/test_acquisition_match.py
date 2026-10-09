import csv
import unicodedata

import numpy as np
import tifffile

from ichnos_image.acquisition_match import (
    estimate_offsets,
    manifest_rows,
    match,
    scan_root,
    tiff_datetime,
)

CAMERA_AHEAD_MIN = 50


def _tif(path, stamp):
    """A tiny TIFF carrying an Image-Pro style DateTime tag."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, np.zeros((4, 4), np.uint8), extratags=[(306, "s", 0, stamp, True)])


def _camera(hh, mm, ss=10):
    """Camera stamp for a wall-clock time, with the camera clock 50 min ahead."""
    total = hh * 60 + mm + CAMERA_AHEAD_MIN
    h, m = divmod(total, 60)
    return f"10/09/2026 {h % 12 or 12:02d}:{m:02d}:{ss:02d}.000 {'PM' if h >= 12 else 'AM'}"


def _field(root, condition, name, red, green, bf):
    base = root / condition / name / "01a"
    _tif(base / f"{name} [01a][Alexa 594].tif", _camera(*red))
    _tif(base / f"{name} [01a][Alexa 488].tif", _camera(*green))
    _tif(root / condition / f"Image000 {name}.tif", _camera(*bf))


def _row(dose, rep, red, green, bf, folder):
    return {"section": "nt_dtt_0910_2h", "date": "2026-10-09", "stressor": "DTT", "strain": "non_transfected",
                "tp_label": "2h", "dose_uM": dose, "rep": rep, "bf": bf, "red": red, "green": green, "exposure_acq_ms": "600",
                "sampling_time_hours": "", "measurement_time_hours": "3.200", "folder": folder, "flags": ""}


def _setup(tmp_path):
    root = tmp_path / unicodedata.normalize("NFD", "igem DTT non transfected")
    # As on 2026-10-09: the 100 uM field was saved in 200/2h under a "50-2Η" name, and vice versa.
    _field(root, "200/2h", "50-2Η", red=(15, 33), green=(15, 33), bf=(15, 34))
    _field(root, "50/2h", "50-2Η", red=(15, 39), green=(15, 39), bf=(15, 40))
    _field(root, "100/2h", "100-2Η", red=(15, 45), green=(15, 45), bf=(15, 46))
    _field(root, "100/2h", "100-2Η extra", red=(16, 30), green=(16, 30), bf=(16, 31))  # nothing in the log
    rows = [
        _row("100", "R1", "15:33", "15:33", "15:34", "D:igem DTT non transfected/200/2h"),
        _row("50", "R1", "15:39", "15:39", "15:40", "D:igem DTT non transfected/50/2h"),
        _row("200", "R1", "15:45", "15:45", "15:46", "D:igem DTT non transfected/100/2h"),
        _row("25", "R3", "15:30", "15:30", "15:32", "D:igem DTT non transfected/25/2h"),  # no image
    ]
    return scan_root("D", root), rows


def test_reads_image_pro_timestamp(tmp_path):
    _tif(tmp_path / "a.tif", "10/06/2026 06:11:10.138 PM")
    assert tiff_datetime(tmp_path / "a.tif").strftime("%Y-%m-%d %H:%M:%S") == "2026-10-06 18:11:10"


def test_offset_is_estimated_from_same_folder_pairs(tmp_path):
    images, rows = _setup(tmp_path)
    offsets = estimate_offsets(images, rows)
    assert abs(offsets["2026-10-09"] - (CAMERA_AHEAD_MIN - 20 / 60)) < 0.05


def test_images_follow_time_not_folder_name(tmp_path):
    images, rows = _setup(tmp_path)
    report = match(images, rows, estimate_offsets(images, rows))
    by_path = {r["red_path"]: r for r in report if r["status"] == "matched"}
    doses = {rows[r["log_index"]]["dose_uM"]: r["red_path"] for r in by_path.values()}
    assert "/200/2h/" in doses["100"] and "/100/2h/" in doses["200"]
    assert all(r["folder_agrees"] for r in by_path.values())


def test_unmatched_images_and_rows_are_reported(tmp_path):
    images, rows = _setup(tmp_path)
    report = match(images, rows, estimate_offsets(images, rows))
    statuses = sorted(r["status"] for r in report)
    assert statuses.count("no_log_row_within_tolerance") == 1
    assert next(r for r in report if r["status"] == "log_row_without_image")["log_index"] == 3


def test_manifest_carries_log_times_exposure_and_bright_field(tmp_path):
    images, rows = _setup(tmp_path)
    offsets = estimate_offsets(images, rows)
    match(images, rows, offsets)
    manifest = manifest_rows(images, rows, offsets, objective="40X")
    assert len(manifest) == 3
    first = manifest[0]
    assert first["sample_id"] == "nt_dtt_0910_2h_100uM_R1"
    assert first["exposure_ms_green"] == first["exposure_ms_red"] == "600"
    assert first["measurement_time_hours"] == "3.200"
    assert first["bright_field_path"].endswith("Image000 50-2Η.tif")
    assert [m["acquisition_order"] for m in manifest] == [0, 1, 2]
    assert len({m["biological_replicate_id"] for m in manifest}) == 1


def test_manifest_is_readable_by_csv(tmp_path):
    from ichnos_image.acquisition_match import MANIFEST_COLUMNS, write_csv
    images, rows = _setup(tmp_path)
    offsets = estimate_offsets(images, rows)
    match(images, rows, offsets)
    out = tmp_path / "images.csv"
    write_csv(manifest_rows(images, rows, offsets, objective="40X"), out, MANIFEST_COLUMNS)
    with out.open(encoding="utf-8") as handle:
        assert next(csv.reader(handle)) == MANIFEST_COLUMNS
