"""Manifest problems found before a run, with the row and column they are in."""
import pandas as pd

from ichnos_vail.manifest_check import check_manifest

COLUMNS = ["session_id", "sample_id", "condition_id", "timepoint", "acquisition_order", "green_path", "red_path",
           "bright_field_path", "exposure_ms_green", "exposure_ms_red", "nd_filter_green", "nd_filter_red",
           "objective", "burner_hours", "lamp_warmup_minutes", "measurement_time_hours"]


def write(tmp_path, rows, drop=()):
    for name in ("g.tif", "r.tif", "bf.tif"):
        (tmp_path / name).write_bytes(b"x")
    base = dict(timepoint="0", acquisition_order="0", green_path="g.tif", red_path="r.tif", bright_field_path="bf.tif",
                exposure_ms_green="400", exposure_ms_red="400", nd_filter_green="", nd_filter_red="",
                objective="40X", burner_hours="", lamp_warmup_minutes="", measurement_time_hours="2")
    table = pd.DataFrame([{**base, **r} for r in rows])[[c for c in COLUMNS if c not in drop]]
    path = tmp_path / "images.csv"
    table.to_csv(path, index=False)
    return path


GOOD = [dict(session_id="20261101_h2o2_tr", condition_id="h2o2_0uM_2h", sample_id="a"),
        dict(session_id="20261101_h2o2_tr", condition_id="h2o2_300uM_2h", sample_id="b"),
        dict(session_id="20261103_h2o2_nt", condition_id="h2o2_0uM_2h", sample_id="c")]


def test_a_complete_manifest_has_no_problems_and_is_summarised(tmp_path):
    _, issues, summary = check_manifest(write(tmp_path, GOOD))
    assert issues == []
    assert summary["rows"] == 3
    assert summary["sessions"]["20261101_h2o2_tr"]["doses"] == [0.0, 300.0]


def test_errors_name_the_row_and_column(tmp_path):
    rows = [dict(GOOD[0]), dict(GOOD[1], exposure_ms_red="0", measurement_time_hours="-1", objective="100X",
                                red_path="missing.tif"),
            dict(GOOD[1]), dict(GOOD[2])]  # row 3 repeats row 2's identity
    _, issues, summary = check_manifest(write(tmp_path, rows))
    errors = {(i["row"], i["column"]) for i in issues if i["level"] == "error"}
    assert {(2, "exposure_ms_red"), (2, "measurement_time_hours"), (2, "objective"), (2, "red_path")} <= errors
    assert (2, "sample_id") in errors and (3, "sample_id") in errors
    assert summary["errors"] == len([i for i in issues if i["level"] == "error"])


def test_missing_required_columns_are_errors_and_missing_analysis_columns_warnings(tmp_path):
    _, issues, _ = check_manifest(write(tmp_path, GOOD, drop=("objective", "measurement_time_hours")))
    levels = {i["column"]: i["level"] for i in issues if i["row"] == ""}
    assert levels["objective"] == "error" and levels["measurement_time_hours"] == "warning"


def test_design_warnings_for_analysis(tmp_path):
    rows = [dict(session_id="20261101_h2o2_tr", condition_id="h2o2_300uM_2h", sample_id="a"),
            dict(session_id="20261101_h2o2_tr", condition_id="mystery", sample_id="b"),
            dict(session_id="my culture", condition_id="h2o2_0uM_x", sample_id="c")]
    _, issues, _ = check_manifest(write(tmp_path, rows))
    text = " | ".join(i["message"] for i in issues if i["level"] == "warning")
    assert "no 0 µM culture" in text and "dose unknown" in text
    assert "no non-transfected h2o2" in text and "not <yyyymmdd>_<stressor>_<tr|nt>" in text
    assert not [i for i in issues if i["level"] == "error"]


def test_an_unreadable_file_is_one_clear_error(tmp_path):
    table, issues, _ = check_manifest(tmp_path / "nope.csv")
    assert table is None and issues[0]["level"] == "error" and "cannot read" in issues[0]["message"]
