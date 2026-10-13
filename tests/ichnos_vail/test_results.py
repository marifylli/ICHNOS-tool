"""Reading script outputs for display."""
import inspect
import json
import re

import pandas as pd

from ichnos_vail import results
from ichnos_image import field_decoder


def test_every_verdict_the_decoder_can_give_has_a_colour_and_an_instruction():
    source = inspect.getsource(field_decoder._verdict) + inspect.getsource(field_decoder._time_verdict)
    produced = set(re.findall(r'return "([^"]+)"', source))
    assert produced and produced == set(results.VERDICTS)


def test_track_record_wording():
    good = {"validated_accuracy": 0.9, "validated_info_gain_bits": 1.2, "chance": 0.2, "reference_days": ["a", "b", "c"]}
    bad = {**good, "validated_accuracy": 0.19, "validated_info_gain_bits": -1.78}
    assert "better than chance" in results.track_record_text(good) and "no better" not in results.track_record_text(good)
    assert "no better than chance" in results.track_record_text(bad)
    assert "unavailable" in results.track_record_text({"validated_accuracy": None, "reference_days": ["a"]})


def test_sessions_are_parsed_where_named_by_the_convention(tmp_path):
    path = tmp_path / "cells.csv"
    pd.DataFrame({"session_id": ["20261005_h2o2_tr", "20261008_h2o2_nt", "user_day", "20261005_h2o2_tr"],
                  "x": range(4)}).to_csv(path, index=False)
    table = results.sessions(path).set_index("session_id")
    assert table.loc["20261008_h2o2_nt", "strain"] == "nt" and pd.isna(table.loc["user_day", "stressor"])
    assert len(table) == 3


def test_analysis_and_decode_folders_are_read(tmp_path):
    (tmp_path / "induction").mkdir()
    pd.DataFrame({"day": ["d"], "stressor": ["h2o2"], "dose_uM": [300], "readout": ["red"], "effect": [2.0],
                  "ci_low": [1.0], "ci_high": [3.0], "mde": [1.0], "q_bh": [0.01], "detected_bh": [True]}
                 ).to_csv(tmp_path / "induction" / "effects.csv", index=False)
    (tmp_path / "model_compare").mkdir()
    (tmp_path / "model_compare" / "compare.png").write_bytes(b"png")
    tables = results.analysis_tables(tmp_path)
    assert len(tables["Detected effects (q < 0.05)"]) == 1
    assert [t for t, _ in results.analysis_figures(tmp_path)] == ["Model against data"]
    pd.DataFrame({"unit": ["field"], "verdict": ["no dose information"]}).to_csv(tmp_path / "decoded.csv", index=False)
    (tmp_path / "summary.json").write_text(json.dumps({"chance": 0.2}))
    summary, decoded = results.decode_result(tmp_path)
    assert summary["chance"] == 0.2 and decoded["verdict"].iloc[0] == "no dose information"


def test_time_record_wording():
    assert "cannot be decoded" in results.time_record_text({"time_available": False, "time_reason": "one window"})
    unscored = {"time_available": True, "time_validated_accuracy": None, "time_windows_in_reference": ["0–2 h"]}
    assert "unavailable" in results.time_record_text(unscored)
    good = {"time_available": True, "time_validated_accuracy": 0.9, "time_validated_info_gain_bits": 1.1,
            "time_chance": 0.25, "time_windows_in_reference": ["0–2 h", "2–4 h"], "time_scored_fields": 30,
            "time_unscorable_fields": 0}
    assert "better than chance" in results.time_record_text(good)
