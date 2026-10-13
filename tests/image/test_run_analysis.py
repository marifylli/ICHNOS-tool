"""One configuration runs the analysis steps on a cell table and records what it did."""
import json
import runpy
from pathlib import Path

import pandas as pd
import pytest

from tests.image.test_model_compare import synthetic

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "run_analysis.py"


def cell_table(tmp_path):
    cells = synthetic()
    second_day = cells[cells["session_id"] == "20261007_h2o2_tr"].assign(session_id="20261006_h2o2_tr")
    path = tmp_path / "cells_per_s.csv"
    pd.concat([cells, second_day]).to_csv(path, index=False)
    return path


def write_config(tmp_path, **entries):
    path = tmp_path / "analysis.json"
    path.write_text(json.dumps(entries))
    return path


def main():
    return runpy.run_path(str(SCRIPT))["main"]


def test_steps_run_into_their_folders_with_a_manifest(tmp_path):
    cell_table(tmp_path)
    config = write_config(tmp_path, cells="cells_per_s.csv", steps=["timecourse", "induction", "decoder_lodo"],
                          n_boot=30, decoders=["absolute", "paired"])
    out = tmp_path / "run"
    assert main()(["--config", str(config), "--out-dir", str(out)]) == 0
    assert (out / "timecourse" / "timecourse.csv").exists()
    assert (out / "induction" / "effects.csv").exists()
    assert (out / "decoder_lodo" / "metrics.csv").exists()
    manifest = json.loads((out / "run-manifest.json").read_text())
    assert manifest["status"] == "ok"
    assert [s["step"] for s in manifest["steps"]] == ["timecourse", "induction", "decoder_lodo"]
    assert manifest["config"]["cells"] == str((tmp_path / "cells_per_s.csv").resolve())
    assert "induction/effects.csv" in manifest["outputs"] and len(manifest["input_sha256"]) == 64


def test_unknown_keys_and_steps_are_refused(tmp_path):
    cell_table(tmp_path)
    with pytest.raises(SystemExit, match="unknown configuration keys"):
        main()(["--config", str(write_config(tmp_path, cells="cells_per_s.csv", nboot=10)),
                "--out-dir", str(tmp_path / "a")])
    with pytest.raises(SystemExit, match="unknown steps"):
        main()(["--config", str(write_config(tmp_path, cells="cells_per_s.csv", steps=["peaks"])),
                "--out-dir", str(tmp_path / "b")])
    assert not (tmp_path / "a").exists() and not (tmp_path / "b").exists()


def test_an_existing_output_directory_is_refused(tmp_path):
    cell_table(tmp_path)
    (tmp_path / "run").mkdir()
    with pytest.raises(SystemExit, match="already exists"):
        main()(["--config", str(write_config(tmp_path, cells="cells_per_s.csv")), "--out-dir", str(tmp_path / "run")])


def test_a_failing_step_stops_the_run_and_is_recorded(tmp_path):
    path = cell_table(tmp_path)
    pd.read_csv(path).drop(columns=["corrected_mean_red_per_s"]).to_csv(path, index=False)
    config = write_config(tmp_path, cells="cells_per_s.csv", steps=["induction", "decoder_lodo"], n_boot=10)
    out = tmp_path / "run"
    assert main()(["--config", str(config), "--out-dir", str(out)]) == 1
    manifest = json.loads((out / "run-manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert [s["step"] for s in manifest["steps"]] == ["induction"]
    assert manifest["steps"][0]["status"].startswith("stopped")
