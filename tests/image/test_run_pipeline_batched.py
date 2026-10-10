import json
import runpy
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
batched = runpy.run_path(str(ROOT / "scripts" / "run_pipeline_batched.py"))


def _manifest(tmp_path, n):
    path = tmp_path / "images.csv"
    pd.DataFrame({"session_id": ["s"] * n, "sample_id": [f"f{i}" for i in range(n)]}).to_csv(path, index=False)
    return path


def _fake_runner(calls, fail_on=None):
    def runner(manifest, out_csv, extra):
        calls.append((manifest.name, list(extra)))
        rows = pd.read_csv(manifest)
        if fail_on is not None and fail_on in set(rows.sample_id):
            return -9  # what the shell reports as "killed"
        pd.DataFrame({"sample_id": rows.sample_id, "cell_id": 1}).to_csv(out_csv, index=False)
        Path(str(out_csv) + ".manifest.json").write_text(json.dumps({"refused_image_sets": []}))
        return 0
    return runner


def test_batches_cover_every_row_once_in_order(tmp_path):
    calls = []
    summary = batched["run_batched"](_manifest(tmp_path, 7), tmp_path / "out", 3, ["--x", "1"], _fake_runner(calls))
    assert [c[0] for c in calls] == ["part_000.csv", "part_001.csv", "part_002.csv"]
    assert all(c[1] == ["--x", "1"] for c in calls)
    merged = pd.read_csv(tmp_path / "out" / "cells.csv")
    assert list(merged.sample_id) == [f"f{i}" for i in range(7)]
    assert [b["rows"] for b in summary["batches"]] == [[0, 3], [3, 6], [6, 7]]


def test_failed_batch_stops_and_rerun_resumes(tmp_path):
    manifest, out = _manifest(tmp_path, 6), tmp_path / "out"
    calls = []
    first = batched["run_batched"](manifest, out, 2, [], _fake_runner(calls, fail_on="f3"))
    assert first["batches"][-1]["status"] == "failed_exit_-9"
    assert not (out / "cells.csv").exists()
    calls.clear()
    batched["run_batched"](manifest, out, 2, [], _fake_runner(calls))
    assert [c[0] for c in calls] == ["part_001.csv", "part_002.csv"]  # part 0 reused
    assert len(pd.read_csv(out / "cells.csv")) == 6


def test_refuses_to_overwrite_merged_output(tmp_path):
    manifest, out = _manifest(tmp_path, 2), tmp_path / "out"
    batched["run_batched"](manifest, out, 5, [], _fake_runner([]))
    with pytest.raises(SystemExit):
        batched["run_batched"](manifest, out, 5, [], _fake_runner([]))


def test_arguments_after_double_dash_go_to_the_pipeline(tmp_path, monkeypatch):
    seen = {}
    monkeypatch.setitem(batched["main"].__globals__, "run_batched",
                        lambda manifest, out_dir, size, extra: seen.update(size=size, extra=extra) or {"batches": []})
    batched["main"](["--manifest", "m.csv", "--out-dir", str(tmp_path), "--batch-size", "5",
                     "--", "--bleed-default", "0.05"])
    assert seen == {"size": 5, "extra": ["--bleed-default", "0.05"]}


def test_batch_without_manifest_is_not_taken_as_done(tmp_path):
    manifest, out = _manifest(tmp_path, 2), tmp_path / "out"
    (out / "parts").mkdir(parents=True)
    (out / "parts" / "part_000.cells.csv").write_text("sample_id,cell_id\nf0,1\n")  # interrupted publish
    with pytest.raises(SystemExit, match="no manifest"):
        batched["run_batched"](manifest, out, 5, [], _fake_runner([]))
