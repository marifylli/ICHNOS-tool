"""The GUI's commands are the scripts' own command lines, checked before they start."""
import json
import runpy
import sys

from ichnos_vail import commands
from ichnos_vail.jobs import ROOT


def test_images_to_cells_runs_the_batched_pipeline_then_normalisation(tmp_path):
    manifest = tmp_path / "images.csv"
    manifest.write_text("session_id\n")
    cmds, errors = commands.images_to_cells(manifest, tmp_path / "run", batch_size=10, gpu=True,
                                            bleed_default=0.05, saturation=255, green_extraction="G")
    assert errors == []
    batched, normalise = (commands.as_shell(c) for c in cmds)
    assert batched.startswith("python3 scripts/run_pipeline_batched.py --manifest")
    assert "-- --segmentation-method cellpose --segmentation-source brightfield --saturation-value 255 " \
           "--cellpose-gpu --bleed-default 0.05 --green-extraction G" in batched
    assert normalise.startswith("python3 scripts/normalize_exposure.py --cells")
    assert cmds[1][-1].endswith("cells_per_s.csv")


def test_images_to_cells_problems_are_reported_before_starting(tmp_path):
    (tmp_path / "run").mkdir()
    (tmp_path / "run" / "cells.csv").write_text("x\n")
    _, errors = commands.images_to_cells(tmp_path / "missing.csv", tmp_path / "run", bleed_default=None)
    text = " ".join(errors)
    assert "manifest not found" in text and "already exists" in text and "bleed-through" in text
    cmds, _ = commands.images_to_cells(tmp_path / "m.csv", tmp_path / "x", segmentation="sparse", gpu=True)
    assert "--cellpose-gpu" not in cmds[0]


def test_controls_replace_the_default_bleed(tmp_path):
    controls = tmp_path / "controls.csv"
    controls.write_text("session_id\n")
    cmds, errors = commands.images_to_cells(tmp_path / "m.csv", tmp_path / "x", controls=controls)
    assert "--controls" in cmds[0] and "--bleed-default" not in cmds[0]


def test_decode_command_and_its_checks(tmp_path):
    ref, new = tmp_path / "ref.csv", tmp_path / "new.csv"
    ref.write_text("x\n")
    new.write_text("x\n")
    cmds, errors = commands.decode(ref, new, "h2o2", tmp_path / "out", exclude=["a_h2o2_tr", "b_h2o2_tr"])
    assert errors == []
    assert cmds[0].count("--exclude-reference-session") == 2
    (tmp_path / "out").mkdir()
    _, errors = commands.decode(ref, tmp_path / "nope.csv", "", tmp_path / "out", decoder="model")
    assert len(errors) == 4


def test_the_analysis_configuration_is_one_run_analysis_accepts(tmp_path):
    cells = tmp_path / "cells.csv"
    cells.write_text("x\n")
    config = tmp_path / "cfg.json"
    commands.write_analysis_config(config, cells, steps=["induction"], n_boot=50,
                                   clearance_per_h=[None, 2.0], decoders=["paired"])
    loaded = runpy.run_path(str(ROOT / "scripts" / "run_analysis.py"))["load_config"](config)
    assert loaded["steps"] == ["induction"] and loaded["clearance_per_h"] == [None, 2.0]
    assert json.loads(config.read_text())["cells"] == str(cells.resolve())
    cmds, errors = commands.analysis(config, cells, tmp_path / "out")
    assert errors == [] and cmds[0][0] == sys.executable


def test_align_writes_a_new_manifest_in_the_output_folder(tmp_path):
    manifest = tmp_path / "images.csv"
    manifest.write_text("x\n")
    cmds, errors = commands.align(manifest, tmp_path / "aligned", search_up=1)
    assert errors == []
    shown = commands.as_shell(cmds[0])
    assert "--write-manifest" in shown and "images_aligned.csv" in shown and "--search-up 1" in shown


def test_time_windows_are_checked_and_passed_on(tmp_path):
    ref = tmp_path / "ref.csv"
    ref.write_text("x\n")
    cmds, errors = commands.decode(ref, ref, "h2o2", tmp_path / "o", time_windows="0,3,6,9")
    assert errors == [] and cmds[0][cmds[0].index("--time-windows") + 1] == "0,3,6,9"
    for bad in ("0,2", "0,4,2", "a,b,c", "-1,2,4"):
        _, errors = commands.decode(ref, ref, "h2o2", tmp_path / "o", time_windows=bad)
        assert any("time windows" in e for e in errors), bad
