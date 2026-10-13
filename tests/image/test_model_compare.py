"""Model prediction against measured induction, as fold over the same day's 0 µM."""
import runpy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ichnos_image import model_compare as mc
from ichnos_image import timecourse as tc

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "compare_model.py"
TIMES = [1.0, 2.0, 3.0, 4.0]


def synthetic(seed=0, cells=60, stressed_red=24.0, extra_late_time=False):
    """Non-transfected cells at 5 grey/s; transfected 0 µM at 12; 300 µM red at ``stressed_red``.

    Plus a dim transfected population (below the non-transfected level) that the
    reporter-positive filter must drop.
    """
    rng = np.random.default_rng(seed)
    rows = []

    def add(session, dose, label, t, green, red, n=cells):
        for k in range(n):
            rows.append(dict(session_id=session, condition_id=f"h2o2_{dose}uM_{label}", sample_id=f"s{dose}_R1",
                             cell_id=k, qc_pass=True, edge_flag=False, measurement_time_hours=t,
                             corrected_mean_green_per_s=green + rng.normal(0, 0.5),
                             corrected_mean_red_per_s=red + rng.normal(0, 0.5)))
    for i, t in enumerate(TIMES):
        add("20261007_h2o2_tr", 0, f"{i}h", t, 12.0, 12.0)
        add("20261007_h2o2_tr", 300, f"{i}h", t + 0.05, 12.0, stressed_red)
        add("20261007_h2o2_tr", 300, f"dim{i}h", t + 0.05, 4.0, 4.0, n=cells // 2)
        add("20261008_h2o2_nt", 0, f"{i}h", t, 5.0, 5.0)
    if extra_late_time:
        add("20261007_h2o2_tr", 300, "late", 9.0, 12.0, stressed_red)
    return pd.DataFrame(rows)


def prepared(**kw):
    cells = tc.annotate(synthetic(**kw))
    return tc.subtract_baseline(cells, tc.autofluorescence_baseline(cells))


def test_reporter_positive_keeps_bright_transfected_and_all_non_transfected_cells():
    cells = prepared()
    kept = mc.reporter_positive(cells)
    assert (kept["strain"] == "nt").sum() == (cells["strain"] == "nt").sum()
    tr = kept[kept["strain"] == "tr"]
    assert not tr["condition_id"].str.contains("dim").any()
    assert len(tr) == 2 * len(TIMES) * 60


def test_data_fold_is_stressed_over_control_after_the_baseline():
    """Baseline 5: control 12 -> 7 above it, stressed red 24 -> 19, fold 19/7; green unchanged."""
    fold = mc.data_fold(prepared(), n_boot=100, seed=1).set_index("readout")
    assert fold.loc["red", "fold"] == pytest.approx(19 / 7, rel=0.03)
    assert fold.loc["red", "fold_ci_low"] > 2.4
    assert fold.loc["green", "fold"] == pytest.approx(1.0, abs=0.03)
    assert fold.loc["green", "fold_ci_low"] < 1 < fold.loc["green", "fold_ci_high"]


def test_imaging_times_stay_inside_the_control_range():
    times = mc.imaging_times(prepared(extra_late_time=True))
    assert times == {("20261007_h2o2_tr", 300.0): [1.05, 2.05, 3.05]}


def fake_run(variant, dose, times, clearance, profile):
    t = np.asarray(times, float)
    level = 2.0 * (1 + dose / 100 / (1 + (clearance or 0)))
    return {"green": np.full(t.shape, level), "red": np.full(t.shape, 2 * level)}


def test_model_fold_per_scenario_and_unknown_stressors_skipped():
    times = {("a_h2o2_tr", 300.0): [1.0, 2.0], ("b_none_tr", 300.0): [1.0]}
    model = mc.model_fold(times, {"a_h2o2_tr": "h2o2", "b_none_tr": "none"}, clearances=(None, 2.0), run=fake_run)
    assert set(model["session_id"]) == {"a_h2o2_tr"}
    constant = model[model["clearance_per_h"].isna()].set_index("readout")["model_fold"]
    cleared = model[model["clearance_per_h"] == 2.0].set_index("readout")["model_fold"]
    assert constant["red"] == pytest.approx(4.0) and constant["green"] == pytest.approx(4.0)
    assert cleared["red"] == pytest.approx(2.0)


def test_compare_marks_whether_the_model_falls_inside_the_interval():
    data = pd.DataFrame({"session_id": ["a", "a"], "day": ["d", "d"], "stressor": ["h2o2"] * 2,
                         "dose_uM": [300.0, 300.0], "readout": ["green", "red"], "n_timepoints": [2, 2],
                         "control_level": [1.0, 1.0], "fold": [1.1, 1.1], "fold_ci_low": [0.8, 0.8],
                         "fold_ci_high": [1.5, 1.5], "fold_mde": [0.4, 0.4], "q_bh": [1.0, 1.0]})
    model = pd.DataFrame({"session_id": ["a", "a"], "stressor": ["h2o2"] * 2, "variant": ["ox"] * 2,
                          "dose_uM": [300.0, 300.0], "readout": ["green", "red"], "clearance_per_h": [np.nan] * 2,
                          "times_h": ["1 2"] * 2, "model_fold": [1.3, 4.0]})
    table = mc.compare(data, model).set_index("readout")
    assert table.loc["green", "model_inside_ci"] and not table.loc["red", "model_inside_ci"]
    assert table.loc["red", "model_over_data"] == pytest.approx(4.0 / 1.1)


def test_the_real_model_rises_with_h2o2_and_less_with_clearance():
    pytest.importorskip("roadrunner")
    times = {("s", 600.0): [2.0, 4.0]}
    model = mc.model_fold(times, {"s": "h2o2"}, clearances=(None, 3.0))
    red = model[model["readout"] == "red"]
    red = red.set_index(red["clearance_per_h"].fillna(0))["model_fold"]
    assert red.loc[0.0] > 2.0
    assert 1.0 < red.loc[3.0] < red.loc[0.0]


def test_script_writes_outputs_and_refuses_an_existing_directory(tmp_path):
    pytest.importorskip("roadrunner")
    path = tmp_path / "cells_per_s.csv"
    synthetic().to_csv(path, index=False)
    main = runpy.run_path(str(SCRIPT))["main"]
    out = tmp_path / "cmp"
    assert main(["--cells", str(path), "--out-dir", str(out), "--n-boot", "50", "--clearance", "none", "3"]) == 0
    for name in ("data_fold.csv", "model_fold.csv", "compare.csv", "compare.png", "run.json"):
        assert (out / name).exists()
    table = pd.read_csv(out / "compare.csv")
    assert len(table) == 2 * 2  # two channels x two scenarios
    with pytest.raises(SystemExit, match="already exists"):
        main(["--cells", str(path), "--out-dir", str(out)])
