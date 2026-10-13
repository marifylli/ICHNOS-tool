"""Within-day induction test: each dose against that day's 0 µM."""
import runpy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ichnos_image import induction
from ichnos_image import timecourse as tc

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "detect_induction.py"
TIMES = [1.0, 2.0, 3.0, 4.0]


def synthetic(seed=0, cells=80, responding_day="20261007", effect_per_log_dose=2.0):
    """Two transfected days and a non-transfected day.

    Each day has its own transformation level (a constant offset), which a
    within-day comparison must ignore. Only ``responding_day`` responds, with
    red rising with log10(dose).
    """
    rng = np.random.default_rng(seed)
    rows = []

    def add(session, dose, label, t, green, red):
        for k in range(cells):
            rows.append(dict(session_id=session, condition_id=f"h2o2_{dose}uM_{label}", sample_id=f"s{dose}_R1",
                             cell_id=k, qc_pass=True, edge_flag=False, measurement_time_hours=t,
                             corrected_mean_green_per_s=green + rng.normal(0, 1.0),
                             corrected_mean_red_per_s=red + rng.normal(0, 1.0)))
    for day, level in (("20261005", 6.0), ("20261007", 9.0)):
        for i, t in enumerate(TIMES):
            for dose in (0, 20, 200, 2000):
                lift = effect_per_log_dose * np.log10(dose) if (dose and day == responding_day) else 0.0
                add(f"{day}_h2o2_tr", dose, f"{i}h", t + 0.02 * (dose > 0), level, level + lift)
    for i, t in enumerate(TIMES):
        add("20261008_h2o2_nt", 0, f"{i}h", t, 5.0, 5.0)
    return pd.DataFrame(rows)


def prepared(**kw):
    cells = tc.annotate(synthetic(**kw))
    return tc.subtract_baseline(cells, tc.autofluorescence_baseline(cells))


def test_a_responding_day_is_found_and_a_flat_day_is_not():
    effects, trend = induction.induction_table(prepared(), n_boot=200, seed=1)
    red = effects[effects["readout"] == "red"].set_index(["day", "dose_uM"])
    assert red.loc[("20261007", 2000), "effect"] == pytest.approx(2.0 * 3.3, abs=0.4)
    assert red.loc[("20261007", 2000), "detected_bh"]
    flat = red.loc["20261005"]
    assert not flat["detected_bh"].any()
    assert ((flat["ci_low"] < 0) & (flat["ci_high"] > 0)).all()
    slope = trend.set_index(["day", "readout"]).loc[("20261007", "red")]
    assert slope["slope_per_log10_dose"] == pytest.approx(2.0, abs=0.3) and slope["ci_low"] > 0


def test_the_transformation_level_of_a_day_does_not_count_as_an_effect():
    """Day 7/10 sits 3 grey/s above day 5/10 at every dose, 0 µM included."""
    effects, _ = induction.induction_table(prepared(responding_day="none"), n_boot=100, seed=2)
    assert effects["effect"].abs().max() < 0.6
    assert not effects["detected_bh"].any()


def test_minimum_detectable_effect_shrinks_with_more_cells():
    few, _ = induction.induction_table(prepared(cells=20, responding_day="none"), n_boot=150, seed=3)
    many, _ = induction.induction_table(prepared(cells=320, responding_day="none"), n_boot=150, seed=3)
    assert many["mde"].median() < 0.7 * few["mde"].median()


def test_benjamini_hochberg_matches_the_textbook_values():
    q = induction._benjamini_hochberg(np.array([0.01, 0.04, 0.03, 0.20]))
    assert q.tolist() == pytest.approx([0.04, 0.0533333, 0.0533333, 0.20], abs=1e-6)


def test_script_writes_outputs_and_refuses_an_existing_directory(tmp_path):
    path = tmp_path / "cells_per_s.csv"
    synthetic().to_csv(path, index=False)
    main = runpy.run_path(str(SCRIPT))["main"]
    out = tmp_path / "ind"
    assert main(["--cells", str(path), "--out-dir", str(out), "--n-boot", "50"]) == 0
    assert {"effects.csv", "trend.csv", "run.json"} <= {p.name for p in out.iterdir()}
    with pytest.raises(SystemExit, match="already exists"):
        main(["--cells", str(path), "--out-dir", str(out)])
