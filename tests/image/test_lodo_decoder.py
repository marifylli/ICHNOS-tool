"""Single-field dose decoding, validated leave-one-day-out."""
import runpy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ichnos_image import lodo_decoder as ld

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "evaluate_decoder_lodo.py"
DOSES = [0.0, 100.0, 1000.0]


def fields(slope=3.0, days=3, day_spread=0.3, noise=0.3, seed=0):
    """Field features for ``days`` days; red rises by ``slope`` per 10-fold dose (0 µM counts as 1 µM)."""
    rng = np.random.default_rng(seed)
    rows = []
    for k in range(days):
        session = f"2026100{k + 1}_h2o2_tr"
        offset = rng.normal(0, day_spread)
        for dose in DOSES:
            lift = slope * np.log10(max(dose, 1.0))
            for t in (1.0, 2.0, 3.0):
                condition = f"h2o2_{dose:g}uM_{t:g}h"
                for r in range(3):
                    red = 5 + offset + lift + rng.normal(0, noise)
                    rows.append(dict(session_id=session, day=session[:8], stressor="h2o2", dose_uM=dose,
                                     condition_id=condition, sample_id=f"{condition}_R{r + 1}", time_h=t + 0.01 * r,
                                     n_cells=50, n_positive=20, green=4 + offset + rng.normal(0, noise), red=red,
                                     positive=0.2 + rng.normal(0, 0.02), green_pos=6 + rng.normal(0, 0.1),
                                     red_pos=red + 3))
    return ld.add_paired(pd.DataFrame(rows))


def scores(predictions, unit="field"):
    s = ld.summarize(predictions)
    return s[(s["unit"] == unit) & (s["held_out"] == "mean over days")].iloc[0]


def test_a_dose_signal_is_decoded_on_held_out_days():
    for decoder in ("absolute", "paired"):
        predictions, folds = ld.leave_one_day_out(fields(), decoder)
        assert folds["between_day_estimated"].all() and (folds["n_training_days"] == 2).all()
        s = scores(predictions)
        assert s["accuracy"] > 0.9 and s["info_gain_bits"] > 1.0
        assert s["credible80_coverage"] >= 0.8


def test_without_signal_the_posterior_stays_near_the_prior():
    """The property step 11 asks for: no dose information in the data -> no confidence in the output."""
    predictions, _ = ld.leave_one_day_out(fields(slope=0.0), "paired")
    s = scores(predictions)
    assert s["entropy_vs_prior"] > 0.9
    assert abs(s["info_gain_bits"]) < 0.3
    assert s["accuracy"] < 0.6


def test_pooling_the_fields_of_a_sample_sharpens_a_real_signal():
    predictions, _ = ld.leave_one_day_out(fields(slope=0.4, noise=0.6), "paired")
    assert scores(predictions, "sample")["mean_p_true"] > scores(predictions, "field")["mean_p_true"]


def test_the_paired_decoder_removes_a_day_offset_the_absolute_one_cannot():
    """Held-out day sits far from the training days at every dose (another transformation)."""
    f = fields(slope=1.0, day_spread=0.0, noise=0.2)
    shifted = f["session_id"] == "20261003_h2o2_tr"
    f.loc[shifted, ["green", "red"]] += 4.0
    f = ld.add_paired(f.drop(columns=[c for c in f.columns if c.startswith(("d_", "log_fold"))]))
    absolute, _ = ld.leave_one_day_out(f, "absolute")
    paired, _ = ld.leave_one_day_out(f, "paired")
    held = lambda p: p[(p["held_out"] == "20261003_h2o2_tr") & (p["unit"] == "field")]  # noqa: E731
    assert held(absolute)["outside_training"].mean() > 0.9
    assert held(paired)["p_true"].mean() > 0.8


def test_model_decoder_uses_model_means_and_is_confidently_wrong_when_the_model_is():
    f = fields(slope=1.0, noise=0.1)
    right = lambda s, d, t: np.array([0.0, np.log((5 + 3 + np.log10(max(d, 1))) / 8)])  # noqa: E731
    wrong = lambda s, d, t: np.array([0.0, 3 * np.log10(max(d, 1))])  # noqa: E731
    good, _ = ld.leave_one_day_out(f, "model", right)
    bad, _ = ld.leave_one_day_out(f, "model", wrong)
    assert scores(good)["accuracy"] > 0.9
    assert scores(bad)["info_gain_bits"] < -1.0
    with pytest.raises(ValueError, match="model_log_fold"):
        ld.leave_one_day_out(f, "model")


def test_one_training_day_leaves_the_between_day_spread_unestimated():
    _, folds = ld.leave_one_day_out(fields(days=2), "absolute")
    assert not folds["between_day_estimated"].any()


def test_a_dose_missing_from_training_skips_the_fold():
    f = fields()
    f = f[~((f["session_id"] != "20261001_h2o2_tr") & (f["dose_uM"] == 100.0))]
    _, folds = ld.leave_one_day_out(f, "absolute")
    skipped = folds.set_index("held_out")["skipped"]
    assert "100" in skipped["20261001_h2o2_tr"] and skipped["20261002_h2o2_tr"] == ""


def test_paired_reference_leaves_out_the_scored_sample():
    f = pd.DataFrame({"session_id": ["s"] * 3, "dose_uM": [0.0, 0.0, 50.0], "condition_id": ["a", "b", "c"],
                      "time_h": [1.0, 3.0, 2.0], "green": [1.0, 3.0, 10.0], "red": [0.0, 0.0, 0.0],
                      "positive": [0.1, 0.1, 0.1], "green_pos": [1.0, 3.0, 4.0], "red_pos": [1.0, 1.0, 1.0]})
    out = ld.add_paired(f).set_index("condition_id")
    assert out.loc["c", "d_green"] == pytest.approx(10 - 2)          # interpolated between 1 and 3
    assert out.loc["a", "d_green"] == pytest.approx(1 - 3)           # its own value left out
    assert out.loc["c", "log_fold_green_pos"] == pytest.approx(np.log(4 / 2))


def test_credible_set_takes_the_most_probable_doses_until_the_mass():
    assert sorted(ld.credible_set(np.array([0.1, 0.6, 0.3]))) == [1, 2]
    assert list(ld.credible_set(np.array([0.05, 0.9, 0.05]))) == [1]


def test_script_writes_outputs_and_refuses_an_existing_directory(tmp_path):
    from tests.image.test_model_compare import synthetic
    path = tmp_path / "cells_per_s.csv"
    cells = synthetic()
    other = cells[cells["session_id"] == "20261007_h2o2_tr"].assign(session_id="20261006_h2o2_tr")
    pd.concat([cells, other]).to_csv(path, index=False)
    main = runpy.run_path(str(SCRIPT))["main"]
    out = tmp_path / "lodo"
    assert main(["--cells", str(path), "--out-dir", str(out), "--decoders", "absolute", "paired"]) == 0
    for name in ("fields.csv", "predictions.csv", "folds.csv", "metrics.csv", "posterior.png", "run.json"):
        assert (out / name).exists()
    assert set(pd.read_csv(out / "predictions.csv")["held_out"]) == {"20261006_h2o2_tr", "20261007_h2o2_tr"}
    with pytest.raises(SystemExit, match="already exists"):
        main(["--cells", str(path), "--out-dir", str(out)])


def test_without_a_time_the_reference_is_the_mean_of_the_days_controls():
    f = pd.DataFrame({"session_id": ["s"] * 3, "dose_uM": [0.0, 0.0, np.nan], "condition_id": ["a", "b", "x"],
                      "time_h": [np.nan, np.nan, np.nan], "green": [1.0, 3.0, 10.0], "red": [0.0, 0.0, 0.0],
                      "positive": [0.1, 0.1, 0.1], "green_pos": [1.0, 3.0, 4.0], "red_pos": [1.0, 1.0, 1.0]})
    out = ld.add_paired(f).set_index("condition_id")
    assert out.loc["x", "d_green"] == pytest.approx(10 - 2)
