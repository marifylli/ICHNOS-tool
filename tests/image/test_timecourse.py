"""Autofluorescence baseline, population time courses and peak timing."""
import runpy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ichnos_image import timecourse as tc

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "fluorescence_timecourse.py"
TIMES = [0.5, 1.5, 2.5, 3.5, 5.0]
LABELS = ["0h", "1h", "2h", "3h", "4h"]


def response(t, peak=3.3, width=1.2, amplitude=40.0):
    return amplitude * np.exp(-((t - peak) ** 2) / (2 * width**2))


def synthetic(seed=0, cells_per_field=60, base=(5.0, 10.0)):
    """Two transfected days and one non-transfected day of one stressor.

    Every sample keeps its sample_id across time points, as the 5/10 series
    does, so a field is only identified together with its condition.
    """
    rng = np.random.default_rng(seed)
    rows = []
    def add(session, dose, label, time, rep, green, red):
        for k in range(cells_per_field):
            rows.append(dict(
                session_id=session, condition_id=f"h2o2_{dose}uM_{label}", sample_id=f"s_{dose}uM_R{rep}",
                cell_id=k, qc_pass=True, edge_flag=False, measurement_time_hours=time,
                corrected_mean_green_per_s=green + rng.normal(0, 2.0),
                corrected_mean_red_per_s=red + rng.normal(0, 2.0)))
    for session in ("20261005_h2o2_tr", "20261007_h2o2_tr"):
        for label, t in zip(LABELS, TIMES):
            for rep in (1, 2):
                g = base[0] + response(t)
                add(session, 75, label, t + 0.05 * rep, rep, g, base[1] + 0.5 * response(t, peak=4.3))
    for label, t in zip(LABELS, TIMES):
        add("20261005_h2o2_tr", 0, label, t + 0.2, 1, base[0] + 1.0, base[1])   # unstressed control
        for rep in (1, 2):
            add("20261008_h2o2_nt", 75, label, t, rep, base[0], base[1])
    cells = pd.DataFrame(rows)
    # an edge cell and a failed one, both very bright: must be ignored
    bad = cells.iloc[[0, 1]].copy()
    bad["corrected_mean_green_per_s"] = 1e4
    bad["edge_flag"] = [True, False]
    bad["qc_pass"] = [True, False]
    return pd.concat([cells, bad], ignore_index=True)


def test_annotate_parses_and_refuses_what_it_cannot_read():
    cells = tc.annotate(synthetic())
    assert set(cells["strain"]) == {"tr", "nt"} and set(cells["dose_uM"]) == {0.0, 75.0}
    medium = pd.DataFrame({"session_id": ["20261005_none_tr"], "condition_id": ["medium_medium"]})
    assert tc.annotate(medium)["dose_uM"].isna().all()
    with pytest.raises(ValueError, match="condition_id"):
        tc.annotate(pd.DataFrame({"session_id": ["20261005_h2o2_tr"], "condition_id": ["dtt_10uM_1h"]}))
    with pytest.raises(ValueError, match="session_id"):
        tc.annotate(pd.DataFrame({"session_id": ["h2o2"], "condition_id": ["h2o2_10uM_1h"]}))


def test_baseline_comes_from_non_transfected_cells_and_is_subtracted():
    cells = tc.annotate(synthetic())
    baseline = tc.autofluorescence_baseline(cells)
    row = baseline.set_index("stressor").loc["h2o2"]
    assert row["baseline_green"] == pytest.approx(5.0, abs=0.3)
    assert row["baseline_red"] == pytest.approx(10.0, abs=0.3)
    out = tc.subtract_baseline(cells, baseline)
    nt = out[(out["strain"] == "nt") & tc.usable(out)]
    assert nt["green_af"].median() == pytest.approx(0.0, abs=0.3)


def test_fields_are_not_merged_across_time_points():
    cells = tc.subtract_baseline(tc.annotate(synthetic()), tc.autofluorescence_baseline(tc.annotate(synthetic())))
    fields = tc.field_summaries(cells)
    one_day = fields[(fields["session_id"] == "20261005_h2o2_tr") & (fields["dose_uM"] == 75)]
    assert len(one_day) == len(TIMES) * 2          # same sample_ids, five time points
    assert (one_day["n_cells"] == 60).all()        # the edge and failed cells were dropped


@pytest.mark.parametrize("vertex", [2.2, 2.9, 3.4])
def test_parabola_recovers_a_peak_between_samples(vertex):
    t = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    v = -(t - vertex) ** 2 + 10
    peak = tc.refine_peak(t, v)
    assert peak["bracketed"] and peak["t_peak"] == pytest.approx(vertex)


def test_a_peak_at_the_edge_is_reported_as_a_bound():
    rising = tc.refine_peak([1, 2, 3], [1, 2, 3])
    falling = tc.refine_peak([1, 2, 3], [3, 2, 1])
    assert (rising["position"], rising["bracketed"]) == ("last_sample", False)
    assert falling["position"] == "first_sample"


def test_peak_time_and_interval_on_a_known_curve():
    cells = tc.annotate(synthetic())
    cells = tc.subtract_baseline(cells, tc.autofluorescence_baseline(cells))
    fields = tc.field_summaries(cells)
    curves = tc.timecourse(fields)
    curve = curves[(curves["session_id"] == "20261005_h2o2_tr") & (curves["dose_uM"] == 75)]
    peaks = tc.bootstrap_peaks(cells, curve, fields, n_boot=200, seed=1).set_index("channel")
    green, red = peaks.loc["green"], peaks.loc["red"]
    assert green["bracketed"] and 2.9 < green["t_peak"] < 3.7
    assert green["t_peak_p05"] <= green["t_peak"] <= green["t_peak_p95"]
    assert green["share_at_sampled_max"] > 0.9
    assert red["t_peak"] > green["t_peak"]          # the slower reporter peaks later
    assert set(peaks.index) == {"green", "red"}       # no "peak" of the timer ratio


def test_excess_over_the_unstressed_control_at_the_peak_time():
    cells = tc.annotate(synthetic())
    cells = tc.subtract_baseline(cells, tc.autofluorescence_baseline(cells))
    curves = tc.timecourse(tc.field_summaries(cells))
    stressed = curves[(curves["session_id"] == "20261005_h2o2_tr") & (curves["dose_uM"] == 75)]
    t_peak = float(stressed.loc[stressed["green"].idxmax(), "time_h"])
    peaks = pd.DataFrame([
        {"session_id": "20261005_h2o2_tr", "dose_uM": 75.0, "channel": "green", "t_max_sampled": t_peak},
        {"session_id": "20261005_h2o2_tr", "dose_uM": 0.0, "channel": "green", "t_max_sampled": 3.7},
        {"session_id": "20261007_h2o2_tr", "dose_uM": 75.0, "channel": "green", "t_max_sampled": 3.55}])
    excess = tc.excess_over_control(peaks, curves)
    assert excess[0] == pytest.approx(response(t_peak) - 1.0, abs=1.5)   # control sits 1 above baseline
    assert np.isnan(excess[1]) and np.isnan(excess[2])              # the control itself; a day without one


def test_script_writes_its_outputs_and_refuses_an_existing_directory(tmp_path):
    path = tmp_path / "cells_per_s.csv"
    synthetic().to_csv(path, index=False)
    main = runpy.run_path(str(SCRIPT))["main"]
    out = tmp_path / "tc"
    assert main(["--cells", str(path), "--out-dir", str(out), "--n-boot", "50"]) == 0
    for name in ("baseline.csv", "baseline_diagnostics.csv", "cells_af.csv", "fields.csv",
                 "timecourse.csv", "peaks.csv", "timecourse_h2o2.png", "run.json"):
        assert (out / name).exists(), name
    peaks = pd.read_csv(out / "peaks.csv")
    assert set(peaks["day"].astype(str)) == {"20261005", "20261007"}   # transfected days only
    with pytest.raises(SystemExit, match="already exists"):
        main(["--cells", str(path), "--out-dir", str(out)])


def test_script_asks_for_exposure_normalisation_first(tmp_path):
    path = tmp_path / "cells.csv"
    synthetic().drop(columns=["corrected_mean_red_per_s"]).to_csv(path, index=False)
    main = runpy.run_path(str(SCRIPT))["main"]
    with pytest.raises(SystemExit, match="normalize_exposure"):
        main(["--cells", str(path), "--out-dir", str(tmp_path / "x")])
