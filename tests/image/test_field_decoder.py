"""Decoding a new culture against a reference experiment, with a verdict in words."""
import json
import runpy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ichnos_image import field_decoder as fd

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "decode_field.py"
DOSES = (0, 100, 1000)
TIMES = (1.0, 2.0, 3.0)


def _cells(rng, session, condition, sample, t, green, red, n=40):
    return [dict(session_id=session, condition_id=condition, sample_id=sample, cell_id=k, qc_pass=True,
                 edge_flag=False, measurement_time_hours=t,
                 corrected_mean_green_per_s=green + rng.normal(0, 0.8),
                 corrected_mean_red_per_s=red + rng.normal(0, 0.8)) for k in range(n)]


def reference(slope=4.0, seed=0, days=3):
    """Transfected days with red rising by ``slope`` per 10-fold dose, plus a non-transfected day."""
    rng = np.random.default_rng(seed)
    rows = []
    for d in range(days):
        session, offset = f"2026100{d + 1}_h2o2_tr", rng.normal(0, 0.3)
        for dose in DOSES:
            lift = slope * np.log10(max(dose, 1))
            for t in TIMES:
                for r in range(3):
                    rows += _cells(rng, session, f"h2o2_{dose}uM_{t:g}h", f"d{dose}_{t:g}_R{r}", t + 0.01 * r,
                                   12 + offset, 12 + offset + lift)
    for t in TIMES:
        rows += _cells(rng, "20261009_h2o2_nt", f"h2o2_0uM_{t:g}h", f"nt_{t:g}", t, 5, 5)
    return pd.DataFrame(rows)


def new_culture(dose=1000, slope=4.0, shift=0.0, with_control=True, seed=7):
    rng = np.random.default_rng(seed)
    rows = []
    for t in (1.5, 2.5):
        if with_control:
            for r in range(2):
                rows += _cells(rng, "my_day", "h2o2_0uM_ctrl", f"ctrl_{t:g}_{r}", t + 0.01 * r, 12 + shift, 12 + shift)
        for r in range(3):
            rows += _cells(rng, "my_day", "mystery", f"m_{t:g}_{r}", t + 0.02 * r, 12 + shift,
                           12 + shift + slope * np.log10(max(dose, 1)))
    return pd.DataFrame(rows)


def mystery(decoded, unit="field"):
    return decoded[(decoded["condition_id"] == "mystery") & (decoded["unit"] == unit)]


def test_a_validated_decoder_names_the_dose_of_a_new_culture():
    trained = fd.train(reference(), "h2o2", "paired")
    assert trained["validated_accuracy"] > 0.9 and trained["validated_gain_bits"] > 1
    out = mystery(fd.decode(new_culture(1000), trained))
    assert (out["verdict"] == "dose estimate").all()
    assert (out["map_dose_uM"] == 1000).all()
    assert out["known_dose_uM"].isna().all()   # 'mystery' carries no dose


def test_without_a_dose_signal_no_dose_is_claimed():
    trained = fd.train(reference(slope=0.0), "h2o2", "paired")
    out = mystery(fd.decode(new_culture(1000, slope=0.0), trained))
    assert not (out["verdict"] == "dose estimate").any()
    assert set(out["verdict"]) <= {"no dose information", "not validated"}


def test_a_culture_unlike_the_reference_is_flagged_not_decoded():
    trained = fd.train(reference(), "h2o2", "absolute")
    out = mystery(fd.decode(new_culture(100, shift=25.0), trained))
    assert (out["verdict"] == "unlike reference").all()


def test_the_paired_decoder_asks_for_a_same_day_control():
    trained = fd.train(reference(), "h2o2", "paired")
    out = mystery(fd.decode(new_culture(100, with_control=False), trained))
    assert (out["verdict"] == "needs a 0 µM image").all()


def test_too_few_cells_is_reported():
    trained = fd.train(reference(), "h2o2", "paired", min_cells=20)
    small = new_culture(1000)
    small = small[~((small["condition_id"] == "mystery") & (small["cell_id"] >= 5))]
    out = mystery(fd.decode(small, trained))
    assert (out["verdict"] == "too few cells").all()


def test_the_model_decoder_is_not_offered_and_a_bad_table_is_explained():
    with pytest.raises(ValueError, match="not offered"):
        fd.train(reference(), "h2o2", "model")
    trained = fd.train(reference(), "h2o2", "paired")
    with pytest.raises(ValueError, match="normalize_exposure"):
        fd.decode(new_culture().drop(columns=["corrected_mean_red_per_s"]), trained)


def test_excluding_a_reference_session_removes_it_from_training():
    trained = fd.train(reference(), "h2o2", "paired", exclude_sessions=("20261003_h2o2_tr",))
    assert trained["reference_days"] == ["20261001_h2o2_tr", "20261002_h2o2_tr"]


def test_script_writes_outputs(tmp_path):
    ref, new = tmp_path / "ref.csv", tmp_path / "new.csv"
    reference().to_csv(ref, index=False)
    new_culture(100).to_csv(new, index=False)
    main = runpy.run_path(str(SCRIPT))["main"]
    out = tmp_path / "decoded"
    assert main(["--reference-cells", str(ref), "--cells", str(new), "--stressor", "h2o2", "--out-dir", str(out)]) == 0
    decoded = pd.read_csv(out / "decoded.csv")
    assert (decoded.loc[decoded["condition_id"] == "mystery", "map_dose_uM"] == 100).all()
    summary = json.loads((out / "summary.json").read_text())
    assert summary["reference_days"] == ["20261001_h2o2_tr", "20261002_h2o2_tr", "20261003_h2o2_tr"]
    assert (out / "validation.csv").exists()
    with pytest.raises(SystemExit, match="already exists"):
        main(["--reference-cells", str(ref), "--cells", str(new), "--stressor", "h2o2", "--out-dir", str(out)])


def test_worked_example_hides_doses_and_decodes_back(tmp_path):
    hide = runpy.run_path(str(SCRIPT.parent / "make_worked_example.py"))["hide_doses"]
    ref = reference()
    cells, key = hide(ref, "20261003_h2o2_tr")
    stressed = ~cells["condition_id"].str.startswith("h2o2_0uM_")
    assert not cells.loc[stressed, "condition_id"].str.contains("uM").any()
    assert not cells.loc[stressed, "sample_id"].str.contains(r"d\d", regex=True).any()
    assert set(key["dose_uM"]) == {100.0, 1000.0} and (cells["session_id"] == "user_day").all()
    trained = fd.train(ref, "h2o2", "paired", exclude_sessions=("20261003_h2o2_tr",))
    decoded = fd.decode(cells, trained)
    fields = decoded[(decoded["unit"] == "field") & decoded["condition_id"].str.startswith("sample")]
    truth = fields["condition_id"].map(key.set_index("sample")["dose_uM"])
    assert (fields["map_dose_uM"] == truth).mean() > 0.9
