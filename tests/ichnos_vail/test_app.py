"""The Streamlit app itself: every page renders, and the decode page runs a real job end to end."""
import json
import time
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

from tests.ichnos_vail.test_manifest_check import GOOD, write  # noqa: E402
from tests.image.test_field_decoder import new_culture, reference, timed_reference, untimed_culture  # noqa: E402

APP = Path(__file__).resolve().parents[2] / "app" / "ichnos_vail.py"
PAGES = ["Start", "1 · Check images", "2 · Images → cells", "3 · Analyse an experiment",
         "4 · Decode a new culture", "Jobs"]


@pytest.fixture
def app(tmp_path, monkeypatch):
    monkeypatch.setenv("ICHNOS_VAIL_JOBS_DIR", str(tmp_path / "jobs"))
    at = AppTest.from_file(str(APP), default_timeout=60)
    at.run()
    assert not at.exception
    return at


def go(at, page):
    at.sidebar.radio(key="page").set_value(page).run()
    assert not at.exception, at.exception


def button(at, label):
    return next(b for b in at.button if b.label == label)


def texts(at):
    return " ".join([m.value for m in at.markdown] + [s.value for s in at.success] + [w.value for w in at.warning]
                    + [e.value for e in at.error] + [c.value for c in at.caption])


def test_every_page_renders(app):
    for page in PAGES:
        go(app, page)


def test_the_manifest_page_checks_a_manifest(app, tmp_path):
    path = write(tmp_path, GOOD)
    go(app, "1 · Check images")
    app.text_input(key="w_manifest_path").input(str(path)).run()
    button(app, "Check").click().run()
    assert not app.exception
    values = {m.label: m.value for m in app.metric}
    assert values["Fields"] == "3" and values["Errors"] == "0"
    assert "No problems found." in texts(app)
    assert "python3 scripts/align_brightfield.py" in " ".join(c.value for c in app.code)


def test_a_start_with_bad_inputs_shows_errors_and_starts_nothing(app, tmp_path):
    go(app, "4 · Decode a new culture")
    ref = tmp_path / "ref.csv"
    reference().to_csv(ref, index=False)
    app.text_input(key="w_decode_reference").input(str(ref)).run()
    app.text_input(key="w_decode_cells").input(str(tmp_path / "missing.csv")).run()
    button(app, "Decode").click().run()
    assert "new cell table not found" in texts(app)
    assert not (tmp_path / "jobs").exists() or not any((tmp_path / "jobs").iterdir())


def test_decoding_a_new_culture_from_the_gui(app, tmp_path):
    ref, new, out = tmp_path / "ref.csv", tmp_path / "new.csv", tmp_path / "decoded"
    reference().to_csv(ref, index=False)
    new_culture(1000).to_csv(new, index=False)
    go(app, "4 · Decode a new culture")
    app.text_input(key="w_decode_reference").input(str(ref)).run()
    app.text_input(key="w_decode_cells").input(str(new)).run()
    app.text_input(key="w_decode_out").input(str(out)).run()
    assert "python3 scripts/decode_field.py" in " ".join(c.value for c in app.code)
    button(app, "Decode").click().run()
    assert not app.exception
    job = next((tmp_path / "jobs").iterdir())
    end = time.time() + 120
    while json.loads((job / "status.json").read_text()).get("state") not in ("finished", "failed"):
        assert time.time() < end, (job / "log.txt").read_text()
        time.sleep(0.5)
    assert json.loads((job / "status.json").read_text())["state"] == "finished", (job / "log.txt").read_text()
    app.run()
    shown = texts(app)
    assert "Track record" in shown and "better than chance" in shown
    assert "dose estimate" in shown
    decoded = pd.read_csv(out / "decoded.csv")
    assert (decoded.loc[decoded["condition_id"] == "mystery", "map_dose_uM"] == 1000).all()


def test_the_decode_page_shows_time_since_stress(app, tmp_path):
    ref, new, out = tmp_path / "ref.csv", tmp_path / "new.csv", tmp_path / "decoded"
    timed_reference().to_csv(ref, index=False)
    untimed_culture(5.0).to_csv(new, index=False)
    go(app, "4 · Decode a new culture")
    app.text_input(key="w_decode_reference").input(str(ref)).run()
    app.text_input(key="w_decode_cells").input(str(new)).run()
    app.text_input(key="w_decode_out").input(str(out)).run()
    button(app, "Decode").click().run()
    job = next((tmp_path / "jobs").iterdir())
    end = time.time() + 180
    while json.loads((job / "status.json").read_text()).get("state") not in ("finished", "failed"):
        assert time.time() < end, (job / "log.txt").read_text()
        time.sleep(0.5)
    assert json.loads((job / "status.json").read_text())["state"] == "finished", (job / "log.txt").read_text()
    app.run()
    shown = texts(app)
    assert "Time since the stress" in shown and "time estimate" in shown
    decoded = pd.read_csv(out / "decoded.csv")
    assert (decoded.loc[decoded["condition_id"] == "mystery", "time_map"] == "4–6 h").all()
