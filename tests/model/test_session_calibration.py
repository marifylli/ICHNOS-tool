import json

import numpy as np
import pytest

from ichnos.calibration import (
    calibrate_image_ratios, fit_session_calibration,
    load_session_calibration, model_reference_id, save_session_calibration,
)


def fit(x=(0.2, 0.5, 0.8), y=(0.3, 0.75, 1.2), **kwargs):
    options = dict(session_id="s1", reference_id="forward-map-1",
                   source_kind="synthetic", source="known scale 1.5")
    options.update(kwargs)
    return fit_session_calibration(x, y, **options)


def test_recovers_known_scale_and_corrects_unseen_ratios():
    calibration = fit()
    assert calibration.c_session == pytest.approx(1.5)
    assert calibration.rmse == pytest.approx(0, abs=1e-15)
    raw = np.array([0.6, 0.9])
    corrected = calibrate_image_ratios(
        raw, calibration, session_id="s1", reference_id="forward-map-1"
    )
    np.testing.assert_allclose(corrected, [0.4, 0.6])
    np.testing.assert_array_equal(raw, [0.6, 0.9])


def test_nonunit_f_is_already_in_reference_ratios():
    # Model f=2: these references are 2 * intrinsic red/green.
    model_ratio = 2 * np.array([0.2, 0.4, 0.6])
    calibration = fit(model_ratio, 1.5 * model_ratio)
    assert calibration.c_session == pytest.approx(1.5)
    corrected = calibrate_image_ratios(
        [1.5], calibration, session_id="s1", reference_id="forward-map-1"
    )
    np.testing.assert_allclose(corrected, [1.0])


def test_seeded_noise_recovery_on_independent_holdout():
    rng = np.random.default_rng(20261005)
    x = np.linspace(0.1, 2, 100)
    y = 1.5 * x * (1 + rng.uniform(-0.02, 0.02, len(x)))
    calibration = fit(x[:70], y[:70])
    assert calibration.c_session == pytest.approx(1.5, rel=0.01)
    corrected = calibrate_image_ratios(
        y[70:], calibration, session_id="s1", reference_id="forward-map-1"
    )
    assert np.max(np.abs(corrected / x[70:] - 1)) < 0.03
    assert calibration.rmse > 0


@pytest.mark.parametrize("x,y", [
    ([0, .5, 1], [.2, .3, .4]),
    ([.1, .5], [.2, .3]),
    ([.1, .5, 1], [.2, .3]),
    ([.1, np.nan, 1], [.2, .3, .4]),
    ([.1, .5, 1], [.2, np.inf, .4]),
    ([.1, .5, 1], [.2, -.3, .4]),
    ([True, .5, 1], [.2, .3, .4]),
    ([.1, .5, 1], [0, 0, 0]),
    ([[.1, .5, 1]], [[.2, .3, .4]]),
])
def test_rejects_invalid_reference_pairs(x, y):
    with pytest.raises(ValueError):
        fit(x, y)


@pytest.mark.parametrize("options", [
    {"session_id": ""}, {"reference_id": ""}, {"source": ""},
    {"source_kind": "measured"},
])
def test_requires_explicit_provenance(options):
    with pytest.raises(ValueError):
        fit(**options)


@pytest.mark.parametrize("session_id,reference_id", [
    ("s2", "forward-map-1"), ("s1", "changed-f-or-protocol"),
])
def test_rejects_wrong_session_or_forward_model(session_id, reference_id):
    with pytest.raises(ValueError):
        calibrate_image_ratios(
            [.5], fit(), session_id=session_id, reference_id=reference_id
        )


@pytest.mark.parametrize("ratios", [[np.nan], [np.inf], [-1], [True], []])
def test_rejects_invalid_measurement_ratios(ratios):
    with pytest.raises(ValueError):
        calibrate_image_ratios(
            ratios, fit(), session_id="s1", reference_id="forward-map-1"
        )


def test_json_roundtrip_preserves_synthetic_status_and_references(tmp_path):
    path = tmp_path / "calibration.json"
    calibration = fit()
    save_session_calibration(calibration, path)
    assert load_session_calibration(path) == calibration
    data = json.loads(path.read_text())
    assert data["source_kind"] == "synthetic"
    assert data["experimental_instrument_validated"] is False
    assert "not microscope calibration" in data["purpose"]
    assert data["f_already_in_model_ratio"] is True
    assert data["model_ratios"] == [.2, .5, .8]
    with pytest.raises(FileExistsError):
        save_session_calibration(calibration, path)


@pytest.mark.parametrize("key,value", [
    ("ratio_direction", "green/red"), ("c_session", 2.0),
    ("experimental_instrument_validated", True),
    ("source_kind", "experimental_reference"),
])
def test_rejects_modified_or_relabelled_artifact(tmp_path, key, value):
    path = tmp_path / "calibration.json"
    data = fit().to_dict()
    data[key] = value
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        load_session_calibration(path)


def test_reference_id_tracks_model_f_and_protocol_but_not_timestamp():
    metadata = dict(
        model_sbml_sha256="abc", parameter_profile={"f": 2},
        protocol={"dose": 75}, observation_times_hours=[0, 1, 2],
        initialization={"method": "finite"}, solver_used={"atol": 1e-10},
        fluorescence={"f": 2, "eps": 1e-9}, created_at_utc="today",
    )
    identifier = model_reference_id(metadata)
    assert model_reference_id({**metadata, "created_at_utc": "tomorrow"}) == identifier
    assert model_reference_id({**metadata, "fluorescence": {"f": 3}}) != identifier
    assert model_reference_id({**metadata, "protocol": {"dose": 50}}) != identifier
