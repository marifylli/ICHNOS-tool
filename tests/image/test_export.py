"""Stage 7/8: build_records folds per-image QC checks (focus, registration
shift, lamp warmup) plus per-cell checks (saturation, edge-touching) into a
single qc_pass column, without ever dropping rows itself.
"""
from ichnos_image import export, extract


def _one_cell_features(sat_flag=False):
    return [
        extract.CellFeatures(
            cell_id=1,
            area_px=500,
            raw_mean_green=1000.0,
            raw_mean_red=200.0,
            corrected_mean_green=900.0,
            corrected_mean_red=150.0,
            integrated_green=450000.0,
            integrated_red=75000.0,
            sat_flag=sat_flag,
            local_background_green=10.0,
            local_background_red=5.0,
        )
    ]


def _build(**overrides):
    kwargs = dict(
        session_id="s1",
        timepoint=0,
        edge_flagged_ids=set(),
        focus_score=8000.0,
        registration_shift_px=0.5,
        exposure_ms_green=800.0,
        exposure_ms_red=2000.0,
        nd_filter_green=0.0,
        nd_filter_red=0.0,
        objective="60x",
        burner_hours=5.0,
        lamp_warmup_minutes=30.0,
        acquisition_order=1,
    )
    kwargs.update(overrides)
    return export.build_records(_one_cell_features(), **kwargs)


def test_qc_pass_when_everything_within_thresholds():
    records = _build()
    assert records[0].qc_pass is True
    assert records[0].lamp_flag is False


def test_qc_fails_on_lamp_warmup():
    records = _build(lamp_warmup_minutes=5.0)  # below the 15-minute protocol threshold
    assert records[0].lamp_flag is True
    assert records[0].qc_pass is False


def test_qc_fails_on_low_focus_score():
    records = _build(focus_score=1e-10)  # far below the calibrated (normalized) threshold
    assert records[0].qc_pass is False


def test_qc_fails_on_large_registration_shift():
    records = _build(registration_shift_px=10.0)  # far above the calibrated 3.5px threshold
    assert records[0].qc_pass is False


def test_qc_fails_on_saturated_cell():
    records = export.build_records(
        _one_cell_features(sat_flag=True),
        session_id="s1", timepoint=0, edge_flagged_ids=set(),
        focus_score=8000.0, registration_shift_px=0.5,
        exposure_ms_green=800.0, exposure_ms_red=2000.0,
        nd_filter_green=0.0, nd_filter_red=0.0, objective="60x",
        burner_hours=5.0, lamp_warmup_minutes=30.0, acquisition_order=1,
    )
    assert records[0].qc_pass is False


def test_qc_fails_on_edge_touching_cell():
    records = _build(edge_flagged_ids={1})
    assert records[0].edge_flag is True
    assert records[0].qc_pass is False


def test_build_records_never_drops_rows_regardless_of_qc():
    records = _build(lamp_warmup_minutes=1.0, focus_score=0.0, registration_shift_px=99.0)
    assert len(records) == 1
    assert records[0].qc_pass is False
