import numpy as np
import pandas as pd
import pytest

from ichnos.calibration import fit_session_calibration
from ichnos.population import CalibrationBinding, summarize_cells


def cell_frame(session="s1", sample="sample-a", order=1):
    return pd.DataFrame({
        "session_id": [session] * 4, "sample_id": [sample] * 4,
        "condition_id": ["dose75"] * 4, "timepoint": [1] * 4,
        "cell_id": [1, 2, 3, 4], "acquisition_order": [order] * 4,
        "qc_pass": [True, True, True, False],
        "corrected_mean_green": [100., 200., 50., 10.],
        "corrected_mean_red": [60., 180., 60., 1000.],
        "ratio_red_green": [.6, .9, 1.2, 100.],
        "sampling_time_hours": [.5] * 4, "measurement_time_hours": [.75] * 4,
        "exposure_ms_green": [100.] * 4, "exposure_ms_red": [200.] * 4,
        "nd_filter_green": [0.] * 4, "nd_filter_red": [0.] * 4,
        "objective": ["40X"] * 4,
    })


def bindings(session="s1", kind="synthetic"):
    calibration = fit_session_calibration(
        [.2, .4, .6], [.3, .6, .9], session_id=session,
        reference_id="reference", source_kind=kind, source="test reference",
    )
    return {session: CalibrationBinding(calibration, "reference")}


def test_qc_filtered_median_and_calibration_preserve_original_data():
    frame = cell_frame()
    original = frame.copy(deep=True)
    row = summarize_cells(frame, data_kind="synthetic", calibrations=bindings()).iloc[0]
    assert row.n_cells_total == 4
    assert row.n_cells_qc_pass == row.n_cells_used == 3
    assert row.ratio_red_green_median == pytest.approx(.9)
    assert row.calibrated_ratio_red_green_median == pytest.approx(.6)
    assert row.ratio_red_green_q25 == pytest.approx(.75)
    assert row.calibrated_ratio_red_green_q25 == pytest.approx(.5)
    assert row.sampling_time_hours == .5
    assert row.measurement_time_hours == .75
    pd.testing.assert_frame_equal(frame, original)


def test_separates_samples_and_sessions_but_combines_fields_of_view():
    frame = pd.concat([
        cell_frame(), cell_frame(order=2), cell_frame(sample="sample-b"),
        cell_frame(session="s2"),
    ], ignore_index=True)
    output = summarize_cells(frame, data_kind="synthetic")
    assert len(output) == 3
    first = output.iloc[0]
    assert first.n_cells_total == 8
    assert first.n_cells_used == 6
    assert first.n_images == 2
    assert output.c_session.isna().all()
    assert output.calibrated_ratio_red_green_median.isna().all()


def test_unknown_times_remain_unknown():
    frame = cell_frame()
    frame["sampling_time_hours"] = None
    frame["measurement_time_hours"] = None
    row = summarize_cells(frame, data_kind="synthetic").iloc[0]
    assert pd.isna(row.sampling_time_hours)
    assert pd.isna(row.measurement_time_hours)


@pytest.mark.parametrize("invalid", [0., np.nan])
def test_low_signal_cells_are_excluded_and_small_groups_retained(invalid):
    frame = cell_frame()
    frame.loc[0, "corrected_mean_green"] = invalid
    row = summarize_cells(frame, data_kind="synthetic", calibrations=bindings()).iloc[0]
    assert row.n_cells_total == 4
    assert row.n_cells_used == 2
    assert pd.isna(row.ratio_red_green_median)
    assert pd.isna(row.calibrated_ratio_red_green_median)


def test_all_qc_failed_group_is_not_dropped():
    frame = cell_frame()
    frame["qc_pass"] = False
    row = summarize_cells(frame, data_kind="synthetic").iloc[0]
    assert row.n_cells_used == 0
    assert pd.isna(row.ratio_red_green_median)


@pytest.mark.parametrize("field", ["session_id", "sample_id"])
def test_requires_explicit_sample_and_session(field):
    frame = cell_frame()
    frame.loc[0, field] = None
    with pytest.raises(ValueError, match=field):
        summarize_cells(frame, data_kind="synthetic")


@pytest.mark.parametrize("field", ["measurement_time_hours", "exposure_ms_red", "objective"])
def test_refuses_inconsistent_acquisition_within_group(field):
    frame = cell_frame()
    frame.loc[0, field] = "60X" if field == "objective" else 999
    with pytest.raises(ValueError, match=field):
        summarize_cells(frame, data_kind="synthetic")


def test_rejects_ratio_with_wrong_direction():
    frame = cell_frame()
    frame["ratio_red_green"] = 1 / frame.ratio_red_green
    with pytest.raises(ValueError, match="differs from red/green"):
        summarize_cells(frame, data_kind="synthetic")


def test_rejects_duplicate_cells():
    with pytest.raises(ValueError, match="duplicate cell"):
        summarize_cells(pd.concat([cell_frame(), cell_frame()]), data_kind="synthetic")


@pytest.mark.parametrize("value", ["False", 1])
def test_rejects_truthy_nonboolean_qc(value):
    frame = cell_frame()
    frame["qc_pass"] = value
    with pytest.raises(ValueError, match="booleans"):
        summarize_cells(frame, data_kind="synthetic")


def test_synthetic_calibration_cannot_be_used_on_real_images():
    with pytest.raises(ValueError, match="synthetic calibration"):
        summarize_cells(cell_frame(), data_kind="experimental", calibrations=bindings())


def test_experimental_references_are_supported_without_claiming_validation():
    row = summarize_cells(
        cell_frame(), data_kind="experimental",
        calibrations=bindings(kind="experimental_reference"),
    ).iloc[0]
    assert row.calibration_source_kind == "experimental_reference"


def test_refuses_partial_session_coverage():
    frame = pd.concat([cell_frame(), cell_frame(session="s2")])
    with pytest.raises(ValueError, match="missing calibration"):
        summarize_cells(frame, data_kind="synthetic", calibrations=bindings())


def test_refuses_mismatched_forward_reference():
    binding = bindings()["s1"]
    with pytest.raises(ValueError, match="mismatch"):
        summarize_cells(cell_frame(), data_kind="synthetic", calibrations={
            "s1": CalibrationBinding(binding.calibration, "different-model")
        })


def test_refuses_changed_settings_between_samples_in_one_calibrated_session():
    second = cell_frame(sample="sample-b")
    second["exposure_ms_red"] = 400
    frame = pd.concat([cell_frame(), second])
    with pytest.raises(ValueError, match="exposure_ms_red"):
        summarize_cells(frame, data_kind="synthetic", calibrations=bindings())
