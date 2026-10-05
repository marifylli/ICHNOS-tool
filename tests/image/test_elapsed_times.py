import numpy as np
import pandas as pd
import pytest

from ichnos_image import ImageSet
from ichnos_image.export import build_records, export_csv
from ichnos_image.extract import CellFeatures


def _image_set(**timing):
    return ImageSet(
        green=np.zeros((8, 8)),
        red=np.zeros((8, 8)),
        session_id="session-1",
        timepoint=1,
        acquisition_order=1,
        exposure_ms_green=100.0,
        exposure_ms_red=100.0,
        nd_filter_green=0.0,
        nd_filter_red=0.0,
        objective="60X",
        burner_hours=5.0,
        lamp_warmup_minutes=30.0,
        **timing,
    )


def _records(image_set):
    feature = CellFeatures(
        cell_id=1,
        area_px=16,
        raw_mean_green=100.0,
        raw_mean_red=50.0,
        corrected_mean_green=90.0,
        corrected_mean_red=40.0,
        integrated_green=1440.0,
        integrated_red=640.0,
        sat_flag=False,
        local_background_green=0.0,
        local_background_red=0.0,
    )
    return build_records(
        [feature],
        session_id=image_set.session_id,
        timepoint=image_set.timepoint,
        edge_flagged_ids=set(),
        focus_score=100.0,
        registration_shift_px=0.0,
        exposure_ms_green=image_set.exposure_ms_green,
        exposure_ms_red=image_set.exposure_ms_red,
        nd_filter_green=image_set.nd_filter_green,
        nd_filter_red=image_set.nd_filter_red,
        objective=image_set.objective,
        burner_hours=image_set.burner_hours,
        lamp_warmup_minutes=image_set.lamp_warmup_minutes,
        acquisition_order=image_set.acquisition_order,
        sampling_time_hours=image_set.sampling_time_hours,
        measurement_time_hours=image_set.measurement_time_hours,
    )


def test_fractional_times_remain_distinct_in_csv(tmp_path):
    image_set = _image_set(
        sampling_time_hours=0.5,
        measurement_time_hours=0.75,
    )
    path = export_csv(_records(image_set), tmp_path / "cells.csv")
    row = pd.read_csv(path).iloc[0]

    assert row["timepoint"] == 1
    assert row["sampling_time_hours"] == 0.5
    assert row["measurement_time_hours"] == 0.75


def test_missing_times_remain_unknown_in_csv(tmp_path):
    image_set = _image_set()
    path = export_csv(_records(image_set), tmp_path / "cells.csv")
    row = pd.read_csv(path).iloc[0]

    assert pd.isna(row["sampling_time_hours"])
    assert pd.isna(row["measurement_time_hours"])


@pytest.mark.parametrize(
    "field",
    ["sampling_time_hours", "measurement_time_hours"],
)
@pytest.mark.parametrize(
    "value",
    [-1.0, float("nan"), float("inf"), True],
)
def test_invalid_elapsed_times_are_rejected(field, value):
    with pytest.raises(ValueError, match=field):
        _image_set(**{field: value})