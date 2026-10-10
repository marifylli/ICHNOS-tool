import numpy as np
import pandas as pd
import pytest
from PIL import Image

from ichnos_image.exposure import (
    add_per_second_columns,
    field_background,
    fit_background_vs_exposure,
)


def _cells(**overrides):
    base = {"session_id": ["s", "s"], "exposure_ms_green": [800.0, 1600.0], "exposure_ms_red": [800.0, 1600.0],
            "corrected_mean_green": [4.0, 8.0], "corrected_mean_red": [10.0, 20.0],
            "integrated_green": [40.0, 80.0], "integrated_red": [100.0, 200.0]}
    base.update(overrides)
    return pd.DataFrame(base)


def test_per_second_columns_divide_by_exposure_in_seconds():
    out = add_per_second_columns(_cells())
    assert list(out.corrected_mean_green_per_s) == [5.0, 5.0]
    assert list(out.integrated_red_per_s) == [125.0, 125.0]
    assert list(out.corrected_mean_green) == [4.0, 8.0]  # originals untouched


def test_unequal_channel_exposures_are_refused():
    with pytest.raises(ValueError, match="differ"):
        add_per_second_columns(_cells(exposure_ms_red=[800.0, 1200.0]))


def test_missing_exposure_is_refused():
    with pytest.raises(ValueError, match="positive"):
        add_per_second_columns(_cells(exposure_ms_green=[800.0, None]))


def test_linear_background_is_recovered():
    t = np.array([600, 800, 1200, 1600, 2400], float)
    fields = pd.DataFrame({"session_id": "s", "exposure_ms": t, "mode": 10 + 20 * t / 1000})
    fit = fit_background_vs_exposure(fields, "mode").iloc[0]
    assert fit.fit == "linear" and fit.offset == pytest.approx(10) and fit.rate_per_s == pytest.approx(20)
    assert fit.r2 == pytest.approx(1.0)


def test_single_exposure_session_is_reported_not_fitted():
    fields = pd.DataFrame({"session_id": ["a", "a", "b", "b"], "exposure_ms": [600, 600, 800, 1600],
                           "mode": [20.0, 22.0, 26.0, 42.0]})
    fits = fit_background_vs_exposure(fields, "mode").set_index("session_id")
    assert fits.loc["a", "fit"] == "single_exposure" and fits.loc["a", "mean"] == 21.0
    assert fits.loc["b", "fit"] == "linear" and fits.loc["all", "n_fields"] == 4


def test_field_background_on_an_8_bit_frame(tmp_path):
    image = np.full((50, 60), 30, np.uint8)
    image[:5, :5] = 255
    Image.fromarray(image).save(tmp_path / "f.tif")
    stats = field_background(tmp_path / "f.tif")
    assert abs(stats["mode"] - 30) < 1 and stats["median"] == 30
    assert stats["saturated_fraction"] == pytest.approx(25 / 3000)
