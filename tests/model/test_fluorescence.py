import numpy as np
import pytest

from ichnos.fluorescence import (
    RATIO_DIRECTION,
    model_fluorescence,
)
from ichnos.simulate import SimulationResult, SolverSettings


def make_result():
    return SimulationResult(
        time=np.array([0.0, 1.0]),
        columns=[
            "Observed_Green",
            "Reporter_red",
            "Measured_Ratio_RG",
            "Total_red_pool",
        ],
        values=np.array([
            [10.0, 5.0, 1.0, 100.0],
            [20.0, 15.0, 1.5, 200.0],
        ]),
        solver=SolverSettings(),
    )


def test_maps_mature_red_and_preserves_scaled_model_ratio():
    signals = model_fluorescence(make_result())

    assert RATIO_DIRECTION == "red/green"
    np.testing.assert_array_equal(signals["green"], [10.0, 20.0])
    np.testing.assert_array_equal(signals["red"], [5.0, 15.0])
    np.testing.assert_array_equal(
        signals["ratio_red_green"], [1.0, 1.5]
    )


def test_missing_mature_red_does_not_fall_back_to_total_pool():
    result = make_result()
    result.columns[1] = "another_observable"

    with pytest.raises(ValueError, match="Reporter_red"):
        model_fluorescence(result)


@pytest.mark.parametrize("value", [np.nan, np.inf, -1.0])
@pytest.mark.parametrize("column", [0, 1, 2])
def test_rejects_invalid_signals(value, column):
    result = make_result()
    result.values[0, column] = value

    with pytest.raises(ValueError, match="finite non-negative"):
        model_fluorescence(result)


def test_returned_signals_do_not_modify_simulation():
    result = make_result()
    signals = model_fluorescence(result)
    signals["green"][0] = 999.0

    assert result["Observed_Green"][0] == 10.0