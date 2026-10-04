import pytest

from ichnos.protocol import ProtocolError, StressProtocol


@pytest.mark.parametrize("variant", ["ox", "er"])
def test_zero_dose_is_valid(variant):
    protocol = StressProtocol(variant=variant, dose=0)

    assert protocol.dose == 0.0
    assert protocol.dose_units == "uM"


@pytest.mark.parametrize(
    "changes",
    [
        {"variant": "copper"},
        {"dose": -1},
        {"dose": float("nan")},
        {"dose": float("inf")},
        {"dose": True},
        {"dose": "invalid"},
        {"dose_units": "mM"},
        {"equilibration_horizon_hours": 0},
        {"equilibration_horizon_hours": -1},
        {"equilibration_horizon_hours": float("nan")},
        {"equilibration_tolerance": 0},
        {"equilibration_tolerance": -1},
        {"equilibration_tolerance": float("inf")},
        {"clears": "yes"},
    ],
)
def test_invalid_protocol_is_rejected(changes):
    arguments = {"variant": "ox", "dose": 75}
    arguments.update(changes)

    with pytest.raises(ProtocolError):
        StressProtocol(**arguments)