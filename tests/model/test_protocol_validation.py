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

@pytest.mark.parametrize(
    "rate",
    [
        None,
        0,
        -1,
        float("nan"),
        float("inf"),
        True,
        "invalid",
    ],
)
def test_clearance_requires_explicit_positive_rate(rate):
    with pytest.raises(ProtocolError):
        StressProtocol(
            variant="ox",
            dose=75,
            clears=True,
            clearance_rate_per_hour=rate,
        )


def test_constant_protocol_rejects_unused_clearance_rate():
    with pytest.raises(ProtocolError):
        StressProtocol(
            variant="ox",
            dose=75,
            clearance_rate_per_hour=0.5,
        )