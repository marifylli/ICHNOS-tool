import libsbml
import numpy as np
import pytest
import roadrunner

from ichnos import simulate as sim


@pytest.fixture
def decay_runner():
    document = libsbml.SBMLDocument(3, 2)
    model = document.createModel()
    model.setId("decay")

    x = model.createParameter()
    x.setId("X")
    x.setValue(10.0)
    x.setConstant(False)

    k = model.createParameter()
    k.setId("k")
    k.setValue(0.5)
    k.setConstant(True)

    rule = model.createRateRule()
    rule.setVariable("X")
    rule.setMath(libsbml.parseL3Formula("-k * X"))

    runner = roadrunner.RoadRunner(
        libsbml.writeSBMLToString(document)
    )
    runner.selections = ["time", "X"]
    return runner


@pytest.mark.parametrize(
    "times",
    [
        [0.5, 1.0, 2.0, 3.0],
        [0.75, 2.0, 4.0],
        [0.0, 0.53, 1.08, 2.12, 3.04],
    ],
)
def test_observations_match_analytical_solution(decay_runner, times):
    settings = sim.SolverSettings()
    result = sim.simulate_at_times(
        decay_runner,
        times_hours=times,
        solver=settings,
    )

    np.testing.assert_allclose(
        result.time, times, rtol=0, atol=1e-12
    )
    np.testing.assert_allclose(
        result["X"],
        10.0 * np.exp(-0.5 * np.asarray(times)),
        rtol=1e-6,
        atol=1e-8,
    )
    assert result.solver == settings
    assert result.n_points == len(times)


def test_continuation_preserves_state_and_uses_onset_time(decay_runner):
    decay_runner.simulate(0, 50, 2)
    decay_runner["X"] = 20.0
    times = [0.0, 0.5, 1.0]

    result = sim.simulate_at_times(
        decay_runner,
        times_hours=times,
        solver=sim.SolverSettings(),
        reset=False,
    )

    np.testing.assert_allclose(result.time, times)
    np.testing.assert_allclose(
        result["X"],
        20.0 * np.exp(-0.5 * np.asarray(times)),
        rtol=1e-6,
        atol=1e-8,
    )


@pytest.mark.parametrize(
    "times",
    [
        [],
        [0],
        [-0.5, 1],
        [1, 0.5],
        [1, 1],
        [float("nan"), 1],
        [1, float("inf")],
        [True, 1],
        [[0.5, 1]],
        ["invalid"],
    ],
)
def test_invalid_times_are_rejected_before_runner_access(times):
    with pytest.raises(ValueError):
        sim.simulate_at_times(object(), times_hours=times)