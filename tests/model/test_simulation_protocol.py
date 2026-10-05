"""Name resolution, the shared simulator, and the stress protocol."""
import io
import contextlib

import libsbml
import pytest
import numpy as np

from ichnos import build, naming, protocol as proto, simulate as sim


VARIANTS = ("ox", "er")


def _loaded(variant):
    with contextlib.redirect_stdout(io.StringIO()):
        sbml = build.build_variant_sbml_string(variant)
        model = libsbml.readSBMLFromString(sbml).getModel()
        runner = sim.load_model(sbml, model=model)
    return sbml, model, runner, naming.id_to_name_map(model)


# --- naming ---------------------------------------------------------------


def test_missing_name_raises_rather_than_returning_none():
    """A lookup that returns None means a parameter is never set and the run
    continues with the old value -- a wrong model, not an error.
    """
    _, model, _, _ = _loaded("ox")
    with pytest.raises(naming.NameNotFoundError):
        naming.resolve_id(model, "no_such_parameter")
    assert naming.try_resolve_id(model, "no_such_parameter") is None


@pytest.mark.parametrize("variant, expected", [("ox", "S_ox"), ("er", "S_er")])
def test_stress_parameter_resolves(variant, expected):
    _, model, _, _ = _loaded(variant)
    assert naming.resolve_id(model, expected, kinds=("parameter",))


# --- simulator ------------------------------------------------------------


@pytest.mark.parametrize("variant", VARIANTS)
def test_assignment_rule_observables_are_in_the_output(variant):
    """RoadRunner's default selections are floating species and rate-rule
    states only. Observed_Green and Measured_Ratio_RG are assignment-rule
    parameters -- between them the whole measured readout -- and are absent
    from a default simulate(). The shared simulator must select them.
    """
    _, model, runner, names = _loaded(variant)
    with contextlib.redirect_stdout(io.StringIO()):
        result = sim.simulate(runner, end=5, n_points=51, id_to_name=names)
    assert result.has("Observed_Green")
    assert result.has("Measured_Ratio_RG")


@pytest.mark.parametrize("variant", VARIANTS)
def test_solver_settings_travel_with_the_result(variant):
    _, _, runner, names = _loaded(variant)
    settings = sim.SolverSettings(
        absolute_tolerance=1e-9,
        relative_tolerance=1e-7,
        stiff=False,
    )

    result = sim.simulate(
        runner,
        end=1,
        n_points=11,
        id_to_name=names,
        solver=settings,
    )

    integrator = runner.getIntegrator()
    assert integrator.getValue("absolute_tolerance") == pytest.approx(
        settings.absolute_tolerance
    )
    assert integrator.getValue("relative_tolerance") == pytest.approx(
        settings.relative_tolerance
    )
    assert bool(integrator.getValue("stiff")) == settings.stiff
    assert result.solver == settings


@pytest.mark.parametrize("variant", VARIANTS)
def test_load_model_applies_all_solver_settings(variant):
    sbml, model, _, _ = _loaded(variant)
    settings = sim.SolverSettings(
        absolute_tolerance=1e-9,
        relative_tolerance=1e-7,
        stiff=False,
    )

    runner = sim.load_model(sbml, solver=settings, model=model)

    assert sim.read_solver_settings(runner) == settings


@pytest.mark.parametrize("variant", VARIANTS)
def test_simulate_without_settings_preserves_runner_settings(variant):
    sbml, model, _, names = _loaded(variant)
    settings = sim.SolverSettings(
        absolute_tolerance=1e-9,
        relative_tolerance=1e-7,
        stiff=False,
    )
    runner = sim.load_model(sbml, solver=settings, model=model)

    result = sim.simulate(
        runner,
        end=1,
        n_points=11,
        id_to_name=names,
    )

    assert sim.read_solver_settings(runner) == settings
    assert result.solver == settings


def test_peak_summary_flags_an_undersampled_peak():
    """A_ox peaks around 4 minutes. On a coarse grid the reported peak is
    whichever sample lands nearest, which was once 0.5342 at 24 minutes
    instead of 0.8622 at 4.2 -- wrong by 38% and plausible-looking.
    """
    _, model, runner, names = _loaded("ox")
    with contextlib.redirect_stdout(io.StringIO()):
        runner[naming.resolve_id(model, "S_ox", kinds=("parameter",))] = 400.0
        coarse = sim.simulate(runner, end=200, n_points=500, id_to_name=names)
    rows = {r["name"]: r for r in sim.peak_summary(coarse, ("A_ox",))}
    assert rows["A_ox"]["undersampled"] is True


# --- protocol -------------------------------------------------------------


@pytest.mark.parametrize("variant, dose", [("ox", 400.0), ("er", 2345.0)])
def test_run_protocol_equilibrates_then_applies_the_bolus(variant, dose):
    _, model, runner, names = _loaded(variant)
    protocol = proto.StressProtocol(variant=variant, dose=dose)
    with contextlib.redirect_stdout(io.StringIO()):
        result, equilibration = proto.run_protocol(
            runner, protocol, duration_hours=12, n_points=2001
        )

    assert equilibration.converged
    assert equilibration.observables_checked
    # A real baseline, not the all-zero initial condition.
    assert result["Observed_Green"][0] > 0
    # Time is reported from stress onset.
    assert result.time[0] == 0.0


def test_equilibration_failure_is_not_silent():
    """A fixed duration is not evidence of convergence. Too short a horizon
    must raise rather than hand back a baseline that is still moving.
    """
    _, model, runner, names = _loaded("ox")
    protocol = proto.StressProtocol(
        variant="ox", dose=400.0, equilibration_horizon_hours=0.01
    )
    with pytest.raises(proto.ProtocolError):
        with contextlib.redirect_stdout(io.StringIO()):
            proto.run_protocol(runner, protocol, duration_hours=1, n_points=101)


@pytest.mark.parametrize("variant", VARIANTS)
def test_stress_does_not_clear_in_the_current_model(variant):
    """S_ox and S_er are constant parameters: there is no clearance rule, so
    stress never decays. A decoder calibrated here estimates an initial bolus,
    not an arbitrary exposure history.
    """
    _, model, _, _ = _loaded(variant)
    assert not proto.has_rule_for(model, proto.STRESS_PARAMETER[variant])
    assert proto.StressProtocol(variant=variant, dose=1.0).clears is False


@pytest.mark.parametrize("variant", VARIANTS)
def test_add_clearance_refuses_a_second_rule(variant):
    """The research helper assumed the stress parameter had no rule. A second
    rule for one variable is invalid SBML and would change the input dynamics
    with no visible sign.
    """
    _, model, _, _ = _loaded(variant)
    proto.add_clearance(model, variant, 0.1)
    assert proto.has_rule_for(model, proto.STRESS_PARAMETER[variant])
    with pytest.raises(proto.DuplicateRuleError):
        proto.add_clearance(model, variant, 0.1)


def test_unsupported_variant_has_no_stress_parameter():
    with pytest.raises(proto.ProtocolError):
        proto.StressProtocol(variant="cu", dose=1.0).stress_parameter()


@pytest.mark.parametrize(
    "variant, dose, times",
    [
        ("ox", 75.0, [0.5, 1.0, 2.0, 3.0]),
        ("er", 500.0, [0.75, 2.0, 4.0]),
    ],
)
def test_protocol_returns_requested_observation_times(
    variant, dose, times
):
    _, model, runner, names = _loaded(variant)
    protocol = proto.StressProtocol(variant=variant, dose=dose)

    result, equilibration = proto.run_protocol_at_times(
        runner,
        protocol,
        times_hours=times,
        id_to_name=names,
    )

    assert equilibration.converged
    assert result.time.tolist() == pytest.approx(times)
    assert result.n_points == len(times)
    assert result.has("Observed_Green")
    assert result.has("Measured_Ratio_RG")

    stress_id = naming.resolve_id(
        model, protocol.stress_parameter(), kinds=("parameter",)
    )
    assert runner[stress_id] == pytest.approx(dose)

@pytest.mark.parametrize("variant", VARIANTS)
def test_clearance_matches_exponential_decay(variant):
    sbml, original_model, _, _ = _loaded(variant)
    protocol = proto.StressProtocol(
        variant=variant,
        dose=75,
        clears=True,
        clearance_rate_per_hour=0.5,
    )
    runner, model = proto.load_protocol_model(sbml, protocol)
    names = naming.id_to_name_map(model)

    result, equilibration = proto.run_protocol_at_times(
        runner,
        protocol,
        times_hours=[0, 0.5, 1, 2, 4],
        id_to_name=names,
    )

    assert equilibration.converged
    np.testing.assert_allclose(
        result[protocol.stress_parameter()],
        75 * np.exp(-0.5 * result.time),
        rtol=1e-6,
        atol=1e-8,
    )
    assert result.exposure == {
        "variant": variant,
        "initial_dose": 75.0,
        "dose_units": "uM",
        "time_units": "hour",
        "model": "first_order_decay",
        "clearance_rate_per_hour": 0.5,
    }

    # Preparing the exposure did not modify the original model.
    assert not proto.has_rule_for(
        original_model, protocol.stress_parameter()
    )


@pytest.mark.parametrize("variant", VARIANTS)
@pytest.mark.parametrize(
    "mismatch",
    ["missing_rule", "different_rate", "constant_requested"],
)
def test_protocol_rejects_mismatched_loaded_exposure(
    variant, mismatch
):
    sbml, _, _, _ = _loaded(variant)

    loaded_protocol = proto.StressProtocol(
        variant=variant,
        dose=75,
        clears=mismatch != "missing_rule",
        clearance_rate_per_hour=(
            None if mismatch == "missing_rule" else 0.5
        ),
    )
    requested_protocol = proto.StressProtocol(
        variant=variant,
        dose=75,
        clears=mismatch != "constant_requested",
        clearance_rate_per_hour=(
            None
            if mismatch == "constant_requested"
            else 0.2 if mismatch == "different_rate" else 0.5
        ),
    )
    runner, model = proto.load_protocol_model(
        sbml, loaded_protocol
    )

    with pytest.raises(proto.ProtocolError):
        proto.run_protocol_at_times(
            runner,
            requested_protocol,
            times_hours=[0.5, 1],
            id_to_name=naming.id_to_name_map(model),
        )