"""Stress exposure and model initialization before stress onset.

Initialization starts from the source SBML initial conditions.
The default workflow checks convergence of selected readouts during
a zero-stress interval. This is a computational equilibrium assumption,
not evidence of an experimentally equilibrated culture.

Finite preincubation is available as a separate assumed scenario.
The model does not explicitly represent glucose-to-galactose switching.

Stress is constant unless an explicit first-order clearance rate is
supplied. Post-stress time is reported in hours from stress onset.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Optional

import libsbml
import numpy as np
import math

from . import naming, simulate as sim


#: The stress input parameter for each variant.
STRESS_PARAMETER = {"ox": "S_ox", "er": "S_er"}


class ProtocolError(ValueError):
    pass


class DuplicateRuleError(ProtocolError):
    """A rule already governs this quantity; adding another would be invalid
    SBML and would silently change the model."""


@dataclass(frozen=True)
class Equilibration:
    """How the baseline was reached, and whether it actually converged."""

    horizon_hours: float
    converged: bool
    criterion: str
    max_relative_drift: float
    observables_checked: tuple[str, ...]

    def require_converged(self) -> None:
        if not self.converged:
            raise ProtocolError(
                f"pre-equilibration did not converge within {self.horizon_hours} h "
                f"(max relative drift {self.max_relative_drift:.3e}, criterion "
                f"{self.criterion}). Extend the horizon or loosen the criterion "
                "deliberately; do not treat a fixed duration as convergence."
            )


@dataclass(frozen=True)
class StressProtocol:
    """Initial stress dose with constant or assumed first-order exposure.

    Dose is in uM and numerical model time is in hours.
    Clearance requires an explicitly supplied positive rate in h^-1.
    This exposure assumption is not an experimental calibration.
    """

    variant: str
    dose: float
    dose_units: str = "uM"
    clears: bool = False
    equilibration_horizon_hours: float = 50.0
    equilibration_tolerance: float = 1e-6
    solver: sim.SolverSettings = field(default_factory=sim.SolverSettings)
    clearance_rate_per_hour: float | None = None

    def __post_init__(self) -> None:
        if self.variant not in STRESS_PARAMETER:
            raise ProtocolError(
                f"unsupported variant: {self.variant!r}"
            )

        if self.dose_units != "uM":
            raise ProtocolError(
                "dose_units must be 'uM'; convert the dose explicitly "
                "before creating the protocol"
            )

        if not isinstance(self.clears, bool):
            raise ProtocolError("clears must be a boolean")

        if self.clears:
            object.__setattr__(
                self,
                "clearance_rate_per_hour",
                _positive_clearance_rate(
                    self.clearance_rate_per_hour
                ),
            )
        elif self.clearance_rate_per_hour is not None:
            raise ProtocolError(
                "clearance_rate_per_hour requires clears=True"
            )

        numeric_fields = (
            ("dose", False),
            ("equilibration_horizon_hours", True),
            ("equilibration_tolerance", True),
        )

        for name, strictly_positive in numeric_fields:
            raw_value = getattr(self, name)

            if isinstance(raw_value, (bool, np.bool_)):
                raise ProtocolError(f"{name} must be a number")

            try:
                value = float(raw_value)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ProtocolError(
                    f"{name} must be a finite number"
                ) from exc

            if not math.isfinite(value):
                raise ProtocolError(
                    f"{name} must be a finite number"
                )

            if strictly_positive and value <= 0:
                raise ProtocolError(f"{name} must be positive")

            if not strictly_positive and value < 0:
                raise ProtocolError(
                    f"{name} must be non-negative"
                )

            object.__setattr__(self, name, value)

    def stress_parameter(self) -> str:
        if self.variant not in STRESS_PARAMETER:
            raise ProtocolError(f"no stress parameter recorded for variant {self.variant!r}")
        return STRESS_PARAMETER[self.variant]

def _positive_clearance_rate(value) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ProtocolError(
            "clearance rate must be a positive finite number in h^-1"
        )

    try:
        rate = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ProtocolError(
            "clearance rate must be explicitly supplied in h^-1"
        ) from exc

    if not math.isfinite(rate) or rate <= 0:
        raise ProtocolError(
            "clearance rate must be a positive finite number in h^-1"
        )

    return rate

def has_rule_for(model: libsbml.Model, name: str) -> bool:
    """Whether a rule already governs the named quantity."""
    target_id = naming.try_resolve_id(model, name)
    return target_id is not None and target_id in naming.rule_targets(model)


def add_clearance(
    model: libsbml.Model,
    variant: str,
    rate_per_hour: float,
) -> None:
    """Add assumed first-order stress decay; numerical time is hours."""
    rate = _positive_clearance_rate(rate_per_hour)
    name = STRESS_PARAMETER.get(variant)

    if name is None:
        raise ProtocolError(f"unsupported variant: {variant!r}")

    if has_rule_for(model, name):
        raise DuplicateRuleError(
            f"{name} already has a rule"
        )

    param_id = naming.resolve_id(
        model, name, kinds=("parameter",)
    )

    if model.getInitialAssignmentBySymbol(param_id) is not None:
        raise ProtocolError(
            f"{name} has an initial assignment; resolve it first"
        )

    expression = libsbml.parseL3Formula(
        f"-{rate!r} * {param_id}"
    )
    if expression is None:
        raise ProtocolError("could not parse the clearance rate law")

    model.getParameter(param_id).setConstant(False)
    rule = model.createRateRule()
    rule.setVariable(param_id)
    rule.setMath(expression)

def load_protocol_model(
    sbml_string: str,
    protocol: StressProtocol,
):
    """Load a fresh model with the requested exposure assumption."""
    document = libsbml.readSBMLFromString(sbml_string)
    model = document.getModel()

    if model is None:
        raise ProtocolError("SBML contains no model")

    stress_id = naming.resolve_id(
        model,
        protocol.stress_parameter(),
        kinds=("parameter",),
    )

    if has_rule_for(model, protocol.stress_parameter()):
        raise DuplicateRuleError(
            "input SBML already governs stress; "
            "supply the baseline model"
        )

    if not model.getParameter(stress_id).getConstant():
        raise ProtocolError(
            "baseline stress parameter must be constant"
        )

    if protocol.clears:
        add_clearance(
            model,
            protocol.variant,
            protocol.clearance_rate_per_hour,
        )

    prepared_sbml = libsbml.writeSBMLToString(document)
    runner = sim.load_model(
        prepared_sbml,
        solver=protocol.solver,
        model=model,
    )
    return runner, model.clone()


def _check_runner_exposure(
    runner,
    protocol: StressProtocol,
) -> None:
    document = libsbml.readSBMLFromString(
        runner.getCurrentSBML()
    )
    model = document.getModel()

    if model is None:
        raise ProtocolError("runner contains no SBML model")

    stress_id = naming.resolve_id(
        model,
        protocol.stress_parameter(),
        kinds=("parameter",),
    )
    parameter = model.getParameter(stress_id)
    rule = model.getRuleByVariable(stress_id)

    if not protocol.clears:
        if rule is not None or not parameter.getConstant():
            raise ProtocolError(
                "constant exposure requested but loaded "
                "stress is not constant"
            )
        return

    if (
        rule is None
        or not rule.isRate()
        or parameter.getConstant()
    ):
        raise ProtocolError(
            "clearance requested but no compatible rate rule is loaded"
        )

    expected = libsbml.parseL3Formula(
        f"-{protocol.clearance_rate_per_hour!r} * {stress_id}"
    )
    if (
        libsbml.formulaToL3String(rule.getMath())
        != libsbml.formulaToL3String(expected)
    ):
        raise ProtocolError(
            "loaded clearance law differs from the protocol; "
            "use load_protocol_model"
        )


def _record_exposure(result, protocol: StressProtocol) -> None:
    result.exposure = {
        "variant": protocol.variant,
        "initial_dose": protocol.dose,
        "dose_units": protocol.dose_units,
        "time_units": "hour",
        "model": (
            "first_order_decay" if protocol.clears else "constant"
        ),
        "clearance_rate_per_hour": (
            protocol.clearance_rate_per_hour
        ),
    }

def _record_initialization(
    result,
    preparation: Equilibration,
    *,
    method: str,
) -> None:
    result.initialization = {
        "method": method,
        "starting_state": "source_sbml",
        "zero_stress_duration_hours": preparation.horizon_hours,
        "observables_converged": preparation.converged,
        "criterion": preparation.criterion,
        "max_relative_drift": preparation.max_relative_drift,
        "observables_checked": list(preparation.observables_checked),
        "experimental_initial_state_validated": False,
    }

def equilibrate(
    runner,
    protocol: StressProtocol,
    *,
    id_to_name: Optional[dict[str, str]] = None,
    observables: tuple[str, ...] = ("Observed_Green", "Measured_Ratio_RG"),
    n_points: int = 2001,
) -> Equilibration:
    """Run at zero stress until the observables stop moving.

    Leaves the runner at the equilibrated state, so a following simulate(...,
    reset=False) continues from the baseline rather than from zero.

    Convergence is judged on the relative change over the last tenth of the
    run, not on having run for a particular length of time.
    """
    _check_runner_exposure(runner, protocol)
    stress_name = protocol.stress_parameter()

    runner.reset()
    _set_runner_value(runner, id_to_name, stress_name, 0.0)

    result = sim.simulate(
        runner,
        end=protocol.equilibration_horizon_hours,
        n_points=n_points,
        id_to_name=id_to_name,
        solver=protocol.solver,
        reset=False,
    )

    checked = tuple(name for name in observables if result.has(name))
    drift = 0.0
    tail_start = max(1, int(0.9 * result.n_points))
    for name in checked:
        series = result[name]
        final = series[-1]
        tail = series[tail_start:]
        scale = max(abs(final), 1e-12)
        drift = max(drift, float(np.max(np.abs(tail - final)) / scale))

    return Equilibration(
        horizon_hours=protocol.equilibration_horizon_hours,
        converged=bool(checked) and drift <= protocol.equilibration_tolerance,
        criterion=(
            f"max |x(t) - x(T)| / |x(T)| over the final 10% of the run "
            f"<= {protocol.equilibration_tolerance:g}"
        ),
        max_relative_drift=drift,
        observables_checked=checked,
    )


def _set_runner_value(runner, id_to_name: Optional[dict[str, str]], name: str, value: float) -> None:
    """Set a quantity on a loaded runner, by readable name."""
    target_id = None
    if id_to_name:
        for element_id, element_name in id_to_name.items():
            if element_name == name:
                target_id = element_id
                break
    target_id = target_id or name
    try:
        runner[target_id] = value
    except Exception as exc:  # roadrunner raises its own error types
        raise ProtocolError(f"could not set {name!r} (id {target_id!r}) on the model") from exc


def apply_stress(runner, protocol: StressProtocol, id_to_name: Optional[dict[str, str]] = None) -> None:
    """Switch the stress input on at the protocol's dose. Call after
    equilibrate(); the baseline state is preserved."""
    _set_runner_value(runner, id_to_name, protocol.stress_parameter(), protocol.dose)


def run_protocol(
    runner,
    protocol: StressProtocol,
    *,
    duration_hours: float,
    n_points: int = 2001,
    id_to_name: Optional[dict[str, str]] = None,
    require_equilibrium: bool = True,
) -> tuple[sim.SimulationResult, Equilibration]:
    """Pre-equilibrate at zero stress, apply the bolus, then run.

    The returned result's time axis starts at 0 **at stress onset**, which is
    the quantity a decoder estimates. The equilibration record travels with
    it so that a run which never reached baseline cannot be mistaken for one
    that did.
    """
    equilibration = equilibrate(runner, protocol, id_to_name=id_to_name)
    if require_equilibrium:
        equilibration.require_converged()

    apply_stress(runner, protocol, id_to_name=id_to_name)
    result = sim.simulate(
        runner,
        end=duration_hours,
        n_points=n_points,
        id_to_name=id_to_name,
        solver=protocol.solver,
        reset=False,
    )
    _record_exposure(result, protocol)
    _record_initialization(
        result,
        equilibration,
        method=(
            "equilibrium_assumption"
            if require_equilibrium
            else "unchecked_baseline"
        ),
    )
    return result, equilibration

def run_protocol_at_times(
    runner,
    protocol: StressProtocol,
    *,
    times_hours,
    id_to_name: Optional[dict[str, str]] = None,
    require_equilibrium: bool = True,
) -> tuple[sim.SimulationResult, Equilibration]:
    """Pre-equilibrate, apply stress and observe at specified onset times."""
    # Validate before changing the runner.
    times = sim.validate_observation_times(times_hours)

    equilibration = equilibrate(
        runner, protocol, id_to_name=id_to_name
    )
    if require_equilibrium:
        equilibration.require_converged()

    apply_stress(runner, protocol, id_to_name=id_to_name)

    result = sim.simulate_at_times(
        runner,
        times_hours=times,
        id_to_name=id_to_name,
        solver=protocol.solver,
        reset=False,
    )
    _record_exposure(result, protocol)
    _record_initialization(
        result,
        equilibration,
        method=(
            "equilibrium_assumption"
            if require_equilibrium
            else "unchecked_baseline"
        ),
    )
    return result, equilibration

def run_protocol_after_preincubation(
    runner,
    protocol: StressProtocol,
    *,
    preincubation_hours,
    times_hours,
    id_to_name: Optional[dict[str, str]] = None,
) -> tuple[sim.SimulationResult, Equilibration]:
    """Run finite zero-stress preparation from source SBML initial conditions.

    This is an assumed scenario, not a validated experimental initial state.
    Convergence is reported but is not required.
    """
    if isinstance(preincubation_hours, (bool, np.bool_)):
        raise ProtocolError(
            "preincubation_hours must be a positive finite number"
        )

    try:
        duration = float(preincubation_hours)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ProtocolError(
            "preincubation_hours must be explicitly supplied in hours"
        ) from exc

    if not math.isfinite(duration) or duration <= 0:
        raise ProtocolError(
            "preincubation_hours must be a positive finite number"
        )

    finite_protocol = replace(
        protocol,
        equilibration_horizon_hours=duration,
    )

    result, preparation = run_protocol_at_times(
        runner,
        finite_protocol,
        times_hours=times_hours,
        id_to_name=id_to_name,
        require_equilibrium=False,
    )

    _record_initialization(
        result,
        preparation,
        method="finite_preincubation_assumption",
    )
    return result, preparation

def _record_initialization(
    result,
    preparation: Equilibration,
    *,
    method: str,
) -> None:
    result.initialization = {
        "method": method,
        "starting_state": "source_sbml",
        "zero_stress_duration_hours": preparation.horizon_hours,
        "observables_converged": preparation.converged,
        "criterion": preparation.criterion,
        "max_relative_drift": preparation.max_relative_drift,
        "observables_checked": list(preparation.observables_checked),
        "experimental_initial_state_validated": False,
    }
