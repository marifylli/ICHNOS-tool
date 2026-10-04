"""What is done to the cells, and when: pre-equilibration, stress onset, and
whether the stress clears.

Two things here exist because getting them wrong is invisible in the output.

First, every state in the source SBML starts at zero. A run from those initial
conditions is a model starting from nothing, not a culture at baseline, and
its early transient is an artifact of the initial condition rather than a
response to stress. So a protocol pre-equilibrates at zero stress first, and
records how long it ran and what convergence criterion it met -- a fixed
duration is not evidence of convergence, and what equilibrates in 50 h under
one parameter profile may not under another.

Second, the merged model as it stands has NO clearance rule: `S_ox` and
`S_er` are constant parameters, so stress never decays. Adding clearance is a
change to the model, not a setting, and `add_clearance()` below refuses to add
a second rule for a quantity that already has one.

After pre-equilibration, time is reported relative to stress onset. That is
the `t` a decoder estimates.
"""
from __future__ import annotations

from dataclasses import dataclass, field
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
    """A bolus applied at onset, held at a constant level thereafter.

    `clears` is False for the model as it stands. A decoder calibrated under
    this protocol estimates the initial bolus dose and the time since onset;
    it does not distinguish an unknown bolus from continuous exposure or from
    repeated pulses without further information.
    """

    variant: str
    dose: float
    dose_units: str = "uM"
    clears: bool = False
    equilibration_horizon_hours: float = 50.0
    equilibration_tolerance: float = 1e-6
    solver: sim.SolverSettings = field(default_factory=sim.SolverSettings)

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


def has_rule_for(model: libsbml.Model, name: str) -> bool:
    """Whether a rule already governs the named quantity."""
    target_id = naming.try_resolve_id(model, name)
    return target_id is not None and target_id in naming.rule_targets(model)


def add_clearance(model: libsbml.Model, variant: str, rate_per_hour: float) -> None:
    """Add dS/dt = -k * S for the variant's stress input.

    Guarded on purpose. The original helper in the research code assumed the
    stress parameter had no rule, which held for the models it ran against and
    stops holding the moment one gains a rule. A second rule for the same
    variable is invalid SBML, and if it were accepted it would change the
    input dynamics with no visible sign.

    This changes the model. It belongs with a parameter profile that was
    determined under clearance -- bolting it onto a profile fitted without it
    produces a configuration nobody has checked.
    """
    name = STRESS_PARAMETER.get(variant)
    if name is None:
        raise ProtocolError(f"no stress parameter recorded for variant {variant!r}")
    if has_rule_for(model, name):
        raise DuplicateRuleError(
            f"{name} already has a rule; refusing to add a second one. "
            "If clearance is meant to replace the existing rule, remove that "
            "rule explicitly."
        )

    param_id = naming.resolve_id(model, name, kinds=("parameter",))
    model.getParameter(param_id).setConstant(False)

    rule = model.createRateRule()
    rule.setVariable(param_id)
    math = libsbml.parseL3Formula(f"-{rate_per_hour} * {param_id}")
    if math is None:
        raise ProtocolError("could not parse the clearance rate law")
    rule.setMath(math)


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
    return result, equilibration
