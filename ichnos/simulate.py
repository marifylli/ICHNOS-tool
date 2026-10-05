"""One simulator, used by everything that runs the model.

Calibration, decoder validation and any Jacobian all have to run the model the
same way, or the forward map a decoder inverts is not the model the decoder
was validated against. Having a single entry point is what makes that true by
construction rather than by discipline.

Loading still goes through tellurium, exactly as in Ichnos_PULSE. Keeping the
same loader means the migration changes no behaviour; moving to libroadrunner
directly is a separate change with its own verification.

Extracted from Ichnos_PULSE python/run_ichnos.py and the sensitivity script
@ e66de65c, with the plotting, printing, implicit export and study-specific
stress grids left behind.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np


DEFAULT_ABSOLUTE_TOLERANCE = 1e-10
DEFAULT_RELATIVE_TOLERANCE = 1e-8


@dataclass(frozen=True)
class SolverSettings:
    """Recorded alongside every result, and into any calibration artifact: a
    forward map is only reproducible if the tolerances that produced it are
    known."""

    absolute_tolerance: float = DEFAULT_ABSOLUTE_TOLERANCE
    relative_tolerance: float = DEFAULT_RELATIVE_TOLERANCE
    stiff: bool = True

def apply_solver_settings(runner, settings: SolverSettings) -> None:
    """Apply the requested settings to the active integrator."""
    integrator = runner.getIntegrator()
    integrator.setValue("absolute_tolerance", settings.absolute_tolerance)
    integrator.setValue("relative_tolerance", settings.relative_tolerance)
    integrator.setValue("stiff", settings.stiff)


def read_solver_settings(runner) -> SolverSettings:
    """Read the settings actually used by the active integrator."""
    integrator = runner.getIntegrator()
    return SolverSettings(
        absolute_tolerance=float(integrator.getValue("absolute_tolerance")),
        relative_tolerance=float(integrator.getValue("relative_tolerance")),
        stiff=bool(integrator.getValue("stiff")),
    )

@dataclass
class SimulationResult:
    time: np.ndarray
    columns: list[str]
    values: np.ndarray  # (n_timepoints, n_columns), excluding time
    solver: SolverSettings
    exposure: dict | None = None
    initialization: dict | None = None
    n_points: int = field(init=False)

    def __post_init__(self):
        self.n_points = len(self.time)

    def __getitem__(self, name: str) -> np.ndarray:
        if name == "time":
            return self.time
        return self.values[:, self.columns.index(name)]

    def final(self, name: str) -> float:
        return float(self[name][-1])

    def has(self, name: str) -> bool:
        return name in self.columns


def load_model(
    sbml_string: str,
    solver: Optional[SolverSettings] = None,
    model: Optional["libsbml.Model"] = None,
):
    """Load an SBML string into a RoadRunner instance via tellurium.

    If `model` is given, the runner's selections are widened to include
    assignment-rule parameters. This is not cosmetic. RoadRunner's default
    output is floating species plus rate-rule states only, so
    `Observed_Green`, `Measured_Ratio_RG`, `Ratio_RG_FRET`, `Total_red_pool`
    and `b_fret` -- all assignment-rule parameters, and between them the
    entire measured readout of this circuit -- are simply absent from a
    default simulate(). Code that looks for them finds nothing and, depending
    on how it checks, either fails or quietly reports on an empty set.
    """
    import tellurium as te

    solver = solver or SolverSettings()
    runner = te.loadSBMLModel(sbml_string)

    apply_solver_settings(runner, solver)

    if model is not None:
        select_observables(runner, model)
    return runner


def observable_ids(model) -> list[str]:
    """Ids of the assignment-rule parameters -- the model's readouts."""
    out = []
    for rule in model.getListOfRules():
        if rule.getElementName() != "assignmentRule" or not rule.isSetVariable():
            continue
        target = model.getParameter(rule.getVariable())
        if target is not None:
            out.append(rule.getVariable())
    return out


def select_observables(runner, model) -> list[str]:
    """Add the assignment-rule readouts to the runner's default selections.

    Returns the full selection list actually set.
    """
    existing = list(runner.selections)
    for element_id in observable_ids(model):
        if element_id not in existing and f"[{element_id}]" not in existing:
            existing.append(element_id)
    runner.selections = existing
    return existing


def _readable_columns(raw_columns, id_to_name: dict[str, str]) -> list[str]:
    out = []
    for col in raw_columns:
        if col == "time":
            out.append(col)
            continue
        inner = col.strip("[]")
        out.append(id_to_name.get(inner, inner))
    return out


def simulate(
    runner,
    *,
    start: float = 0.0,
    end: float,
    n_points: int,
    id_to_name: Optional[dict[str, str]] = None,
    solver: Optional[SolverSettings] = None,
    reset: bool = True,
) -> SimulationResult:
    """Run the model over [start, end].

    `reset=False` continues from the runner's current state, which is how a
    protocol switches stress on without discarding the pre-equilibrated
    baseline.

    Column names come back readable (`Observed_Green`) rather than as GUIDs,
    but the underlying lookups elsewhere always use real ids.
    """
    if reset:
        runner.reset()

    if solver is not None:
        apply_solver_settings(runner, solver)

    raw = runner.simulate(start, end, n_points)
    columns = _readable_columns(raw.colnames, id_to_name or {})
    array = np.asarray(raw)

    return SimulationResult(
        time=array[:, 0],
        columns=columns[1:],
        values=array[:, 1:],
        solver=read_solver_settings(runner),
    )

def validate_observation_times(times_hours) -> np.ndarray:
    """Finite, increasing observation times in hours after onset."""
    try:
        raw = np.asarray(times_hours, dtype=object)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "times_hours must be a one-dimensional sequence"
        ) from exc

    if raw.ndim != 1 or raw.size == 0:
        raise ValueError(
            "times_hours must be a non-empty one-dimensional sequence"
        )

    if any(isinstance(value, (bool, np.bool_)) for value in raw):
        raise ValueError("times_hours must not contain booleans")

    try:
        times = raw.astype(float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(
            "times_hours must contain finite numbers"
        ) from exc

    if not np.isfinite(times).all():
        raise ValueError("times_hours must contain finite numbers")

    if (times < 0).any():
        raise ValueError("times_hours must be non-negative")

    if (np.diff(times) <= 0).any():
        raise ValueError("times_hours must be strictly increasing")

    if times[-1] <= 0:
        raise ValueError(
            "times_hours must include at least one positive time"
        )

    return times


def simulate_at_times(
    runner,
    *,
    times_hours,
    id_to_name: Optional[dict[str, str]] = None,
    solver: Optional[SolverSettings] = None,
    reset: bool = True,
) -> SimulationResult:
    """Simulate from onset and return only the requested observation times.

    With reset=False, preserve the current state and treat it as the state
    at onset. This allows continuation from a pre-equilibrated baseline.
    """
    times = validate_observation_times(times_hours)
    includes_onset = times[0] == 0.0
    output_times = (
        times if includes_onset else np.concatenate(([0.0], times))
    )

    if reset:
        runner.reset()

    if solver is not None:
        apply_solver_settings(runner, solver)

    raw = runner.simulate(times=output_times.tolist())
    columns = _readable_columns(raw.colnames, id_to_name or {})
    array = np.asarray(raw)

    if not includes_onset:
        array = array[1:]

    return SimulationResult(
        time=array[:, 0].copy(),
        columns=columns[1:],
        values=array[:, 1:].copy(),
        solver=read_solver_settings(runner),
    )

def peak_summary(result: SimulationResult, names: tuple[str, ...]) -> list[dict]:
    """Peak, time-to-peak and plateau for the named observables.

    Carries over the undersampling guard from run_ichnos.py, and makes it a
    returned flag rather than a printed warning. It matters: an adaptive
    sensor can peak within minutes, and A_ox once came out as 0.5342 at 24
    minutes instead of 0.8622 at 4.2 minutes on a 0.4 h grid -- wrong by 38%
    and entirely plausible-looking.
    """
    rows = []
    dt = float(result.time[1] - result.time[0]) if result.n_points > 1 else float("nan")
    for name in names:
        if not result.has(name):
            continue
        series = result[name]
        index = int(np.argmax(series))
        peak = float(series[index])
        t_peak = float(result.time[index])
        plateau = float(series[-1])
        undersampled = bool(dt > 0 and t_peak <= 5 * dt and peak > plateau)
        rows.append(
            {
                "name": name,
                "peak": peak,
                "t_peak": t_peak,
                "plateau": plateau,
                "peak_over_plateau": peak / plateau if plateau > 0 else float("nan"),
                "undersampled": undersampled,
                "grid_hours": dt,
            }
        )
    return rows
