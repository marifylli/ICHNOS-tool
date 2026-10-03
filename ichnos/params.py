"""Parameter profiles: complete, executable configurations of one variant.

A profile is all-or-nothing. It names every free parameter the variant has,
and it is applied whole. Values are never taken from two profiles at once,
because a parameter set that runs is not the same thing as a parameter set
that was jointly estimated, and mixing two of them produces a model whose
behaviour nobody has checked.

Each parameter carries its own `status` and `source`, which is a statement
about where that number came from, not a grouping key. The `default` profile
in particular mixes fitted, swept and inherited values and records
`joint_fit: false` to say so explicitly.

Rule-governed quantities cannot be set by a profile. `P` is the case that
matters: the builder wires the reporter onto the dynamic, assignment-rule `P`
from TIP-TetR and sets aside the reporter's own static `P = 0.5`. A profile
that assigned `P` would either be silently overwritten by the rule or, worse,
turn a dynamic quantity into a constant.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import libsbml
import math

from . import naming
from .config import model_path

PARAMETERS_FILE = "parameters.yaml"

#: Values that `status` may take. See parameters.yaml for what each means.
VALID_STATUSES = ("fitted", "swept", "non_identifiable", "inherited", "measured")


class ProfileError(ValueError):
    """The requested profile does not exist, or is not applicable as given."""


class IncompleteProfileError(ProfileError):
    """The profile does not cover every free parameter of the variant."""


class GovernedParameterError(ProfileError):
    """The profile tries to set a quantity that a rule governs."""

class InvalidParameterValueError(ProfileError):
    """A parameter value is outside the supported domain."""


class UnitMismatchError(ProfileError):
    """Profile units do not match the SBML unit identifier."""


STRICTLY_POSITIVE_PARAMETERS = {
    "K_R",
    "Kd_TIP_TetR",
    "K_act_ox",
    "K_act_er",
    "n",
    "n_ox",
    "n_er",
    "f",
    "eps",
}

@dataclass(frozen=True)
class ParameterEntry:
    name: str
    value: float
    units: str
    status: str
    source: Optional[str]
    note: Optional[str] = None

    def __post_init__(self):
        if isinstance(self.value, bool):
            raise InvalidParameterValueError(
                f"{self.name}: booleans are not parameter values"
            )

        try:
            value = float(self.value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise InvalidParameterValueError(
                f"{self.name}: value must be numeric"
            ) from exc

        if not math.isfinite(value) or value < 0:
            raise InvalidParameterValueError(
                f"{self.name}: value must be finite and non-negative"
            )

        if self.name in STRICTLY_POSITIVE_PARAMETERS and value <= 0:
            raise InvalidParameterValueError(
                f"{self.name}: value must be strictly positive"
            )

        if self.name in {"E", "P_min"} and value > 1:
            raise InvalidParameterValueError(
                f"{self.name}: value must be between 0 and 1"
            )

        if not isinstance(self.units, str) or not self.units.strip():
            raise UnitMismatchError(
                f"{self.name}: units must be explicitly specified"
            )

        object.__setattr__(self, "value", value)

    @property
    def has_recorded_origin(self) -> bool:
        return self.source is not None


@dataclass(frozen=True)
class Profile:
    variant: str
    name: str
    description: str
    joint_fit: bool
    parameters: dict[str, ParameterEntry]
    governed: dict[str, str]
    initial_state: dict[str, float]

    def values(self) -> dict[str, float]:
        return {name: entry.value for name, entry in self.parameters.items()}

    def by_status(self, status: str) -> list[str]:
        return sorted(n for n, e in self.parameters.items() if e.status == status)

    def provenance_summary(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for entry in self.parameters.values():
            counts[entry.status] = counts.get(entry.status, 0) + 1
        return counts


def _load_yaml() -> dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover - environment-dependent
        raise ImportError(
            "PyYAML is required to read parameter profiles; "
            "install this package with the [model] extra"
        ) from exc
    return yaml.safe_load(model_path(PARAMETERS_FILE).read_text())


def available_profiles(variant: str) -> list[str]:
    data = _load_yaml()
    if variant not in data:
        raise ProfileError(f"no profiles recorded for variant {variant!r}")
    return sorted(k for k in data[variant] if not k.startswith("_"))


def load_profile(variant: str, name: str = "default") -> Profile:
    """Read one profile. Does not touch a model."""
    data = _load_yaml()
    if variant not in data:
        raise ProfileError(f"no profiles recorded for variant {variant!r}")
    if name not in data[variant]:
        raise ProfileError(
            f"variant {variant!r} has no profile {name!r}; "
            f"available: {available_profiles(variant)}"
        )
    raw = data[variant][name]

    entries = {}
    for param_name, spec in raw["parameters"].items():
        status = spec.get("status", "inherited")
        if status not in VALID_STATUSES:
            raise ProfileError(
                f"{variant}/{name}/{param_name}: status {status!r} is not one of "
                f"{VALID_STATUSES}"
            )
        entries[param_name] = ParameterEntry(
            name=param_name,
            value=spec["value"],
            units=spec.get("units"),
            status=status,
            source=spec.get("source"),
            note=spec.get("note"),
        )

    return Profile(
        variant=variant,
        name=name,
        description=raw.get("description", ""),
        joint_fit=bool(raw.get("joint_fit", False)),
        parameters=entries,
        governed=dict(raw.get("governed", {})),
        initial_state=dict(raw.get("initial_state", {})),
    )


def free_parameter_names(model: libsbml.Model) -> set[str]:
    """Names of the parameters a profile is expected to cover: every parameter
    that is not set by a rule or an initial assignment.
    """
    governed = set(naming.rule_targets(model)) | naming.initial_assignment_targets(model)
    return {
        (p.getName() or p.getId())
        for p in model.getListOfParameters()
        if p.getId() not in governed
    }

def _checked_parameter(
    model: libsbml.Model,
    name: str,
    entry: ParameterEntry,
):
    if name != entry.name:
        raise ProfileError(
            f"profile key {name!r} differs from "
            f"entry name {entry.name!r}"
        )

    param_id = naming.resolve_id(
        model, name, kinds=("parameter",)
    )

    governed = (
        set(naming.rule_targets(model))
        | naming.initial_assignment_targets(model)
    )

    if param_id in governed:
        raise GovernedParameterError(
            f"{name} is governed by a rule or initial assignment"
        )

    parameter = model.getParameter(param_id)

    if not parameter.isSetUnits():
        raise UnitMismatchError(
            f"{name}: SBML units are unspecified"
        )

    expected = parameter.getUnits()

    if entry.units != expected:
        raise UnitMismatchError(
            f"{name}: profile units {entry.units!r} differ from "
            f"SBML units {expected!r}; "
            "no automatic conversion is performed"
        )

    return parameter

def check_profile_against_model(profile: Profile, model: libsbml.Model) -> None:
    """Raise unless the profile is a complete, applicable configuration.

    Three failure modes, all of which would otherwise surface as a run that
    looks fine and is wrong:

      - the profile omits a free parameter, so that parameter silently keeps
        whatever the SBML happened to carry;
      - the profile names a parameter the model does not have, usually a
        typo or a leftover from another variant;
      - the profile sets a rule-governed quantity.
    """
    free = free_parameter_names(model)
    covered = set(profile.parameters)

    missing = free - covered
    if missing:
        raise IncompleteProfileError(
            f"profile {profile.variant}/{profile.name} does not cover "
            f"{sorted(missing)}. A profile must be a complete configuration; "
            "a partial one leaves parameters at whatever the SBML carries."
        )

    unknown = covered - free
    governed_names = {
        (model.getElementBySId(i).getName() or i) if model.getElementBySId(i) else i
        for i in set(naming.rule_targets(model)) | naming.initial_assignment_targets(model)
    }
    governed_overlap = unknown & governed_names
    if governed_overlap:
        raise GovernedParameterError(
            f"profile {profile.variant}/{profile.name} tries to set "
            f"{sorted(governed_overlap)}, which a rule or initial assignment "
            "governs. Setting these either does nothing or changes the "
            "circuit; P in particular must stay the dynamic assignment-rule "
            "parameter the builder wires the reporter onto."
        )

    truly_unknown = unknown - governed_names
    if truly_unknown:
        raise ProfileError(
            f"profile {profile.variant}/{profile.name} names "
            f"{sorted(truly_unknown)}, which this model does not have"
        )
        for name, entry in profile.parameters.items():
            _checked_parameter(model, name, entry)


def apply_profile(
    model: libsbml.Model,
    profile: Profile,
    *,
    check: bool = True,
) -> None:
    """Validate all supplied entries before writing parameter values.

    check=False skips completeness checks only. Units and governed
    targets are still checked for every supplied entry.
    """
    if check:
        check_profile_against_model(profile, model)

    updates = [
        (_checked_parameter(model, name, entry), entry.value)
        for name, entry in profile.parameters.items()
    ]

    for parameter, value in updates:
        parameter.setValue(value)


def describe(profile: Profile) -> str:
    """A short, honest summary of what a profile is, for logs and manifests."""
    counts = profile.provenance_summary()
    breakdown = ", ".join(f"{n} {status}" for status, n in sorted(counts.items()))
    joint = (
        "jointly fitted"
        if profile.joint_fit
        else "NOT a joint fit -- values reached the model by different routes"
    )
    return (
        f"{profile.variant}/{profile.name}: {len(profile.parameters)} parameters "
        f"({breakdown}); {joint}"
    )
