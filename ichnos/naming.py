"""Resolving a human-readable name to an SBML id.

The submodels were exported by different COPASI/SimBiology sessions, so their
ids are GUIDs that do not line up across files. Everything in this project
therefore matches by *name*. That makes name resolution load-bearing, and a
silent failure here is a wrong model rather than an error: looking up a name
that does not exist and getting `None` back means a parameter is never set and
the run proceeds with the old value.

So this module raises. Missing is an error; ambiguous is an error. The caller
never gets `None` and never gets an arbitrary pick between two matches.

Extracted from Ichnos_PULSE python/ichnos_io.py and ichnos_core.py @ e66de65c.
The lookups there returned None on a miss and took the first match on a
collision; the resolution logic is otherwise unchanged.
"""
from __future__ import annotations

import libsbml


class NameNotFoundError(KeyError):
    """No element of the requested kind carries this name."""


class AmbiguousNameError(KeyError):
    """More than one element carries this name, so the reference is unsafe."""


def _display_name(element) -> str:
    return element.getName() or element.getId()


def _collect(model: libsbml.Model, kinds: tuple[str, ...]) -> list:
    elements = []
    if "parameter" in kinds:
        elements.extend(model.getListOfParameters())
    if "species" in kinds:
        elements.extend(model.getListOfSpecies())
    if "compartment" in kinds:
        elements.extend(model.getListOfCompartments())
    return elements


def resolve_id(
    model: libsbml.Model,
    name: str,
    kinds: tuple[str, ...] = ("parameter", "species"),
) -> str:
    """The SBML id of the single element named `name`.

    Searches parameters then species by default. Raises NameNotFoundError if
    nothing matches and AmbiguousNameError if more than one does -- including
    a name that exists as both a parameter and a species, which is exactly the
    case where picking one silently would be worst.
    """
    matches = [e for e in _collect(model, kinds) if _display_name(e) == name]
    if not matches:
        available = sorted(_display_name(e) for e in _collect(model, kinds))
        raise NameNotFoundError(
            f"no {'/'.join(kinds)} named {name!r} in this model. "
            f"Available: {available}"
        )
    if len(matches) > 1:
        raise AmbiguousNameError(
            f"{name!r} matches {len(matches)} elements "
            f"({[e.getId() for e in matches]}); the reference is unsafe"
        )
    return matches[0].getId()


def try_resolve_id(model: libsbml.Model, name: str, **kwargs) -> str | None:
    """resolve_id() returning None when the name is absent.

    An ambiguous name still raises -- absence can be a legitimate question
    ("does this variant have A_ox?"), but a collision is always a defect.
    """
    try:
        return resolve_id(model, name, **kwargs)
    except NameNotFoundError:
        return None


def id_to_name_map(model: libsbml.Model) -> dict[str, str]:
    """Every species and parameter id mapped to its readable name.

    For display only -- plot legends, printed values. Merge and simulation
    logic always works on real ids.
    """
    mapping = {}
    for element in _collect(model, ("species", "parameter")):
        mapping[element.getId()] = _display_name(element)
    return mapping


def rule_targets(model: libsbml.Model) -> dict[str, str]:
    """Element id -> rule kind, for every id a rule assigns to.

    Used to refuse to set a quantity that a rule governs: writing to it either
    does nothing or changes the circuit's structure.
    """
    return {
        rule.getVariable(): rule.getElementName()
        for rule in model.getListOfRules()
        if rule.isSetVariable()
    }


def initial_assignment_targets(model: libsbml.Model) -> set[str]:
    """Ids whose initial value comes from an initial assignment."""
    return {a.getSymbol() for a in model.getListOfInitialAssignments()}
