"""A profile is a complete, executable configuration, and `default` must
reproduce today's behaviour exactly.
"""
import io
import contextlib

import libsbml
import pytest

from ichnos import build, params
from dataclasses import replace



VARIANTS = ("ox", "er")


def _model(variant):
    with contextlib.redirect_stdout(io.StringIO()):
        sbml = build.build_variant_sbml_string(variant)
    return sbml, libsbml.readSBMLFromString(sbml).getModel()


@pytest.mark.parametrize("variant", VARIANTS)
def test_default_profile_is_a_no_op(variant):
    """The strongest statement available about `default`: applying it changes
    nothing, so every result produced before profiles existed is still the
    result the profile produces.
    """
    sbml, model = _model(variant)
    profile = params.load_profile(variant, "default")
    params.apply_profile(model, profile)
    assert libsbml.writeSBMLToString(model.getSBMLDocument()) == sbml


@pytest.mark.parametrize("variant", VARIANTS)
def test_default_profile_covers_every_free_parameter(variant):
    _, model = _model(variant)
    profile = params.load_profile(variant, "default")
    params.check_profile_against_model(profile, model)
    assert set(profile.parameters) == params.free_parameter_names(model)


@pytest.mark.parametrize("variant", VARIANTS)
def test_default_profile_is_not_claimed_to_be_a_joint_fit(variant):
    """The set runs and is internally consistent. It was not jointly
    estimated, and must not be presented as if it had been.
    """
    profile = params.load_profile(variant, "default")
    assert profile.joint_fit is False
    assert "NOT a joint fit" in params.describe(profile)


def test_provenance_is_recorded_per_parameter():
    """Different origins do not make different profiles, but they are
    recorded. The oxidative set mixes values fitted to Delaunay 2000, values
    chosen by sweeping for a target shape, and a parameter that is not
    identifiable on its own.
    """
    profile = params.load_profile("ox", "default")

    assert profile.parameters["K_act_ox"].status == "fit_summary"
    assert "Delaunay" in profile.parameters["K_act_ox"].source
    assert profile.parameters["n_ox"].status == "fit_summary"

    assert profile.parameters["k_off_ox"].status == "swept"
    assert profile.parameters["d_x_ox"].status == "swept"

    assert profile.parameters["k_on_ox"].status == "model_convention"

    # Values without a recorded derivation say so rather than borrowing a
    # neighbouring parameter's citation.
    inherited = profile.by_status("inherited")
    assert inherited
    for name in inherited:
        assert profile.parameters[name].source is None


@pytest.mark.parametrize("variant", VARIANTS)
def test_profile_may_not_set_a_rule_governed_quantity(variant):
    """P is the case that matters: it is an assignment-rule parameter, and
    the builder wires the reporter onto it rather than onto the reporter's own
    static P = 0.5. A profile that set P would change the circuit.
    """
    _, model = _model(variant)
    profile = params.load_profile(variant, "default")
    assert "P" in profile.governed
    assert profile.governed["P"] == "assignmentRule"
    assert "P" not in profile.parameters

    tampered = dict(profile.parameters)
    tampered["P"] = params.ParameterEntry("P", 0.5, "dimensionless", "inherited", None)
    broken = params.Profile(
        variant=variant, name="broken", description="", joint_fit=False,
        parameters=tampered, governed=profile.governed,
        initial_state=profile.initial_state,
    )
    with pytest.raises(params.GovernedParameterError):
        params.check_profile_against_model(broken, model)


@pytest.mark.parametrize("variant", VARIANTS)
def test_incomplete_profile_is_rejected(variant):
    """Half a profile would leave the rest of the parameters at whatever the
    SBML carries, producing a mixture nobody chose.
    """
    _, model = _model(variant)
    profile = params.load_profile(variant, "default")
    trimmed = dict(profile.parameters)
    trimmed.pop(sorted(trimmed)[0])
    partial = params.Profile(
        variant=variant, name="partial", description="", joint_fit=False,
        parameters=trimmed, governed=profile.governed,
        initial_state=profile.initial_state,
    )
    with pytest.raises(params.IncompleteProfileError):
        params.check_profile_against_model(partial, model)


def test_unknown_profile_and_variant_raise():
    with pytest.raises(params.ProfileError):
        params.load_profile("ox", "does-not-exist")
    with pytest.raises(params.ProfileError):
        params.load_profile("cu", "default")


@pytest.mark.parametrize("variant", VARIANTS)
def test_initial_state_is_all_zero(variant):
    """Recorded because it is the reason pre-equilibration exists: the source
    SBML starts every state at zero, which is not a culture at baseline.
    """
    profile = params.load_profile(variant, "default")
    assert profile.initial_state
    assert all(value == 0 for value in profile.initial_state.values())

def _changed_profile(profile, name, **changes):
    entries = dict(profile.parameters)
    entries[name] = replace(entries[name], **changes)
    return replace(profile, parameters=entries)


@pytest.mark.parametrize(
    "name, value",
    [
        ("S_ox", -1),
        ("S_ox", float("nan")),
        ("S_ox", float("inf")),
        ("S_ox", True),
        ("S_ox", "invalid"),
        ("K_act_ox", 0),
        ("n_ox", 0),
        ("eps", 0),
        ("f", 0),
        ("E", 1.1),
        ("P_min", 1.1),
    ],
)
def test_invalid_parameter_values_are_rejected(name, value):
    profile = params.load_profile("ox")

    with pytest.raises(params.InvalidParameterValueError):
        _changed_profile(profile, name, value=value)


@pytest.mark.parametrize("variant", VARIANTS)
@pytest.mark.parametrize("check", [True, False])
def test_invalid_units_do_not_partially_modify_model(variant, check):
    _, model = _model(variant)
    profile = params.load_profile(variant)

    # An earlier valid change must not be written if a later entry fails.
    profile = _changed_profile(profile, "b", value=8)
    profile = _changed_profile(
        profile, "eps", units="dimensionless"
    )

    before = libsbml.writeSBMLToString(model.getSBMLDocument())

    with pytest.raises(params.UnitMismatchError):
        params.apply_profile(model, profile, check=check)

    after = libsbml.writeSBMLToString(model.getSBMLDocument())
    assert after == before


@pytest.mark.parametrize("units", [None, "", " "])
def test_unspecified_profile_units_are_rejected(units):
    profile = params.load_profile("ox")

    with pytest.raises(params.UnitMismatchError):
        _changed_profile(profile, "K_act_ox", units=units)


def test_zero_stress_is_allowed():
    _, model = _model("ox")
    profile = _changed_profile(
        params.load_profile("ox"),
        "S_ox",
        value=0,
    )

    params.apply_profile(model, profile)