import libsbml
import numpy as np
import pytest

from ichnos import protocol as p


def _baseline_sbml(variant):
    document = libsbml.SBMLDocument(3, 2)
    model = document.createModel()
    model.setId("finite_preincubation")

    parameters = (
        (p.STRESS_PARAMETER[variant], 0.0, True),
        ("X", 0.0, False),
        ("Observed_Green", 0.0, False),
        ("Measured_Ratio_RG", 1.0, False),
    )

    for name, value, constant in parameters:
        parameter = model.createParameter()
        parameter.setId(name)
        parameter.setName(name)
        parameter.setValue(value)
        parameter.setConstant(constant)

    rule = model.createRateRule()
    rule.setVariable("X")
    rule.setMath(libsbml.parseL3Formula("1 - 0.5 * X"))

    for name, formula in (
        ("Observed_Green", "X"),
        ("Measured_Ratio_RG", "1"),
    ):
        rule = model.createAssignmentRule()
        rule.setVariable(name)
        rule.setMath(libsbml.parseL3Formula(formula))

    return libsbml.writeSBMLToString(document)


@pytest.mark.parametrize("variant", ["ox", "er"])
@pytest.mark.parametrize("duration", [2.0, 4.0])
def test_finite_preincubation_preserves_state(variant, duration):
    protocol = p.StressProtocol(variant=variant, dose=75.0)
    runner, _ = p.load_protocol_model(
        _baseline_sbml(variant),
        protocol,
    )

    runner.timeCourseSelections = [
        *runner.timeCourseSelections,
        p.STRESS_PARAMETER[variant],
    ]

    result, preparation = p.run_protocol_after_preincubation(
        runner,
        protocol,
        preincubation_hours=duration,
        times_hours=[0.0, 0.5, 1.0],
    )

    expected = 2.0 * (
        1.0 - np.exp(-0.5 * (duration + result.time))
    )
    np.testing.assert_allclose(
        result["Observed_Green"],
        expected,
        rtol=1e-6,
        atol=1e-8,
    )
    np.testing.assert_allclose(result.time, [0.0, 0.5, 1.0])
    np.testing.assert_allclose(
        result[p.STRESS_PARAMETER[variant]],
        75.0,
    )

    assert preparation.horizon_hours == duration
    assert preparation.converged is False

    metadata = result.initialization
    assert metadata["method"] == "finite_preincubation_assumption"
    assert metadata["starting_state"] == "source_sbml"
    assert metadata["zero_stress_duration_hours"] == duration
    assert metadata["observables_converged"] is False
    assert metadata["experimental_initial_state_validated"] is False
    assert metadata["observables_checked"] == [
        "Observed_Green",
        "Measured_Ratio_RG",
    ]

    assert result.exposure["initial_dose"] == 75.0
    assert protocol.equilibration_horizon_hours == 50.0


@pytest.mark.parametrize(
    "duration",
    [None, 0, -1, float("nan"), float("inf"), True, "invalid"],
)
def test_invalid_preincubation_rejected_before_runner_access(duration):
    protocol = p.StressProtocol(variant="ox", dose=75.0)

    with pytest.raises(p.ProtocolError, match="preincubation_hours"):
        p.run_protocol_after_preincubation(
            object(),
            protocol,
            preincubation_hours=duration,
            times_hours=[0.5, 1.0],
        )