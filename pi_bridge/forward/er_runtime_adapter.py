"""Prepare ER scenario SBML in memory, without changing source models."""

from __future__ import annotations

import libsbml

from ichnos import build, naming, params


EXPECTED_IDS = {"frozen", "er_m2_n4"}
EXPECTED_OVERRIDES = {
    "K_act_er",
    "n_er",
    "k_on_er",
    "k_off_er",
    "d_x_er",
}


def prepare_er_scenario_sbml(scenario: dict) -> tuple[str, dict]:
    """Return runtime SBML and verified applied parameter values."""

    scenario_id = scenario["scenario_id"]
    overrides = scenario["parameter_overrides"]

    if scenario_id not in EXPECTED_IDS:
        raise ValueError(f"Unknown ER scenario: {scenario_id}")

    if scenario_id == "frozen" and overrides:
        raise ValueError("Frozen ER scenario must have no overrides")

    if scenario_id == "er_m2_n4":
        if set(overrides) != EXPECTED_OVERRIDES:
            raise ValueError("Incorrect ER PI override parameter set")

    sbml = build.build_variant_sbml_string(
        "er", save_sbml=False
    )

    document = libsbml.readSBMLFromString(sbml)
    model = document.getModel()

    if model is None:
        raise RuntimeError("ER SBML contains no model")

    # Apply the Tool's canonical parameter profile.
    baseline_profile = params.load_profile("er", "default")
    params.apply_profile(model, baseline_profile)

    if overrides:
        entries = {}

        for name, value in overrides.items():
            parameter_id = naming.resolve_id(
                model, name, kinds=("parameter",)
            )
            parameter = model.getParameter(parameter_id)

            entries[name] = params.ParameterEntry(
                name=name,
                value=float(value),
                units=parameter.getUnits(),
                status="physics_informed",
                source="Ichnos_PULSE/python/pi_scenarios.py",
                note=f"ER scenario {scenario_id}",
            )

        pi_profile = params.Profile(
            variant="er",
            name=scenario_id,
            description="Frozen ER PI runtime overrides",
            joint_fit=False,
            parameters=entries,
            governed={},
            initial_state={},
        )

        params.apply_profile(model, pi_profile, check=False)

    applied = {}

    for name in EXPECTED_OVERRIDES:
        parameter_id = naming.resolve_id(
            model, name, kinds=("parameter",)
        )
        parameter = model.getParameter(parameter_id)
        applied[name] = parameter.getValue()

    for name, expected in overrides.items():
        if abs(applied[name] - float(expected)) > 1e-10:
            raise RuntimeError(f"ER override mismatch: {name}")

    return libsbml.writeSBMLToString(document), applied