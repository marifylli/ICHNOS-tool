"""Check ER PI overrides on an in-memory SBML model."""

import json
from pathlib import Path

import libsbml

from ichnos import build, naming, params


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "scenarios" / "er_pi_scenarios.json"


def main():
    data = json.loads(CONTRACT.read_text(encoding="utf-8"))

    scenario = next(
        s for s in data["scenarios"]
        if s["scenario_id"] == "er_m2_n4"
    )

    overrides = scenario["parameter_overrides"]

    sbml = build.build_variant_sbml_string("er", save_sbml=False)
    document = libsbml.readSBMLFromString(sbml)
    model = document.getModel()

    # Apply the Tool's existing baseline profile first.
    baseline_profile = params.load_profile("er", "default")
    params.apply_profile(model, baseline_profile)

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
            note="Frozen ER M2 scenario override",
        )

    pi_profile = params.Profile(
        variant="er",
        name="er_m2_n4",
        description="Frozen ER physics-informed runtime overrides",
        joint_fit=False,
        parameters=entries,
        governed={},
        initial_state={},
    )

    # Partial profile: retain Tool validation of individual entries.
    params.apply_profile(model, pi_profile, check=False)

    print("\nER runtime override diagnostic")
    print("------------------------------")

    for name, expected in overrides.items():
        parameter_id = naming.resolve_id(
            model, name, kinds=("parameter",)
        )
        actual = model.getParameter(parameter_id).getValue()

        print(f"{name:12s} expected={expected} actual={actual}")

        if abs(actual - float(expected)) > 1e-10:
            raise RuntimeError(f"Override failed: {name}")

    print("\nRESULT: ER OVERRIDES APPLIED SUCCESSFULLY.")


if __name__ == "__main__":
    main()
    