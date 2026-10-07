"""Validate the frozen OX physics-informed scenario contract."""

from __future__ import annotations

import json
from pathlib import Path


CONTRACT = Path(__file__).with_name("ox_pi_scenarios.json")

EXPECTED_K_CLEAR = {
    "frozen": 0.0,
    "ox_pi_reference": 1.95,
    "ox_pi_systematic": 3.30,
}


def main() -> None:
    data = json.loads(CONTRACT.read_text(encoding="utf-8"))

    assert data["schema_version"] == "1.0"
    assert data["variant"] == "ox"

    input_model = data["input_model"]
    assert input_model["type"] == "exponential_effective_input_decay"
    assert input_model["equation"] == "dS/dt = -k_clear * S"
    assert input_model["k_clear_units"] == "h^-1"

    scenarios = {
        item["scenario_id"]: item
        for item in data["scenarios"]
    }

    assert set(scenarios) == set(EXPECTED_K_CLEAR)

    for scenario_id, expected_k_clear in EXPECTED_K_CLEAR.items():
        scenario = scenarios[scenario_id]

        actual = float(scenario["k_clear"])
        assert abs(actual - expected_k_clear) < 1e-12, (
            f"{scenario_id}: expected k_clear={expected_k_clear}, "
            f"got {actual}"
        )

        assert scenario["parameter_overrides"] == {}, (
            f"{scenario_id}: OX PI scenarios must not override "
            "canonical OX kinetic parameters"
        )

    # Scientific interpretation guard.
    interpretation = input_model["interpretation"].lower()

    assert "effective sensor drive" in interpretation
    assert "must not automatically be interpreted" in interpretation
    assert "extracellular h2o2 clearance" in interpretation

    print("OX PI scenario contract validation")
    print("----------------------------------")
    print(f"Contract: {CONTRACT}")
    print(f"Variant:  {data['variant']}")
    print()

    for scenario_id in EXPECTED_K_CLEAR:
        scenario = scenarios[scenario_id]
        print(
            f"{scenario_id:20s} "
            f"k_clear={float(scenario['k_clear']):.2f} h^-1  "
            f"kinetic_overrides={len(scenario['parameter_overrides'])}"
        )

    print()
    print("RESULT: CONTRACT VALID.")
    print(
        "All frozen OX PI scenarios preserve canonical kinetic parameters "
        "and vary only the effective-input decay assumption."
    )


if __name__ == "__main__":
    main()