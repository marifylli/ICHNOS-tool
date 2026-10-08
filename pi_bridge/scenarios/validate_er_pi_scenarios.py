"""Validate the frozen ER physics-informed scenario contract."""

from __future__ import annotations

import json
from pathlib import Path


CONTRACT = Path(__file__).with_name("er_pi_scenarios.json")

EXPECTED_OVERRIDES = {
    "K_act_er": 929.5,
    "n_er": 4.0,
    "k_on_er": 3.474,
    "k_off_er": 13.674,
    "d_x_er": 1.755,
}


def main() -> None:
    data = json.loads(CONTRACT.read_text(encoding="utf-8"))

    assert data["schema_version"] == "1.0"
    assert data["bundle_id"] == "ichnos_er_pi_scenarios_v1"
    assert data["variant"] == "er"

    source = data["scientific_source"]
    assert source["repository"] == "Ichnos_PULSE"
    assert source["registry"] == "python/pi_scenarios.py"
    assert source["registry_source_commit"] == (
        "f4f07176ddf8fce8f0af3c8374b4811b09c9e297"
    )

    input_model = data["input_model"]
    assert input_model["equation"] == "dS/dt = -k_clear * S"
    assert input_model["k_clear_units"] == "h^-1"

    scenarios_list = data["scenarios"]
    assert len(scenarios_list) == 2

    scenarios = {
        scenario["scenario_id"]: scenario
        for scenario in scenarios_list
    }

    assert set(scenarios) == {"frozen", "er_m2_n4"}

    frozen = scenarios["frozen"]
    assert frozen["parameter_overrides"] == {}
    assert float(frozen["k_clear"]) == 0.0

    er_pi = scenarios["er_m2_n4"]
    assert float(er_pi["k_clear"]) == 0.5032
    assert er_pi["parameter_overrides"] == EXPECTED_OVERRIDES

    interpretation = input_model["interpretation"].lower()
    assert "effective er sensor drive" in interpretation
    assert "physical dtt clearance" in interpretation

    print("ER PI scenario contract validation")
    print("----------------------------------")
    print(f"Contract: {CONTRACT}")
    print(f"Variant:  {data['variant']}")
    print()
    print("frozen     k_clear=0.0000 h^-1  kinetic_overrides=0")
    print("er_m2_n4   k_clear=0.5032 h^-1  kinetic_overrides=5")
    print()
    print("RESULT: CONTRACT VALID.")


if __name__ == "__main__":
    main()