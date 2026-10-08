"""Generate frozen OX/ER PI scenario contracts reproducibly."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent

SOURCE = {
    "repository": "Ichnos_PULSE",
    "registry": "python/pi_scenarios.py",
    "registry_source_commit": "f4f07176ddf8fce8f0af3c8374b4811b09c9e297",
    "validated_against_repository_head": "56d064f",
}

CONTRACTS = {
    "ox": {
        "schema_version": "1.0",
        "bundle_id": "ichnos_ox_pi_scenarios_v1",
        "variant": "ox",
        "scientific_source": SOURCE,
        "input_model": {
            "type": "exponential_effective_input_decay",
            "equation": "dS/dt = -k_clear * S",
            "k_clear_units": "h^-1",
            "interpretation": (
                "k_clear parameterises decay of effective sensor drive "
                "and must not automatically be interpreted as physical "
                "extracellular H2O2 clearance."
            ),
        },
        "scenarios": [
            {
                "scenario_id": "frozen",
                "parameter_overrides": {},
                "k_clear": 0.0,
                "role": "constant-effective-input baseline",
            },
            {
                "scenario_id": "ox_pi_reference",
                "parameter_overrides": {},
                "k_clear": 1.95,
                "role": "PI reference effective-input decay scenario",
            },
            {
                "scenario_id": "ox_pi_systematic",
                "parameter_overrides": {},
                "k_clear": 3.30,
                "role": (
                    "systematic alternative representing "
                    "quantification-method uncertainty"
                ),
            },
        ],
    },
    "er": {
        "schema_version": "1.0",
        "bundle_id": "ichnos_er_pi_scenarios_v1",
        "variant": "er",
        "scientific_source": SOURCE,
        "input_model": {
            "type": "scenario_specific_effective_input",
            "equation": "dS/dt = -k_clear * S",
            "k_clear_units": "h^-1",
            "interpretation": (
                "k_clear describes decay of effective ER sensor drive "
                "and does not by itself establish physical DTT clearance."
            ),
        },
        "scenarios": [
            {
                "scenario_id": "frozen",
                "parameter_overrides": {},
                "k_clear": 0.0,
                "role": "canonical ER baseline with constant effective input",
            },
            {
                "scenario_id": "er_m2_n4",
                "parameter_overrides": {
                    "K_act_er": 929.5,
                    "n_er": 4.0,
                    "k_on_er": 3.474,
                    "k_off_er": 13.674,
                    "d_x_er": 1.755,
                },
                "k_clear": 0.5032,
                "role": "validated ER M2 physics-informed scenario",
            },
        ],
    },
}

EXPECTED_HASHES = {
    "ox": "0e8ec503b0bdc84e4addab4f26e93a266d4e9b746a3fbb47701e7eb649c5b261",
    "er": "6e54921440af799f178275fc0d856f52d6fa7425e5bdded15fb24c0cc4e1f163",
}


def canonical_bytes(data: dict) -> bytes:
    return json.dumps(
        data,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def main():
    print("PI contract generation")
    print("----------------------")

    # Verify all generated content before writing any file.
    for variant, contract in CONTRACTS.items():
        actual = hashlib.sha256(canonical_bytes(contract)).hexdigest()
        expected = EXPECTED_HASHES[variant]

        print(f"{variant.upper()} hash: {actual}")

        if actual != expected:
            raise RuntimeError(
                f"{variant.upper()} contract hash mismatch; "
                "no files written"
            )

    for variant, contract in CONTRACTS.items():
        destination = HERE / f"{variant}_pi_scenarios.json"

        if destination.exists():
            existing = json.loads(
                destination.read_text(encoding="utf-8")
            )

            if canonical_bytes(existing) != canonical_bytes(contract):
                raise RuntimeError(
                    f"Existing contract differs: {destination}. "
                    "Refusing to overwrite."
                )

            print(f"{variant.upper()}: existing contract verified")
            continue

        with destination.open("x", encoding="utf-8") as file:
            json.dump(contract, file, indent=2, ensure_ascii=False)
            file.write("\n")

        print(f"{variant.upper()}: contract created")

    print("\nRESULT: CONTRACTS VERIFIED OR CREATED.")


if __name__ == "__main__":
    main()