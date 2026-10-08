"""Build provenance-tracked ER PI forward tables."""

from __future__ import annotations

import hashlib
import json
import platform
from dataclasses import asdict
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import libsbml
import numpy as np

from ichnos import naming, protocol
from ichnos.fluorescence import model_fluorescence

from er_runtime_adapter import prepare_er_scenario_sbml


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "scenarios" / "er_pi_scenarios.json"
MANIFEST = ROOT / "scenarios" / "er_pi_bundle_manifest.json"
OUTPUT = ROOT / "results" / "er_pi_forward_tables.json"

DOSES = [25, 30, 35, 40, 45, 50]
TIMES = [0.5, 1.0]
PREINCUBATION_HOURS = 2.0

EXPECTED_SCENARIOS = {
    "frozen": {"overrides": {}, "k_clear": 0.0},
    "er_m2_n4": {
        "overrides": {
            "K_act_er": 929.5,
            "n_er": 4.0,
            "k_on_er": 3.474,
            "k_off_er": 13.674,
            "d_x_er": 1.755,
        },
        "k_clear": 0.5032,
    },
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_bytes(data: dict) -> bytes:
    return json.dumps(
        data,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def load_verified_contract():
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))

    actual_hash = sha256_bytes(canonical_bytes(contract))

    if actual_hash != manifest["bundle_hash"]:
        raise RuntimeError("ER contract hash does not match manifest")

    if contract["variant"] != "er":
        raise RuntimeError("Incorrect contract variant")

    if manifest["bundle_id"] != contract["bundle_id"]:
        raise RuntimeError("ER bundle ID mismatch")

    if manifest["schema_version"] != contract["schema_version"]:
        raise RuntimeError("ER schema version mismatch")

    scenarios = {
        s["scenario_id"]: s
        for s in contract["scenarios"]
    }

    if len(contract["scenarios"]) != len(scenarios):
        raise RuntimeError("Duplicate ER scenario IDs")

    if set(scenarios) != set(EXPECTED_SCENARIOS):
        raise RuntimeError("Unexpected ER scenario set")

    for scenario_id, expected in EXPECTED_SCENARIOS.items():
        scenario = scenarios[scenario_id]

        if scenario["parameter_overrides"] != expected["overrides"]:
            raise RuntimeError(
                f"Unexpected kinetic overrides: {scenario_id}"
            )

        if scenario["k_clear"] != expected["k_clear"]:
            raise RuntimeError(
                f"Unexpected clearance value: {scenario_id}"
            )

    return contract, manifest, actual_hash


def build_scenario(scenario: dict, bundle_hash: str) -> dict:
    scenario_id = scenario["scenario_id"]

    baseline_sbml, applied = prepare_er_scenario_sbml(scenario)
    baseline_hash = sha256_bytes(baseline_sbml.encode("utf-8"))

    k_clear = (
        None if scenario_id == "frozen"
        else float(scenario["k_clear"])
    )

    ratio_rows = []
    green_rows = []
    run_records = []

    for dose in DOSES:
        exposure = protocol.StressProtocol(
            variant="er",
            dose=float(dose),
            clears=k_clear is not None,
            clearance_rate_per_hour=k_clear,
        )

        runner, loaded_model = protocol.load_protocol_model(
            baseline_sbml, exposure
        )

        prepared_hash = sha256_bytes(
            runner.getCurrentSBML().encode("utf-8")
        )

        names = naming.id_to_name_map(loaded_model)

        result, _ = protocol.run_protocol_after_preincubation(
            runner,
            exposure,
            times_hours=np.asarray(TIMES, dtype=float),
            id_to_name=names,
            preincubation_hours=PREINCUBATION_HOURS,
        )

        fluorescence = model_fluorescence(result)

        ratios = np.asarray(
            fluorescence["ratio_red_green"], dtype=float
        )
        greens = np.asarray(
            fluorescence["green"], dtype=float
        )

        if ratios.shape != (len(TIMES),):
            raise RuntimeError("Unexpected ER ratio output shape")

        if greens.shape != (len(TIMES),):
            raise RuntimeError("Unexpected ER green output shape")

        if not np.isfinite(ratios).all():
            raise RuntimeError("Non-finite ER ratio prediction")

        if not np.isfinite(greens).all():
            raise RuntimeError("Non-finite ER green prediction")

        ratio_rows.append(ratios.tolist())
        green_rows.append(greens.tolist())

        run_records.append({
            "dose_uM": float(dose),
            "exposure": result.exposure,
            "initialization": result.initialization,
            "solver_used": asdict(result.solver),
            "prepared_model_sbml_sha256": prepared_hash,
        })

    return {
        "scenario_id": scenario_id,
        "variant": "er",
        "source_kind": "model_generated",
        "experimentally_validated": False,
        "bundle_hash": bundle_hash,
        "parameter_overrides": scenario["parameter_overrides"],
        "applied_kinetic_parameters": applied,
        "k_clear_per_hour": scenario["k_clear"],
        "effective_input": (
            "constant" if k_clear is None
            else "exponential_effective_input_decay"
        ),
        "doses_uM": DOSES,
        "times_hours": TIMES,
        "initialization": "finite-preincubation",
        "preincubation_hours": PREINCUBATION_HOURS,
        "ratio_direction": "red/green",
        "ratio_observable": "Measured_Ratio_RG",
        "green_observable": "Observed_Green",
        "ratios": ratio_rows,
        "greens": green_rows,
        "provenance": {
            "model_sbml_sha256": baseline_hash,
            "runs": run_records,
            "runtime": {
                "python": platform.python_version(),
                "packages": {
                    name: package_version(name)
                    for name in (
                        "numpy",
                        "python-libsbml",
                        "libroadrunner",
                        "tellurium",
                        "PyYAML",
                    )
                },
            },
            "code_sha256": {
                name: sha256_bytes(
                    (ROOT / "forward" / name).read_bytes()
                )
                for name in (
                    "build_er_pi_forward_tables.py",
                    "er_runtime_adapter.py",
                )
            },
            "scope": (
                "Discrete model-generated dose/time predictions; "
                "not calibrated experimental measurements"
            ),
        },
    }


def main():
    if OUTPUT.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing output: {OUTPUT}"
        )

    contract, manifest, bundle_hash = load_verified_contract()

    print("ER PI forward simulations")
    print("-------------------------")
    print(f"Bundle hash: {bundle_hash}")

    results = []

    for scenario in contract["scenarios"]:
        scenario_id = scenario["scenario_id"]
        print(f"\nRunning {scenario_id}...")

        result = build_scenario(scenario, bundle_hash)
        results.append(result)

        print(
            f"  Ratios: {len(result['ratios'])} x "
            f"{len(result['ratios'][0])}"
        )
        print(
            f"  Greens: {len(result['greens'])} x "
            f"{len(result['greens'][0])}"
        )

    payload = {
        "schema_version": "1.0",
        "bundle_id": contract["bundle_id"],
        "bundle_hash": bundle_hash,
        "scientific_source": contract["scientific_source"],
        "input_model": contract["input_model"],
        "source_kind": "model_generated",
        "experimentally_validated": False,
        "scenarios": results,
    }

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT.open("x", encoding="utf-8") as file:
        json.dump(payload, file, indent=2, ensure_ascii=False)
        file.write("\n")

    print(f"\nOutput: {OUTPUT}")
    print("RESULT: ER FORWARD TABLES WRITTEN.")


if __name__ == "__main__":
    main()