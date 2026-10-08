"""Build comparable forward predictions for frozen OX PI scenarios."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from ichnos.decoder import _build_response_table


ROOT = Path(__file__).resolve().parents[1]
SCENARIOS = ROOT / "scenarios"
OUTPUT = ROOT / "results" / "ox_pi_forward_tables.json"

DOSES_UM = [25, 30, 35, 40, 45, 50]
TIMES_HOURS = [0.5, 1.0]


def canonical_bytes(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def main():
    contract = json.loads(
        (SCENARIOS / "ox_pi_scenarios.json").read_text(encoding="utf-8")
    )
    manifest = json.loads(
        (SCENARIOS / "ox_pi_bundle_manifest.json").read_text(encoding="utf-8")
    )

    actual_hash = hashlib.sha256(canonical_bytes(contract)).hexdigest()

    if actual_hash != manifest["bundle_hash"]:
        raise ValueError("PI contract hash does not match manifest")

    if contract["variant"] != "ox":
        raise ValueError("Expected OX scenario bundle")

    expected = {
        "frozen": 0.0,
        "ox_pi_reference": 1.95,
        "ox_pi_systematic": 3.30,
    }
    scenarios = contract["scenarios"]

    if len(scenarios) != len(expected):
        raise ValueError("Unexpected scenario count")

    for scenario in scenarios:
        sid = scenario["scenario_id"]
        if sid not in expected:
            raise ValueError(f"Unexpected scenario: {sid}")
        if float(scenario["k_clear"]) != expected[sid]:
            raise ValueError(f"Unexpected k_clear for {sid}")
        if scenario["parameter_overrides"] != {}:
            raise ValueError(f"Kinetic overrides not supported: {sid}")

    if {s["scenario_id"] for s in scenarios} != set(expected):
        raise ValueError("Missing or duplicate PI scenarios")

    outputs = []

    for scenario in scenarios:
        sid = scenario["scenario_id"]
        k_clear = float(scenario["k_clear"])

        # None means constant input in the existing Tool.
        clearance = None if sid == "frozen" else k_clear

        table, greens = _build_response_table(
            variant="ox",
            doses_uM=DOSES_UM,
            times_hours=TIMES_HOURS,
            initialization="finite-preincubation",
            preincubation_hours=2.0,
            clearance_rate_per_hour=clearance,
        )

        ratios = np.asarray(table.ratios, dtype=float)
        greens = np.asarray(greens, dtype=float)

        expected_shape = (len(DOSES_UM), len(TIMES_HOURS))
        if ratios.shape != expected_shape or greens.shape != expected_shape:
            raise RuntimeError(f"{sid}: unexpected forward table shape")

        if not np.isfinite(ratios).all() or not np.isfinite(greens).all():
            raise RuntimeError(f"{sid}: non-finite model predictions")

        outputs.append({
            "scenario_id": sid,
            "k_clear_per_hour": k_clear,
            "tool_clearance_rate_per_hour": clearance,
            "ratio_red_green": ratios.tolist(),
            "observed_green": greens.tolist(),
            "tool_provenance": table.provenance,
        })

        print(f"{sid:20s} k_clear={k_clear:.2f} shape={ratios.shape}")

    result = {
        "artifact_type": "ox_pi_forward_tables",
        "bundle_id": contract["bundle_id"],
        "bundle_hash": actual_hash,
        "variant": "ox",
        "initialization": "finite-preincubation",
        "preincubation_hours": 2.0,
        "dose_grid_uM": DOSES_UM,
        "time_grid_hours": TIMES_HOURS,
        "observables": {
            "ratio": "model red/green",
            "green": "model Observed_Green",
            "camera_units": False,
        },
        "scenarios": outputs,
    }

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    if OUTPUT.exists():
        raise FileExistsError(
            f"Output already exists: {OUTPUT}. "
            "Preserve the previous result before running again."
        )
    OUTPUT.write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )

    print(f"\nOutput: {OUTPUT}")
    print("RESULT: FORWARD TABLES WRITTEN.")


if __name__ == "__main__":
    main()
    