"""Build deterministic provenance manifest for the ER PI contract."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
CONTRACT = HERE / "er_pi_scenarios.json"
MANIFEST = HERE / "er_pi_bundle_manifest.json"


def canonical_json_bytes(data: dict) -> bytes:
    return json.dumps(
        data,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def main() -> None:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))

    bundle_hash = hashlib.sha256(
        canonical_json_bytes(contract)
    ).hexdigest()

    manifest = {
        "bundle_id": contract["bundle_id"],
        "schema_version": contract["schema_version"],
        "hash_algorithm": "sha256",
        "hash_scope": "canonical JSON content of er_pi_scenarios.json",
        "bundle_hash": bundle_hash,
        "scientific_source": contract["scientific_source"],
    }

    MANIFEST.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print("ER PI bundle manifest")
    print("---------------------")
    print(f"Contract:    {CONTRACT}")
    print(f"Manifest:    {MANIFEST}")
    print(f"Bundle hash: {bundle_hash}")
    print()
    print("RESULT: MANIFEST WRITTEN.")


if __name__ == "__main__":
    main()