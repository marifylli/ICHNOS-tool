"""Generate and validate both frozen PI scenario bundles."""

from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent

SCRIPTS = [
    "generate_pi_contracts.py",
    "build_bundle_manifest.py",
    "build_er_bundle_manifest.py",
    "validate_ox_pi_scenarios.py",
    "validate_er_pi_scenarios.py",
]


def main():
    for script in SCRIPTS:
        print(f"\n>>> {script}", flush=True)

        subprocess.run(
            [sys.executable, str(HERE / script)],
            check=True,
        )

    print("\nRESULT: BOTH PI BUNDLES PREPARED AND VALIDATED.")


if __name__ == "__main__":
    main()