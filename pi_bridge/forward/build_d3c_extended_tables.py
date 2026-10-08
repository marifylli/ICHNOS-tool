"""Build D3c extended-time forward tables using existing PI builders."""

from __future__ import annotations

import importlib.util
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

TIMES_HOURS = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0]

BUILDERS = {
    "ox": HERE / "build_ox_pi_forward_tables.py",
    "er": HERE / "build_er_pi_forward_tables.py",
}

OUTPUTS = {
    "ox": ROOT / "results" / "ox_pi_forward_tables_d3c.json",
    "er": ROOT / "results" / "er_pi_forward_tables_d3c.json",
}


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)

    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load builder: {path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


def main():
    # Preflight: do not start expensive simulations if outputs exist.
    for variant, output in OUTPUTS.items():
        if output.exists():
            raise FileExistsError(
                f"{variant.upper()} D3c output already exists: {output}"
            )

    modules = {
        variant: load_module(f"d3c_{variant}_builder", path)
        for variant, path in BUILDERS.items()
    }

    # Set the new grid and output paths only in these loaded modules.
    # The original source files remain unchanged.
    modules["ox"].TIMES_HOURS = list(TIMES_HOURS)
    modules["ox"].OUTPUT = OUTPUTS["ox"]

    modules["er"].TIMES = list(TIMES_HOURS)
    modules["er"].OUTPUT = OUTPUTS["er"]

    print("D3c EXTENDED TIME-GRID SIMULATIONS")
    print("----------------------------------")
    print(f"Time grid (hours after stress): {TIMES_HOURS}")
    print("Preincubation: 2 hours")
    print("Dose grid: 25, 30, 35, 40, 45, 50 uM")
    print("Model-generated outputs; not experimentally validated.")

    for variant in ("ox", "er"):
        print(f"\n>>> Running {variant.upper()}", flush=True)
        modules[variant].main()

    print("\nRESULT: D3c EXTENDED FORWARD TABLES WRITTEN.")


if __name__ == "__main__":
    main()
    