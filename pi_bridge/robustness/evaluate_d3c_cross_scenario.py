"""D3c: extended-time synthetic cross-scenario comparison.

Reuses the validated D3b nearest-grid algorithm.
No new simulations, experimental tolerance, or noise model.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

D3B_SCRIPT = (
    Path(__file__).resolve().parent
    / "evaluate_cross_scenario.py"
)

OUTPUT = RESULTS / "d3c_cross_scenario_v1.json"


def load_d3b():
    spec = importlib.util.spec_from_file_location(
        "d3b_cross_scenario_core",
        D3B_SCRIPT,
    )

    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load D3b evaluator")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


def main():
    if OUTPUT.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing result: {OUTPUT}"
        )

    for variant in ("ox", "er"):
        source = RESULTS / f"{variant}_pi_forward_tables_d3c.json"
        if not source.is_file():
            raise FileNotFoundError(source)

    d3b = load_d3b()

    # Redirect the existing D3b loader to the D3c source files.
    original_load_variant = d3b.load_variant

    def load_d3c_variant(variant):
        original_results = d3b.RESULTS

        try:
            # D3b's loader constructs the filename internally.
            # Reuse its implementation through a temporary Path adapter.
            class D3cResultsPath:
                def __truediv__(self, filename):
                    return RESULTS / filename.replace(
                        "_pi_forward_tables.json",
                        "_pi_forward_tables_d3c.json",
                    )

            d3b.RESULTS = D3cResultsPath()
            return original_load_variant(variant)
        finally:
            d3b.RESULTS = original_results

    d3b.load_variant = load_d3c_variant

    all_rows = []
    all_summary = []

    for variant in ("ox", "er"):
        rows, summary = d3b.evaluate_variant(variant)
        all_rows.extend(rows)
        all_summary.extend(summary)

    output = {
        "analysis_id": "d3c_synthetic_cross_scenario_v1",
        "analysis_type": "deterministic_nearest_grid_comparison",
        "observation": "model_generated_red_green",
        "experimentally_validated": False,
        "time_grid_hours": [
            0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0
        ],
        "dose_grid_uM": [25, 30, 35, 40, 45, 50],
        "source_files": {
            "ox": "ox_pi_forward_tables_d3c.json",
            "er": "er_pi_forward_tables_d3c.json",
        },
        "time_grid_interpretation": (
            "Discrete comparison at stored time points only; "
            "not continuous time inference."
        ),
        "matching_rule": (
            "Minimum absolute Red/Green residual; exact floating-point "
            "ties only; no measurement tolerance."
        ),
        "status_interpretation": {
            "exact_match": "Exactly equal stored numerical ratio.",
            "nearest_only": (
                "Unique nearest candidate; not an accepted "
                "measurement-level match."
            ),
            "tied_nearest": (
                "Multiple equally nearest candidates; "
                "no unique dose/time assignment."
            ),
        },
        "summary": all_summary,
        "comparisons": all_rows,
    }

    with OUTPUT.open("x", encoding="utf-8") as file:
        json.dump(output, file, indent=2, allow_nan=False)
        file.write("\n")

    print("D3c EXTENDED CROSS-SCENARIO EVALUATION")
    print("--------------------------------------")

    for row in all_summary:
        print(
            f"{row['variant'].upper():2} "
            f"{row['truth_scenario']:18} -> "
            f"{row['assumed_scenario']:18} | "
            f"exact={row['exact_match_count']:2} "
            f"nearest={row['nearest_only_count']:2} "
            f"ties={row['tied_nearest_count']:2} "
            f"time_changed={row['nearest_time_changed_count']:2} "
            f"dose_changed={row['nearest_dose_changed_count']:2}"
        )

    print(f"\nRESULT: {OUTPUT}")
    print("NOTE: Nearest-grid comparisons are not validated decoding.")


if __name__ == "__main__":
    main()
    