"""D3b: deterministic synthetic cross-scenario grid comparison.

This is a descriptive nearest-grid analysis, not a validated decoder.
No measurement tolerance, noise model, or experimental calibration is assumed.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

# Only exact equality is accepted. No arbitrary measurement tolerance.
# This is not a claim of physical or experimental indistinguishability.


def load_variant(variant):
    path = RESULTS / f"{variant}_pi_forward_tables.json"
    data = json.loads(path.read_text(encoding="utf-8"))

    scenarios = {}

    for scenario in data["scenarios"]:
        scenario_id = scenario["scenario_id"]

        if variant == "ox":
            doses = data["dose_grid_uM"]
            times = data["time_grid_hours"]
            ratios = scenario["ratio_red_green"]
        else:
            doses = scenario["doses_uM"]
            times = scenario["times_hours"]
            ratios = scenario["ratios"]

        if len(ratios) != len(doses):
            raise ValueError(f"{variant}/{scenario_id}: dose shape mismatch")

        points = []

        for i, dose in enumerate(doses):
            if len(ratios[i]) != len(times):
                raise ValueError(
                    f"{variant}/{scenario_id}: time shape mismatch"
                )

            for j, time in enumerate(times):
                ratio = float(ratios[i][j])

                if not math.isfinite(ratio):
                    raise ValueError(
                        f"{variant}/{scenario_id}: nonfinite ratio"
                    )

                points.append({
                    "dose_uM": float(dose),
                    "time_hours": float(time),
                    "ratio_red_green": ratio,
                })

        if scenario_id in scenarios:
            raise ValueError(f"Duplicate scenario ID: {scenario_id}")

        scenarios[scenario_id] = points

    return scenarios


def nearest_match(observation, candidates):
    distances = [
        abs(
            observation["ratio_red_green"]
            - candidate["ratio_red_green"]
        )
        for candidate in candidates
    ]

    minimum = min(distances)

    nearest = [
        candidate
        for candidate, distance in zip(candidates, distances)
        if distance == minimum
    ]

    if len(nearest) > 1:
        status = "tied_nearest"
    elif minimum == 0.0:
        status = "exact_match"
    else:
        status = "nearest_only"

    result = {
        "status": status,
        "minimum_absolute_ratio_residual": minimum,
        "nearest_count": len(nearest),
        "nearest_candidates": nearest,
    }

    if len(nearest) == 1:
        candidate = nearest[0]
        result["nearest_dose_difference_uM"] = (
            candidate["dose_uM"] - observation["dose_uM"]
        )
        result["nearest_time_difference_hours"] = (
            candidate["time_hours"] - observation["time_hours"]
        )
    else:
        result["nearest_dose_difference_uM"] = None
        result["nearest_time_difference_hours"] = None

    return result


def evaluate_variant(variant):
    scenarios = load_variant(variant)
    rows = []
    summary = []

    for truth_id, observations in scenarios.items():
        for assumed_id, candidates in scenarios.items():
            pair_results = []

            for observation in observations:
                match = nearest_match(observation, candidates)

                row = {
                    "variant": variant,
                    "truth_scenario": truth_id,
                    "assumed_scenario": assumed_id,
                    "truth_dose_uM": observation["dose_uM"],
                    "truth_time_hours": observation["time_hours"],
                    "truth_ratio_red_green": observation["ratio_red_green"],
                    **match,
                }

                rows.append(row)
                pair_results.append(row)

            summary.append({
                "variant": variant,
                "truth_scenario": truth_id,
                "assumed_scenario": assumed_id,
                "n_observations": len(pair_results),
                "exact_match_count": sum(
                    r["status"] == "exact_match"
                    for r in pair_results
                ),
                "nearest_only_count": sum(
                    r["status"] == "nearest_only"
                    for r in pair_results
                ),
                "tied_nearest_count": sum(
                    r["status"] == "tied_nearest"
                    for r in pair_results
                ),
                "maximum_nearest_ratio_residual": max(
                    r["minimum_absolute_ratio_residual"]
                    for r in pair_results
                ),
                "nearest_time_changed_count": sum(
                    r["nearest_time_difference_hours"] not in (None, 0.0)
                    for r in pair_results
                ),
                "nearest_dose_changed_count": sum(
                    r["nearest_dose_difference_uM"] not in (None, 0.0)
                    for r in pair_results
                ),
            })

    return rows, summary


def main():
    all_rows = []
    all_summary = []

    for variant in ("ox", "er"):
        rows, summary = evaluate_variant(variant)
        all_rows.extend(rows)
        all_summary.extend(summary)

    output = {
        "analysis_id": "d3b_synthetic_cross_scenario_v1",
        "analysis_type": "deterministic_nearest_grid_comparison",
        "observation": "model_generated_red_green",
        "experimentally_validated": False,
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

    destination = RESULTS / "d3b_cross_scenario_v1.json"

    # Preserve previous results; never silently overwrite.
    with destination.open("x", encoding="utf-8") as file:
        json.dump(output, file, indent=2, allow_nan=False)
        file.write("\n")

    print("D3b SYNTHETIC CROSS-SCENARIO EVALUATION")
    print("---------------------------------------")

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

    print(f"\nRESULT: {destination}")
    print("NOTE: Nearest-grid comparisons are not validated decoding.")


if __name__ == "__main__":
    main()
    