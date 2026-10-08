
"""D5: compare initial-dose and integrated-effective-exposure robustness.

Reads existing D3c synthetic nearest-grid comparisons.
No simulations, no experimental calibration, no file overwrites.
"""

import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "d3c_cross_scenario_v1.json"
OUTPUT = ROOT / "results" / "d5_effective_exposure_v1.json"

K_CLEAR = {
    "ox": {
        "frozen": 0.0,
        "ox_pi_reference": 1.95,
        "ox_pi_systematic": 3.30,
    },
    "er": {
        "frozen": 0.0,
        "er_m2_n4": 0.5032,
    },
}

EXPECTED_SCENARIO_PAIRS = {"ox": 9, "er": 4}
EXPECTED_COUNTS = {"ox": 432, "er": 192}


def exposure(dose, time, k):
    """Integral of dose * exp(-k*tau) from 0 to time."""
    assert dose > 0 and math.isfinite(dose)
    assert time >= 0 and math.isfinite(time)
    assert k >= 0 and math.isfinite(k)

    if k == 0:
        return dose * time

    return dose * (-math.expm1(-k * time)) / k


def main():
    if OUTPUT.exists():
        raise SystemExit(f"Refusing to overwrite: {OUTPUT}")

    data = json.loads(SOURCE.read_text(encoding="utf-8"))

    assert data["experimentally_validated"] is False
    assert len(data["comparisons"]) == 624

    doses = [float(x) for x in data["dose_grid_uM"]]
    times = [float(x) for x in data["time_grid_hours"]]

    assert len(doses) == 6 and len(set(doses)) == 6
    assert len(times) == 8 and len(set(times)) == 8
    assert all(x > 0 and math.isfinite(x) for x in doses + times)

    records = []
    grouped = defaultdict(list)
    seen = set()

    for row in data["comparisons"]:
        variant = row["variant"]
        truth = row["truth_scenario"]
        assumed = row["assumed_scenario"]

        assert variant in K_CLEAR
        assert truth in K_CLEAR[variant]
        assert assumed in K_CLEAR[variant]

        d_true = float(row["truth_dose_uM"])
        t_true = float(row["truth_time_hours"])

        assert d_true in doses and t_true in times

        key = (variant, truth, assumed, d_true, t_true)
        assert key not in seen, f"Duplicate comparison: {key}"
        seen.add(key)

        candidates = row["nearest_candidates"]
        assert len(candidates) == row["nearest_count"]
        assert len(candidates) == 1, (
            "D5 requires unique nearest candidates; "
            f"found {len(candidates)} for {key}"
        )

        candidate = candidates[0]
        d_est = float(candidate["dose_uM"])
        t_est = float(candidate["time_hours"])

        assert d_est in doses and t_est in times

        k_true = K_CLEAR[variant][truth]
        k_assumed = K_CLEAR[variant][assumed]

        e_true = exposure(d_true, t_true, k_true)
        e_est = exposure(d_est, t_est, k_assumed)

        assert e_true > 0 and e_est > 0

        # Absolute log ratios are dimensionless and symmetric
        # under exchanging estimated and true quantities.
        dose_log_error = abs(math.log(d_est / d_true))
        exposure_log_error = abs(math.log(e_est / e_true))

        improvement = dose_log_error - exposure_log_error

        record = {
            "variant": variant,
            "truth_scenario": truth,
            "assumed_scenario": assumed,
            "truth_dose_uM": d_true,
            "truth_time_hours": t_true,
            "nearest_dose_uM": d_est,
            "nearest_time_hours": t_est,
            "truth_exposure_uM_hours": e_true,
            "estimated_exposure_uM_hours": e_est,
            "absolute_dose_log_ratio": dose_log_error,
            "absolute_exposure_log_ratio": exposure_log_error,
            "exposure_improvement_log_units": improvement,
            "minimum_absolute_ratio_residual": float(
                row["minimum_absolute_ratio_residual"]
            ),
        }

        assert all(
            math.isfinite(v)
            for v in (
                e_true,
                e_est,
                dose_log_error,
                exposure_log_error,
                improvement,
                record["minimum_absolute_ratio_residual"],
            )
        )

        records.append(record)
        grouped[(variant, truth, assumed)].append(record)

    for variant, expected in EXPECTED_COUNTS.items():
        assert sum(r["variant"] == variant for r in records) == expected

    summary = []

    for (variant, truth, assumed), rows in sorted(grouped.items()):
        assert len(rows) == 48

        # A strict comparison without arbitrary tolerance.
        better = sum(
            r["exposure_improvement_log_units"] > 0 for r in rows
        )
        worse = sum(
            r["exposure_improvement_log_units"] < 0 for r in rows
        )
        equal = len(rows) - better - worse

        summary.append({
            "variant": variant,
            "truth_scenario": truth,
            "assumed_scenario": assumed,
            "n": len(rows),
            "exposure_better_count": better,
            "dose_better_count": worse,
            "equal_count": equal,
            "mean_absolute_dose_log_ratio": sum(
                r["absolute_dose_log_ratio"] for r in rows
            ) / len(rows),
            "mean_absolute_exposure_log_ratio": sum(
                r["absolute_exposure_log_ratio"] for r in rows
            ) / len(rows),
        })

    for variant, expected in EXPECTED_SCENARIO_PAIRS.items():
        assert sum(s["variant"] == variant for s in summary) == expected

    payload = {
        "analysis_id": "d5_integrated_effective_exposure_v1",
        "source_file": SOURCE.name,
        "experimentally_validated": False,
        "posterior_inference": False,
        "exposure_model": "S_eff(t) = D0 * exp(-k_clear*t)",
        "exposure_definition": "Integral of S_eff from stress onset to observation",
        "exposure_units": "uM*hours",
        "comparison_metric": (
            "Absolute natural log of estimated/true quantity; "
            "smaller is better"
        ),
        "limitations": [
            "Synthetic nearest-grid analysis only",
            "No experimental noise or calibration",
            "No posterior credible regions or coverage",
            "Assumes effective drive amplitude equals input dose",
            "Scenario k_clear values are frozen contract constants",
            "No uncertainty threshold or statistical significance",
        ],
        "k_clear_per_hour": K_CLEAR,
        "summary": summary,
        "comparisons": records,
    }

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT.open("x", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, allow_nan=False)
        f.write("\n")

    print("PASS | D5 input integrity and exposure calculations")
    print(f"PASS | {len(records)} comparisons, {len(summary)} scenario pairs")

    for variant in ("ox", "er"):
        cross = [
            s for s in summary
            if s["variant"] == variant
            and s["truth_scenario"] != s["assumed_scenario"]
        ]

        better = sum(s["exposure_better_count"] for s in cross)
        worse = sum(s["dose_better_count"] for s in cross)
        equal = sum(s["equal_count"] for s in cross)

        print(
            f"{variant.upper()} cross-scenario | "
            f"exposure_better={better} | "
            f"dose_better={worse} | equal={equal}"
        )

    print(f"WROTE | {OUTPUT}")


if __name__ == "__main__":
    main()

