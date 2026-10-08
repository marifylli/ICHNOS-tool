
"""D4 deterministic nearest-grid scenario sensitivity maps; no simulations."""

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "d3c_cross_scenario_v1.json"
OUT = ROOT / "results" / "d4_robustness_maps_v1.json"


def main():
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite: {OUT}")

    data = json.loads(SOURCE.read_text(encoding="utf-8"))

    assert data.get("experimentally_validated") is False

    doses = [float(x) for x in data["dose_grid_uM"]]
    times = [float(x) for x in data["time_grid_hours"]]

    assert len(doses) == 6
    assert len(times) == 8
    assert len(set(doses)) == 6
    assert len(set(times)) == 8
    assert all(math.isfinite(x) for x in doses + times)

    comparisons = data["comparisons"]
    assert len(comparisons) == 624

    grouped = {}

    for row in comparisons:
        variant, truth, assumed = (
            row[k]
            for k in (
                "variant",
                "truth_scenario",
                "assumed_scenario",
            )
        )

        key = (variant, truth, assumed)

        assert row["truth_dose_uM"] in doses
        assert row["truth_time_hours"] in times

        point = (
            row["truth_dose_uM"],
            row["truth_time_hours"],
        )

        grid = grouped.setdefault(key, {})

        assert point not in grid, f"Duplicate: {key}, {point}"

        candidates = row["nearest_candidates"]

        assert len(candidates) == row["nearest_count"]
        assert candidates

        residual = float(
            row["minimum_absolute_ratio_residual"]
        )

        assert math.isfinite(residual)
        assert residual >= 0

        assert all(
            c["dose_uM"] in doses
            and c["time_hours"] in times
            for c in candidates
        )

        assert all(
            math.isclose(
                abs(
                    c["ratio_red_green"]
                    - row["truth_ratio_red_green"]
                ),
                residual,
                rel_tol=1e-9,
                abs_tol=1e-12,
            )
            for c in candidates
        )

        grid[point] = row

    assert len(grouped) == 13, (
        f"Unexpected scenario pairs: {len(grouped)}"
    )

    maps = []

    for (variant, truth, assumed), grid in sorted(
        grouped.items()
    ):
        assert len(grid) == 48

        dose_map = []
        time_map = []
        residual_map = []
        ambiguity_map = []

        changed_dose = 0
        changed_time = 0
        ambiguous = 0

        for dose in doses:
            drow = []
            trow = []
            rrow = []
            arow = []

            for time in times:
                row = grid[(dose, time)]
                candidates = row["nearest_candidates"]

                # Multiple equally nearest candidates:
                # do not invent a unique inference.
                if len(candidates) == 1:
                    c = candidates[0]

                    delta_d = (
                        float(c["dose_uM"]) - dose
                    )

                    delta_t = (
                        float(c["time_hours"]) - time
                    )

                    assert math.isclose(
                        delta_d,
                        float(row["nearest_dose_difference_uM"]),
                            rel_tol=0.0,
                            abs_tol=1e-10,
                    )

                    assert math.isclose(
                        delta_t,
                        float(row["nearest_time_difference_hours"]),
                            rel_tol=0.0,
                            abs_tol=1e-10,
                    )

                    changed_dose += delta_d != 0
                    changed_time += delta_t != 0

                else:
                    delta_d = None
                    delta_t = None
                    ambiguous += 1

                drow.append(delta_d)
                trow.append(delta_t)

                rrow.append(
                    float(
                        row["minimum_absolute_ratio_residual"]
                    )
                )

                arow.append(len(candidates))

            dose_map.append(drow)
            time_map.append(trow)
            residual_map.append(rrow)
            ambiguity_map.append(arow)

        maps.append(
            {
                "variant": variant,
                "truth_scenario": truth,
                "assumed_scenario": assumed,
                "signed_delta_dose_uM": dose_map,
                "signed_delta_time_hours": time_map,
                "minimum_absolute_ratio_residual": residual_map,
                "nearest_candidate_count": ambiguity_map,
                "summary": {
                    "points": 48,
                    "changed_dose_count": changed_dose,
                    "changed_time_count": changed_time,
                    "ambiguous_count": ambiguous,
                },
            }
        )

    payload = {
        "analysis_id": "d4_nearest_grid_robustness_v1",
        "source_file": str(SOURCE.relative_to(ROOT)),
        "experimentally_validated": False,
        "posterior_inference": False,
        "interpretation": (
            "Signed nearest-grid changes "
            "(assumed candidate minus synthetic truth); "
            "no noise model, coverage, credible regions, "
            "or acceptance threshold. "
            "Null denotes tied nearest candidates."
        ),
        "dose_grid_uM": doses,
        "time_grid_hours": times,
        "maps": maps,
    }

    OUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Exclusive creation: never overwrite existing results.
    with OUT.open("x", encoding="utf-8") as f:
        json.dump(
            payload,
            f,
            indent=2,
            allow_nan=False,
        )
        f.write("\n")

    print(
        "PASS | 13 maps x 48 points = 624 comparisons"
    )

    for variant in ("ox", "er"):
        cross = [
            m
            for m in maps
            if m["variant"] == variant
            and m["truth_scenario"] != m["assumed_scenario"]
        ]

        time_changed = sum(
            m["summary"]["changed_time_count"]
            for m in cross
        )

        dose_changed = sum(
            m["summary"]["changed_dose_count"]
            for m in cross
        )

        ambiguous_count = sum(
            m["summary"]["ambiguous_count"]
            for m in cross
        )

        print(
            f"{variant.upper()} cross-scenario | "
            f"time_changed={time_changed} | "
            f"dose_changed={dose_changed} | "
            f"ambiguous={ambiguous_count}"
        )

    print(f"WROTE | {OUT}")


if __name__ == "__main__":
    main()

    