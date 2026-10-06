"""Build a model-generated dose table or decode model-scale ratios."""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from ichnos.decoder import build_dose_table, decode_dose, load_dose_table, save_dose_table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    create = commands.add_parser("build")
    create.add_argument("--variant", choices=["ox", "er"], required=True)
    create.add_argument("--profile", default="default")
    create.add_argument("--doses-uM", type=float, nargs="+", required=True)
    create.add_argument("--times-hours", type=float, nargs="+", required=True)
    create.add_argument("--initialization", choices=["equilibrium", "finite-preincubation"], required=True)
    create.add_argument("--preincubation-hours", type=float)
    create.add_argument("--clearance-rate-per-hour", type=float)
    create.add_argument("--out", type=Path, required=True)
    decode = commands.add_parser("decode")
    decode.add_argument("--table", type=Path, required=True)
    decode.add_argument("--observations", type=Path, required=True)
    decode.add_argument("--ratio-tolerance", type=float, required=True)
    decode.add_argument("--green-floor", type=float)
    decode.add_argument("--data-kind", choices=["synthetic", "experimental"], required=True)
    decode.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.out.exists():
            raise FileExistsError(f"output path already exists: {args.out}")
        if args.command == "build":
            table = build_dose_table(
                variant=args.variant, profile_name=args.profile, doses_uM=args.doses_uM,
                times_hours=args.times_hours, initialization=args.initialization,
                preincubation_hours=args.preincubation_hours,
                clearance_rate_per_hour=args.clearance_rate_per_hour,
            )
            save_dose_table(table, args.out)
        else:
            table = load_dose_table(args.table)
            frame = pd.read_csv(args.observations)
            result = decode_dose(
                table, times_hours=frame["time_hours"].to_numpy(),
                calibrated_ratios=frame["calibrated_ratio_red_green"].to_numpy(),
                ratio_tolerance=args.ratio_tolerance,
                green_measurements=(frame["corrected_mean_green"].to_numpy()
                                    if args.green_floor is not None else None),
                green_floor=args.green_floor,
            )
            result["data_kind"] = args.data_kind
            result["observations_csv_sha256"] = hashlib.sha256(
                args.observations.read_bytes()
            ).hexdigest()
            result["warning"] = "model-generated discrete compatibility; not validated experimental dose recovery"
            with args.out.open("x", encoding="utf-8") as handle:
                handle.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
        print(f"Saved {args.command} result to {args.out}")
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
