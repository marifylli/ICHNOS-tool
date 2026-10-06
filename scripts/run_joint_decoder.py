"""Decode discrete dose and elapsed time from already calibrated ratios."""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd
from ichnos.decoder import load_dose_table
from ichnos.joint_decoder import decode_joint_dose_time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--elapsed-times-hours", type=float, nargs="+", required=True)
    parser.add_argument("--ratio-tolerance", type=float, required=True)
    parser.add_argument("--green-floor", type=float)
    parser.add_argument("--data-kind", choices=["synthetic", "experimental"], required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.out.exists():
            raise FileExistsError(f"output path already exists: {args.out}")
        frame = pd.read_csv(args.observations)
        result = decode_joint_dose_time(
            load_dose_table(args.table), relative_times_hours=frame["relative_time_hours"].to_numpy(),
            elapsed_time_grid_hours=args.elapsed_times_hours,
            calibrated_ratios=frame["calibrated_ratio_red_green"].to_numpy(),
            ratio_tolerance=args.ratio_tolerance,
            green_measurements=(frame["corrected_mean_green"].to_numpy() if args.green_floor is not None else None),
            green_floor=args.green_floor,
        )
        result.update(data_kind=args.data_kind, observations_csv_sha256=hashlib.sha256(args.observations.read_bytes()).hexdigest())
        encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
        with args.out.open("x", encoding="utf-8") as handle:
            handle.write(encoded)
        print(f"Saved joint decoder result to {args.out}")
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
