"""Check dose interpolation or decode ratios with the checked interpolant."""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from ichnos.decoder import load_dose_table
from ichnos.continuous_decoder import (
    decode_continuous_dose, load_interpolation_validation,
    save_interpolation_validation, validate_interpolation,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("check")
    check.add_argument("--allowed-error", type=float, required=True)
    decode = commands.add_parser("decode")
    decode.add_argument("--observations", type=Path, required=True)
    decode.add_argument("--validation", type=Path, required=True)
    decode.add_argument("--ratio-tolerance", type=float, required=True)
    decode.add_argument("--green-floor", type=float)
    decode.add_argument("--data-kind", choices=["synthetic", "experimental"], required=True)
    for command in (check, decode):
        command.add_argument("--table", type=Path, required=True)
        command.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.out.exists():
            raise FileExistsError(f"output path already exists: {args.out}")
        table = load_dose_table(args.table)
        if args.command == "check":
            validation = validate_interpolation(table, allowed_error=args.allowed_error)
            save_interpolation_validation(validation, args.out)
            print(f"Interpolation check passed={validation['passed']} -> {args.out}")
        else:
            frame = pd.read_csv(args.observations)
            result = decode_continuous_dose(
                table, times_hours=frame.time_hours.to_numpy(),
                calibrated_ratios=frame.calibrated_ratio_red_green.to_numpy(),
                ratio_tolerance=args.ratio_tolerance,
                interpolation_validation=load_interpolation_validation(args.validation),
                green_measurements=(frame.corrected_mean_green.to_numpy() if args.green_floor is not None else None),
                green_floor=args.green_floor,
            )
            result["data_kind"] = args.data_kind
            result["observations_csv_sha256"] = hashlib.sha256(args.observations.read_bytes()).hexdigest()
            encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
            with args.out.open("x", encoding="utf-8") as handle:
                handle.write(encoded)
            print(f"{result['status']} -> {args.out}")
    except (ValueError, KeyError, OSError, AttributeError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
