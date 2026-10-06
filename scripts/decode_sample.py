"""Decode one sample from population summaries with provenance checks."""
import argparse
import json
from pathlib import Path

from ichnos.sample_decoder import TIME_FIELDS, decode_sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("samples-csv", "summary-metadata", "dose-table", "calibration-reference", "out"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--session-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--condition-id")
    parser.add_argument("--time-field", choices=TIME_FIELDS, required=True)
    parser.add_argument("--ratio-tolerance", type=float, required=True)
    parser.add_argument("--mode", choices=["discrete", "continuous", "joint"], default="discrete")
    parser.add_argument("--interpolation-validation", type=Path)
    parser.add_argument("--elapsed-time-grid-hours", type=float, nargs="+")
    args = parser.parse_args()
    try:
        if args.out.exists():
            raise FileExistsError(f"output path already exists: {args.out}")
        options = vars(args).copy()
        output = options.pop("out")
        result = decode_sample(**options)
        encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
        with output.open("x", encoding="utf-8") as handle:
            handle.write(encoded)
        print(f"{result['status']} -> {output}")
    except (ValueError, KeyError, OSError, AttributeError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
