"""Summarize image-pipeline cell CSVs and optionally apply session calibration."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from ichnos.calibration import load_session_calibration
from ichnos.population import CalibrationBinding, summarize_cells


def run(args):
    if args.out_dir.exists():
        raise FileExistsError(f"output path already exists: {args.out_dir}")
    cells = pd.read_csv(
        args.cells,
        dtype={"session_id": str, "sample_id": str, "condition_id": str},
    )
    calibrations = None
    if args.calibrations is not None:
        manifest = pd.read_csv(args.calibrations, dtype=str)
        required = {"session_id", "reference_id", "calibration_path"}
        if not required.issubset(manifest.columns):
            raise ValueError("calibration manifest requires session_id, reference_id, calibration_path")
        if manifest[list(required)].isna().any().any():
            raise ValueError("calibration manifest values must not be missing")
        if manifest.session_id.duplicated().any():
            raise ValueError("duplicate calibration session_id")
        calibrations = {}
        for row in manifest.itertuples():
            path = Path(row.calibration_path)
            if not path.is_absolute():
                path = args.calibrations.parent / path
            calibration = load_session_calibration(path)
            calibrations[row.session_id] = CalibrationBinding(calibration, row.reference_id)
    summaries = summarize_cells(
        cells, data_kind=args.data_kind, calibrations=calibrations,
        min_cells=args.min_cells, green_floor=args.green_floor,
    )
    metadata = {
        "schema_version": 1,
        "input_cell_csv_sha256": hashlib.sha256(args.cells.read_bytes()).hexdigest(),
        "data_kind": args.data_kind,
        "group_by": ["session_id", "sample_id", "condition_id", "timepoint"],
        "method": "median of QC-passing finite per-cell red/green ratios",
        "quartiles": "descriptive cell distribution; not confidence intervals",
        "cells_are_independent_biological_replicates": False,
        "min_cells": args.min_cells,
        "green_floor_corrected_image_units": args.green_floor,
        "green_floor_is_experimentally_calibrated": False,
        "calibration_equation": "calibrated_ratio = original_image_ratio / c_session",
        "f_applied_again": False,
        "experimental_instrument_validated": False,
        "calibrations": {
            session: binding.calibration.to_dict()
            for session, binding in (calibrations or {}).items()
        },
        "limitations": [
            "does not establish camera linearity, model validity or experimental accuracy",
            "acquisition settings must match those represented by session references",
            "unknown elapsed times remain unknown; timepoint is not substituted for hours",
        ],
    }
    encoded = json.dumps(metadata, indent=2, allow_nan=False) + "\n"
    args.out_dir.mkdir(parents=True)
    summaries.to_csv(args.out_dir / "samples.csv", index=False)
    (args.out_dir / "metadata.json").write_text(encoded, encoding="utf-8")
    print(f"{len(cells)} cell rows -> {len(summaries)} sample summaries -> {args.out_dir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--data-kind", choices=["synthetic", "experimental"], required=True)
    parser.add_argument("--calibrations", type=Path)
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--green-floor", type=float, default=0.0)
    args = parser.parse_args()
    try:
        run(args)
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
