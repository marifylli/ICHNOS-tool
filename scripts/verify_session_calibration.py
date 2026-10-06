"""Verify relative session calibration with synthetic ratios, not lab images."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from ichnos.calibration import (
    calibrate_image_ratios, fit_session_calibration,
    model_reference_id, save_session_calibration,
)


def run(args):
    if args.out_dir.exists():
        raise FileExistsError(f"output path already exists: {args.out_dir}")
    if not np.isfinite(args.scale) or args.scale <= 0:
        raise ValueError("scale must be positive and finite")
    if not np.isfinite(args.relative_noise) or not 0 <= args.relative_noise < 1:
        raise ValueError("relative-noise must be in [0, 1)")
    if args.seed < 0:
        raise ValueError("seed must be non-negative")
    metadata = json.loads((args.protocol_dir / "metadata.json").read_text())
    reference_id = model_reference_id(metadata)
    with (args.protocol_dir / "fluorescence.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    model_ratios = np.array([float(row["ratio_red_green"]) for row in rows])
    times = np.array([float(row["time_hours"]) for row in rows])
    if len(rows) < 6:
        raise ValueError("need at least six observation times for fitting and holdout")
    if times.tolist() != metadata["observation_times_hours"]:
        raise ValueError("CSV times differ from protocol metadata")
    if not np.isfinite(model_ratios).all() or (model_ratios <= 0).any():
        raise ValueError("synthetic verification requires positive finite model ratios")

    # Bounded zero-mean multiplicative noise; no clipping or altered samples.
    rng = np.random.default_rng(args.seed)
    noise = rng.uniform(-args.relative_noise, args.relative_noise, len(rows))
    image_ratios = args.scale * model_ratios * (1 + noise)
    fit_mask = np.arange(len(rows)) % 2 == 0
    session_id = "synthetic-verification"
    calibration = fit_session_calibration(
        model_ratios[fit_mask], image_ratios[fit_mask],
        session_id=session_id, reference_id=reference_id,
        source_kind="synthetic",
        source=f"model-generated ratios; seed={args.seed}; alternating fit/holdout times",
    )
    corrected = calibrate_image_ratios(
        image_ratios, calibration,
        session_id=session_id, reference_id=reference_id,
    )
    holdout_rmse = float(np.sqrt(np.mean(
        (corrected[~fit_mask] - model_ratios[~fit_mask]) ** 2
    )))
    verification = {
        "purpose": "verify scale recovery and correction in software with synthetic ratios",
        "why_synthetic": "known ground truth is available; experimental session references are not supplied",
        "source_kind": "synthetic",
        "experimental_instrument_validated": False,
        "generation_equation": "image_ratio = known_scale * model_ratio * (1 + noise)",
        "noise_distribution": "uniform[-relative_noise, +relative_noise]",
        "seed": args.seed,
        "known_scale": args.scale,
        "relative_noise": args.relative_noise,
        "estimated_scale": calibration.c_session,
        "scale_relative_error": abs(calibration.c_session / args.scale - 1),
        "holdout_rmse_model_ratio_units": holdout_rmse,
        "n_fit": int(fit_mask.sum()),
        "n_holdout": int((~fit_mask).sum()),
        "split": "alternating observation times; holdout is not used in fitting",
        "model_reference_id": reference_id,
        "protocol_metadata": metadata,
        "limitations": [
            "synthetic ratios, not synthetic microscope images",
            "not biological validation or experimental uncertainty coverage",
            "assumes one multiplicative ratio scale with zero intercept",
            "does not establish instrument linearity, detection floor or absolute intensity calibration",
            "f is already included in Measured_Ratio_RG; not independently estimated",
        ],
    }
    encoded = json.dumps(verification, indent=2, allow_nan=False) + "\n"
    args.out_dir.mkdir(parents=True)
    save_session_calibration(calibration, args.out_dir / "calibration.json")
    (args.out_dir / "verification.json").write_text(encoded, encoding="utf-8")
    with (args.out_dir / "synthetic_ratios.csv").open("x", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "time_hours", "split", "model_ratio_red_green",
            "image_ratio_red_green", "calibrated_ratio_red_green",
        ])
        for index, time in enumerate(times):
            writer.writerow([
                float(time), "fit" if fit_mask[index] else "holdout",
                float(model_ratios[index]), float(image_ratios[index]),
                float(corrected[index]),
            ])
    print(f"Synthetic verification saved to {args.out_dir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--scale", type=float, default=1.5)
    parser.add_argument("--relative-noise", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=20261005)
    args = parser.parse_args()
    try:
        run(args)
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
