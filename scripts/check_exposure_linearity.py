"""Step 2a: is the cell-free background linear in exposure?

    python3 scripts/check_exposure_linearity.py \\
        --manifest data/acquisition/manifest/images.csv \\
        --out-dir outputs/linearity_20261010

Reads every green and red frame of the manifest once (no segmentation),
records its background (histogram mode, 1st percentile, median) and fits
background = offset + rate * exposure per session and pooled.

Writes into a new --out-dir:
  fields.csv   one row per field and channel with its background statistics
  fits.csv     the fits, per session and pooled, per channel and statistic
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ichnos_image.exposure import field_background, fit_background_vs_exposure


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    if args.out_dir.exists():
        raise SystemExit(f"output path already exists: {args.out_dir}")

    manifest = pd.read_csv(args.manifest)
    records = []
    for number, row in enumerate(manifest.itertuples(), start=1):
        for channel in ("green", "red"):
            stats = field_background(getattr(row, f"{channel}_path"))
            records.append({"session_id": row.session_id, "sample_id": row.sample_id,
                            "acquisition_order": row.acquisition_order, "channel": channel,
                            "exposure_ms": float(getattr(row, f"exposure_ms_{channel}")), **stats})
        if number % 50 == 0:
            print(f"{number}/{len(manifest)} fields read", flush=True)
    fields = pd.DataFrame(records)

    fits = pd.concat([fit_background_vs_exposure(fields[fields.channel == channel], value).assign(channel=channel)
                      for channel in ("green", "red") for value in ("mode", "p01")], ignore_index=True)

    args.out_dir.mkdir(parents=True)
    fields.to_csv(args.out_dir / "fields.csv", index=False)
    fits.to_csv(args.out_dir / "fits.csv", index=False)

    shown = fits[fits.value == "mode"]
    with pd.option_context("display.width", 160, "display.max_columns", 20):
        print(shown.reindex(columns=["channel", "session_id", "n_fields", "exposures_ms", "fit", "mean",
                                     "offset", "rate_per_s", "r2"])
              .round(3).to_string(index=False))
    print(f"{len(manifest)} fields -> {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
