"""Measure, for each row of an image manifest, whether its bright-field frame
shows the same field as the fluorescence pair -- and by how much it is shifted.

Writes the brightfield_shift_dy / brightfield_shift_dx columns that
run_pipeline.py reads when --segmentation-source brightfield is used, and a
z column saying how much to believe each one.

Every bright-field frame is scored twice: against the fluorescence it is
paired with, and against the fluorescence of a different row. The second is
the null -- what this very mask scores on a field it certainly does not
belong to -- so the verdict is a comparison of two measured columns rather
than a threshold chosen in advance. That matters: the first version of this
check used a fixed threshold borrowed from a comparison of a different kind
and called 24 real matches misses.

Usage
  python scripts/align_brightfield.py --manifest images.csv --out aligned.csv
  python scripts/align_brightfield.py --manifest images.csv --out aligned.csv --write-back
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ichnos_image import manifest as manifest_io  # noqa: E402
from ichnos_image import segment  # noqa: E402
from ichnos_image.align import align_mask  # noqa: E402
from ichnos_image.image_io import load_image  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--max-shift", type=int, default=200,
                        help="largest displacement considered, in pixels")
    parser.add_argument("--min-size", type=int, default=200,
                        help="smallest bright-field object kept, in pixels")
    parser.add_argument("--write-back", action="store_true",
                        help="also write the shifts into the manifest itself")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    image_sets = manifest_io._build_image_sets(args.manifest, image_loader=load_image)
    usable = [(i, s) for i, s in enumerate(image_sets) if s.bright_field is not None]
    if not usable:
        raise SystemExit("no row of the manifest has a bright_field_path")

    rng = np.random.default_rng(args.seed)
    rows = []
    for position, (index, image_set) in enumerate(usable):
        labels = segment.segment_cells(
            image_set.bright_field, method="transmitted", min_size=args.min_size
        )
        mask = labels > 0
        matched = align_mask(mask, image_set.red, max_shift=args.max_shift)
        on_green = align_mask(mask, image_set.green, max_shift=args.max_shift)
        others = [p for p in range(len(usable)) if p != position]
        unrelated = (align_mask(mask, usable[int(rng.choice(others))][1].red,
                                max_shift=args.max_shift)
                     if others else None)
        rows.append(dict(
            row=index, session_id=image_set.session_id, sample_id=image_set.sample_id,
            acquisition_order=image_set.acquisition_order,
            cells_fraction=round(matched.foreground_fraction, 4),
            z_matched=round(matched.z, 1),
            z_unrelated=round(unrelated.z, 1) if unrelated else np.nan,
            z_green=round(on_green.z, 1),
            contrast_red=round(matched.contrast, 3),
            contrast_green=round(on_green.contrast, 3),
            brightfield_shift_dy=matched.dy, brightfield_shift_dx=matched.dx,
            shift_px=round(matched.shift_px, 1),
            at_search_limit=bool(max(abs(matched.dy), abs(matched.dx)) >= args.max_shift),
        ))
        print(f"row {index:3d} {image_set.sample_id or '':>16}: cells "
              f"{matched.foreground_fraction:5.1%}  z own {matched.z:6.1f}  "
              f"vs unrelated {rows[-1]['z_unrelated']:6.1f}  "
              f"contrast {matched.contrast:5.2f}  shift {matched.shift_px:6.1f}px"
              f"{'  [at search limit]' if rows[-1]['at_search_limit'] else ''}")

    result = pd.DataFrame(rows)
    null = result.z_unrelated.dropna()
    cut = float(np.percentile(null, 95)) if len(null) else np.nan
    matched_rows = result.z_matched > cut

    print("\n--- calibration, measured on this manifest ---")
    print(f"z on an unrelated field: median {null.median():.1f}, 95th percentile {cut:.1f}")
    print(f"z on its own field:      median {result.z_matched.median():.1f}, "
          f"max {result.z_matched.max():.1f}")
    print(f"\n{int(matched_rows.sum())} of {len(result)} bright-field frames score above "
          f"anything the unrelated comparisons reach.")
    if matched_rows.any():
        print(f"their median shift: {result.loc[matched_rows, 'shift_px'].median():.1f} px; "
              f"median contrast inside/outside: "
              f"{result.loc[matched_rows, 'contrast_red'].median():.2f}")
    if (~matched_rows).any():
        print("\nrows below the cut -- check at_search_limit (a shift pinned to the limit "
              "means the true one is larger, which usually means a different field) and "
              "z_green (if both channels say no, it is the field; if only red does, the "
              "red frame may simply be empty):")
        print(result.loc[~matched_rows, ["row", "sample_id", "z_matched", "z_unrelated",
                                         "z_green", "shift_px", "at_search_limit"]]
              .to_string(index=False))

    result.to_csv(args.out, index=False)
    print(f"\nwrote {args.out}")

    if args.write_back:
        manifest = pd.read_csv(args.manifest)
        for column in ("brightfield_shift_dy", "brightfield_shift_dx"):
            if column not in manifest:
                manifest[column] = np.nan
        for _, item in result[matched_rows].iterrows():
            manifest.loc[int(item.row), "brightfield_shift_dy"] = item.brightfield_shift_dy
            manifest.loc[int(item.row), "brightfield_shift_dx"] = item.brightfield_shift_dx
        # Rows that did not clear the null keep no shift AND no bright field:
        # segmenting them off a fluorescence channel instead would put two
        # mask sources in one run, which is the confusion mask_source exists
        # to prevent. They are refused by name when the run asks for
        # bright-field masks.
        if "bright_field_path" in manifest:
            manifest.loc[result.loc[~matched_rows, "row"].astype(int), "bright_field_path"] = ""
        manifest.to_csv(args.manifest, index=False)
        print(f"wrote the shifts back into {args.manifest}; "
              f"cleared bright_field_path on {int((~matched_rows).sum())} unmatched row(s)")


if __name__ == "__main__":
    main()
