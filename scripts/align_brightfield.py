"""Find each row's bright-field frame, say whether it shows the same field as
the fluorescence pair, and measure by how much the two are displaced.

Writes bright_field_path and brightfield_shift_dy / brightfield_shift_dx back
into the manifest, which is what run_pipeline.py reads for
--segmentation-source brightfield.

Candidates are discovered rather than declared. A manifest built before
anyone noticed the bright-field frames has no bright_field_path at all, and
a folder may hold several bright-field files with nothing in their names to
say which field each belongs to. So every file matching --brightfield-glob
beside the fluorescence pair is tried, and the one that fits best is kept --
if any of them fits at all.

Fit is decided against a null measured here, not against a threshold. Each
bright-field frame is also scored on the fluorescence of a different row:
what this very mask reaches on a field it certainly does not belong to. The
verdict compares the two measured columns. An earlier version of this check
used a fixed threshold taken from a comparison of a different kind and
called 24 real matches misses.

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

from ichnos_image import segment  # noqa: E402
from ichnos_image.align import align_mask  # noqa: E402
from ichnos_image.image_io import load_image  # noqa: E402
from ichnos_image.manifest import _load_bright_field  # noqa: E402


def resolve(manifest_path: Path, value: str) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else manifest_path.parent / path


def fluorescence_plane(path: Path) -> np.ndarray:
    image = load_image(path)
    return image.mean(axis=-1) if image.ndim == 3 else image


def candidates_for(row_paths: list[Path], pattern: str) -> list[Path]:
    """Bright-field files sitting beside the fluorescence pair.

    The pair's own files are excluded by path, not by name pattern, so a
    glob that happens to match them cannot make an image its own mask
    source.
    """
    folder = row_paths[0].parent
    taken = {p.resolve() for p in row_paths}
    return [p for p in sorted(folder.glob(pattern)) if p.resolve() not in taken]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--brightfield-glob", default="Image*.tif",
                        help="which files beside the fluorescence pair to try as "
                             "bright-field (default: Image*.tif); ignored for rows that "
                             "already name a bright_field_path")
    parser.add_argument("--max-shift", type=int, default=200,
                        help="largest displacement considered, in pixels")
    parser.add_argument("--min-size", type=int, default=200,
                        help="smallest bright-field object kept, in pixels")
    parser.add_argument("--write-back", action="store_true",
                        help="write bright_field_path and the shifts into the manifest")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest)
    for column in ("green_path", "red_path"):
        if column not in manifest:
            raise SystemExit(f"the manifest has no {column} column")

    work = []
    for index, row in manifest.iterrows():
        pair = [resolve(args.manifest, row.green_path), resolve(args.manifest, row.red_path)]
        declared = row.get("bright_field_path", "")
        if isinstance(declared, str) and declared.strip():
            found = [resolve(args.manifest, declared)]
        else:
            found = candidates_for(pair, args.brightfield_glob)
        if found:
            work.append((index, pair, found))
    if not work:
        raise SystemExit(
            f"no row has a bright_field_path, and nothing matching "
            f"{args.brightfield_glob!r} sits beside the fluorescence files. "
            f"Pass --brightfield-glob with the pattern your files use."
        )
    print(f"{len(work)} of {len(manifest)} rows have at least one bright-field candidate")

    # The null has to have enough draws for its 95th percentile to mean
    # anything. One unrelated comparison per row is plenty for a session of
    # a couple of dozen; a manifest of three rows would be deciding against
    # three numbers, so short manifests take several each.
    draws = max(1, min(5, -(-20 // max(len(work), 1))))
    print(f"scoring each against {draws} unrelated field(s) as the null")

    rng = np.random.default_rng(args.seed)
    rows = []
    for position, (index, pair, found) in enumerate(work):
        red = fluorescence_plane(pair[1])
        green = fluorescence_plane(pair[0])

        best = None
        for candidate in found:
            labels = segment.segment_cells(
                _load_bright_field(candidate), method="transmitted", min_size=args.min_size
            )
            placed = align_mask(labels > 0, red, max_shift=args.max_shift)
            if best is None or np.nan_to_num(placed.z, nan=-np.inf) > best[1].z:
                best = (candidate, placed, labels > 0)
        candidate, placed, mask = best

        others = [p for p in range(len(work)) if p != position]
        nulls = []
        for other in rng.choice(others, size=min(draws, len(others)), replace=False) if others else []:
            nulls.append(align_mask(mask, fluorescence_plane(work[int(other)][1][1]),
                                    max_shift=args.max_shift).z)
        on_green = align_mask(mask, green, max_shift=args.max_shift)

        rows.append(dict(
            row=index, brightfield=str(candidate.relative_to(args.manifest.parent)
                                       if candidate.is_relative_to(args.manifest.parent)
                                       else candidate),
            n_candidates=len(found),
            cells_fraction=round(placed.foreground_fraction, 4),
            z_matched=round(placed.z, 1),
            z_unrelated=round(float(np.nanmax(nulls)), 1) if nulls else np.nan,
            z_green=round(on_green.z, 1),
            contrast_red=round(placed.contrast, 3),
            contrast_green=round(on_green.contrast, 3),
            brightfield_shift_dy=placed.dy, brightfield_shift_dx=placed.dx,
            shift_px=round(placed.shift_px, 1),
            at_search_limit=bool(max(abs(placed.dy), abs(placed.dx)) >= args.max_shift),
        ))
        label = f"{pair[0].parent.name}/{pair[0].stem}"[-34:]
        print(f"{label:>34}: {len(found)} candidate(s)  cells {placed.foreground_fraction:5.1%}  "
              f"z own {placed.z:6.1f}  vs unrelated {rows[-1]['z_unrelated']:6.1f}  "
              f"contrast {placed.contrast:5.2f}  shift {placed.shift_px:6.1f}px"
              f"{'  [at search limit]' if rows[-1]['at_search_limit'] else ''}")

    result = pd.DataFrame(rows)
    # Each row's null is the worst (highest) score its mask reached on a
    # field it does not belong to, so the cut is drawn against the best
    # that chance managed, not the average.
    null = result.z_unrelated.dropna()
    cut = float(np.percentile(null, 95)) if len(null) else np.nan
    matched = result.z_matched > cut

    print("\n--- calibration, measured on this manifest ---")
    print(f"z on an unrelated field: median {null.median():.1f}, 95th percentile {cut:.1f}")
    print(f"z on its own field:      median {result.z_matched.median():.1f}, "
          f"max {result.z_matched.max():.1f}")
    print(f"\n{int(matched.sum())} of {len(result)} rows have a bright-field frame that "
          f"scores above anything the unrelated comparisons reach.")
    if matched.any():
        print(f"their median shift: {result.loc[matched, 'shift_px'].median():.1f} px; "
              f"median contrast inside/outside: "
              f"{result.loc[matched, 'contrast_red'].median():.2f}")
    if (~matched).any():
        print("\nrows below the cut -- at_search_limit means the best displacement is pinned "
              "to the limit, so the true one is larger, which usually means a different "
              "field; z_green says whether the green channel agrees:")
        print(result.loc[~matched, ["row", "brightfield", "z_matched", "z_unrelated",
                                    "z_green", "shift_px", "at_search_limit"]]
              .to_string(index=False))

    result.to_csv(args.out, index=False)
    print(f"\nwrote {args.out}")

    if args.write_back:
        for column in ("bright_field_path", "brightfield_shift_dy", "brightfield_shift_dx"):
            if column not in manifest:
                manifest[column] = "" if column == "bright_field_path" else np.nan
        manifest["bright_field_path"] = manifest["bright_field_path"].astype(object)
        # Only the rows that cleared the null get a bright field. The rest
        # get none, so a bright-field run refuses them by name rather than
        # quietly segmenting them off a fluorescence channel -- two mask
        # sources in one run is the confusion mask_source exists to prevent.
        for _, item in result.iterrows():
            keep = bool(item.z_matched > cut)
            manifest.loc[int(item.row), "bright_field_path"] = item.brightfield if keep else ""
            manifest.loc[int(item.row), "brightfield_shift_dy"] = (
                item.brightfield_shift_dy if keep else np.nan)
            manifest.loc[int(item.row), "brightfield_shift_dx"] = (
                item.brightfield_shift_dx if keep else np.nan)
        manifest.to_csv(args.manifest, index=False)
        print(f"wrote {int(matched.sum())} bright-field path(s) and their shifts into "
              f"{args.manifest}; left {int((~matched).sum())} row(s) without one")


if __name__ == "__main__":
    main()
