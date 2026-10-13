"""Turn one reference day into a "new culture" with hidden doses, to try decode_field.py on known answers.

    python3 scripts/make_worked_example.py --cells outputs/full_20261012_cellpose/cells_all_per_s.csv \\
        --session 20261007_h2o2_tr --out-dir outputs/example_user

Writes into --out-dir (which must not exist yet):
  cells_per_s.csv  that session's cells; every stressed condition renamed sampleNN (sample ids too),
                   the 0 µM conditions kept, so they act as the same-day control; session renamed user_day
  key.csv          sampleNN -> original condition and dose, to check the decoded answers

Then decode it with that session left out of the reference:
    python3 scripts/decode_field.py --reference-cells <same cells> --cells outputs/example_user/cells_per_s.csv \\
        --stressor h2o2 --exclude-reference-session 20261007_h2o2_tr --out-dir outputs/example_user/decoded
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def hide_doses(cells: pd.DataFrame, session: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    day = cells[cells["session_id"].astype(str) == session].copy()
    if day.empty:
        raise ValueError(f"no cells for session {session!r}")
    stressor = session.split("_")[1]
    control = re.compile(rf"^{re.escape(stressor)}_0(?:\.0+)?uM(?:_|$)")
    names, rows = {}, []
    for condition in sorted(day["condition_id"].astype(str).unique()):
        if not control.match(condition):
            names[condition] = f"sample{len(names) + 1:02d}"
            dose = re.match(rf"^{re.escape(stressor)}_(\d+(?:\.\d+)?)uM", condition)
            rows.append({"sample": names[condition], "condition_id": condition,
                         "dose_uM": float(dose.group(1)) if dose else None})
    stressed = day["condition_id"].isin(names)
    day.loc[stressed, "sample_id"] = (day.loc[stressed, "condition_id"].map(names) + "_"
                                      + day.loc[stressed].groupby("condition_id")["sample_id"]
                                      .transform(lambda s: s.astype("category").cat.codes + 1).astype(str))
    day.loc[stressed, "condition_id"] = day.loc[stressed, "condition_id"].map(names)
    day["session_id"] = "user_day"
    return day, pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--session", required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    if args.out_dir.exists():
        raise SystemExit(f"output directory already exists: {args.out_dir}")
    try:
        cells, key = hide_doses(pd.read_csv(args.cells, low_memory=False), args.session)
    except ValueError as exc:
        raise SystemExit(f"error: {exc}")
    args.out_dir.mkdir(parents=True)
    cells.to_csv(args.out_dir / "cells_per_s.csv", index=False)
    key.to_csv(args.out_dir / "key.csv", index=False)
    print(f"{len(cells)} cells, {len(key)} hidden samples; wrote {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
