"""Which dose was my culture exposed to? Decode new images against a reference data set.

    python3 scripts/decode_field.py \\
        --reference-cells outputs/full_20261012_cellpose/cells_all_per_s.csv \\
        --cells my_run/cells_per_s.csv --stressor h2o2 --out-dir my_run/decoded

--reference-cells  exposure-normalised cells of the reference experiment: transfected days
                   at known doses plus non-transfected cells (autofluorescence baseline)
--cells            the new images, through scripts/run_pipeline.py and scripts/normalize_exposure.py.
                   Name the same-day unstressed culture's condition <stressor>_0uM_<label>
                   (e.g. h2o2_0uM_ctrl); give unknown samples any other condition name.
--decoder          paired (default; needs that 0 µM culture) or absolute (one image, nothing else)
--time-windows     edges of the time-since-stress windows in hours (default 0,2,4,6,8)

Writes into --out-dir (which must not exist yet):
  decoded.csv     per field and per sample: posterior over the reference doses, most probable dose,
                  80 % credible set, verdict; and the same for the time since the stress
                  (time_map, time_credible80, pt_<window>, time_verdict)
  validation.csv  how the dose decoder scored on held-out reference days (leave-one-day-out)
  validation_time.csv  the same for the time decoder
  summary.json    settings, reference days, input fingerprints, verdict counts

Read the verdict before the dose: see ichnos_image/field_decoder.py and docs/user_guide.md.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ichnos_image import field_decoder as fd


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference-cells", type=Path, required=True)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--stressor", required=True, help="as in the reference session ids, e.g. h2o2 or dtt")
    ap.add_argument("--decoder", default="paired", choices=list(fd.OFFERED))
    ap.add_argument("--exclude-reference-session", action="append", default=[],
                    help="leave a reference session out of training (repeatable)")
    ap.add_argument("--time-windows", default="0,2,4,6,8",
                    help="edges of the time-since-stress windows, hours, comma-separated")
    ap.add_argument("--min-cells", type=int, default=20)
    ap.add_argument("--positive-quantile", type=float, default=0.9)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)

    if args.out_dir.exists():
        raise SystemExit(f"output directory already exists: {args.out_dir}")
    try:
        edges = tuple(float(x) for x in args.time_windows.split(",") if x.strip())
    except ValueError:
        raise SystemExit("error: --time-windows must be numbers separated by commas, e.g. 0,2,4,6,8")
    try:
        trained = fd.train(pd.read_csv(args.reference_cells, low_memory=False), args.stressor, args.decoder,
                           quantile=args.positive_quantile, min_cells=args.min_cells,
                           exclude_sessions=tuple(args.exclude_reference_session), time_edges=edges)
        decoded = fd.decode(pd.read_csv(args.cells, low_memory=False), trained)
    except ValueError as exc:
        raise SystemExit(f"error: {exc}")

    args.out_dir.mkdir(parents=True)
    decoded.to_csv(args.out_dir / "decoded.csv", index=False)
    trained["validation"].to_csv(args.out_dir / "validation.csv", index=False)
    time = trained["time"]
    if time["available"]:
        time["validation"].to_csv(args.out_dir / "validation_time.csv", index=False)
    counts = decoded[decoded["unit"] == "field"]["verdict"].value_counts().to_dict()
    time_counts = decoded[decoded["unit"] == "field"]["time_verdict"].value_counts().to_dict()
    record = time.get("record", {})
    (args.out_dir / "summary.json").write_text(json.dumps({
        "stressor": args.stressor, "decoder": args.decoder, "reference_days": trained["reference_days"],
        "reference_doses_uM": trained["doses"], "excluded_reference_sessions": args.exclude_reference_session,
        "validated_accuracy": trained["validated_accuracy"], "validated_info_gain_bits": trained["validated_gain_bits"],
        "chance": 1 / len(trained["doses"]), "field_verdicts": counts,
        "time_windows": time["labels"], "time_available": time["available"], "time_reason": time.get("reason", ""),
        "time_windows_in_reference": [time["labels"][w] for w in time.get("windows", [])],
        "time_validated_accuracy": record.get("accuracy"), "time_validated_info_gain_bits": record.get("info_gain_bits"),
        "time_chance": record.get("chance"), "time_scored_fields": record.get("n_scored"),
        "time_unscorable_fields": record.get("n_unscorable"), "time_field_verdicts": time_counts,
        "reference_cells": str(args.reference_cells), "reference_sha256": sha256(args.reference_cells),
        "cells": str(args.cells), "cells_sha256": sha256(args.cells),
        "min_cells": args.min_cells, "positive_quantile": args.positive_quantile}, indent=2))

    pd.set_option("display.width", 200)
    print(f"reference: {len(trained['reference_days'])} day(s) of {args.stressor}, doses "
          f"{', '.join(f'{d:g}' for d in trained['doses'])} µM; decoder '{args.decoder}'")
    acc, gain = trained["validated_accuracy"], trained["validated_gain_bits"]
    if pd.notna(acc):
        print(f"track record on held-out reference days: accuracy {acc:.2f} (chance {1 / len(trained['doses']):.2f}), "
              f"information gain {gain:+.2f} bits per field")
    else:
        print("track record: not available (a single reference day cannot be held out)")
    if not time["available"]:
        print(f"time since stress: not available — {time['reason']}")
    elif record.get("n_scored"):
        print(f"time since stress, windows {', '.join(time['labels'][w] for w in time['windows'])}: on held-out "
              f"reference days accuracy {record['accuracy']:.2f} (chance {record['chance']:.2f}), information gain "
              f"{record['info_gain_bits']:+.2f} bits per field ({record['n_scored']} scored, "
              f"{record['n_unscorable']} in windows no other day covers)")
    else:
        print("time since stress: track record not available (no held-out field in a window other days cover)")
    shown = decoded[["unit", "condition_id", "sample_id", "known_dose_uM", "known_time_h", "n_cells",
                     "map_dose_uM", "credible80_uM", "verdict", "time_map", "time_credible80", "time_verdict"]]
    print(shown.round(2).to_string(index=False))
    print(f"\nverdicts (fields): dose {counts}; time {time_counts}\nwrote {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
