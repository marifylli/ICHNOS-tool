"""Run scripts/run_pipeline.py over a large manifest in batches, then merge.

    python3 scripts/run_pipeline_batched.py --manifest data/acquisition/manifest/images.csv \\
        --out-dir outputs/full_20261010 --batch-size 20 -- \\
        --bleed-default 0.05 --segmentation-method sparse --saturation-value 255

Why: run_pipeline.py loads every image of the manifest before processing
the first one (about 85 MB per field: green, red and bright-field as
float64 plus three masks), so 294 fields need ~25 GB and the operating
system kills the process. Fields are processed independently of each other
(only the bleed coefficient is per session), so running the manifest in
batches gives the same per-cell records with bounded memory.

Layout of --out-dir:
  parts/part_NNN.csv              manifest rows of each batch
  parts/part_NNN.cells.csv        run_pipeline.py output, with its .manifest.json
  cells.csv                       all batches concatenated, in manifest order
  batches.json                    per-batch status, refused fields and timing

Batches whose output already exists are skipped, so an interrupted run is
resumed by repeating the same command. Everything after ``--`` is passed to
run_pipeline.py unchanged.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def run_one(manifest: Path, out_csv: Path, extra: list[str]) -> int:
    """Run run_pipeline.py on one batch; return its exit status."""
    command = [sys.executable, str(ROOT / "scripts" / "run_pipeline.py"),
               "--manifest", str(manifest), "--out", str(out_csv), *extra]
    return subprocess.run(command, check=False).returncode


def _concatenate(parts: list[Path], final: Path) -> None:
    """Join the batch CSVs as text, so values are byte-identical to a single run."""
    header = None
    with final.open("x", encoding="utf-8", newline="") as out:
        for part in parts:
            with part.open(encoding="utf-8", newline="") as handle:
                first = handle.readline()
                if header is None:
                    header = first
                    out.write(first)
                elif first != header:
                    raise SystemExit(f"{part} has different columns from the first batch")
                out.writelines(handle)


def run_batched(manifest: Path, out_dir: Path, batch_size: int, extra: list[str], runner=run_one) -> dict:
    rows = pd.read_csv(manifest, dtype=str, keep_default_na=False)
    if rows.empty:
        raise SystemExit("manifest has no rows")
    parts = out_dir / "parts"
    parts.mkdir(parents=True, exist_ok=True)
    final = out_dir / "cells.csv"
    if final.exists():
        raise SystemExit(f"output path already exists: {final}")

    batches = []
    for number, start in enumerate(range(0, len(rows), batch_size)):
        part_manifest = parts / f"part_{number:03d}.csv"
        part_cells = parts / f"part_{number:03d}.cells.csv"
        chunk = rows.iloc[start:start + batch_size]
        record = {"part": number, "rows": [start, start + len(chunk)], "cells_csv": str(part_cells)}
        sidecar = Path(str(part_cells) + ".manifest.json")
        if part_cells.exists() and sidecar.exists():
            # The manifest is written last and marks a complete batch.
            record["status"] = "skipped_existing"
        elif part_cells.exists():
            raise SystemExit(f"{part_cells} has no manifest, so that batch did not finish; "
                             f"delete it and rerun to redo the batch")
        else:
            chunk.to_csv(part_manifest, index=False)
            began = time.time()
            print(f"batch {number}: manifest rows {start}-{start + len(chunk) - 1}", flush=True)
            code = runner(part_manifest, part_cells, extra)
            record["seconds"] = round(time.time() - began, 1)
            record["status"] = "ok" if code == 0 and part_cells.exists() and sidecar.exists() else f"failed_exit_{code}"
        if sidecar.exists():
            record["refused"] = json.loads(sidecar.read_text()).get("refused_image_sets", [])
        batches.append(record)
        if record["status"].startswith("failed"):
            break

    summary = {"manifest": str(manifest), "batch_size": batch_size, "pipeline_args": extra, "batches": batches}
    (out_dir / "batches.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    failed = [b for b in batches if b["status"].startswith("failed")]
    if failed:
        print(f"batch {failed[0]['part']} failed; fix it and rerun the same command to resume", file=sys.stderr)
        return summary
    _concatenate([Path(b["cells_csv"]) for b in batches], final)
    refused = sum(len(b.get("refused", [])) for b in batches)
    print(f"{len(batches)} batches, {len(rows)} fields, {refused} refused -> {final}")
    return summary


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    extra = argv[argv.index("--") + 1:] if "--" in argv else []
    own = argv[:argv.index("--")] if "--" in argv else argv
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--batch-size", type=int, default=20)
    args = ap.parse_args(own)
    if args.batch_size < 1:
        raise SystemExit("--batch-size must be at least 1")
    summary = run_batched(args.manifest, args.out_dir, args.batch_size, extra)
    return 1 if any(b["status"].startswith("failed") for b in summary["batches"]) else 0


if __name__ == "__main__":
    raise SystemExit(main())
