"""Run the analysis of an exposure-normalised cell table from one configuration file.

    python3 scripts/run_analysis.py --config examples/analysis_202610.json \\
        --out-dir outputs/analysis_20261013

Steps, each into its own sub-folder of --out-dir (which must not exist yet):
  timecourse      scripts/fluorescence_timecourse.py  baseline, time courses, peak times
  induction       scripts/detect_induction.py         each dose vs the same day's 0 µM, detection limits
  model_compare   scripts/compare_model.py            model-predicted vs measured fold over 0 µM
  decoder_lodo    scripts/evaluate_decoder_lodo.py    single-field dose decoding, leave-one-day-out

Configuration (JSON; paths relative to the configuration file; unknown keys are refused):
  cells                    exposure-normalised cell table (scripts/normalize_exposure.py output)
  steps                    which of the steps above to run, in this order (default: all)
  n_boot, seed             resampling for timecourse, induction and model_compare
  positive_quantile        reporter-positive threshold (quantile of non-transfected cells)
  clearance_per_h          model_compare scenarios, e.g. [null, 1, 3]; null = constant stress
  decoder_clearance_per_h  clearance for the model decoder (null = constant stress)
  profile                  model parameter profile
  decoders                 decoder_lodo decoders
  min_cells                decoder_lodo: smallest field scored

Writes run-manifest.json: the configuration, the input's sha256, the git
commit (and whether the working tree had changes), and the sha256 of every
output file. The image steps that produce the cell table (manifest,
bright-field alignment, segmentation, exposure normalisation) take hours and
are run separately; docs/architecture.md lists them.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import runpy
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEPS = {
    "timecourse": "fluorescence_timecourse.py",
    "induction": "detect_induction.py",
    "model_compare": "compare_model.py",
    "decoder_lodo": "evaluate_decoder_lodo.py",
}
DEFAULTS = {
    "steps": list(STEPS), "n_boot": 1000, "seed": 0, "positive_quantile": 0.9,
    "clearance_per_h": [None, 1.0, 3.0], "decoder_clearance_per_h": None, "profile": "default",
    "decoders": ["absolute", "paired", "model"], "min_cells": 20,
}


def load_config(path: Path) -> dict:
    raw = json.loads(path.read_text())
    unknown = set(raw) - set(DEFAULTS) - {"cells"}
    if unknown:
        raise ValueError(f"unknown configuration keys: {sorted(unknown)}")
    if "cells" not in raw:
        raise ValueError("configuration needs 'cells'")
    config = {**DEFAULTS, **raw}
    bad = [s for s in config["steps"] if s not in STEPS]
    if bad:
        raise ValueError(f"unknown steps {bad}; choose from {list(STEPS)}")
    config["cells"] = str((path.parent / raw["cells"]).resolve())
    return config


def step_argv(step: str, config: dict, out: Path) -> list[str]:
    common = ["--cells", config["cells"], "--out-dir", str(out)]
    boot = ["--n-boot", str(config["n_boot"]), "--seed", str(config["seed"])]
    quantile = ["--positive-quantile", str(config["positive_quantile"])]
    if step == "timecourse":
        return common + boot
    if step == "induction":
        return common + boot + quantile
    if step == "model_compare":
        scenarios = ["none" if c is None else str(c) for c in config["clearance_per_h"]]
        return common + boot + quantile + ["--profile", config["profile"], "--clearance", *scenarios]
    if step == "decoder_lodo":
        clearance = config["decoder_clearance_per_h"]
        return (common + quantile + ["--decoders", *config["decoders"], "--profile", config["profile"],
                                     "--min-cells", str(config["min_cells"]),
                                     "--clearance", "none" if clearance is None else str(clearance)])
    raise ValueError(step)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def git_state() -> dict:
    def git(*args):
        result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)
        return result.stdout.strip() if result.returncode == 0 else None
    commit = git("rev-parse", "HEAD")
    status = git("status", "--porcelain", "--untracked-files=no")
    return {"commit": commit, "tracked_changes": None if status is None else bool(status)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)

    if args.out_dir.exists():
        raise SystemExit(f"output directory already exists: {args.out_dir}")
    try:
        config = load_config(args.config)
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"configuration error: {exc}")
    cells = Path(config["cells"])
    if not cells.exists():
        raise SystemExit(f"cell table not found: {cells}")

    args.out_dir.mkdir(parents=True)
    manifest = {"schema_version": 1, "config_file": str(args.config.resolve()), "config": config,
                "input_sha256": sha256(cells), "code": git_state(), "steps": []}
    status = 0
    for step in config["steps"]:
        out = args.out_dir / step
        print(f"\n=== {step} ===", flush=True)
        started = time.time()
        entry = {"step": step, "script": f"scripts/{STEPS[step]}", "argv": step_argv(step, config, out)}
        try:
            code = runpy.run_path(str(ROOT / "scripts" / STEPS[step]), run_name="ichnos_step")["main"](entry["argv"])
            entry["status"] = "ok" if code in (0, None) else f"exit {code}"
        except SystemExit as exc:
            entry["status"] = f"stopped: {exc}"
        entry["seconds"] = round(time.time() - started, 1)
        manifest["steps"].append(entry)
        if entry["status"] != "ok":
            status = 1
            print(f"{step} failed ({entry['status']}); later steps not run", file=sys.stderr)
            break
    manifest["outputs"] = {str(p.relative_to(args.out_dir)): sha256(p)
                           for p in sorted(args.out_dir.rglob("*")) if p.is_file()}
    manifest["status"] = "ok" if status == 0 else "failed"
    (args.out_dir / "run-manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\n{manifest['status']}: {len(manifest['steps'])} step(s); wrote {args.out_dir}")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
