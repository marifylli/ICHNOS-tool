"""Build the exact command lines ichnos_vail runs, and check the inputs before anything starts.

Every builder returns argument lists (no shell), so a path with spaces or
quotes is passed as one argument and nothing typed into ichnos_vail is ever
interpreted as a command. ``as_shell`` renders the same lists for display,
so a user can copy them and run them without ichnos_vail.
"""
from __future__ import annotations

import json
import shlex
import sys
from pathlib import Path

from .jobs import ROOT

SCRIPTS = ROOT / "scripts"
SEGMENTATION = ("cellpose", "brightfield_rim", "transmitted", "sparse", "otsu")
EXTRACTION = ("", "R", "G", "B", "sum", "luminance")


def _script(name: str) -> list[str]:
    return [sys.executable, "-u", str(SCRIPTS / name)]


def as_shell(command: list[str]) -> str:
    """Readable form: 'python3 scripts/x.py …' with paths inside the repository shortened."""
    shown = []
    for i, arg in enumerate(command):
        if i == 0 and arg == sys.executable:
            shown.append("python3")
            continue
        if arg == "-u" and i == 1:
            continue
        try:
            arg = str(Path(arg).relative_to(ROOT)) if Path(arg).is_absolute() else arg
        except ValueError:
            pass
        shown.append(arg)
    return shlex.join(shown)


def _new_folder(path: Path, what: str) -> list[str]:
    return [f"{what} already exists: {path} — choose a new folder"] if path.exists() else []


def _existing_file(path: Path, what: str) -> list[str]:
    return [] if path.is_file() else [f"{what} not found: {path}"]


def align(manifest: Path, out_dir: Path, search_up: int = 2) -> tuple[list[list[str]], list[str]]:
    """Find and align bright-field frames; writes out_dir/align.csv and out_dir/images_aligned.csv."""
    errors = _existing_file(manifest, "manifest") + _new_folder(out_dir, "output folder")
    command = _script("align_brightfield.py") + ["--manifest", str(manifest), "--out", str(out_dir / "align.csv"),
                                                  "--write-manifest", str(out_dir / "images_aligned.csv"),
                                                  "--search-up", str(int(search_up))]
    return [command], errors


def images_to_cells(manifest: Path, out_dir: Path, *, batch_size: int = 20, segmentation: str = "cellpose",
                    source: str = "brightfield", gpu: bool = True, bleed_default: float | None = 0.05,
                    controls: Path | None = None, saturation: float = 255, green_extraction: str = "",
                    red_extraction: str = "", max_foreground: float | None = None
                    ) -> tuple[list[list[str]], list[str]]:
    """Segmentation and extraction in batches (resumable), then per-second normalisation."""
    errors = _existing_file(manifest, "manifest")
    if (out_dir / "cells.csv").exists():
        errors.append(f"{out_dir / 'cells.csv'} already exists: this run is complete; choose a new folder")
    if segmentation not in SEGMENTATION:
        errors.append(f"unknown segmentation method {segmentation!r}")
    if controls is None and bleed_default is None:
        errors.append("give either a bleed-through value or a controls file")
    if controls is not None:
        errors += _existing_file(controls, "controls file")
    if not batch_size or batch_size < 1:
        errors.append("batch size must be at least 1")
    pipeline = ["--segmentation-method", segmentation, "--segmentation-source", source,
                "--saturation-value", f"{saturation:g}"]
    if segmentation == "cellpose" and gpu:
        pipeline.append("--cellpose-gpu")
    if controls is not None:
        pipeline += ["--controls", str(controls)]
    elif bleed_default is not None:
        pipeline += ["--bleed-default", f"{bleed_default:g}"]
    if green_extraction:
        pipeline += ["--green-extraction", green_extraction]
    if red_extraction:
        pipeline += ["--red-extraction", red_extraction]
    if max_foreground is not None:
        pipeline += ["--max-foreground-fraction", f"{max_foreground:g}"]
    batched = _script("run_pipeline_batched.py") + ["--manifest", str(manifest), "--out-dir", str(out_dir),
                                                     "--batch-size", str(int(batch_size)), "--", *pipeline]
    normalise = _script("normalize_exposure.py") + ["--cells", str(out_dir / "cells.csv"),
                                                     "--out", str(out_dir / "cells_per_s.csv")]
    return [batched, normalise], errors


def write_analysis_config(path: Path, cells: Path, *, steps: list[str], n_boot: int = 1000, seed: int = 0,
                          positive_quantile: float = 0.9, clearance_per_h: list[float | None] | None = None,
                          decoders: list[str] | None = None, min_cells: int = 20) -> dict:
    """The JSON run_analysis.py reads; ``cells`` is stored as an absolute path."""
    config = {"cells": str(Path(cells).resolve()), "steps": steps, "n_boot": int(n_boot), "seed": int(seed),
              "positive_quantile": float(positive_quantile),
              "clearance_per_h": clearance_per_h if clearance_per_h is not None else [None, 1.0, 3.0],
              "decoder_clearance_per_h": None, "profile": "default",
              "decoders": decoders or ["absolute", "paired", "model"], "min_cells": int(min_cells)}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(config, indent=2))
    return config


def analysis(config: Path, cells: Path, out_dir: Path) -> tuple[list[list[str]], list[str]]:
    errors = _existing_file(cells, "cell table") + _new_folder(out_dir, "output folder")
    return [_script("run_analysis.py") + ["--config", str(config), "--out-dir", str(out_dir)]], errors


def decode(reference: Path, cells: Path, stressor: str, out_dir: Path, *, decoder: str = "paired",
           exclude: list[str] | None = None, min_cells: int = 20,
           time_windows: str = "0,2,4,6,8") -> tuple[list[list[str]], list[str]]:
    errors = (_existing_file(reference, "reference cell table") + _existing_file(cells, "new cell table")
              + _new_folder(out_dir, "output folder"))
    if not stressor:
        errors.append("choose a stressor")
    if decoder not in ("paired", "absolute"):
        errors.append("decoder must be paired or absolute")
    try:
        edges = [float(x) for x in str(time_windows).split(",") if x.strip()]
        if len(edges) < 3 or any(b <= a for a, b in zip(edges, edges[1:])) or edges[0] < 0:
            raise ValueError
    except ValueError:
        errors.append("time windows: at least 3 increasing hours ≥ 0, separated by commas, e.g. 0,2,4,6,8")
    command = _script("decode_field.py") + ["--reference-cells", str(reference), "--cells", str(cells),
                                             "--stressor", stressor, "--decoder", decoder,
                                             "--min-cells", str(int(min_cells)), "--time-windows", str(time_windows),
                                             "--out-dir", str(out_dir)]
    for session in exclude or []:
        command += ["--exclude-reference-session", session]
    return [command], errors
