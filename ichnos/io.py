"""Saving a merged model as an artifact: a latest copy, a timestamped archive
copy that is never overwritten, and a JSON manifest recording exactly which
sources and parameter values went into it.

Migrated from Ichnos_PULSE python/ichnos_io.py @ e66de65c. Function bodies
are unchanged except for one thing: the output directory is now an argument
rather than a module-level constant pointing inside the source checkout.
Artifacts belong in a work directory the caller chooses, not inside an
installed package.

The display helpers that were in this module are in display_names.py.
"""

import os
import json
from datetime import datetime

import libsbml

from .config import SHARED_PARAM_NAMES
from .merge_checks import find_param_value_anywhere


def _ensure_dir(path):
    """Best-effort directory creation. On some machines os.path.isdir() has
    proven unreliable for a folder inside OneDrive (Files-On-Demand can
    present it as a cloud placeholder / reparse point that Python's stat
    calls don't always resolve consistently) — so we no longer trust isdir()
    to decide anything. Just attempt makedirs and swallow FileExistsError
    unconditionally: if something is there, proceed; if it turns out not to
    be a real usable directory, the subsequent open(..., 'w') call will fail
    with its own clear error instead."""
    try:
        os.makedirs(path, exist_ok=True)
    except FileExistsError:
        pass


def _collect_shared_param_snapshot(model):
    """Records the value+scope of every watched parameter, for the manifest.
    Same watch list as print_shared_parameter_summary — this is the written-
    to-disk counterpart of what that prints to the console."""
    watch_names = sorted(SHARED_PARAM_NAMES | {
        "b", "u_w", "delta_w", "a_TetR", "k_deg_TetR", "K_R", "n",
        "k_deg_TIP", "P_min",
    })
    snapshot = {}
    for name in watch_names:
        value, scope = find_param_value_anywhere(model, name)
        snapshot[name] = {"value": value, "scope": scope}
    return snapshot


def save_merged_sbml(sbml_str, variant, m_tetr, source_paths, out_dir):
    """Writes the merged SBML to two places:
      - exportsbml/merged_<variant>.sbml — 'latest' pointer, overwritten every run
      - exportsbml/archive/merged_<variant>_<timestamp>.sbml — never overwritten,
        plus a sidecar .json manifest recording exactly which source .sbml files
        (and their mtimes) and which key parameter values went into THIS merge.

    The archive copy exists so that downstream work (sensitivity analysis,
    the no-feedback comparison circuit, etc.) can always point back at the
    exact merged model a given result came from, instead of relying on
    whatever 'merged_<variant>.sbml' happens to contain right now.

    Returns (latest_path, archive_path, manifest_path).
    """
    out_dir = str(out_dir)
    _ensure_dir(out_dir)
    archive_dir = os.path.join(out_dir, "archive")
    _ensure_dir(archive_dir)

    latest_path = os.path.join(out_dir, f"merged_{variant}.sbml")
    with open(latest_path, "w", encoding="utf-8") as f:
        f.write(sbml_str)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_path = os.path.join(archive_dir, f"merged_{variant}_{timestamp}.sbml")
    with open(archive_path, "w", encoding="utf-8") as f:
        f.write(sbml_str)

    manifest = {
        "variant": variant,
        "built_at": datetime.now().isoformat(timespec="seconds"),
        "source_files": {
            label: {
                "path": path,
                "modified": datetime.fromtimestamp(os.path.getmtime(path)).isoformat(timespec="seconds")
                if os.path.isfile(path) else None,
            }
            for label, path in source_paths.items()
        },
        "key_parameters_in_merged_model": _collect_shared_param_snapshot(m_tetr),
    }
    manifest_path = archive_path.replace(".sbml", ".manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)

    print(f"  Merged SBML (latest):  {latest_path}")
    print(f"  Merged SBML (archive): {archive_path}")
    print(f"  Manifest:              {manifest_path}")
    return latest_path, archive_path, manifest_path


