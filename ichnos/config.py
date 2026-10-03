"""Model configuration: which SBML files make up each variant, and the rules
the merge has to respect.

Migrated from Ichnos_PULSE python/ichnos_config.py @ e66de65c. The merge
rules and the calibration-only exclusion are preserved exactly; what changed
is how the SBML files are located. They are packaged resources inside
ichnos/models/ now, reached through importlib.resources, rather than paths
pointing at a sibling clone's integration/ directory. A wheel installed into
a clean environment therefore carries its own models.

This module previously held the image pipeline's QC thresholds. Those moved
to ichnos_image/instrument.py, so that model configuration and instrument
configuration cannot collide.
"""
from __future__ import annotations

from importlib import resources
from pathlib import Path


MODELS_PACKAGE = "ichnos.models"


def model_path(filename: str) -> Path:
    """Filesystem path to one of the packaged SBML resources.

    Goes through importlib.resources so this works from an installed wheel,
    not only from a source checkout.
    """
    with resources.as_file(resources.files(MODELS_PACKAGE) / filename) as path:
        if not path.exists():
            raise FileNotFoundError(f"packaged model {filename!r} is missing from {MODELS_PACKAGE}")
        return Path(path)


TIP_TETR_MODEL_FILE = "TIP_TetR_binding.sbml"
REPORTER_MODEL_FILE = "reporter_module_v2.sbml"

# Two variants, each a merge of three SBML files: TIP-TetR, the reporter, and
# one sensing module. These four files are a library of two variants, not
# four modules that all combine at once.
#
# 2026-09-10 (carried over from the source): the "ox" sensing file changed
# from oxidative_module_v3.sbml (static Hill, ox_step) to ox_adaptive.sbml
# (adaptive sensor with A_ox/X_ox buffer-node states). The ER variant is
# UNCHANGED -- ERModule.sbml is still the static Hill version.
VARIANTS: dict[str, dict] = {
    "er": {
        "sensing_file": "ERModule.sbml",
        "tip_name": "TIP_er",
        # 2026-09-11 (carried over): ERModule gained an adaptive sensor
        # (A_er/X_er) AND a simulated reporter pair (R_imm/R_mat) used only
        # to fit the module against the Pincus 2010 time course. R_imm/R_mat
        # are NOT part of the circuit -- the real readout is the tandem-timer
        # reporter module -- so they are excluded from the merge (confirmed
        # with the team). Listed BY NAME rather than id, consistent with
        # every other matching decision in this codebase.
        "calibration_only_species": ("R_imm", "R_mat"),
    },
    "ox": {
        "sensing_file": "ox_adaptive.sbml",
        "tip_name": "TIP_ox",
    },
}

# There is no working copper/CuSO4 model. The software must reject that
# choice rather than silently fall back to another variant.
UNSUPPORTED_VARIANTS = {"cu": "no copper/CuSO4 sensing model or calibration exists"}

SHARED_PARAM_NAMES = {"mu", "P"}
SHARED_PARAM_TOLERANCE = 1e-9

CROSS_VARIANT_SHARED_NAMES = {
    "k_deg_TIP",  # post-production degradation of TIP; same TIP coding sequence
                  # in ox/er/copper (confirmed wet lab 2026-08-17) -> must match.
}
CROSS_VARIANT_ABS_TOL = 1e-9


def sensing_path(variant: str) -> Path:
    """Packaged SBML for a variant's sensing module."""
    if variant in UNSUPPORTED_VARIANTS:
        raise ValueError(f"variant {variant!r} is not supported: {UNSUPPORTED_VARIANTS[variant]}")
    if variant not in VARIANTS:
        raise KeyError(f"unknown variant {variant!r}; available: {sorted(VARIANTS)}")
    return model_path(VARIANTS[variant]["sensing_file"])
