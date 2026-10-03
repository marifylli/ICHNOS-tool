"""The migrated builder must produce exactly the model the source repository
produces, from packaged resources, with no filesystem side effects.

The strongest check available at this stage is reproduction: Step 3 moved
code, it did not change the model, so any difference in the output is a
migration bug. test_golden_hashes below pins what the builder currently
emits, and the equivalence against Ichnos_PULSE @ e66de65c was verified
directly at migration time (byte-identical SBML, bit-identical trajectories
for both variants).
"""
import hashlib

import pytest

from ichnos import build, config


VARIANTS = ("ox", "er")

# SHA256 of the merged SBML string, verified byte-identical against
# Ichnos_PULSE python/ichnos_core.build_variant_sbml_string() @ e66de65c on
# 2026-10-03. A change here means the merged model changed; that may be
# intended, but it is never a side effect of moving files.
GOLDEN_SBML_SHA256 = {
    "ox": "49684d48907e",
    "er": "af7132b2e89b",
}


def _sha12(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()[:12]


@pytest.mark.parametrize("variant", VARIANTS)
def test_golden_hashes(variant):
    sbml = build.build_variant_sbml_string(variant)
    assert _sha12(sbml) == GOLDEN_SBML_SHA256[variant]


@pytest.mark.parametrize("variant", VARIANTS)
def test_build_is_deterministic(variant):
    assert build.build_variant_sbml_string(variant) == build.build_variant_sbml_string(variant)


def test_the_two_variants_differ():
    """Four SBML files are a library of two variants, not four modules that
    all merge at once. ox and er share TIP-TetR and the reporter and differ
    only in the sensing module.
    """
    assert build.build_variant_sbml_string("ox") != build.build_variant_sbml_string("er")


@pytest.mark.parametrize("variant", VARIANTS)
def test_build_writes_nothing_by_default(variant, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    build.build_variant_sbml_string(variant)
    assert list(tmp_path.iterdir()) == []


def test_saving_requires_an_explicit_directory():
    with pytest.raises(ValueError, match="out_dir"):
        build.build_variant_sbml_string("ox", save_sbml=True)


def test_save_writes_latest_archive_and_manifest(tmp_path):
    build.build_variant_sbml_string("ox", save_sbml=True, out_dir=tmp_path)
    assert (tmp_path / "merged_ox.sbml").is_file()
    archive = tmp_path / "archive"
    assert archive.is_dir()
    assert any(p.suffix == ".sbml" for p in archive.iterdir())
    assert any(p.suffix == ".json" for p in archive.iterdir())


def test_copper_variant_is_rejected_not_substituted():
    """There is no working copper/CuSO4 model. The tool must refuse rather
    than quietly build something else.
    """
    with pytest.raises(ValueError, match="not supported"):
        build.build_variant_sbml_string("cu")
    with pytest.raises(ValueError, match="not supported"):
        config.sensing_path("cu")


def test_unknown_variant_raises():
    with pytest.raises(KeyError):
        build.build_variant_sbml_string("hydrogen-peroxide")
