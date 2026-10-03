"""The SBML files must travel with the package.

The source repository found them at a path relative to a sibling checkout's
integration/ directory. Here they are package data, reached through
importlib.resources, so an installed wheel carries its own models. These
tests fail if that wiring breaks, which is otherwise only discovered by a
user whose install cannot find a model.
"""
import hashlib

import pytest

from ichnos import config


# The four executable sources, copied unedited from Ichnos_PULSE
# integration/ @ e66de65c. Hashes pin that they were not edited in transit.
SOURCE_SBML_SHA256 = {
    "TIP_TetR_binding.sbml": "217e162a61f5",
    "reporter_module_v2.sbml": "738009d231b7",
    "ox_adaptive.sbml": "843a41b1f465",
    "ERModule.sbml": "b49e9c0e2847",
}


def test_all_four_models_resolve():
    for filename in SOURCE_SBML_SHA256:
        path = config.model_path(filename)
        assert path.is_file()
        assert path.stat().st_size > 0

@pytest.mark.parametrize("variant", ["ox", "er"])
def test_packaged_default_profile_loads(variant):
    from ichnos.params import load_profile

    assert config.model_path("parameters.yaml").is_file()

    profile = load_profile(variant, "default")
    assert profile.variant == variant
    assert profile.name == "default"
    assert profile.parameters

def test_missing_resource_raises_clearly():
    with pytest.raises(FileNotFoundError):
        config.model_path("no_such_module.sbml")


@pytest.mark.parametrize("variant, expected_file", [("ox", "ox_adaptive.sbml"), ("er", "ERModule.sbml")])
def test_sensing_path_per_variant(variant, expected_file):
    assert config.sensing_path(variant).name == expected_file


def test_models_are_parseable_sbml():
    import libsbml

    for filename in SOURCE_SBML_SHA256:
        doc = libsbml.readSBMLFromFile(str(config.model_path(filename)))
        assert doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR) == 0, filename
        assert doc.getModel() is not None, filename


def test_source_files_unedited_in_transit():
    for filename, expected in SOURCE_SBML_SHA256.items():
        digest = hashlib.sha256(config.model_path(filename).read_bytes()).hexdigest()[:12]
        assert digest == expected, f"{filename} differs from the Ichnos_PULSE source"
