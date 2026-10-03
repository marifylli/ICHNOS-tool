"""ichnos_image.metadata: reads back what tifffile can find on a TIFF.
Uses a self-created synthetic TIFF (no real Olympus file available) --
verifies the read mechanics work, not that it matches real Olympus output.
"""
import numpy as np
import pytest

tifffile = pytest.importorskip("tifffile")

from ichnos_image.metadata import read_tiff_metadata  # noqa: E402


def test_reads_image_description_and_shape(tmp_path):
    arr = (np.random.default_rng(0).random((32, 40)) * 1000).astype("uint16")
    path = tmp_path / "synthetic.tiff"
    tifffile.imwrite(path, arr, description="Exposure=800ms;Objective=60x;PixelSize=0.1076um")

    meta = read_tiff_metadata(path)

    assert meta["shape"] == (32, 40)
    assert meta["dtype"] == "uint16"
    assert "Exposure=800ms" in meta["image_description"]
    assert "tags" in meta and "ImageWidth" in meta["tags"]


def test_missing_tifffile_raises_clear_error(monkeypatch, tmp_path):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "tifffile":
            raise ImportError("simulated missing tifffile")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="tifffile"):
        read_tiff_metadata(tmp_path / "does_not_matter.tiff")
