"""segment_cells(method="cellpose") against a real cellpose install, when
available. Skipped in the default test environment (cellpose+torch are
heavy optional deps, not in requirements.txt) -- run with cellpose installed
to exercise this. This test exists because the cellpose API changed between
major versions (models.Cellpose(model_type=...) was removed in cellpose 4,
replaced by models.CellposeModel(...).eval(..., channels=None)) and the
previous version of _segment_cellpose here was silently broken against a
real cellpose>=4 install until this was checked by hand.
"""
import numpy as np
import pytest

cellpose = pytest.importorskip("cellpose")

from ichnos_image import segment  # noqa: E402


def _small_synthetic_bf(size=128, n_cells=6, seed=0):
    rng = np.random.default_rng(seed)
    bf = np.full((size, size), 0.5)
    yy, xx = np.mgrid[0:size, 0:size]
    for cy, cx in rng.integers(15, size - 15, size=(n_cells, 2)):
        blob = ((yy - cy) ** 2 + (xx - cx) ** 2) < 8**2
        bf[blob] -= 0.3
    return bf


def test_cellpose_segmentation_runs_and_finds_cells():
    bf = _small_synthetic_bf()
    labels = segment.segment_cells(bf, method="cellpose")
    assert labels.shape == bf.shape
    assert labels.max() >= 1


def test_cellpose_model_can_be_reused_across_calls():
    model = segment.load_cellpose_model()
    bf1 = _small_synthetic_bf(seed=1)
    bf2 = _small_synthetic_bf(seed=2)
    labels1 = segment.segment_cells(bf1, method="cellpose", cellpose_model=model)
    labels2 = segment.segment_cells(bf2, method="cellpose", cellpose_model=model)
    assert labels1.max() >= 1
    assert labels2.max() >= 1
