"""segment_cells(resize_factor=...) mechanics: output shape matches the
original image and labels stay spatially sane after the downsample/upsample
round-trip. Does NOT assert accuracy for method="otsu" at low resize_factor
-- scripts/evaluate_resize_accuracy_tradeoff.py found that breaks down hard
(F1 0.83 -> 0.0 between resize_factor 1.0 and 0.5) because otsu's internal
edge-detection kernels are fixed in pixels, unlike cellpose which degrades
gracefully (see segment.segment_cells's docstring for the measured numbers).
"""
import numpy as np

from ichnos_image import segment


def _bf_scene(size=200, n_cells=15, seed=0):
    rng = np.random.default_rng(seed)
    bf = np.full((size, size), 0.5)
    yy, xx = np.mgrid[0:size, 0:size]
    for cy, cx in rng.integers(20, size - 20, size=(n_cells, 2)):
        blob = ((yy - cy) ** 2 + (xx - cx) ** 2) < 8**2
        bf[blob] -= 0.3
    return bf


def test_resize_factor_preserves_output_shape():
    bf = _bf_scene()
    for resize_factor in (1.0, 0.75, 0.5, 0.25):
        labels = segment.segment_cells(bf, method="otsu", resize_factor=resize_factor)
        assert labels.shape == bf.shape
        assert labels.dtype == np.int32


def test_resize_factor_close_to_one_finds_similar_cell_count():
    bf = _bf_scene()
    full = segment.segment_cells(bf, method="otsu", resize_factor=1.0)
    moderate = segment.segment_cells(bf, method="otsu", resize_factor=0.75)
    assert abs(int(full.max()) - int(moderate.max())) <= 3
