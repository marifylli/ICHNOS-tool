import numpy as np
import pytest

from ichnos_image import segment


@pytest.mark.parametrize("min_size", [30, 60])
def test_sparse_filters_final_watershed_labels(monkeypatch, min_size):
    rng = np.random.default_rng(42)
    image = rng.normal(20, 1, (128, 128))
    image[40:80, 40:80] += 100

    # A retained component can split into labels below the size threshold.
    split_labels = np.zeros(image.shape, dtype=np.int32)
    split_labels[10, 10:10 + min_size - 1] = 7
    split_labels[20, 10:10 + min_size] = 12
    split_labels[30, 10:10 + min_size + 1] = 19

    def controlled_watershed(binary):
        assert binary.any()
        return split_labels.copy()

    monkeypatch.setattr(segment, "_split_touching", controlled_watershed)

    result = segment.segment_cells(
        image, method="sparse", min_size=min_size
    )

    # Below threshold is rejected; at/above threshold retain IDs and pixels.
    assert not np.any(result == 7)
    for label in (12, 19):
        np.testing.assert_array_equal(
            result == label, split_labels == label
        )
    assert set(np.unique(result)) == {0, 12, 19}
