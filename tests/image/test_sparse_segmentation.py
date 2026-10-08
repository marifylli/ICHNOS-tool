"""segment_cells(method="sparse") must hold up where "otsu" does not: on a
frame whose cells are a few percent of the pixels, and as that frame gets
dim. And foreground_fraction must separate both failure modes from a real
field, since neither shows up in the object count."""
import numpy as np
import pytest

from ichnos_image import segment


def sparse_field(
    seed=0, n_cells=25, shape=(320, 320), radius=8, amplitude=40.0, noise=3.0, background=30.0
):
    """Disc cells on a noisy, unevenly illuminated background.

    Deliberately shaped like the team's real frames rather than like a
    textbook test image: one histogram mode, cells at ~2% of the pixels, and
    a slow illumination gradient across the frame.
    """
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    vignette = 1.0 - 0.3 * (((yy - shape[0] / 2) ** 2 + (xx - shape[1] / 2) ** 2) / (shape[0] / 2) ** 2)
    image = background * vignette
    truth = np.zeros(shape, dtype=bool)
    step = shape[0] // int(np.ceil(np.sqrt(n_cells)))
    centres = [(step // 2 + step * (i // 5), step // 2 + step * (i % 5)) for i in range(n_cells)]
    for cy, cx in centres:
        inside = (yy - cy) ** 2 + (xx - cx) ** 2 <= radius**2
        truth |= inside
        image[inside] += amplitude
    return image + rng.normal(0, noise, shape), truth


def test_cells_are_a_small_minority_of_the_test_frame():
    """Guards the fixture itself: if the cells ever stop being sparse, the
    tests below would stop testing the case they exist for."""
    _, truth = sparse_field()
    assert 0.01 < truth.mean() < 0.06


def test_finds_roughly_the_right_number_of_cells():
    image, _ = sparse_field(n_cells=25)
    labels = segment.segment_cells(image, method="sparse", min_size=60)
    assert 20 <= labels.max() <= 30


def test_foreground_stays_plausible_as_the_frame_dims():
    """The real session's frames ranged over roughly this much contrast, and
    'otsu' swung from 0.1% to 87% foreground across them."""
    for amplitude in (40.0, 16.0, 8.0):
        image, _ = sparse_field(amplitude=amplitude)
        labels = segment.segment_cells(image, method="sparse", min_size=60)
        assert 0.005 < segment.foreground_fraction(labels) < 0.1, amplitude


def test_every_cell_is_found_and_the_masks_sit_on_them():
    """Recall is what the method is for; precision is limited by the
    smoothing that buys the robustness, and the docstring says so."""
    image, truth = sparse_field()
    labels = segment.segment_cells(image, method="sparse", min_size=60)
    found = labels > 0
    assert (found & truth).sum() / truth.sum() > 0.95
    assert (found & truth).sum() / found.sum() > 0.6


def test_masks_are_inflated_by_the_smoothing():
    """Documented, not incidental: an inflated mask pulls a per-cell mean
    toward the local background, so anyone reading absolute intensities
    needs to know it is there."""
    image, truth = sparse_field()
    labels = segment.segment_cells(image, method="sparse", min_size=60)
    assert 1.2 < (labels > 0).sum() / truth.sum() < 2.0


def test_boundary_erosion_tightens_the_masks():
    image, truth = sparse_field()
    loose = segment.segment_cells(image, method="sparse", min_size=60) > 0
    tight = segment.segment_cells(image, method="sparse", min_size=60, boundary_erosion_px=2) > 0
    assert (tight & truth).sum() / tight.sum() > (loose & truth).sum() / loose.sum()
    assert (tight & truth).sum() / truth.sum() > 0.9


def test_touching_cells_are_separated():
    image, _ = sparse_field(n_cells=4, shape=(120, 120), radius=12)
    image[50:70, 40:80] += 40.0  # bridge two neighbours
    labels = segment.segment_cells(image, method="sparse", min_size=60)
    assert labels.max() >= 2


def test_empty_frame_yields_no_cells_rather_than_noise():
    rng = np.random.default_rng(0)
    image = rng.normal(30, 3, (200, 200))
    labels = segment.segment_cells(image, method="sparse", min_size=60)
    assert segment.foreground_fraction(labels) < 0.01


def test_a_flat_frame_returns_an_empty_mask_instead_of_dividing_by_zero():
    labels = segment.segment_cells(np.full((64, 64), 42.0), method="sparse")
    assert labels.max() == 0


def test_raising_the_threshold_admits_fewer_cells():
    image, _ = sparse_field(amplitude=12.0)
    loose = segment.segment_cells(image, method="sparse", min_size=60, noise_sigmas=2.0)
    strict = segment.segment_cells(image, method="sparse", min_size=60, noise_sigmas=6.0)
    assert segment.foreground_fraction(strict) < segment.foreground_fraction(loose)


@pytest.mark.parametrize(
    "bad", [dict(noise_sigmas=0), dict(background_sigma_px=-1), dict(boundary_erosion_px=-1)]
)
def test_invalid_parameters_are_rejected(bad):
    image, _ = sparse_field()
    with pytest.raises(ValueError):
        segment.segment_cells(image, method="sparse", **bad)


def test_unknown_method_still_rejected():
    with pytest.raises(ValueError):
        segment.segment_cells(np.zeros((20, 20)), method="sparsey")


def test_foreground_fraction_accepts_labels_and_masks():
    labels = np.zeros((10, 10), dtype=int)
    labels[:2, :] = 7  # one object, 20% of the frame
    assert segment.foreground_fraction(labels) == pytest.approx(0.2)
    assert segment.foreground_fraction(labels > 0) == pytest.approx(0.2)


def test_foreground_fraction_separates_the_two_failure_modes():
    """Object count cannot: a shattered background and a real field both
    produce many cell-sized objects, and a missed field and a genuinely
    sparse one both produce few."""
    shattered = np.zeros((100, 100), dtype=int)
    shattered[:60, :] = np.arange(6000).reshape(60, 100) // 50 + 1
    missed = np.zeros((100, 100), dtype=int)
    missed[:2, :2] = 1
    assert segment.foreground_fraction(shattered) > 0.2
    assert segment.foreground_fraction(missed) < 0.001
