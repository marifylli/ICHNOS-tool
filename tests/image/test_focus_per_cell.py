"""focus.score_cells must fall with defocus, ignore exposure/gain, and judge
each cell against its own background; focus.channel_agreement must catch a
channel acquired at a different focal plane."""
import numpy as np
import pytest
from scipy.ndimage import gaussian_filter

from ichnos_image import focus


def _field(seed=0, n_cells=12, shape=(220, 220), radius=7, amplitude=40.0, noise=3.0):
    """Disc 'cells' on a noisy background, with their label mask."""
    rng = np.random.default_rng(seed)
    labels = np.zeros(shape, dtype=int)
    image = np.full(shape, 30.0)
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    centres = [(30 + 50 * (i // 4), 30 + 50 * (i % 4)) for i in range(n_cells)]
    for cell_id, (cy, cx) in enumerate(centres, start=1):
        inside = (yy - cy) ** 2 + (xx - cx) ** 2 <= radius**2
        labels[inside] = cell_id
        image[inside] += amplitude
    return labels, image + rng.normal(0, noise, shape)


def _median_score(labels, image, **kwargs):
    scores = [s.focus_score for s in focus.score_cells(labels, image, **kwargs) if s.scored]
    assert scores, "no cell could be scored"
    return float(np.median(scores))


def test_score_falls_monotonically_with_defocus():
    labels, image = _field()
    scores = [_median_score(labels, gaussian_filter(image, sigma)) for sigma in (0, 2, 4, 6)]
    assert scores == sorted(scores, reverse=True), scores
    assert scores[0] > 1.5 * scores[-1]


def test_score_is_invariant_to_exposure_and_gain():
    labels, image = _field()
    baseline = _median_score(labels, image)
    for gain in (0.25, 4.0):
        assert _median_score(labels, image * gain) == pytest.approx(baseline, rel=1e-9)


def test_dim_cell_is_not_penalised_for_being_dim():
    """A faint but sharp cell must not score like a blurred one: the metric
    is normalized per cell, so low amplitude alone is not defocus."""
    labels, bright = _field(amplitude=80.0, noise=1.0)
    _, faint = _field(amplitude=8.0, noise=1.0)
    assert _median_score(labels, faint) > 0.6 * _median_score(labels, bright)


def test_blurred_cell_approaches_the_no_edge_floor():
    labels, image = _field()
    assert _median_score(labels, gaussian_filter(image, 12)) == pytest.approx(
        focus.NO_EDGE_SCORE, abs=0.15
    )


def test_score_never_falls_below_the_floor():
    labels, image = _field()
    for sigma in (0, 3, 10):
        scored = [s for s in focus.score_cells(labels, gaussian_filter(image, sigma)) if s.scored]
        assert all(s.focus_score >= focus.NO_EDGE_SCORE for s in scored)


def test_unscorable_cell_is_flagged_not_defaulted():
    """A cell with no background ring must come back NaN and scored=False,
    never a number that downstream code would compare against a threshold."""
    labels = np.ones((40, 40), dtype=int)  # every pixel is one cell: no ring
    image = np.random.default_rng(0).normal(100, 5, (40, 40))
    (result,) = focus.score_cells(labels, image)
    assert not result.scored and np.isnan(result.focus_score)


def test_shape_mismatch_is_rejected():
    with pytest.raises(ValueError):
        focus.score_cells(np.zeros((10, 10), dtype=int), np.zeros((10, 11)))


def test_agreement_high_when_both_channels_share_a_focal_plane():
    labels, image = _field()
    _, other = _field(seed=1)
    agreement = focus.channel_agreement(
        focus.score_cells(labels, image), focus.score_cells(labels, other)
    )
    assert agreement.agrees and agreement.agreement > 0.6


def test_agreement_low_when_one_channel_is_defocused():
    labels, image = _field()
    agreement = focus.channel_agreement(
        focus.score_cells(labels, gaussian_filter(image, 6)),
        focus.score_cells(labels, image),
    )
    assert not agreement.agrees
    assert agreement.median_green < agreement.median_red


def test_agreement_reports_not_informative_with_too_few_shared_cells():
    labels, image = _field(n_cells=2)
    agreement = focus.channel_agreement(
        focus.score_cells(labels, image), focus.score_cells(labels, image)
    )
    assert not agreement.agrees and np.isnan(agreement.agreement)
