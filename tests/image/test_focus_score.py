"""segment.focus_score() must be comparable across bit depths/exposure
levels, not just within one camera's raw pixel-value scale -- this matters
concretely: the team's real Olympus camera (SC30) is 8-bit/channel (0-255),
while the public reference images FOCUS_SCORE_THRESHOLD was calibrated on
are 16-bit (0-65535). A non-normalized variance-of-Laplacian would differ by
~(65535/255)^2 =~ 66000x between otherwise-similar images at those two
depths, making a threshold calibrated on one meaningless on the other.
"""
import numpy as np

from ichnos_image import segment


def _textured_field(size=128, seed=0):
    rng = np.random.default_rng(seed)
    field = rng.normal(0.5, 0.12, size=(size, size))
    yy, xx = np.mgrid[0:size, 0:size]
    for cy, cx in rng.integers(15, size - 15, size=(10, 2)):
        blob = ((yy - cy) ** 2 + (xx - cx) ** 2) < 8**2
        field[blob] -= 0.3
    return np.clip(field, 0, 1)


def test_focus_score_invariant_to_bit_depth_rescaling():
    field = _textured_field()

    score_8bit = segment.focus_score(field * 255)  # simulated 8-bit (SC30) scale
    score_16bit = segment.focus_score(field * 65535)  # simulated 16-bit (public dataset) scale

    assert score_8bit == score_16bit  # same relative content -> same score, regardless of scale


def test_focus_score_still_separates_sharp_from_blurred_at_both_depths():
    from scipy.ndimage import gaussian_filter

    sharp = _textured_field()
    blurred = gaussian_filter(sharp, sigma=3.0)

    for scale in (255, 65535):
        sharp_score = segment.focus_score(sharp * scale)
        blurred_score = segment.focus_score(blurred * scale)
        assert sharp_score > blurred_score * 5  # blur still clearly reduces the score
