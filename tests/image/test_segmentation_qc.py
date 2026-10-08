"""The foreground-fraction guard must stop a frame whose segmentation failed,
before anything runs per cell, and must say which way it failed."""
import numpy as np
import pytest

from ichnos_image import pipeline
from ichnos_image.pipeline import ImageSet, SegmentationQCError


def image_set(**overrides):
    rng = np.random.default_rng(0)
    yy, xx = np.mgrid[:160, :160]
    green = 30 + rng.normal(0, 2, (160, 160))
    for cy, cx in [(40, 40), (40, 110), (110, 40), (110, 110)]:
        green[(yy - cy) ** 2 + (xx - cx) ** 2 <= 64] += 60
    fields = dict(
        session_id="s", timepoint=0, acquisition_order=0, green=green, red=green * 0.8 + 5,
        exposure_ms_green=100, exposure_ms_red=100, nd_filter_green=0, nd_filter_red=0,
        objective="40X", burner_hours=1, lamp_warmup_minutes=30,
    )
    fields.update(overrides)
    return ImageSet(**fields)


def run(labels, image_set_overrides=None, **kwargs):
    """Drive process_image_set with a segmentation of our choosing."""
    import ichnos_image.segment as segment

    original = segment.segment_cells
    segment.segment_cells = lambda *a, **k: labels
    try:
        return pipeline.process_image_set(
            image_set(**(image_set_overrides or {})), bleed_green_to_red=0.05, **kwargs
        )
    finally:
        segment.segment_cells = original


def shattered_labels(shape=(160, 160), share=0.5):
    """What a threshold inside the background produces: many cell-sized
    objects covering a large share of the frame."""
    labels = np.zeros(shape, dtype=np.int32)
    rows = int(shape[0] * share)
    block = np.arange(rows * shape[1]).reshape(rows, shape[1])
    labels[:rows, :] = block // 100 + 1
    return labels


def plausible_labels(shape=(160, 160)):
    labels = np.zeros(shape, dtype=np.int32)
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    for i, (cy, cx) in enumerate([(40, 40), (40, 110), (110, 40), (110, 110)], start=1):
        labels[(yy - cy) ** 2 + (xx - cx) ** 2 <= 64] = i
    return labels


def test_a_plausible_field_passes():
    records = run(plausible_labels())
    assert len(records) == 4


def test_shattered_background_is_refused():
    with pytest.raises(SegmentationQCError) as excinfo:
        run(shattered_labels())
    assert "50.0%" in str(excinfo.value)


def test_a_frame_with_almost_nothing_found_is_refused():
    """The quieter failure: the threshold landed above the signal, so the
    few survivors are the brightest cells rather than a sample of the field.
    Nothing downstream would look wrong."""
    labels = np.zeros((160, 160), dtype=np.int32)
    labels[:2, :4] = 1
    with pytest.raises(SegmentationQCError):
        run(labels)


def test_the_message_names_the_frame_and_suggests_a_method():
    with pytest.raises(SegmentationQCError) as excinfo:
        run(shattered_labels())
    message = str(excinfo.value)
    assert "'s'" in message and "sparse" in message


def test_bounds_are_adjustable():
    labels = shattered_labels(share=0.3)
    with pytest.raises(SegmentationQCError):
        run(labels)
    assert run(labels, max_foreground_fraction=0.5)


def test_the_guard_runs_before_the_per_cell_work():
    """The point of the check is cost: an unguarded failure on the team's
    real session spent 15 hours extracting features from noise. If this ever
    moves after extraction it stops being a guard."""
    import ichnos_image.extract as extract

    original = extract.extract_per_cell
    extract.extract_per_cell = lambda *a, **k: pytest.fail("extraction ran on a refused frame")
    try:
        with pytest.raises(SegmentationQCError):
            run(shattered_labels())
    finally:
        extract.extract_per_cell = original


def test_a_declared_blank_frame_is_allowed_to_be_empty():
    """A medium-only control has no cells by design. Judging it by the rule
    for sample frames reports a segmentation failure where the segmentation
    was right -- which is what aborted the team's first guarded run."""
    labels = np.zeros((160, 160), dtype=np.int32)
    labels[:2, :4] = 1
    records = run(labels, image_set_overrides=dict(expect_cells=False))
    assert len(records) == 1  # whatever was found is still reported, not refused


def test_cells_in_a_declared_blank_frame_are_refused():
    """The inverted failure: if a frame that should be empty is full, either
    it is contaminated or the threshold is reading noise as cells on every
    frame in the session."""
    with pytest.raises(SegmentationQCError) as excinfo:
        run(plausible_labels(), image_set_overrides=dict(expect_cells=False))
    assert "declared cell-free" in str(excinfo.value)


def test_blank_ceiling_is_adjustable():
    labels = plausible_labels()
    assert run(labels, image_set_overrides=dict(expect_cells=False),
               blank_max_foreground_fraction=0.9)
