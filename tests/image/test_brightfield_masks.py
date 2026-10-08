"""Bright-field masks: the only mask source that cannot favour a channel.

Three things have to hold before they can be used. The segmentation has to
find cells in transmitted light, where they are neither bright nor dark; the
alignment has to say whether the bright-field frame shows the same field at
all, and by how much it is displaced; and the pipeline has to apply that
displacement and refuse the frames that have no bright-field rather than
quietly falling back to a fluorescence channel.
"""
import numpy as np
import pytest

from ichnos_image import pipeline, segment
from ichnos_image.align import MaskAlignment, align_mask, shift_mask
from ichnos_image.pipeline import ImageSet, MissingBrightfieldError, SegmentationQCError
from scipy import ndimage as ndi
from skimage import morphology


# ------------------------------------------------------------- fixtures --

def cell_field(seed=1, n=40, shape=(320, 400), radius=9):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    cells = np.zeros(shape, dtype=bool)
    for _ in range(n):
        cy = int(rng.integers(radius * 3, shape[0] - radius * 3))
        cx = int(rng.integers(radius * 3, shape[1] - radius * 3))
        cells |= (yy - cy) ** 2 + (xx - cx) ** 2 <= radius**2
    return cells


def as_brightfield(cells, seed=0):
    """Dark rim, bright halo, smooth uneven background -- the thing that
    makes transmitted light different from a fluorescence frame."""
    rng = np.random.default_rng(seed)
    rim = (ndi.binary_dilation(cells, morphology.disk(2))
           & ~ndi.binary_erosion(cells, morphology.disk(1)))
    halo = ndi.binary_dilation(cells, morphology.disk(4)) & ~cells
    yy, xx = np.mgrid[: cells.shape[0], : cells.shape[1]]
    image = 140 + 12 * np.sin(xx / 160) + 8 * np.cos(yy / 120)
    image[halo] += 25
    image[rim] -= 55
    return image + rng.normal(0, 1.5, cells.shape)


def as_fluorescence(cells, lit=0.3, seed=0, amplitude=70.0):
    """Only some of the cells fluoresce, as in the real session."""
    rng = np.random.default_rng(seed)
    labels, count = ndi.label(cells)
    keep = rng.random(count + 1) < lit
    keep[0] = False
    image = 40 + rng.normal(0, 2.5, cells.shape)
    image[keep[labels]] += amplitude
    return image


# -------------------------------------------------- segmentation itself --

def test_transmitted_masks_the_cells_where_sparse_masks_the_halos():
    """The failure this method exists for. Both find the cells -- recall is
    not the difference -- but 'sparse' looks for bright things, and in
    transmitted light the bright thing is the halo around the cell. It
    returns twice the area it should, so every per-cell mean it produces is
    pulled toward the background it has swallowed."""
    cells = cell_field()
    image = as_brightfield(cells)
    transmitted = segment.segment_cells(image, method="transmitted", min_size=100) > 0
    sparse = segment.segment_cells(image, method="sparse", min_size=100) > 0

    precision = lambda mask: (mask & cells).sum() / max(mask.sum(), 1)  # noqa: E731
    assert (transmitted & cells).sum() / cells.sum() > 0.9
    assert precision(transmitted) > 1.8 * precision(sparse)
    assert sparse.mean() > 1.8 * cells.mean()        # halos included
    assert transmitted.mean() < 1.4 * cells.mean()   # roughly the cells


def test_transmitted_ignores_the_illumination_gradient():
    """Local structure, not brightness: an empty frame that is twice as
    bright on one side must still come back empty."""
    yy, xx = np.mgrid[:200, :260]
    rng = np.random.default_rng(0)
    gradient = 80 + 80 * xx / 260 + rng.normal(0, 1.0, (200, 260))
    labels = segment.segment_cells(gradient, method="transmitted", min_size=100)
    assert segment.foreground_fraction(labels) < 0.01


def test_transmitted_masks_are_tightened_toward_the_cells():
    cells = cell_field()
    image = as_brightfield(cells)
    loose = segment.segment_cells(image, method="transmitted", min_size=100, erosion_px=0) > 0
    tight = segment.segment_cells(image, method="transmitted", min_size=100) > 0
    assert tight.sum() < loose.sum()
    assert (tight & cells).sum() / tight.sum() > (loose & cells).sum() / loose.sum()


def test_erosion_is_backed_off_rather_than_erasing_small_cells():
    """A cell smaller than the erosion radius would vanish, taking the row
    with it and leaving no sign in the output that it ever existed."""
    small = cell_field(radius=4, n=20)
    labels = segment.segment_cells(as_brightfield(small), method="transmitted",
                                   min_size=40, erosion_px=12)
    assert labels.max() > 10


@pytest.mark.parametrize("bad", [dict(texture_window_px=1), dict(noise_sigmas=0),
                                 dict(erosion_px=-1), dict(closing_px=-1)])
def test_invalid_parameters_are_rejected(bad):
    with pytest.raises(ValueError):
        segment.segment_cells(np.zeros((40, 40)), method="transmitted", **bad)


def test_a_flat_frame_yields_no_cells():
    labels = segment.segment_cells(np.full((64, 64), 120.0), method="transmitted")
    assert labels.max() == 0


# --------------------------------------------------------- the alignment --

@pytest.mark.parametrize("dy,dx", [(0, 0), (12, 0), (0, -9), (15, -21), (-30, 25)])
def test_the_displacement_convention_is_exact(dy, dx):
    """On a case with no ambiguity. Anything but an exact answer here is a
    sign or axis error, which would shift every mask the wrong way and show
    up only as a quietly wrong measurement."""
    yy, xx = np.mgrid[:200, :260]
    mask = (yy - 80) ** 2 + (xx - 100) ** 2 <= 100
    image = np.full((200, 260), 40.0)
    image[(yy - 80 - dy) ** 2 + (xx - 100 - dx) ** 2 <= 100] += 80
    found = align_mask(mask, image, max_shift=80)
    assert (found.dy, found.dx) == (dy, dx)


def test_a_real_shift_is_recovered_when_only_some_cells_fluoresce():
    cells = cell_field()
    shift = (17, -23)
    lit = shift_mask(cells, *shift)
    found = align_mask(segment.segment_cells(as_brightfield(cells), method="transmitted",
                                             min_size=100) > 0,
                       as_fluorescence(lit))
    assert abs(found.dy - shift[0]) <= 8 and abs(found.dx - shift[1]) <= 8
    assert found.contrast > 1.05


def test_the_matching_field_stands_clear_of_an_unrelated_one():
    """The whole verdict rests on this separation, and on nothing absolute:
    z has no meaning except against the same number on a field the mask does
    not belong to."""
    cells = cell_field(seed=1)
    mask = segment.segment_cells(as_brightfield(cells), method="transmitted", min_size=100) > 0
    same = align_mask(mask, as_fluorescence(shift_mask(cells, 17, -23)))
    other = align_mask(mask, as_fluorescence(cell_field(seed=2)))
    assert same.z > 2 * other.z


def test_overlap_is_normalised_so_large_shifts_are_not_penalised():
    """Without dividing by the part of the mask still inside the frame, a
    large displacement scores low merely because the mask has run off the
    edge, and the answer is dragged toward zero."""
    yy, xx = np.mgrid[:200, :260]
    mask = (yy - 100) ** 2 + (xx - 130) ** 2 <= 400
    image = np.full((200, 260), 40.0)
    image[(yy - 20) ** 2 + (xx - 30) ** 2 <= 400] += 80  # shifted by (-80, -100)
    found = align_mask(mask, image, max_shift=120)
    assert (found.dy, found.dx) == (-80, -100)


def test_shift_mask_does_not_wrap():
    mask = np.zeros((10, 10), dtype=bool)
    mask[0:2, 0:2] = True
    assert not shift_mask(mask, -3, 0).any()


def test_an_empty_mask_is_reported_rather_than_divided_by():
    found = align_mask(np.zeros((40, 40), dtype=bool), np.ones((40, 40)))
    assert found.foreground_fraction == 0.0 and np.isnan(found.z)


def test_mismatched_shapes_are_rejected():
    with pytest.raises(ValueError, match="same shape"):
        align_mask(np.zeros((10, 10), dtype=bool), np.zeros((10, 12)))


# ------------------------------------------------------------- pipeline --

def image_set(**overrides):
    cells = cell_field(seed=3, n=25, shape=(200, 240), radius=8)
    fields = dict(
        session_id="s", timepoint=0, acquisition_order=0,
        green=as_fluorescence(cells, lit=1.0, seed=4),
        red=as_fluorescence(cells, lit=1.0, seed=5, amplitude=90.0),
        bright_field=as_brightfield(cells),
        exposure_ms_green=100, exposure_ms_red=100, nd_filter_green=0, nd_filter_red=0,
        objective="40X", burner_hours=1, lamp_warmup_minutes=30, sample_id="a",
    )
    fields.update(overrides)
    return ImageSet(**fields)


def run(**kwargs):
    overrides = kwargs.pop("image_set_overrides", {})
    return pipeline.process_image_set(
        image_set(**overrides), bleed_green_to_red=0.0,
        segmentation_method="transmitted", segmentation_source="brightfield",
        segmentation_kwargs=dict(min_size=100), **kwargs,
    )


def test_brightfield_masks_reach_the_records():
    records = run()
    assert records and {r.mask_source for r in records} == {"brightfield"}


def test_the_recorded_shift_is_applied_to_the_masks():
    """A mask displaced by the stated offset must end up where the
    unshifted one was -- otherwise the column is decorative."""
    cells = cell_field(seed=3, n=25, shape=(200, 240), radius=8)
    moved = as_brightfield(shift_mask(cells, -14, 11))
    corrected = run(image_set_overrides=dict(
        bright_field=moved, brightfield_shift_dy=14, brightfield_shift_dx=-11))
    uncorrected = run(image_set_overrides=dict(bright_field=moved))
    assert len(corrected) >= 0.7 * len(run())

    # The corrected masks land on the cells; the uncorrected ones land
    # beside them, on background.
    brightness = lambda rows: float(np.median([r.corrected_mean_green for r in rows]))  # noqa: E731
    assert brightness(corrected) > 2 * brightness(uncorrected)


def test_a_frame_without_a_brightfield_is_refused_not_quietly_switched():
    with pytest.raises(MissingBrightfieldError, match="bright_field"):
        run(image_set_overrides=dict(bright_field=None))


def test_that_refusal_is_both_a_value_error_and_a_qc_error():
    """Called directly it is a caller's mistake; inside a batch it is a
    property of the data and must be survivable."""
    assert issubclass(MissingBrightfieldError, ValueError)
    assert issubclass(MissingBrightfieldError, SegmentationQCError)


def test_a_batch_survives_a_row_with_no_brightfield(tmp_path):
    """The real case: the session's bright-field frames cover one field per
    folder, so folders imaged at two fields have one row without one."""
    out = pipeline.process_experiment(
        [image_set(sample_id="has-one"),
         image_set(sample_id="has-none", bright_field=None, acquisition_order=1)],
        tmp_path / "out.csv", bleed_green_to_red=0.0,
        segmentation_method="transmitted", segmentation_source="brightfield",
        segmentation_kwargs=dict(min_size=100),
    )
    import json
    refused = json.loads((tmp_path / "out.csv.manifest.json").read_text())["refused_image_sets"]
    assert [item["sample_id"] for item in refused] == ["has-none"]
    assert out.exists() and sum(1 for _ in out.open()) > 1


def test_shifting_labels_keeps_them_whole():
    labels = np.zeros((20, 20), dtype=np.int32)
    labels[5:9, 5:9] = 7
    moved = pipeline._shift_labels(labels, 3, -2)
    assert (moved == 7).sum() == 16
    assert moved[8:12, 3:7].all()
