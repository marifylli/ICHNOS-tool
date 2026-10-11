"""segment_cells(method="brightfield_rim"): bright-field masks from each
cell's dark rim, with every size in micrometres."""
import numpy as np
import pytest
from scipy import ndimage as ndi
from skimage import morphology
from skimage.transform import rescale

from ichnos_image import pipeline, segment

PX = 0.1725  # µm per pixel, 40X with the 0.5X adapter


def separated_cells(seed=1, n=25, shape=(400, 480), radius=12):
    """Non-overlapping discs, so one mask per disc is the right answer."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    labels = np.zeros(shape, dtype=np.int32)
    placed = 0
    for _ in range(2000):
        cy = int(rng.integers(3 * radius, shape[0] - 3 * radius))
        cx = int(rng.integers(3 * radius, shape[1] - 3 * radius))
        disc = (yy - cy) ** 2 + (xx - cx) ** 2 <= radius**2
        if (ndi.binary_dilation(disc, morphology.disk(radius)) & (labels > 0)).any():
            continue
        placed += 1
        labels[disc] = placed
        if placed == n:
            break
    return labels


def brightfield(cells, seed=0, rim=(2, 1), halo=4):
    """Dark rim on the cell edge, bright halo outside, uneven background."""
    rng = np.random.default_rng(seed)
    inside = cells > 0
    ring = (ndi.binary_dilation(inside, morphology.disk(rim[0]))
            & ~ndi.binary_erosion(inside, morphology.disk(rim[1])))
    glow = ndi.binary_dilation(inside, morphology.disk(halo)) & ~inside
    yy, xx = np.mgrid[: cells.shape[0], : cells.shape[1]]
    image = 140 + 15 * np.sin(xx / 200) + 10 * np.cos(yy / 150)
    image[glow] += 25
    image[ring] -= 55
    return image + rng.normal(0, 1.5, cells.shape)


def rim(image, pixel_size_um=PX, **kwargs):
    return segment.segment_cells(image, method="brightfield_rim", pixel_size_um=pixel_size_um, **kwargs)


def test_one_mask_per_cell_and_close_to_the_cell():
    cells = separated_cells()
    labels = rim(brightfield(cells))
    assert labels.max() == cells.max()
    found, truth = labels > 0, cells > 0
    iou = (found & truth).sum() / (found | truth).sum()
    assert iou > 0.75


def test_the_halo_is_left_out():
    """The reason for the method: on the real frames 'transmitted' keeps the
    bright halo, which is textured too. Masks must stop at the rim."""
    cells = separated_cells()
    by_rim = rim(brightfield(cells, halo=6)) > 0
    beyond_rim = ~ndi.binary_dilation(cells > 0, morphology.disk(3))
    assert (by_rim & beyond_rim).sum() / by_rim.sum() < 0.02
    assert by_rim.sum() < 1.4 * (cells > 0).sum()


def test_sizes_are_in_micrometres_not_pixels():
    """The same scene sampled twice as finely must give the same cells with
    the same area in µm²; with pixel-based settings it would not."""
    cells = separated_cells(seed=4, n=12, shape=(300, 300))
    image = brightfield(cells)
    coarse = rim(image, PX)
    fine = rim(rescale(image, 2, order=1, preserve_range=True), PX / 2)
    assert coarse.max() == fine.max() == cells.max()
    def area(labels, px):
        return (labels > 0).sum() * px**2

    assert area(fine, PX / 2) == pytest.approx(area(coarse, PX), rel=0.1)


def test_an_empty_frame_with_an_illumination_gradient_stays_empty():
    xx = np.mgrid[:300, :360][1]
    rng = np.random.default_rng(0)
    frame = 80 + 80 * xx / 360 + rng.normal(0, 1.0, (300, 360))
    assert rim(frame).max() == 0


def test_objects_outside_the_area_range_are_dropped():
    cells = separated_cells(seed=2, n=6, shape=(300, 300))
    image = brightfield(cells)
    yy, xx = np.mgrid[:300, :300]
    speck = (yy - 15) ** 2 + (xx - 15) ** 2 <= 9      # ~0.8 µm², debris
    image[speck] -= 55
    labels = rim(image)
    assert labels.max() == cells.max()
    assert labels[15, 15] == 0
    assert rim(image, max_area_um2=5.0).max() == 0  # cells (~13 µm²) above the cap


def test_touching_cells_are_split():
    yy, xx = np.mgrid[:200, :260]
    cells = np.zeros((200, 260), dtype=np.int32)
    cells[(yy - 100) ** 2 + (xx - 110) ** 2 <= 144] = 1
    cells[((yy - 100) ** 2 + (xx - 133) ** 2 <= 144) & (cells == 0)] = 2
    labels = rim(brightfield(cells, rim=(1, 1)), min_solidity=0.0)
    assert len(np.unique(labels[labels > 0])) == 2


def test_a_flat_frame_yields_no_cells():
    assert rim(np.full((64, 64), 120.0)).max() == 0


def test_pixel_size_is_required():
    with pytest.raises(TypeError):
        segment.segment_cells(np.zeros((40, 40)), method="brightfield_rim")


@pytest.mark.parametrize("bad", [{"pixel_size_um": 0}, {"pixel_size_um": float("nan")},
                                 {"rim_sigmas": 0}, {"min_area_um2": 0},
                                 {"max_area_um2": 1.0, "min_area_um2": 4.0}, {"min_solidity": 1.5}])
def test_invalid_parameters_are_rejected(bad):
    kwargs = {"pixel_size_um": PX} | bad
    with pytest.raises(ValueError):
        segment.segment_cells(np.zeros((40, 40)), method="brightfield_rim", **kwargs)


def test_the_pipeline_takes_the_pixel_size_from_the_objective(monkeypatch):
    seen = {}

    def spy(image, method, **kwargs):
        seen.update(method=method, **kwargs)
        raise RuntimeError("stop after segmentation")

    monkeypatch.setattr(segment, "segment_cells", spy)
    cells = separated_cells(seed=3, n=8, shape=(200, 240), radius=10)
    image_set = pipeline.ImageSet(
        session_id="s", timepoint=0, acquisition_order=0,
        green=np.full(cells.shape, 40.0), red=np.full(cells.shape, 40.0),
        bright_field=brightfield(cells),
        exposure_ms_green=100, exposure_ms_red=100, nd_filter_green=0, nd_filter_red=0,
        objective="40X", burner_hours=1, lamp_warmup_minutes=30, sample_id="a",
    )
    with pytest.raises(RuntimeError, match="stop after segmentation"):
        pipeline.process_image_set(image_set, bleed_green_to_red=0.0, segmentation_method="brightfield_rim",
                                   segmentation_source="brightfield")
    assert seen["method"] == "brightfield_rim"
    assert seen["pixel_size_um"] == pytest.approx(PX)
