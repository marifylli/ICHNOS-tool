"""Which image the masks are cut from must be a stated choice, recorded on
every row, because it moves the ratio those rows exist to report."""
import numpy as np
import pytest

from ichnos_image import focus, pipeline
from ichnos_image.pipeline import ImageSet


def two_channel_field(shape=(160, 160), seed=0):
    """Two populations: one bright in green only, one bright in red only.

    This is the situation the choice of mask source decides between. A real
    field is not this stark, but the bias it produces is the same one, and
    here it is unambiguous rather than a matter of degree.
    """
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    green = 30 + rng.normal(0, 1.5, shape)
    red = 30 + rng.normal(0, 1.5, shape)
    for cy, cx in [(40, 40), (40, 80), (40, 120)]:
        green[(yy - cy) ** 2 + (xx - cx) ** 2 <= 81] += 60
    for cy, cx in [(120, 40), (120, 80), (120, 120)]:
        red[(yy - cy) ** 2 + (xx - cx) ** 2 <= 81] += 60
    return green, red


def image_set(**overrides):
    green, red = two_channel_field()
    fields = dict(
        session_id="s", timepoint=0, acquisition_order=0, green=green, red=red,
        exposure_ms_green=100, exposure_ms_red=100, nd_filter_green=0, nd_filter_red=0,
        objective="40X", burner_hours=1, lamp_warmup_minutes=30,
    )
    fields.update(overrides)
    return ImageSet(**fields)


def run(**kwargs):
    return pipeline.process_image_set(
        image_set(**kwargs.pop("image_set_overrides", {})),
        bleed_green_to_red=0.0, segmentation_method="sparse",
        segmentation_kwargs=dict(min_size=60), **kwargs,
    )


def median_ratio(records):
    return float(np.median([r.ratio_red_green for r in records if r.ratio_red_green is not None]))


def test_the_mask_source_moves_the_ratio():
    """The finding this option exists for. If this ever stops holding, the
    option is no longer load-bearing and the default can be revisited."""
    from_green = median_ratio(run(segmentation_source="green"))
    from_red = median_ratio(run(segmentation_source="red"))
    assert from_red > 3 * from_green


def test_the_default_sum_sits_between_the_two_single_channel_choices():
    from_green = median_ratio(run(segmentation_source="green"))
    from_red = median_ratio(run(segmentation_source="red"))
    from_sum = median_ratio(run())  # default
    assert from_green < from_sum < from_red


def test_the_default_finds_cells_of_both_populations():
    """Summing is symmetric, so neither channel's cells are passed over."""
    records = run()
    assert len(records) > len(run(segmentation_source="green"))
    assert len(records) > len(run(segmentation_source="red"))


def test_every_row_records_the_mask_source():
    assert {r.mask_source for r in run()} == {"green+red"}
    assert {r.mask_source for r in run(segmentation_source="red")} == {"red"}


def test_brightfield_is_used_when_asked_for():
    green, red = two_channel_field()
    records = run(image_set_overrides=dict(bright_field=green + red),
                  segmentation_source="brightfield")
    assert {r.mask_source for r in records} == {"brightfield"}


def test_asking_for_brightfield_without_one_is_an_error_not_a_silent_fallback():
    """Falling back to a fluorescence channel would reintroduce the bias the
    caller asked to avoid, and nothing in the output would say so."""
    with pytest.raises(ValueError, match="bright_field"):
        run(segmentation_source="brightfield")


def test_an_unknown_source_is_rejected():
    with pytest.raises(ValueError, match="segmentation_source"):
        run(segmentation_source="gren")


@pytest.mark.parametrize(
    "channel,mask_source,biased",
    [("green", "green", True), ("red", "red", True), ("green", "red", False),
     ("green", "green+red", False), ("red", "brightfield", False)],
)
def test_mask_source_biases_flags_only_the_self_selected_channel(channel, mask_source, biased):
    assert focus.mask_source_biases(channel, mask_source) is biased
