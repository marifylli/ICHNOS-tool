"""A batch must survive one unusable frame. The frames already processed are
work that has been paid for, and discarding them is how the team lost half an
hour to a refusal on the last of 25 frames."""
import json
from pathlib import Path

import numpy as np
import pytest

from ichnos_image import pipeline
from ichnos_image.pipeline import ImageSet, SegmentationQCError, process_experiment


def image_set(order, *, cells=True, **overrides):
    rng = np.random.default_rng(order)
    yy, xx = np.mgrid[:120, :120]
    green = 30 + rng.normal(0, 2, (120, 120))
    if cells:
        for cy, cx in [(30, 30), (30, 90), (90, 30), (90, 90)]:
            green[(yy - cy) ** 2 + (xx - cx) ** 2 <= 64] += 60
    fields = dict(
        session_id="s", timepoint=order, acquisition_order=order, green=green,
        red=green * 0.8 + 5, exposure_ms_green=100, exposure_ms_red=100,
        nd_filter_green=0, nd_filter_red=0, objective="40X", burner_hours=1,
        lamp_warmup_minutes=30, sample_id=f"sample-{order}",
    )
    fields.update(overrides)
    return ImageSet(**fields)


def segmentation(per_sample):
    """Return a fake segment_cells that answers by the frame's sample_id."""
    state = {"calls": 0}

    def fake(image, **kwargs):
        labels = per_sample[state["calls"]]
        state["calls"] += 1
        return labels
    return fake


def good_labels():
    labels = np.zeros((120, 120), dtype=np.int32)
    yy, xx = np.mgrid[:120, :120]
    for i, (cy, cx) in enumerate([(30, 30), (30, 90), (90, 30), (90, 90)], start=1):
        labels[(yy - cy) ** 2 + (xx - cx) ** 2 <= 64] = i
    return labels


def shattered_labels():
    labels = np.zeros((120, 120), dtype=np.int32)
    labels[:72, :] = np.arange(72 * 120).reshape(72, 120) // 60 + 1
    return labels


@pytest.fixture
def fake_segmentation(monkeypatch):
    def install(sequence):
        monkeypatch.setattr(pipeline.segment, "segment_cells", segmentation(sequence))
    return install


def run(tmp_path, image_sets, name="out.csv", **options):
    return process_experiment(
        image_sets, tmp_path / name, bleed_green_to_red=0.05, **options
    )


def manifest_of(out_path):
    return json.loads(Path(str(out_path) + ".manifest.json").read_text())


def test_one_refused_frame_does_not_discard_the_others(tmp_path, fake_segmentation):
    fake_segmentation([good_labels(), shattered_labels(), good_labels()])
    out = run(tmp_path, [image_set(i) for i in range(3)])
    assert out.exists()
    assert sum(1 for _ in out.open()) - 1 == 8  # two good frames, four cells each


def test_the_refusal_is_recorded_in_the_manifest(tmp_path, fake_segmentation):
    """A run that quietly has fewer frames than its manifest must not be
    mistakable for a complete one."""
    fake_segmentation([good_labels(), shattered_labels(), good_labels()])
    out = run(tmp_path, [image_set(i) for i in range(3)])
    refused = manifest_of(out)["refused_image_sets"]
    assert [item["sample_id"] for item in refused] == ["sample-1"]
    assert "60.0%" in refused[0]["reason"]


def test_a_clean_run_records_no_refusals(tmp_path, fake_segmentation):
    fake_segmentation([good_labels(), good_labels()])
    out = run(tmp_path, [image_set(i) for i in range(2)])
    assert manifest_of(out)["refused_image_sets"] == []


def test_a_batch_refused_outright_raises_rather_than_publishing_nothing(tmp_path, fake_segmentation):
    """Carrying on past every frame would publish an empty CSV that looks
    like a successful run of a session with no cells in it."""
    fake_segmentation([shattered_labels(), shattered_labels()])
    with pytest.raises(SegmentationQCError, match="every one of the 2 frames"):
        run(tmp_path, [image_set(i) for i in range(2)])


def test_a_declared_blank_control_does_not_count_as_a_refusal(tmp_path, fake_segmentation):
    """The real case: a medium-only control alongside sample frames."""
    empty = np.zeros((120, 120), dtype=np.int32)
    fake_segmentation([good_labels(), empty])
    out = run(tmp_path, [image_set(0), image_set(1, cells=False, expect_cells=False)])
    assert manifest_of(out)["refused_image_sets"] == []
