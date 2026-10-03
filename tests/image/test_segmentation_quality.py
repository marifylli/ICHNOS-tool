"""Regression guard for the quantitative segmentation-validation finding:
segmentation quality (F1 against synthetic ground truth) shouldn't collapse
under stress-like morphology, and shouldn't silently regress below what was
measured when this was calibrated. Scene generation and scoring come from
tests/helpers/segmentation_scenes.py -- see its docstring for why synthetic,
not real, ground truth is used here.
"""
from tests.helpers import segmentation_scenes as vsq
from ichnos_image import segment


def test_baseline_and_stressed_f1_above_floor():
    for stressed in (False, True):
        scores = []
        for seed in range(3):
            bf, gt = vsq._make_scene(seed=seed, stressed=stressed, n_cells=15, size=200)
            pred = segment.segment_cells(bf, method="otsu", min_size=30)
            scores.append(vsq._match_and_score(pred, gt))
        mean_f1 = sum(s["f1"] for s in scores) / len(scores)
        assert mean_f1 > 0.6, f"stressed={stressed}: F1 {mean_f1:.3f} fell below floor"


def test_stressed_f1_not_much_worse_than_baseline():
    def mean_f1(stressed):
        scores = [
            vsq._match_and_score(
                segment.segment_cells(bf, method="otsu", min_size=30), gt
            )
            for bf, gt in (vsq._make_scene(seed=s, stressed=stressed, n_cells=15, size=200) for s in range(3))
        ]
        return sum(s["f1"] for s in scores) / len(scores)

    baseline_f1 = mean_f1(False)
    stressed_f1 = mean_f1(True)
    assert stressed_f1 > baseline_f1 - 0.15
