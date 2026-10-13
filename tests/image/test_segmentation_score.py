"""Scoring masks against hand-marked cell centres."""
import numpy as np
import pandas as pd
import pytest

from ichnos_image.segmentation_score import read_points, score_labels


def frame():
    labels = np.zeros((100, 100), dtype=np.int32)
    labels[10:20, 10:20] = 1   # one cell, one point -> separated
    labels[10:20, 40:60] = 2   # two points -> merged
    labels[60:70, 10:20] = 3   # no point -> false
    labels[60:70, 60:70] = 4   # unsure point only -> left out
    return labels


def test_each_outcome_is_counted():
    points = np.array([[15, 15], [15, 45], [15, 55], [85, 85]])  # last one is in no mask
    score, detail = score_labels(frame(), points, unsure=np.array([[65, 65]]))
    assert (score.marked, score.detected, score.separated) == (4, 3, 1)
    assert (score.merged_masks, score.false_masks, score.masks_scored) == (1, 1, 3)
    assert score.recall == pytest.approx(1 / 4)
    assert score.precision == pytest.approx(1 / 3)
    assert detail["masks"] == {1: "separated", 2: "merged", 3: "false", 4: "unsure"}
    assert detail["missed"].tolist() == [[85, 85]]


def test_a_click_just_outside_the_mask_still_counts():
    labels = frame()
    near = np.array([[15, 22]])                 # 2 px right of mask 1
    assert score_labels(labels, near, tolerance_px=3)[0].detected == 1
    assert score_labels(labels, near, tolerance_px=1)[0].detected == 0


def test_masks_cut_by_the_crop_edge_are_judged_by_their_whole_centroid():
    """A mask whose centre lies outside the crop belongs to a cell nobody was
    asked to mark; it must not count as false. One whose centre is inside
    must."""
    labels = np.zeros((100, 100), dtype=np.int32)
    labels[45:55, 30:70] = 1     # centroid x = 49.5, outside a crop starting at x0=50
    labels[45:55, 72:80] = 2     # centroid inside the crop
    score, detail = score_labels(labels, np.zeros((0, 2)), window=(0, 50, 100, 50))
    assert detail["masks"] == {2: "false"}
    assert score.false_masks == 1


def test_points_are_in_crop_coordinates():
    labels = np.zeros((100, 100), dtype=np.int32)
    labels[60:70, 60:70] = 1
    score, _ = score_labels(labels, np.array([[5, 5]]), window=(60, 60, 30, 30))
    assert score.separated == 1


def test_a_window_outside_the_frame_is_rejected():
    with pytest.raises(ValueError, match="outside"):
        score_labels(frame(), np.zeros((0, 2)), window=(90, 90, 20, 20))


def test_fiji_results_table_is_read_as_row_column(tmp_path):
    path = tmp_path / "a.points.csv"
    pd.DataFrame({" ": [1, 2], "X": [12.5, 40.0], "Y": [3.0, 77.5]}).to_csv(path, index=False)
    assert read_points(path).tolist() == [[3.0, 12.5], [77.5, 40.0]]


def test_a_table_without_coordinates_is_rejected(tmp_path):
    path = tmp_path / "a.points.csv"
    pd.DataFrame({"Area": [1.0]}).to_csv(path, index=False)
    with pytest.raises(ValueError, match="X and Y"):
        read_points(path)
