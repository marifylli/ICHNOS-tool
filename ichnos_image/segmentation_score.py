"""Score a label mask against cells marked by hand with one point each.

The reference is a set of points, one per cell, clicked on a crop of the
bright-field frame (data/segmentation_reference/README.md). A point marks a
cell's position, not its outline, so this measures detection and
separation, not boundary accuracy:

- a mask holding exactly one marked point is a correctly separated cell;
- a mask holding two or more points has merged cells;
- a mask holding no point is a false mask (debris, noise, the second half
  of a cell cut in two), unless it holds an "unsure" point, in which case
  it is left out of the score;
- a point inside no mask is a missed cell (after looking `tolerance_px`
  around it, since a click lands near the centre, not on it).

recall = separated cells / marked cells
precision = masks with exactly one point / masks scored

Segment the whole frame and pass the crop's window: only masks whose
centroid (over the whole mask, not its visible part) lies inside the crop
are scored as false, because the crop cuts through masks at its edge and
those cells were marked only if their centre was inside it.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi


def read_points(path: str | Path) -> np.ndarray:
    """(row, column) of each point in a Fiji Results table (columns X, Y in pixels)."""
    table = pd.read_csv(path)
    columns = {c.strip().lower(): c for c in table.columns}
    if "x" not in columns or "y" not in columns:
        raise ValueError(f"{path}: needs X and Y columns (Fiji: Set Measurements -> Centroid); "
                         f"found {list(table.columns)}")
    xy = table[[columns["y"], columns["x"]]].to_numpy(dtype=float)
    if not np.isfinite(xy).all():
        raise ValueError(f"{path}: X/Y contain empty or non-numeric values")
    return xy


def _label_at(labels: np.ndarray, point: np.ndarray, tolerance_px: float) -> int:
    """The mask under a point, or the nearest one within tolerance_px; 0 if none."""
    row, col = (int(np.floor(v)) for v in point)
    height, width = labels.shape
    if not (0 <= row < height and 0 <= col < width):
        return 0
    if labels[row, col]:
        return int(labels[row, col])
    reach = int(np.ceil(tolerance_px))
    r0, r1 = max(row - reach, 0), min(row + reach + 1, height)
    c0, c1 = max(col - reach, 0), min(col + reach + 1, width)
    window = labels[r0:r1, c0:c1]
    rows, cols = np.nonzero(window)
    if rows.size == 0:
        return 0
    distance = np.hypot(rows + r0 - row, cols + c0 - col)
    best = int(np.argmin(distance))
    return int(window[rows[best], cols[best]]) if distance[best] <= tolerance_px else 0


@dataclass(frozen=True)
class CropScore:
    marked: int
    unsure: int
    masks_scored: int
    detected: int
    separated: int
    merged_masks: int
    false_masks: int
    recall: float
    precision: float
    f1: float

    def as_dict(self) -> dict:
        return asdict(self)


def score_labels(labels: np.ndarray, points: np.ndarray, unsure: np.ndarray | None = None,
                 tolerance_px: float = 3.0, window: tuple[int, int, int, int] | None = None,
                 ) -> tuple[CropScore, dict]:
    """Score a label mask against points marked on a crop of it.

    labels: the whole frame's label mask. window: (y0, x0, height, width) of
    the crop the points were marked on, in frame pixels; None means the
    whole frame. points / unsure: (row, column) in crop coordinates.

    Returns the score and, per mask id, its status ("separated", "merged",
    "false", "unsure"), plus the "missed" points (crop coordinates).
    """
    if labels.ndim != 2:
        raise ValueError("labels must be 2D")
    y0, x0, height, width = window if window is not None else (0, 0, *labels.shape)
    if y0 < 0 or x0 < 0 or y0 + height > labels.shape[0] or x0 + width > labels.shape[1]:
        raise ValueError(f"window {window} lies outside the {labels.shape} frame")
    offset = np.array([y0, x0], dtype=float)
    points = np.asarray(points, dtype=float).reshape(-1, 2)
    unsure = np.zeros((0, 2)) if unsure is None else np.asarray(unsure, dtype=float).reshape(-1, 2)
    hit = np.array([_label_at(labels, p + offset, tolerance_px) for p in points], dtype=int)
    hit_unsure = {_label_at(labels, p + offset, tolerance_px) for p in unsure} - {0}

    crop = labels[y0:y0 + height, x0:x0 + width]
    ids = np.unique(crop)
    ids = ids[ids != 0]
    centroids = ndi.center_of_mass(np.ones_like(labels), labels, ids) if ids.size else []
    inside = {int(i) for i, (r, c) in zip(ids, centroids)
              if y0 <= r < y0 + height and x0 <= c < x0 + width}
    counts = pd.Series(hit[hit > 0]).value_counts().to_dict()
    scored = (inside | set(counts)) - (hit_unsure - set(counts))

    status = {}
    for i in scored:
        n = counts.get(i, 0)
        status[i] = "separated" if n == 1 else "merged" if n > 1 else "false"
    for i in hit_unsure - set(counts):
        status[i] = "unsure"

    marked = len(points)
    separated = sum(1 for p in hit if p and counts.get(p) == 1)
    masks_scored = len(scored)
    recall = separated / marked if marked else float("nan")
    # a mask with exactly one point and a point alone in its mask are the same pairs
    precision = separated / masks_scored if masks_scored else float("nan")
    f1 = (2 * recall * precision / (recall + precision)) if recall + precision > 0 else 0.0
    score = CropScore(
        marked=marked, unsure=len(unsure), masks_scored=masks_scored,
        detected=int((hit > 0).sum()), separated=separated,
        merged_masks=sum(1 for s in status.values() if s == "merged"),
        false_masks=sum(1 for s in status.values() if s == "false"),
        recall=recall, precision=precision, f1=f1,
    )
    return score, {"masks": status, "missed": points[hit == 0]}
