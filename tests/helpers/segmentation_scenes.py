"""Synthetic bright-field scene generation and segmentation scoring helpers.

EXTRACTED verbatim (function bodies unchanged) from DryLabTool
scripts/validate_segmentation_quantitative.py @ 2692ea89, so that the
segmentation-quality regression tests no longer import a script via
sys.path manipulation. The original script remains in the source repo as a
standalone validation report; only the reusable parts live here.

No real expert-annotated ground truth exists for this project yet; these
scenes are a synthetic proxy, not a replacement for annotating real images.
"""
import numpy as np

from ichnos_image import segment  # noqa: F401  (kept for parity with the source module)

def _make_scene(n_cells=25, size=300, seed=0, stressed=False):
    """Synthetic bright-field-like scene with known ground-truth labels.

    baseline: regular circles, slight phase-contrast-like dark rim.
    stressed: irregular (multi-harmonic radius perturbation) boundary, plus
    an internal low-contrast "vacuole" hole -- a stand-in for the
    swollen-vacuole / deformed morphology the team's protocol flagged as a
    stress-response concern.
    """
    rng = np.random.default_rng(seed)
    bf = np.full((size, size), 0.5)
    ground_truth = np.zeros((size, size), dtype=np.int32)

    centers = rng.integers(20, size - 20, size=(n_cells, 2))
    yy, xx = np.mgrid[0:size, 0:size]
    for idx, (cy, cx) in enumerate(centers, start=1):
        base_r = rng.uniform(7, 10)
        theta = np.arctan2(yy - cy, xx - cx)
        r = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)

        if stressed:
            n_harmonics = 3
            perturb = sum(
                rng.uniform(0.08, 0.18) * np.sin(k * theta + rng.uniform(0, 2 * np.pi))
                for k in range(2, 2 + n_harmonics)
            )
            boundary = base_r * (1 + perturb)
        else:
            boundary = base_r

        cell = r < boundary
        if not cell.any():
            continue

        ground_truth[cell & (ground_truth == 0)] = idx
        bf[cell] -= 0.28  # phase-contrast-like darker interior
        bf[cell & (r < boundary * 0.35)] = 0.5  # bright interior halo (real phase-contrast optics, not an artifact)

        if stressed:
            vacuole = r < (base_r * 0.4)
            bf[cell & vacuole] += 0.15  # vacuole: less contrast than cytoplasm

        bf[cell] += rng.normal(0, 0.01, size=cell.sum())  # texture, avoids flat interior

    return bf, ground_truth


def _match_and_score(pred_labels: np.ndarray, gt_labels: np.ndarray, iou_threshold: float = 0.5):
    gt_ids = [i for i in np.unique(gt_labels) if i != 0]
    pred_ids = [i for i in np.unique(pred_labels) if i != 0]

    ious = []
    matched_gt, matched_pred = set(), set()
    pairs = []
    for gt_id in gt_ids:
        gt_mask = gt_labels == gt_id
        overlapping_pred_ids = np.unique(pred_labels[gt_mask])
        for pred_id in overlapping_pred_ids:
            if pred_id == 0:
                continue
            pred_mask = pred_labels == pred_id
            inter = (gt_mask & pred_mask).sum()
            union = (gt_mask | pred_mask).sum()
            iou = inter / union if union else 0.0
            pairs.append((iou, gt_id, pred_id))

    for iou, gt_id, pred_id in sorted(pairs, key=lambda p: -p[0]):
        if gt_id in matched_gt or pred_id in matched_pred:
            continue
        if iou >= iou_threshold:
            matched_gt.add(gt_id)
            matched_pred.add(pred_id)
            ious.append(iou)

    tp = len(matched_gt)
    fp = len(pred_ids) - len(matched_pred)
    fn = len(gt_ids) - len(matched_gt)
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    mean_iou = float(np.mean(ious)) if ious else 0.0
    return dict(n_gt=len(gt_ids), n_pred=len(pred_ids), tp=tp, fp=fp, fn=fn,
                precision=precision, recall=recall, f1=f1, mean_iou=mean_iou)
