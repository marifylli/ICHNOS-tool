"""Stage 7 QC calibration: pick a defensible focus_score_threshold by tying
segment.focus_score() to something it's actually meant to predict --
segmentation quality -- rather than guessing a number.

Method: take the 3 real DIC images (presumably in reasonably good focus, as
published data), apply increasing synthetic Gaussian blur to simulate
progressive defocus, and at each blur level record focus_score plus how much
segment.segment_cells()'s output has drifted from the sharp-image baseline
segmentation (mask IoU, not just cell count -- an earlier version of this
script used cell count alone and found it a poor proxy: this edge-activity
segmenter turned out to be blur-tolerant enough that count barely moves even
under heavy blur, while the actual mask boundaries visibly degrade). The
threshold is chosen where IoU vs. the sharp baseline starts dropping
noticeably -- i.e. where defocus has actually started changing what gets
measured, not an arbitrary score cutoff.

Usage: python scripts/calibrate_focus_threshold.py
"""
from pathlib import Path

import numpy as np
from skimage.filters import gaussian

from ichnos_image import segment
from ichnos_image.synthesize import load_png16

RAW_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"
DIC_IMAGES = {
    "TRX2_179997": RAW_DIR / "179997_TRX2_dic.png",
    "YAP1_165478": RAW_DIR / "165478_YAP1_dic.png",
    "HAC1_182391": RAW_DIR / "182391_HAC1_dic.png",
}
BLUR_SIGMAS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0]


def mask_iou(mask_a: np.ndarray, mask_b: np.ndarray) -> float:
    a, b = mask_a > 0, mask_b > 0
    union = (a | b).sum()
    return float((a & b).sum() / union) if union else 1.0


def main():
    print(f"{'image':15s} {'sigma':>6s} {'focus_score':>14s} {'n_cells':>8s} {'iou_vs_sharp':>13s}")
    rows = []
    for name, path in DIC_IMAGES.items():
        dic = load_png16(path)
        sharp_labels = segment.segment_cells(dic, method="otsu", min_size=30)
        for sigma in BLUR_SIGMAS:
            blurred = gaussian(dic, sigma=sigma, preserve_range=True) if sigma > 0 else dic
            score = segment.focus_score(blurred)
            labels = segment.segment_cells(blurred, method="otsu", min_size=30)
            n_cells = int(labels.max())
            iou = mask_iou(sharp_labels, labels)
            rows.append(dict(image=name, sigma=sigma, focus_score=score, n_cells=n_cells, iou_vs_sharp=iou))
            print(f"{name:15s} {sigma:6.1f} {score:14.6e} {n_cells:8d} {iou:13.3f}")
        print()

    import pandas as pd

    df = pd.DataFrame(rows)
    out_path = RAW_DIR.parent / "segmentation_test" / "focus_threshold_calibration.csv"
    df.to_csv(out_path, index=False)
    print(f"written to {out_path}")

    # suggest a threshold: the focus_score at the sigma where mask IoU vs. the
    # sharp baseline first drops below 0.8, averaged across images
    suggested = []
    for name in DIC_IMAGES:
        sub = df[df["image"] == name].sort_values("sigma")
        broken = sub[sub["iou_vs_sharp"] < 0.8]
        if len(broken):
            suggested.append(broken.iloc[0]["focus_score"])
    if suggested:
        print(f"\nsuggested focus_score_threshold (mean of per-image break points): {np.mean(suggested):.6e}")
    else:
        print("\nno image dropped below IoU 0.8 vs. its sharp baseline within tested blur range")


if __name__ == "__main__":
    main()
