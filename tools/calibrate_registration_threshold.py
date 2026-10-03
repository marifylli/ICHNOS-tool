"""Stage 7 QC calibration: pick a defensible registration_shift_threshold_px
by tying a known green/red misalignment to the ratio error it actually
causes, rather than guessing a pixel count.

Method: build a synthetic dual-channel pair with a known true ratio and NO
misalignment, then apply a range of known sub-pixel/pixel shifts to the red
channel only (simulating imperfect channel co-registration) *before* any
registration correction. At each shift magnitude, measure the resulting
per-cell ratio's relative error vs. the true ratio. The threshold is the
shift magnitude at which relative ratio error first exceeds 10%.

Usage: python scripts/calibrate_registration_threshold.py
"""
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from scipy.ndimage import shift as ndi_shift

from ichnos_image.synthesize import load_png16

RAW_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"
BASE_IMAGES = {
    "TRX2_179997": RAW_DIR / "179997_TRX2_green.png",
    "YAP1_165478": RAW_DIR / "165478_YAP1_green.png",
}
SHIFTS_PX = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0]
TRUE_RATIO = 0.4
ERROR_TOLERANCE = 0.10


def main():
    print(f"{'image':15s} {'shift_px':>9s} {'mean_rel_error':>15s}")
    rows = []
    for name, path in BASE_IMAGES.items():
        green = load_png16(path)
        true_red = TRUE_RATIO * green
        labels, _ = ndi.label(green > np.percentile(green, 95))

        for shift_px in SHIFTS_PX:
            shifted_red = ndi_shift(true_red, shift=(shift_px, 0.0), order=1, mode="nearest")

            rel_errors = []
            for cell_id in range(1, int(labels.max()) + 1):
                mask = labels == cell_id
                if mask.sum() < 10:
                    continue
                g = green[mask].mean()
                r = shifted_red[mask].mean()
                if g <= 0:
                    continue
                measured_ratio = r / g
                rel_errors.append(abs(measured_ratio - TRUE_RATIO) / TRUE_RATIO)

            mean_rel_error = float(np.mean(rel_errors)) if rel_errors else float("nan")
            rows.append(dict(image=name, shift_px=shift_px, mean_rel_error=mean_rel_error))
            print(f"{name:15s} {shift_px:9.2f} {mean_rel_error:15.4f}")
        print()

    import pandas as pd

    df = pd.DataFrame(rows)
    out_path = RAW_DIR.parent / "segmentation_test" / "registration_threshold_calibration.csv"
    df.to_csv(out_path, index=False)
    print(f"written to {out_path}")

    suggested = []
    for name in BASE_IMAGES:
        sub = df[df["image"] == name].sort_values("shift_px")
        broken = sub[sub["mean_rel_error"] > ERROR_TOLERANCE]
        if len(broken):
            suggested.append(broken.iloc[0]["shift_px"])
    if suggested:
        print(f"\nsuggested registration_shift_threshold_px (min break point, conservative): {min(suggested):.2f}")
    else:
        print(f"\nratio error never exceeded {ERROR_TOLERANCE:.0%} within tested shift range")


if __name__ == "__main__":
    main()
