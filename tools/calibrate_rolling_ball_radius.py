"""Empirically verify correct.suggest_rolling_ball_radius() on real images,
per objective: sweep rolling_ball radius around the formula-suggested value
and confirm retained cell signal has plateaued there (i.e. the ball is no
longer dipping into cells and eating real signal at that radius) rather than
trusting the formula alone.

Uses the real 100X (179997/TRX2) and 60X (210071/IRE1) green-channel images
already in data/raw/ -- a different microscope than the team's Olympus, so
these confirm the *mechanism* works, not real Olympus radii. Re-run against
real Olympus images once available.

Usage: python scripts/calibrate_rolling_ball_radius.py
"""
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

from ichnos_image import correct
from ichnos_image.synthesize import load_png16

RAW_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"

IMAGES = [
    ("179997_TRX2_green.png", "100X", 0.13000),
    ("210071_IRE1_green.png", "60X", 0.10760),
]

# radii to sweep, spanning well below to well above the formula suggestion
RADII_FACTORS = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]


def main():
    for filename, objective, pixel_size_um in IMAGES:
        image = load_png16(RAW_DIR / filename)
        suggested = correct.suggest_rolling_ball_radius(pixel_size_um)

        # cell mask: simple Otsu-style bright-region threshold, just to locate
        # cell peaks for the "is signal being eaten" check below
        threshold = np.percentile(image, 95)
        cell_mask = image > threshold
        labels, n_cells = ndi.label(cell_mask)

        print(f"\n{filename}  objective={objective}  pixel_size={pixel_size_um}um/px  "
              f"suggested_radius={suggested:.1f}px  n_cells~{n_cells}")
        print(f"  {'radius':>8s} {'factor':>7s} {'mean_cell_signal':>17s}")

        for factor in RADII_FACTORS:
            radius = suggested * factor
            corrected, _ = correct.subtract_background(image, method="rolling_ball", rolling_ball_radius=radius)
            mean_cell_signal = float(corrected[cell_mask].mean())
            print(f"  {radius:8.1f} {factor:7.2f} {mean_cell_signal:17.1f}")


if __name__ == "__main__":
    main()
