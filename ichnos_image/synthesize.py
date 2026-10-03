"""Synthetic dual-channel (GFP+mCherry) data generation, for training/validating
the Stage 4 crosstalk-correction algorithm (and, via true_ratio_red_green, for
sanity-checking Stage 6 ratio recovery) before real wet-lab dual-channel data
exists.

Why synthetic, and why from real single-GFP images: the public reference
dataset has no genuine "GFP-only line imaged through the mCherry filter"
control, so there is no real ground truth for the crosstalk coefficient
alone. A naive per-pixel regression of red vs. green on the real dual-channel
entries (IRE1/NUP49, KAR2/NUP49, HSP104/NUP49, TRX2/NUP49) gave inconsistent,
noisy slopes (0.0-0.21) -- because real GFP and real mCherry (NUP49, nuclear
rim) signal spatially co-occur for some proteins (e.g. IRE1 is ER/nuclear-
envelope localized), so a plain linear fit conflates real red signal with
green leak. That is itself evidence for a segmentation-aware or NMF/SSASU
approach in the real Stage 4 algorithm, not something to fix here. Synthetic
data sidesteps the problem: we choose the ground truth, using a real GFP
image only for realistic cell shapes/background/shot noise.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu


# Moved to image_io.py: the runtime pipeline must not import its loader from
# this synthetic-data module. Re-exported here so existing callers keep
# working while they are migrated.
from .image_io import load_png16, save_png16  # noqa: F401


def _label_cells_from_green(green: np.ndarray, min_size: int = 20) -> np.ndarray:
    """Rough per-cell labeling straight from a GFP channel (Otsu + connected
    components) -- only good enough to give each cell its own synthetic
    red/green ratio below. Use segment.segment_cells on a bright-field/DIC
    channel for real segmentation.
    """
    binary = green > threshold_otsu(green)
    labels, n = ndi.label(binary)
    sizes = ndi.sum(binary, labels, index=np.arange(1, n + 1))
    for label_id, size in enumerate(sizes, start=1):
        if size < min_size:
            labels[labels == label_id] = 0
    return labels


def make_synthetic_pair(
    real_green: np.ndarray,
    *,
    bleed_green_to_red: float,
    true_ratio_red_green: float = 0.0,
    ratio_jitter_frac: float = 0.3,
    noise_std: float = 20.0,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a synthetic (observed_green, observed_red, true_red_signal) triple
    from a real GFP-only image.

    - `real_green` supplies realistic background, cell shapes, and shot noise.
    - `true_red_signal` is a synthetic "real" mCherry signal: each connected
      cell in `real_green` (Otsu-thresholded) gets its OWN ratio, drawn from
      Normal(true_ratio_red_green, ratio_jitter_frac * true_ratio_red_green)
      clipped at 0 -- not one ratio applied uniformly to every pixel.
      `true_ratio_red_green` is therefore the *mean* per-cell ratio, and
      per-cell ground truth varies around it.

      This matters: if every cell used the exact same ratio, true_red would
      be a fixed multiple of green everywhere -- the same shape as the
      crosstalk term (also proportional to green) -- making the two provably
      inseparable from a single image, regardless of estimation method. (An
      earlier version of this function did exactly that, and
      correct.estimate_crosstalk_coefficient() correctly recovered
      bleed + true_ratio instead of bleed alone every time -- see
      scripts/validate_crosstalk.py. Per-cell jitter breaks that
      degeneracy, and reflects biology anyway: cells differ in how long ago
      their stress event started, so a real field never has one ratio.)
    - `observed_red = true_red_signal + bleed_green_to_red * green + noise`,
      the same additive model correct.unmix_crosstalk() assumes, so a correct
      correction implementation should recover close to `true_red_signal`.
    """
    rng = rng or np.random.default_rng()
    green = real_green.astype(np.float64)

    if true_ratio_red_green > 0:
        cell_labels = _label_cells_from_green(green)
        n_cells = int(cell_labels.max())
        per_cell_ratio = np.zeros(n_cells + 1)
        per_cell_ratio[1:] = np.clip(
            rng.normal(true_ratio_red_green, ratio_jitter_frac * true_ratio_red_green, size=n_cells),
            0, None,
        )
        ratio_map = per_cell_ratio[cell_labels]
    else:
        ratio_map = np.zeros_like(green)

    true_red = ratio_map * green
    bleed = bleed_green_to_red * green
    observed_green = np.clip(green + rng.normal(0, noise_std, size=green.shape), 0, 65535)
    observed_red = np.clip(true_red + bleed + rng.normal(0, noise_std, size=green.shape), 0, 65535)
    return observed_green, observed_red, true_red


def build_dataset(
    base_images: dict[str, Path],
    out_dir: str | Path,
    *,
    bleed_values: tuple[float, ...] = (0.02, 0.05, 0.08, 0.12, 0.18),
    true_ratio_values: tuple[float, ...] = (0.0, 0.1, 0.3, 0.6, 1.0),
    ratio_jitter_frac: float = 0.3,
    noise_std: float = 20.0,
    seed: int = 0,
) -> Path:
    """Generate one synthetic sample per (base image x bleed x true_ratio)
    combination, write green/red PNGs, and return the path to a labels.csv
    with the ground truth needed to score a crosstalk-correction algorithm
    (bleed_green_to_red) or a ratio decoder (true_ratio_red_green -- the
    *mean* per-cell ratio for that sample; see make_synthetic_pair for why
    individual cells vary around it rather than sharing one exact value).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    rows = []
    for base_name, green_path in base_images.items():
        real_green = load_png16(green_path)
        for bleed in bleed_values:
            for true_ratio in true_ratio_values:
                green, red, _true_red = make_synthetic_pair(
                    real_green,
                    bleed_green_to_red=bleed,
                    true_ratio_red_green=true_ratio,
                    ratio_jitter_frac=ratio_jitter_frac,
                    noise_std=noise_std,
                    rng=rng,
                )
                sample_id = f"{base_name}_bleed{bleed:.2f}_ratio{true_ratio:.2f}"
                save_png16(green, out_dir / f"{sample_id}_green.png")
                save_png16(red, out_dir / f"{sample_id}_red.png")
                rows.append(
                    dict(
                        sample_id=sample_id,
                        base_image=base_name,
                        bleed_green_to_red=bleed,
                        true_ratio_red_green=true_ratio,
                        noise_std=noise_std,
                        green_file=f"{sample_id}_green.png",
                        red_file=f"{sample_id}_red.png",
                    )
                )

    labels_path = out_dir / "labels.csv"
    with open(labels_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return labels_path
