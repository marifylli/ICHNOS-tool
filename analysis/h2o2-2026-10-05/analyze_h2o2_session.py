"""Exploratory per-cell analysis of the 2026-10-05 H2O2 microscopy session.

Scope: QC and descriptive trends only. Not calibrated, not decoder input.
What it does NOT decide: the exposure history between sampling and imaging
(open question 1). It reports elapsed time from H2O2 addition to image
acquisition; whether that equals exposure time depends on that answer.

Steps
  1. Inventory images, parse dose / time label / field from the folder tree.
  2. Read acquisition timestamps from TIFF tags; estimate the computer-clock
     offset against the handwritten red-channel times; flag disagreements.
  3. Flat field per channel = smoothed median of all images (each divided by
     its own median after subtracting a dark level). Approximate; reported.
  4. Linearity check: cell-free background vs recorded exposure.
  5. Segment on mCherry (flat-corrected, high-pass, MAD threshold).
  6. Per cell: annulus background, flat correction, ratio R/G, green/exposure.
  7. Per-cell focus score and contrast-to-noise in each channel
     (ichnos_image.focus); cells below the focus percentile cutoff are
     dropped. Images are flagged either as focus-not-comparable (a channel
     is within the noise, so its focus score cannot be read) or, when both
     channels carry enough contrast, as channels-disagree-on-focus.
  8. Per-condition medians with bootstrap 95% CI, plots, overlays.

Usage
  python analyze_h2o2_session.py --root "/Users/.../igem h2o2" \
      --notes h2o2_session_notes.csv --out-dir h2o2_analysis
"""
from __future__ import annotations

import argparse
import json
import re
import unicodedata
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage import measure, morphology, segmentation

from ichnos_image import focus

CHANNEL_RE = re.compile(r"\[(\d+[a-z])\]\[Alexa (488|594)\]\.tif$")
TIME_ORDER = ["0", "30", "1h", "2h", "3h", "medium"]
SATURATION = 255


def nfc(text):
    return unicodedata.normalize("NFC", text)


def inventory(root: Path) -> pd.DataFrame:
    rows = {}
    for path in sorted(root.rglob("*.tif")):
        match = CHANNEL_RE.search(nfc(path.name))
        if not match:
            continue
        parts = [nfc(p) for p in path.relative_to(root).parts]
        if parts[0].startswith("θρεπτικ"):
            dose, label = "medium", "medium"
        else:
            dose, label = parts[0], parts[1]
        field, wavelength = match.groups()
        key = (dose, label, field)
        rows.setdefault(key, {"dose_label": dose, "time_label": label, "field": field})
        rows[key]["green_path" if wavelength == "488" else "red_path"] = str(path)
    frame = pd.DataFrame(rows.values())
    missing = frame[frame[["green_path", "red_path"]].isna().any(axis=1)]
    if len(missing):
        raise SystemExit(f"incomplete channel pairs:\n{missing}")
    return frame


def tiff_time(path) -> datetime:
    with tifffile.TiffFile(path) as tif:
        stamp = tif.pages[0].tags["DateTime"].value
    return datetime.strptime(stamp, "%m/%d/%Y %I:%M:%S.%f %p")


def notes_time(value: str, day: datetime) -> datetime | None:
    if not isinstance(value, str) or not value.strip():
        return None
    hour, minute = map(int, value.split(":"))
    return day.replace(hour=hour + 12, minute=minute, second=0, microsecond=0)  # PM


def add_times(frame, notes, addition: str, offset_tolerance_min: float):
    frame = frame.merge(notes, on=["dose_label", "time_label", "field"], how="left",
                        validate="one_to_one", indicator=True)
    frame["in_notes"] = frame.pop("_merge") == "both"
    frame["red_time"] = [tiff_time(p) for p in frame.red_path]
    frame["green_time"] = [tiff_time(p) for p in frame.green_path]
    day = frame.red_time.iloc[0]
    noted = [notes_time(v, day) for v in frame.notes_red_time]
    frame["clock_offset_min"] = [
        (t - n).total_seconds() / 60 if n is not None else np.nan
        for t, n in zip(frame.red_time, noted)
    ]
    offset = float(np.nanmedian(frame.clock_offset_min))
    frame["offset_residual_min"] = frame.clock_offset_min - offset
    frame["notes_time_disagrees"] = frame.offset_residual_min.abs() > offset_tolerance_min
    added = notes_time(addition, day) + timedelta(minutes=offset)
    frame["elapsed_h_at_red"] = [(t - added).total_seconds() / 3600 for t in frame.red_time]
    return frame, offset


def load(path):
    return tifffile.imread(path).astype(np.float64)


def flat_field(paths, dark, sigma):
    stack = []
    for path in paths:
        img = load(path) - dark
        stack.append(img / np.median(img))
    flat = ndi.gaussian_filter(np.median(stack, axis=0), sigma)
    return flat / flat.mean()


def segment_red(red, dark, flat, *, k, min_area, max_area, smooth):
    image = ndi.gaussian_filter(ndi.median_filter((red - dark) / flat, 3), smooth)
    high = image - ndi.gaussian_filter(image, 30)
    noise = 1.4826 * np.median(np.abs(high - np.median(high)))
    mask = high > k * noise
    labels = measure.label(ndi.binary_fill_holes(mask))
    sizes = np.bincount(labels.ravel())
    small = sizes < min_area
    small[0] = True
    labels[small[labels]] = 0
    mask = labels > 0
    labels = measure.label(mask)
    labels = segmentation.clear_border(labels)
    areas = np.bincount(labels.ravel())
    too_large = np.flatnonzero(areas > max_area)
    too_large = too_large[too_large > 0]
    return labels, set(too_large.tolist()), mask


def add_focus(cells, labels, green, red, *, sigma):
    """Per-cell focus in both channels, plus this image's channel agreement."""
    green_scores = focus.score_cells(labels, green, smoothing_sigma_px=sigma)
    red_scores = focus.score_cells(labels, red, smoothing_sigma_px=sigma)
    agreement = focus.channel_agreement(green_scores, red_scores)
    for column, scores in (("focus_green", green_scores), ("focus_red", red_scores)):
        by_id = {s.cell_id: s.focus_score if s.scored else np.nan for s in scores}
        cells[column] = cells.cell_id.map(by_id)
    return agreement


def per_cell(labels, oversize, any_mask, green, red, flat_g, flat_r, *, ring=(3, 9)):
    near = morphology.disk(ring[0])
    far = morphology.disk(ring[1])
    foreground = ndi.binary_dilation(any_mask, near)
    rows = []
    for region in measure.regionprops(labels):
        r0, c0, r1, c1 = region.bbox
        pad = ring[1] + 1
        sl = (slice(max(r0 - pad, 0), r1 + pad), slice(max(c0 - pad, 0), c1 + pad))
        cell = labels[sl] == region.label
        annulus = ndi.binary_dilation(cell, far) & ~foreground[sl]
        if annulus.sum() < 20:
            continue
        g, r = green[sl], red[sl]
        bg_g, bg_r = np.median(g[annulus]), np.median(r[annulus])
        noise_g = 1.4826 * np.median(np.abs(g[annulus] - bg_g))
        sig_g = (g[cell].mean() - bg_g) / flat_g[sl][cell].mean()
        sig_r = (r[cell].mean() - bg_r) / flat_r[sl][cell].mean()
        rows.append(dict(
            cell_id=region.label, area_px=region.area,
            centroid_r=region.centroid[0], centroid_c=region.centroid[1],
            green_bg=bg_g, red_bg=bg_r, green_signal=sig_g, red_signal=sig_r,
            green_snr=(g[cell].mean() - bg_g) / max(noise_g, 1e-9) * np.sqrt(region.area),
            saturated=bool((g[cell] >= SATURATION).any() or (r[cell] >= SATURATION).any()),
            oversize_cluster=region.label in oversize,
        ))
    return pd.DataFrame(rows)


def central_background(img, mask, frac=0.25):
    h, w = img.shape
    sl = (slice(int(h * frac), int(h * (1 - frac))), slice(int(w * frac), int(w * (1 - frac))))
    free = ~ndi.binary_dilation(mask, morphology.disk(9))[sl]
    return float(np.median(img[sl][free]))


def bootstrap_ci(values, rng, n=2000):
    values = np.asarray(values)
    if len(values) < 5:
        return np.nan, np.nan
    meds = np.median(rng.choice(values, (n, len(values))), axis=1)
    return tuple(np.percentile(meds, [2.5, 97.5]))


def overlay(red, labels, kept_ids, path, dark, flat):
    disp = (red - dark) / flat
    lo, hi = np.percentile(disp, [1, 99.8])
    disp = np.clip((disp - lo) / (hi - lo), 0, 1)
    rgb = np.dstack([disp] * 3)
    bounds = segmentation.find_boundaries(np.isin(labels, list(kept_ids)) * labels)
    other = segmentation.find_boundaries(labels) & ~bounds
    rgb[bounds] = [0, 1, 0]
    rgb[other] = [1, 0.3, 0]
    plt.imsave(path, rgb[::2, ::2])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--notes", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--addition-time", default="1:41", help="H2O2 addition, notes clock (PM)")
    ap.add_argument("--dark-level", type=float, default=None,
                    help="camera offset in grey levels; default = min of per-image 0.1st percentiles")
    ap.add_argument("--flat-sigma", type=float, default=60)
    ap.add_argument("--k-mad", type=float, default=3.0)
    ap.add_argument("--smooth-sigma", type=float, default=2.5)
    ap.add_argument("--min-area", type=int, default=60)
    ap.add_argument("--max-area", type=int, default=1500)
    ap.add_argument("--min-green-snr", type=float, default=5.0)
    ap.add_argument("--focus-sigma", type=float, default=1.5,
                    help="smoothing before the focus gradient; raise for noisier sensors")
    ap.add_argument("--focus-drop-percentile", type=float, default=25.0,
                    help="drop cells below this percentile of focus_red across the session; 0 keeps all")
    ap.add_argument("--offset-tolerance-min", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=20261005)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=False)
    (args.out_dir / "overlays").mkdir()
    rng = np.random.default_rng(args.seed)

    notes = pd.read_csv(args.notes, dtype=str)
    notes["exposure_acq_ms"] = pd.to_numeric(notes.exposure_acq_ms)
    images, offset = add_times(inventory(args.root), notes, args.addition_time, args.offset_tolerance_min)

    paths_g, paths_r = list(images.green_path), list(images.red_path)
    if args.dark_level is None:
        dark = min(np.percentile(load(p), 0.1) for p in paths_g + paths_r)
    else:
        dark = args.dark_level
    flat_g = flat_field(paths_g, dark, args.flat_sigma)
    flat_r = flat_field(paths_r, dark, args.flat_sigma)
    np.save(args.out_dir / "flat_green.npy", flat_g)
    np.save(args.out_dir / "flat_red.npy", flat_r)

    all_cells, bg_rows = [], []
    for idx, row in images.iterrows():
        green, red = load(row.green_path), load(row.red_path)
        labels, oversize, mask = segment_red(red, dark, flat_r, k=args.k_mad,
                                             min_area=args.min_area, max_area=args.max_area,
                                             smooth=args.smooth_sigma)
        cells = per_cell(labels, oversize, mask, green, red, flat_g, flat_r)
        bg_rows.append(dict(image=idx, green_bg_center=central_background(green, mask),
                            red_bg_center=central_background(red, mask)))
        agreement = None
        if len(cells):
            agreement = add_focus(cells, labels, green, red, sigma=args.focus_sigma)
            cells["usable"] = (~cells.saturated & ~cells.oversize_cluster
                               & (cells.green_snr >= args.min_green_snr)
                               & (cells.green_signal > 0) & (cells.red_signal > 0))
            for key in ("dose_label", "time_label", "field", "exposure_acq_ms", "elapsed_h_at_red"):
                cells[key] = row[key]
            cells["image"] = idx
            all_cells.append(cells)
            kept = set(cells.loc[cells.usable, "cell_id"])
        else:
            kept = set()
        bg_rows[-1].update(
            focus_green=getattr(agreement, "median_green", np.nan),
            focus_red=getattr(agreement, "median_red", np.nan),
            focus_agreement=getattr(agreement, "agreement", np.nan),
            cnr_green=getattr(agreement, "median_cnr_green", np.nan),
            cnr_red=getattr(agreement, "median_cnr_red", np.nan),
            focus_comparable=getattr(agreement, "comparable", False),
            focus_agrees=getattr(agreement, "agrees", False),
            n_focus_cells=getattr(agreement, "n_cells", 0),
        )
        name = f"{row.dose_label}_{row.time_label}_{row.field}.png"
        overlay(red, labels, kept, args.out_dir / "overlays" / name, dark, flat_r)
        focus_note = ""
        if agreement is not None and not agreement.comparable:
            focus_note = (f"  [?] focus not comparable: contrast-to-noise green "
                          f"{agreement.median_cnr_green:.2f} vs red {agreement.median_cnr_red:.2f}")
        elif agreement is not None and not agreement.agrees:
            focus_note = (f"  [!] channels disagree on focus: green {agreement.median_green:.2f} "
                          f"vs red {agreement.median_red:.2f}")
        print(f"{row.dose_label:>6} {row.time_label:>6} {row.field}: "
              f"{len(cells)} objects, {len(kept)} usable{focus_note}")

    images = images.join(pd.DataFrame(bg_rows).set_index("image"))
    cells = pd.concat(all_cells, ignore_index=True)
    cells["ratio_rg"] = cells.red_signal / cells.green_signal
    cells["green_per_s"] = cells.green_signal / (cells.exposure_acq_ms / 1000)

    # One session-wide focus cutoff, not one per image: a per-image cutoff
    # would keep the sharpest quarter of a badly focused image and discard
    # good cells from a sharp one, hiding exactly the variation we are
    # looking for.
    focus_cutoff = None
    if args.focus_drop_percentile > 0:
        scored = cells.focus_red.dropna()
        if len(scored):
            focus_cutoff = float(np.percentile(scored, args.focus_drop_percentile))
            cells["in_focus"] = cells.focus_red >= focus_cutoff
            cells["usable"] &= cells.in_focus.fillna(False)
    if "in_focus" not in cells:
        cells["in_focus"] = True

    summary = []
    for idx, group in cells[cells.usable].groupby("image"):
        lo_r, hi_r = bootstrap_ci(group.ratio_rg, rng)
        good = group.green_per_s.dropna()
        lo_g, hi_g = bootstrap_ci(good, rng) if len(good) else (np.nan, np.nan)
        summary.append(dict(image=idx, n_usable=len(group),
                            ratio_median=group.ratio_rg.median(), ratio_ci_low=lo_r, ratio_ci_high=hi_r,
                            green_per_s_median=good.median() if len(good) else np.nan,
                            green_per_s_ci_low=lo_g, green_per_s_ci_high=hi_g))
    summary = images.drop(columns=["green_path", "red_path"]).join(
        pd.DataFrame(summary).set_index("image"))
    summary.to_csv(args.out_dir / "conditions.csv", index_label="image")
    cells.to_csv(args.out_dir / "cells.csv", index=False)

    # Linearity: background at image centre vs recorded exposure
    lin = summary.dropna(subset=["exposure_acq_ms"])
    fit = {}
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
    for a, col, colour in [(ax[0], "green_bg_center", "g"), (ax[1], "red_bg_center", "r")]:
        x, y = lin.exposure_acq_ms.to_numpy(), lin[col].to_numpy()
        if len(np.unique(x)) >= 2:
            slope, intercept = np.polyfit(x, y, 1)
            r2 = 1 - np.sum((y - (slope * x + intercept)) ** 2) / np.sum((y - y.mean()) ** 2)
            fit[col] = dict(slope_per_ms=slope, intercept=intercept, r2=r2, n=len(x))
            xs = np.linspace(0, x.max() * 1.05, 50)
            a.plot(xs, slope * xs + intercept, colour + "--", lw=1)
        a.scatter(x, y, c=colour, s=18)
        a.axhline(dark, color="k", lw=0.8, ls=":")
        a.set(xlabel="recorded exposure (ms)", ylabel="cell-free background (grey)", title=col)
    fig.tight_layout()
    fig.savefig(args.out_dir / "linearity.png", dpi=150)

    # Trends vs elapsed time
    data = summary[summary.dose_label != "medium"].dropna(subset=["ratio_median"])
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for dose, group in data.groupby("dose_label", sort=False):
        group = group.sort_values("elapsed_h_at_red")
        ax[0].errorbar(group.elapsed_h_at_red, group.ratio_median,
                       yerr=[group.ratio_median - group.ratio_ci_low, group.ratio_ci_high - group.ratio_median],
                       marker="o", capsize=3, label=f"{dose}")
        ax[1].errorbar(group.elapsed_h_at_red, group.green_per_s_median,
                       yerr=[group.green_per_s_median - group.green_per_s_ci_low,
                             group.green_per_s_ci_high - group.green_per_s_median],
                       marker="o", capsize=3, label=f"{dose}")
    ax[0].set(xlabel="hours from H2O2 addition to image", ylabel="median per-cell red/green")
    ax[1].set(xlabel="hours from H2O2 addition to image", ylabel="median per-cell green / exposure (grey/s)")
    ax[0].legend(title="dose label", fontsize=8)
    fig.suptitle("Exploratory, uncalibrated; one field per condition; error bars = cell bootstrap (not biological)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out_dir / "trends.png", dpi=150)

    report = dict(
        scope="exploratory descriptive QC; not calibrated; not decoder input",
        n_images=len(images), dark_level=float(dark), clock_offset_min=offset,
        notes_time_disagreements=summary.loc[summary.notes_time_disagrees,
            ["dose_label", "time_label", "field", "notes_red_time", "offset_residual_min"]].to_dict("records"),
        images_not_in_notes=summary.loc[~summary.in_notes, ["dose_label", "time_label", "field"]].to_dict("records"),
        linearity_fit=fit, parameters={k: str(v) for k, v in vars(args).items()},
        focus=dict(
            cutoff_on_focus_red=focus_cutoff,
            cells_dropped_out_of_focus=int((~cells.in_focus.fillna(True)).sum()),
            images_where_focus_is_not_comparable=summary.loc[
                ~summary.focus_comparable, ["dose_label", "time_label", "field",
                                            "cnr_green", "cnr_red"]
            ].to_dict("records"),
            images_with_channel_focus_disagreement=summary.loc[
                summary.focus_comparable & ~summary.focus_agrees,
                ["dose_label", "time_label", "field",
                 "focus_green", "focus_red", "focus_agreement"]
            ].to_dict("records"),
        ),
        open_questions=["exposure history between sampling and imaging (Q1)",
                        "dose units / final concentration (Q2)",
                        "one flask per dose? (Q3)"],
        caveats=["error bars resample cells within one field, not biological replicates",
                 "a ratio from an image whose channels disagree on focus is a focus artefact",
                 "where focus is not comparable, a low green focus score means low contrast, "
                 "not necessarily defocus",
                 "flat field estimated from the data with an approximate dark level",
                 "mCherry-based segmentation; clusters above max_area excluded"],
    )
    (args.out_dir / "report.json").write_text(json.dumps(report, indent=2, default=float))
    print(json.dumps({k: report[k] for k in ("dark_level", "clock_offset_min",
                                             "notes_time_disagreements", "images_not_in_notes",
                                             "linearity_fit", "focus")}, indent=2, default=float))


if __name__ == "__main__":
    main()
