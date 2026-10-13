# Data contract

What each stage of the real-data route expects and produces
([architecture.md](architecture.md) lists the stages). Columns not listed
here are carried through unchanged. A file that breaks a rule below is
refused with a message; nothing is dropped silently.

## Identifiers

| Field | Form | Example | Meaning |
| --- | --- | --- | --- |
| `session_id` | `<yyyymmdd>_<stressor>_<tr\|nt>` | `20261007_h2o2_tr` | One imaging day of one stressor; `tr` transfected, `nt` non-transfected. `none` as stressor marks medium-only controls. |
| `condition_id` | `<stressor>_<dose>uM_<label>` | `h2o2_300uM_strosi_14:29` | One sample at one time point. The dose is read from here; the stressor must match the session. |
| `sample_id` | free text, unique within a condition | `h2o2_0710_s1429_300uM_R2` | One field (one image). R1–R3 / F1… are fields of the same sample, **not** independent cultures. |
| field key | (`session_id`, `condition_id`, `sample_id`) | | `sample_id` alone repeats across time points of a series. |
| `measurement_time_hours` | ≥ 0 | `3.95` | Hours from stress addition to imaging. The time axis of every analysis. `timepoint` is only an acquisition index. |

Doses are in µM. Stressors: `h2o2` (model variant `ox`), `dtt` (model variant `er`).

## Image manifest (`images.csv`, `images_aligned.csv`)

One row per field. Columns the reader requires (values may be empty where
noted): `session_id`, `timepoint`, `acquisition_order`, `green_path`,
`red_path`, `bright_field_path` (empty if none), `exposure_ms_green`,
`exposure_ms_red`, `nd_filter_green`, `nd_filter_red` (empty if none),
`objective`, `burner_hours`, `lamp_warmup_minutes`. Optional for the
pipeline but required by the analysis: `sample_id`, `condition_id`,
`measurement_time_hours`. Also used when present: `sampling_time_hours`,
`specimen_id`, `biological_replicate_id`, `gain_setting`, `expect_cells`
(false for medium-only frames), `log_flags`.

`align_brightfield.py --write-manifest` adds `brightfield_shift_dy`,
`brightfield_shift_dx` (pixels, bright-field → fluorescence) and fills
`bright_field_path`; it writes a new file and never edits the input.
Paths are resolved relative to the manifest's folder.

Imaging constants of the October data: 40X objective with a 0.5X camera
adapter, 0.1725 µm per pixel; 8-bit frames (saturation at 255);
green and red of one field share one exposure.

## Cell table (`cells_all.csv`, then `cells_all_per_s.csv`)

One row per segmented cell, schema in `ichnos/schema.py`. The columns the
analysis relies on:

| Column | Unit | Note |
| --- | --- | --- |
| `cell_id`, `area_px` | – | Area in pixels of the bright-field mask. |
| `corrected_mean_green`, `corrected_mean_red` | grey levels | Background-subtracted; red corrected for green bleed-through (0.05, indicative). |
| `corrected_mean_green_per_s`, `corrected_mean_red_per_s` | grey levels / s | Added by `normalize_exposure.py`; the intensities every analysis uses. |
| `ratio_red_green` | – | Red/green; exposure-free. |
| `qc_pass`, `edge_flag`, `qc_reasons` | – | A **usable** cell has `qc_pass` true and `edge_flag` false. |
| `mask_source` | – | `brightfield` for the October run. |
| `exposure_ms_*`, `objective`, `acquisition_json` | – | Acquisition settings carried per cell. |

## Analysis outputs (`run_analysis.py`)

Every intensity below is per second and has the autofluorescence baseline
removed: baseline = median of usable non-transfected cells of that stressor
(`green_af`, `red_af`). A **reporter-positive** cell has green and red both
above the 90th percentile of non-transfected cells of its stressor.

| Folder | File | One row per | Key columns |
| --- | --- | --- | --- |
| `timecourse/` | `baseline.csv` | stressor | baseline per channel, n cells |
| | `fields.csv`, `timecourse.csv` | field; transfected curve × time | medians of `green_af`, `red_af`, ratio |
| | `peaks.csv` | curve × channel | peak time, 5–95 % interval, whether it is only a bound |
| `induction/` | `effects.csv` | day × dose × readout (green, red, positive share) | `effect` vs same-day 0 µM, `ci_low/high`, `mde`, `p_two_sided`, `q_bh` |
| | `trend.csv` | day × readout | slope per 10-fold dose, interval, q |
| `model_compare/` | `compare.csv` | day × dose × channel × clearance scenario | measured `fold` (positive cells) and interval, `model_fold`, `model_inside_ci` |
| `decoder_lodo/` | `predictions.csv` | held-out field or sample × decoder | posterior per dose `p_<dose>uM`, `map_dose_uM`, `p_true`, 80 % credible set, `outside_training` |
| | `metrics.csv` | stressor × decoder × unit × held-out day | accuracy vs chance, information gain (bits), entropy vs prior, credible coverage |
| (root) | `run-manifest.json` | run | configuration, input sha256, git commit, output sha256 |

## Rules every stage keeps

- An existing output path is never overwritten; give a new folder.
- Non-transfected sessions feed only the baseline and the positive threshold;
  they are never scored.
- Effects and folds are computed within one day; days are never pooled
  into one control.
- Resampling: cells within fields, fields within time points; seeds are
  recorded so reruns are identical.
