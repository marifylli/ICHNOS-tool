# Run log

Every step on real data is recorded here, whether it passed or failed:
what was run, what came out, and what was decided. Newest entry last.
Failed runs stay in the log; a later entry records the fix.

Entry format: `## YYYY-MM-DD HH:MM — step — PASS | FAIL | PARTIAL`,
then command, result (numbers), cause (if failed) and next action.

## Open questions for the wet lab

| # | Session | Question | Status |
|---|---|---|---|
| Q1 | DTT 6/10, 0 min | Images of 0 and 10 µM appear swapped between folders `0/0` and `10/0` (camera times 13:16 and 13:06). Confirm. | open |
| Q2 | H2O2 5/10, 0h, 0 µM | Log says red at 14:32; the only image in the folder was taken at ~14:24. Which is right? | open |
| Q3 | H2O2 7/10, "0 min", 0 µM | No image near 11:59. Folder `0/1h` holds an image from ~13:00. What is it? | open |
| Q4 | NT DTT 9/10, 45' | 50 µM R3 missing from the log (image at ~14:00); the duplicated "0 R1" line is probably it. | open |
| Q5 | NT H2O2 8/10 | Log covers 0h, 30', 1h only; folders 2h, 3h, 4h (45 fields) have no times. | open |
| Q6 | H2O2 5/10, 0/2h | Second field (01b, ~16:27) is not in the log. | open |
| Q7 | Medium 8/10 | Two medium images (folder `ΘΡΕΠΤΙΚΌ` and `600/0h/ΘΡΕΠΤΙΚΟ`), one log row. | open |
| Q8 | NT DTT 9/10, 2h | 100/200 µM folders swapped. | **resolved 2026-10-10**: log is authoritative (order 100 → 50 → 200) |
| Q9 | Camera | Please image one slide, same field, at Exp 150, 200, 300, 400 and 600, both channels, plus one frame with the shutter closed at each exposure (~5 min). This is the direct linearity test; the data so far cannot provide it. | open |
| Q10 | Camera | The log lists "Offset 10 / 50". Is that 10 for one channel and 50 for the other (which)? | **answered 2026-10-10**: 10 for preview (pvw), 50 for acquisition (acq), i.e. 50 for every saved image |
| Q11 | Camera | Does the capture software apply per-channel colour gains or white balance to the fluorescence frames? The red background sits at ~50–60 at every exposure, which the acquisition offset alone does not explain. | open |

## 2026-10-10 00:44 — Build acquisition log — PASS

- Command: `python3 scripts/build_acquisition_log.py --pdf data/acquisition/wetlab_log.pdf --out data/acquisition/acquisition_log.csv`
- Result: 297 rows, no warnings. Every row has an exposure and an elapsed time except the 2 medium controls.
- Merged in PR #5.

## 2026-10-10 01:11 — Build image manifest, matcher v1 — FAIL

- Command: `python3 scripts/build_image_manifest.py ... --objective 40X --out-dir data/acquisition/manifest`
- Result: 344 red/green pairs found (A 76, B 92, C 104, D 72). 271 matched, 73 images and 26 log rows unmatched, 60 matches in a folder other than the log's.
- Cause: clock offset estimated at ~47 min for 7–9/10 instead of ~49.4. Each image was paired with its nearest log row; with R1–R3 about 2 min apart this biased the offset and shifted matches by one row (R3 of one dose onto R1 of the next).
- Fix: offset from ordered pairing within folders; match within the same folder first, then across folders. Branch `fix/acquisition-match-offset`.

## 2026-10-10 01:20 — Matcher v2, dry run on the v1 timestamps — PASS (preview)

- Run in the cloud workspace on the camera times from the v1 match report (red channel only).
- Offsets: 49.39 (5/10), 49.38 (6/10), 49.46 (7/10), 49.91 (8/10), 49.25 (9/10) min.
- Result: 294 of 297 log rows matched (292 in their own folder, 2 in another folder: Q1). Residual median 0.43 min, max 4.66 min. 50 images without a log row (Q2–Q7).
- Next: rerun on the device with red and green times.

## 2026-10-10 01:35 — Build image manifest, matcher v2 — PASS

- Command: as above, after `rm -r data/acquisition/manifest`.
- Offsets: identical to the dry run (49.25–49.91 min).
- Result: 294 of 297 log rows matched (292 own folder, 2 other folder: Q1). Residual median 0.44 min, max 4.99 min. 3 log rows without image (Q2, Q3, Q4); 50 images without a log row (Q2–Q7).
- Manifest: 294 rows in 9 sessions. 7 rows have no bright-field frame within 3 min; 2 rows (medium) have no elapsed time, as expected.
- Next: step 1, pipeline smoke test on a few fields.

## 2026-10-10 01:31 — Step 1a: pipeline smoke test, 4 fields — PASS

- Command: `python3 scripts/run_pipeline.py --manifest data/acquisition/manifest/smoke.csv --bleed-default 0.05 --segmentation-method sparse --saturation-value 255 --out outputs/smoke_20261010/cells.csv`
- Input: first 4 manifest rows, H2O2 5/10 "0h" at 20, 75, 300, 600 µM (the 0 µM field is Q2). Measurement times 0.95–1.30 h, exposures read from the log.
- Result: 4/4 fields processed, none refused; 319 objects, 96.2 % pass QC (edge 1.3 %, saturation 2.5 %). Mask source: sum of channels.
- Observation: median background-corrected green is 4.3 grey levels (red 10.1) on the 8-bit camera, against raw levels of ~49 (green) and ~65 (red). Per-cell ratios are therefore noisy at these exposures. To be quantified in step 2 (exposure normalisation and linearity).
- Focus (report mode, not enforced): 99 objects focus_pass, 126 low contrast-to-noise, 94 channels disagree. Consistent with the weak green signal; revisit before enforcing focus QC.
- Caveat: bleed-through 0.05 is illustrative, not calibrated (no single-fluorophore controls yet). Photobleaching not corrected.
- Next: step 1b, all 294 fields.

## 2026-10-10 01:54 — Step 1b: pipeline on all 294 fields — FAIL

- Command: `python3 scripts/run_pipeline.py --manifest data/acquisition/manifest/images.csv --bleed-default 0.05 --segmentation-method sparse --saturation-value 255 --out outputs/full_20261010/cells.csv`
- Result: `zsh: killed` after ~20 min, no output written.
- Cause: out of memory. `run_pipeline.py` loads every image of the manifest before processing the first one. Measured on real images: ~85 MB held per field (green, red and bright-field as float64 = 25 MB each, plus three 3 MB masks), so 294 fields need ~25 GB, more than the laptop has.
- Fix: `scripts/run_pipeline_batched.py` runs the manifest in batches (default 20 fields ≈ 1.7 GB of images), resumes after an interruption and joins the batch CSVs as text. Fields are processed independently (only the bleed coefficient is per session), so the result does not depend on batching. Checked on 4 real fields in 2 batches: merged CSV byte-identical to a single run. 4 tests.
- Follow-up (not done): load images lazily inside the pipeline so a single run needs memory for one field only.

## 2026-10-10 09:43 — Step 1b: batched pipeline run — PARTIAL (slow)

- Command: `python3 scripts/run_pipeline_batched.py --manifest data/acquisition/manifest/images.csv --out-dir outputs/full_20261010 --batch-size 20 -- --bleed-default 0.05 --segmentation-method sparse --saturation-value 255`
- Result so far: batches 0–10 done (220 of 294 fields), no field refused. Batches 0–9: 58,267 objects, QC pass 93–98 % per session. No memory failure.
- Problem: ~1–1.5 h per batch of 20 (batch 9: 07:11→08:09, batch 10: 08:09→09:42), i.e. 3–5 min per field. The full run would take ~10 h.
- Observation for later: non-transfected H2O2 fields give ~480 objects per field against ~160–360 in transfected fields; segmentation may be picking up noise or debris where there is no reporter. Check before using them for autofluorescence.

## 2026-10-10 10:00 — Profile one field — cause of the slowness found

- Profiled `run_pipeline.py` on one real field (DTT 6/10, 200 µM, 311 cells) in the cloud workspace: 88 s in total, 65 s of it in 933 calls to a min/max filter inside `extract_per_cell`.
- Cause: for every cell, the erosion and the two dilations used for the measurement region and the background annulus ran on the whole 2048×1536 frame, so the cost grew as cells × frame size.
- Fix: run those steps in each cell's bounding box padded by `erosion_px + annulus_width_px + 1`. That covers everything they can reach, so the values do not change.
- Checks: one field 88 s → 6.8 s; output CSV byte-identical before and after on 1 and on 4 real fields; new test compares against the full-frame computation, including cells on the image border and neighbouring cells; full suite 612 passed, 1 skipped.
- Next: stop the run, copy in the new `extract.py`, rerun the same command. It resumes at batch 11; batches already done are kept, since the output is identical.

## 2026-10-10 15:10 — Install the extraction fix on the laptop — PASS

- State found: batches 0–10 complete (each with its manifest); batch 11 not finished and no lock left behind, so the run had stopped. Batch 11 will be redone.
- Placed on the laptop: `ichnos_image/extract.py` (bounding-box fix), `tests/image/test_extract_bbox_equivalence.py`, `scripts/run_pipeline_batched.py` (a batch counts as done only when its manifest exists), `tests/image/test_run_pipeline_batched.py` (5 tests). The laptop's previous `extract.py` was confirmed identical to `main` before it was replaced.
- Next: run the two test files, then resume the batched run.

## 2026-10-10 15:17 — Step 1b: full pipeline run, resumed with the extraction fix — PASS

- Tests on the laptop: `test_extract_bbox_equivalence.py` + `test_run_pipeline_batched.py`, 6 passed.
- Command: `caffeinate -i python3 scripts/run_pipeline_batched.py --manifest data/acquisition/manifest/images.csv --out-dir outputs/full_20261010 --batch-size 20 -- --bleed-default 0.05 --segmentation-method sparse --saturation-value 255`
- Batches 0–10 reused; batches 11–14 ran in 70, 70, 72 and 50 s (~3.5 s per field, against 3–5 min per field before the fix).
- Result: 15 batches, 294 fields, 0 refused → `outputs/full_20261010/cells.csv` (90,388 objects, 65 MB). Checked in the cloud workspace: every biological field of the manifest is present (292), no duplicated cell identity, no repeated header from the merge.
- Medium-only controls (5/10 and 8/10): no cell records, as declared; 25 and 8 small fluorescent objects found, both within the foreground ceiling ("unresolved fluorescent objects; not confirmed cells").

| Session | Fields | Objects | Objects/field | QC pass | Saturated | Focus pass | Median corr. green | Median corr. red |
|---|---|---|---|---|---|---|---|---|
| 20261005_h2o2_tr | 25 | 3,945 | 158 | 97.1 % | 1.4 % | 28 % | 4.1 | 11.5 |
| 20261006_dtt_tr | 38 | 11,647 | 306 | 94.3 % | 3.4 % | 48 % | 8.5 | 15.1 |
| 20261006_h2o2_tr | 10 | 1,981 | 198 | 96.2 % | 1.7 % | 31 % | 3.9 | 12.8 |
| 20261007_dtt_tr | 54 | 18,886 | 350 | 96.2 % | 1.9 % | 48 % | 7.7 | 17.4 |
| 20261007_h2o2_tr | 49 | 11,167 | 228 | 97.9 % | 0.5 % | 39 % | 6.4 | 13.7 |
| 20261008_h2o2_nt | 45 | 23,284 | 517 | 97.6 % | 1.2 % | 50 % | 7.4 | 17.0 |
| 20261009_dtt_nt | 71 | 19,478 | 274 | 97.8 % | 0.6 % | 42 % | 5.2 | 11.8 |

- Observations to carry forward (not yet interpreted):
  - Corrected intensities are not yet exposure-normalised, so sessions cannot be compared from this table (step 2).
  - Non-transfected cells are as bright as or brighter than transfected ones in corrected green (7.4 vs 3.9–6.4 for H2O2), at an exposure no longer than theirs (800 ms vs 800–2400 ms). This suggests autofluorescence is a large part of the measured signal; to be confirmed after exposure normalisation (steps 2–3).
  - Non-transfected H2O2 fields give 517 objects per field, about twice the transfected ones; check what segmentation picks up there.
  - Focus pass is 28–50 %; focus QC stays in report mode.
- Housekeeping: `outputs/full_20261010/cells.csv.lock` and `.ichnos-1uftwptc/` are left over from the out-of-memory run at 01:54; safe to delete.
- Next: step 2, exposure normalisation and camera linearity.

## 2026-10-10 15:24 — Merge the matcher and pipeline fixes — PASS

- PR #7 (`fix/acquisition-match-offset`) and PR #8 (`fix/pipeline-memory-and-speed`) merged into `main`.

## 2026-10-10 15:45 — Step 2: exposure normalisation and linearity check — code ready

- New: `ichnos_image/exposure.py`, `scripts/check_exposure_linearity.py` (2a), `scripts/normalize_exposure.py` (2b), `tests/image/test_exposure.py` (6 tests). Full suite in the cloud workspace: 625 passed, 1 skipped.
- 2a reads every green and red frame once and fits background = offset + rate × exposure per session and pooled. It is a necessary condition for linearity only: the background is medium autofluorescence plus dark offset, and no field was imaged at several exposures.
- 2b adds `*_per_s` columns (grey levels per second) to the cells CSV and refuses rows whose two channels have different exposures.
- Next: run 2a and 2b on the laptop.

## 2026-10-10 15:37 — Step 2a: background vs exposure — PARTIAL

- Tests on the laptop: `test_exposure.py`, 6 passed.
- Command: `python3 scripts/check_exposure_linearity.py --manifest data/acquisition/manifest/images.csv --out-dir outputs/linearity_20261010`
- Green background (histogram mode) rises with exposure: pooled fit 13.4 + 18.3 grey/s × t, R² 0.70 over 294 fields; per session R² 0.15–0.58 with intercepts 11–20. Sessions imaged at one exposure fall on the pooled line (non-transfected H2O2 8/10 at 800 ms: 29.4 measured, 28.0 predicted; non-transfected DTT 9/10 at 600 ms: 23.6 vs 24.4; medium 8/10: 26.9 vs 28.0). Exception: medium 5/10 at 1500 ms, 81.9 against 40.8 predicted.
- Red background does not follow exposure: pooled R² 0.06, per-session slopes from −34 to +21 grey/s, level ~50–60 at any exposure. Red background is dominated by a constant (camera offset), and within a session exposure is confounded with slide and time (e.g. 7/10 used longer exposures on the earlier slides).
- Reading: consistent with a linear green channel with a dark offset of ~13; nothing can be concluded about the red channel. The direct test is Q9.
- Update 15:47 (Q10 answered): offset 50 applies to every saved (acq) image, both channels. If that 50 is on the camera's 10-bit scale it is 12.5 on the 8-bit images, which matches the green intercept (13.4). This is a hypothesis, not a measurement. The constant red background (~50–60) is not explained by the offset; one possibility is a colour gain applied to the red plane (Q11).
- Consequence: background subtraction removes the offset in both channels, so per-second normalisation of green is supported by these data; for red it is assumed, not shown.

## 2026-10-10 15:37 — Step 2b: exposure-normalised cells CSV — PASS

- Command: `python3 scripts/normalize_exposure.py --cells outputs/full_20261010/cells.csv --out outputs/full_20261010/cells_per_s.csv`
- Result: 90,388 cells with `*_per_s` columns; every row had equal green and red exposures.
- Median corrected green per second, QC-pass objects (computed in the cloud workspace from the same CSV):

| Session | Objects | Green/s median | Green/s 90th pct | Red/s median | Ratio R/G median | Object area median (px) | ≈ diameter (µm) |
|---|---|---|---|---|---|---|---|
| 20261005_h2o2_tr | 3,832 | 2.8 | 7.2 | 8.1 | 2.81 | 143 | 2.3 |
| 20261006_dtt_tr | 10,982 | 6.5 | 21.4 | 11.9 | 1.78 | 263 | 3.2 |
| 20261006_h2o2_tr | 1,905 | 3.0 | 8.8 | 10.1 | 3.30 | 244 | 3.1 |
| 20261007_dtt_tr | 18,164 | 8.1 | 23.5 | 19.1 | 2.08 | 150 | 2.4 |
| 20261007_h2o2_tr | 10,931 | 7.1 | 18.9 | 15.6 | 2.13 | 117 | 2.1 |
| 20261008_h2o2_nt | 22,735 | 9.2 | 19.7 | 21.2 | 2.24 | 79 | 1.7 |
| 20261009_dtt_nt | 19,051 | 8.6 | 27.3 | 19.5 | 2.16 | 99 | 1.9 |

- Finding that changes the plan: the segmented objects are smaller than yeast cells. At 40X with the 0.5X adapter a pixel is 0.1725 µm, so a 4–5 µm cell covers ~420–660 px; the median objects are 79–263 px (1.7–3.2 µm). Segmentation is probably picking up bright parts of cells or debris, most of all in the non-transfected sessions (79–99 px). Per-object means are therefore not comparable between transfected and non-transfected samples, and the higher non-transfected green above should not be read as "reporter below autofluorescence" yet.
- Next: step 1c, check segmentation against the bright-field before step 3 (autofluorescence).
