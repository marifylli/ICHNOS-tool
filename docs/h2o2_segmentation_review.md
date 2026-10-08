# H2O2 exploratory image review

## Processing configuration

The trial processed 24 image pairs, excluding the declared medium-only
control from biological extraction.

Configuration: sparse segmentation, source=sum, min_size=30,
green extraction=G, red extraction=R, saturation_value=255,
background_method=mode, focus_mode=report, bleed_default=0.05.

The bleed coefficient is illustrative, not experimentally calibrated.
Timepoint values 0–4 are acquisition labels, not elapsed hours.
Sample names indicate nominal times of 0, 0.5, 1, 2 and 3 hours;
the acquisition record must confirm these before modelling.

## Computational checks

The trial produced 3,853 object records:
3,746 passed QC and 107 failed.
Border flags affected 59 objects and saturation flags affected 50,
with two objects carrying both flags.

The final sparse-label size filter removes labels smaller than min_size
after watershed while preserving retained label IDs and pixels.

The complete test suite passed: 593 passed, 1 skipped.

## Interpretation rules

- A detected object is not automatically a confirmed individual cell.
- Declared cell-free frames produce no biological cell records.
  Their diagnostics are stored under cell_free_control_reviews in the
  processing manifest, including foreground fraction and object count.
  Exceeding the foreground ceiling is reported as a control warning.
- Focus remains in report mode. Restricting summaries to focus_pass
  lowered median Red/Green and raised median corrected green in all
  24 samples; this does not establish which objects are out of focus.
- Manual flags describe uncertainty, not confirmed segmentation errors.
  Baseline summaries retain the original QC selection.
  Sensitivity summaries additionally exclude flagged identities.
- Manual review requires session_id, sample_id, timepoint,
  acquisition_order, cell_id and review_reason.
  Duplicate or unmatched identities are rejected.
- Cells and fields of view are not assumed to be independent biological
  replicates.

## Selected manual-review sensitivity

Nine objects were flagged across the two 300-3h fields.

| Sample | Objects used before/after | Median Red/Green change | Median corrected green change |
|---|---|---|---|
| 300-3h-01a | 150 / 147 | +1.49% | -9.65% |
| 300-3h-01b | 155 / 149 | +2.12% | -2.72% |

The difference between fields remains after these exclusions.
This review does not assess every detected object.

## Reproduction

Run scripts/run_pipeline.py with the configuration above and the
session's image manifest. Include expect_cells=false for declared
cell-free controls.

Run scripts/summarize_manual_review.py with --cells, --review,
--out-dir and --data-kind experimental to publish baseline,
without-flagged, comparison and annotated-cell CSVs plus metadata.

Use the original image-processing manifest and source images to
reproduce extraction. Review annotations apply only to the corresponding
cell identities and processing run; reassess them after resegmentation.

## Remaining validation

Segmentation accuracy, camera linearity, exposure comparability,
bleed-through calibration and biological/model validity remain
experimentally unvalidated. These exploratory summaries are not
validated dose estimates.
