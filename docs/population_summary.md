# Image measurements, sample summaries and session calibration

Implemented on 2026-10-05. The original cell CSV is preserved. A separate
summary command produces `samples.csv` and a provenance `metadata.json`.

## Sample identity

Add `sample_id` and optional `condition_id` to the image manifest. They are
propagated through `ImageSet`, `CellRecord` and the cell CSV. Numeric-looking
identifiers such as `001` remain strings. Existing image manifests still
run with missing IDs; summarization requires an explicit sample ID for every
cell and never invents one from session or acquisition order.

Assign distinct sample IDs to biological replicates. Multiple fields of
view from the same biological sample share its sample ID and have distinct
`acquisition_order` values. A segmented cell's identity is scoped to its
image; cells are not tracked across fields of view or timepoints.

## Grouping and statistics

Group by `(session_id, sample_id, condition_id, timepoint)`. Never pool
different samples or imaging sessions. Sampling/measurement elapsed times,
exposures, ND settings and objective must agree within a group. Unknown
elapsed times remain unknown and are never replaced by `timepoint`.

Compute the median of individual corrected red/green cell ratios. This is
not the ratio of the sample's mean red and mean green intensities. Only
QC-passing cells with finite, non-negative red and ratios and green above
the specified floor are used. The ratio is checked against red/green to
detect wrong direction or an already modified ratio column.

Export total cells, QC-passing cells, usable cells and number of images.
Q25/Q75 describe the distribution across cells; they are not confidence
intervals or biological replicate counts. At least three usable cells are
required by default; otherwise retain the group and counts with missing
summary values. This threshold is a software setting, not an experimentally
established adequacy criterion. The default green floor of zero only
excludes undefined ratios; it is not a measured detection limit.

## Calibration linkage

Without a calibration, report original sample statistics and leave
calibrated values missing. Never silently assume `c_session = 1`.

With a calibration, divide each usable original cell ratio by its session's
`c_session` and summarize the result in separate calibrated columns.
Preserve the original median and quartiles. FRET and model `f` are not
applied again. Store coefficient, reference ID, source kind and source text
in summary rows, and full calibration artifacts in the metadata.

An explicit calibration manifest supplies:

```csv
session_id,reference_id,calibration_path
synthetic-verification,ACTUAL_REFERENCE_ID,calibration.json
```

Use the actual identifier saved in the calibration artifact; the placeholder
above is not executable. The expected reference identifies the forward model
being compared. Paths are relative to the calibration manifest's directory
unless absolute. Each artifact's session ID must match the corresponding
image session. An artifact from a different session must not be renamed to
make it applicable.

If calibrations are supplied, every input session needs one. Each calibrated
session must have one known exposure per channel, ND setting per channel and
objective. Changed settings require a separate acquisition session and
appropriate references. The acquisition settings must also correspond to
the references used for calibration; this is a provenance obligation, not
an experimentally verified fact inferred from a JSON file.

`--data-kind synthetic` accepts synthetic calibration artifacts only.
`--data-kind experimental` accepts experimental-reference artifacts only;
using a synthetic calibration on declared experimental images is rejected.
Neither path automatically establishes instrument or biological validation.

## Commands

First run the existing image pipeline with the extended image manifest:

```bash
python scripts/run_pipeline.py \
  --manifest images.csv --controls controls.csv --out cells.csv
```

Summarize real measurements without claiming a calibration:

```bash
python scripts/summarize_cells.py \
  --cells cells.csv --data-kind experimental \
  --out-dir outputs/sample-summary
```

For synthetic image measurements with matching synthetic session references:

```bash
python scripts/summarize_cells.py \
  --cells synthetic-cells.csv --data-kind synthetic \
  --calibrations calibrations.csv \
  --out-dir outputs/synthetic-calibrated-summary
```

Output directories must be new. No summary, calibrated ratio, or quartile
is a decoder estimate, confidence interval or experimental accuracy claim.

## Verification

The full test suite passed on 2026-10-05: 312 passed, 1 skipped, with
FutureWarning treated as an error. The whitespace check was clean.

An integration test starts from three separated synthetic cell shapes,
passes them through segmentation, corrections, feature extraction and cell
CSV export, then runs the summary entry point with a known synthetic scale.
It verifies a raw ratio median of 0.6 and a calibrated median of 0.4 for
`c_session = 1.5`, preserves the original CSV bytes and rejects switching
the declared data kind to experimental. These simplified synthetic images
check software integration; they do not represent the microscope's noise,
optics or response.

Other tests cover QC exclusion, insufficient/no usable cells, multiple
fields of view, sample/session separation, unknown times, duplicate cells,
acquisition mismatches, wrong ratio direction and incomplete calibration
coverage. The synthetic ratio calibration procedure and its motivation
remain documented in `session_calibration.md`.
