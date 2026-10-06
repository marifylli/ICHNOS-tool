# Population summaries to dose decoder

Implemented on 2026-10-06. This connects one calibrated sample trajectory
to the dose decoder. Discrete mode is the default; checked continuous mode
is described in `continuous_dose_decoder.md`. Neither estimates onset time.

## What was implemented, how and why

`ichnos/sample_decoder.py` selects a sample by explicit `session_id` and
`sample_id`, with optional `condition_id`. Multiple conditions require an
explicit choice rather than silently combining them. Both numeric-looking
and textual identifiers are preserved as strings.

The caller must choose `sampling_time_hours` or `measurement_time_hours`.
Selected rows are sorted by that field and must match the dose artifact's
observation times exactly. Missing or duplicate times are rejected.
`timepoint` is an identifier and is never substituted for elapsed hours.
This explicit choice is needed because acquisition and sampling time can
differ, and the appropriate biological readout time depends on the protocol.

The adapter uses `calibrated_ratio_red_green_median` directly. It verifies
that this equals `ratio_red_green_median / c_session` and that the coefficient
and source match the recorded session artifact. It never divides a second
time and never applies model `f` again. It requires adequate usable-cell
counts and finite summaries, retaining the population QC rules.

## Provenance and compatibility

The population summary writer now records `samples_csv_sha256`. The adapter
requires this hash to match the supplied CSV before selecting a sample.
The embedded session-calibration artifact is validated with the same rules
as standalone file loading, including its recorded reference fit.

Supply the original protocol-run `metadata.json` that was used to construct
the calibration reference identifier. Its `model_reference_id` must match
the session artifact. The adapter then compares that reference with the
dose-table records:

- variant and complete parameter profile;
- prepared model SBML hash, including the loaded exposure law;
- actual solver settings;
- exposure assumption, clearance rate and units;
- initialization method, zero-stress duration and convergence criterion;
- ratio direction and model `f`/`eps`.

The dose-table builder now saves the prepared-model hash before preparation
or dose application. This makes it directly comparable with the snapshot
recorded by the protocol runner. The source baseline SBML remains included
in the table as before.

Reference dose is deliberately excluded from dose equality checks: it is
the known reference used to calibrate the session scale, not the unknown
sample dose. Reference fitting times may differ from the sample observation
times because the fitted scale is assumed constant over that session. The
sample times themselves must match the chosen dose table.

The summary data kind must agree with the calibration source: synthetic
calibration is rejected for declared experimental measurements. None of
these software checks proves instrument linearity or biological validity.
Acquisition conditions still must match the real calibration references.

## Existing artifacts

Older `samples.csv`/metadata pairs do not have the summary CSV hash. Regenerate
them with the updated `summarize_cells.py`; do not add a hash manually to make
an unchecked file appear verified.

Older dose tables do not have prepared-model identity records. Rebuild them
with the updated `run_decoder.py build` command. They still work with the
original low-level ratio decoder, but this stricter adapter refuses them.
No existing artifact is silently altered or overwritten.

## Command

Run on one selected sample after generating the new summary and table:

```bash
python scripts/decode_sample.py \
  --samples-csv outputs/summary/samples.csv \
  --summary-metadata outputs/summary/metadata.json \
  --dose-table dose-table.json \
  --calibration-reference outputs/reference-protocol/metadata.json \
  --session-id session-1 --sample-id sample-1 \
  --time-field measurement_time_hours \
  --ratio-tolerance 0.01 \
  --out decoded-sample.json
```

The IDs, time choice and tolerance above are illustrative. Use the correct
ones for the actual input. Add `--condition-id` if needed. Output paths must
be new. Ratios, counts and provenance checks happen before output creation.

## Result

The existing decoder's statuses remain unchanged: unique grid match,
ambiguity, no joint grid match or outside the response envelope. The result
also records the sample/session/condition, chosen time field, timepoint IDs,
CSV hash, calibration reference, coefficient, original ratio medians and
usable-cell counts. It explicitly records that calibration was not applied
again and that the result is not experimentally validated.

The adapter uses the sample summary's upstream green-floor filtering and
cell counts. It does not infer a new experimental detection floor from a
median ratio; the lower-level decoder's optional direct green-floor check
is not automatically claimed to have run.

## Synthetic end-to-end verification

Verification on 2026-10-06: 365 passed, 1 skipped in the full suite,
with FutureWarning treated as an error; the 23 new linkage tests passed.
The whitespace check was clean.

The integration fixture runs a known 25 uM reference protocol and fits a
synthetic session coefficient of 1.5. A separate 75 uM trajectory supplies
ratios for three simple synthetic cells per observation time. Images pass
through segmentation, corrections, extraction, cell CSV export and the
population summary command before reaching the new adapter.

The decoder recovers the 75 uM grid dose while preserving the separate
25 uM calibration reference. The synthetic image intensity scale is arbitrary
and the images do not model real microscope noise or optics. This checks
software integration and scale handling, not biological dose recovery.

Additional tests reject wrong sample/session/condition, wrong time fields,
missing summaries, inadequate counts, changed CSVs, inconsistent coefficient
or calibration identity, mismatched forward models/profiles/solvers/exposure/
initialization, and output overwrite. Documentation of synthetic calibration
and its motivation remains in `session_calibration.md`.

## Unknown exposure onset

Use `--mode joint --elapsed-time-grid-hours ...` for a discrete dose/time search.
In this mode the selected time field supplies measurement spacing only; candidate
elapsed times define exposure onset relative to the first selected measurement.
See [sample joint decoding](sample_joint_decoder_linkage.md) for the procedure,
required table coverage and synthetic verification.
