# Integrated image-to-posterior workflow

## Scope and compatibility

The new path connects raw images, per-cell QC, paired ratio/green population
medians and the existing snapshot/posterior APIs. The older ratio-only CLIs
remain available. Calibration, fitted covariance and model parameters are still
held fixed; no experimental accuracy or full uncertainty propagation is claimed.

### 1. Registered saturation and historical impact

The red-channel clipping mask now follows the exact linear interpolation and
nearest-edge rule used to register red to green. Every nonzero interpolation
contribution from a clipped pixel is flagged; nearest-neighbour interpolation
of a binary mask would miss fractional-pixel contamination. Green's mask stays
in the segmentation frame. Both masks are computed before RGB extraction.

API callers should pass `green_saturation_mask` and `red_saturation_mask`.
A legacy combined `raw_saturation_mask` alone is conservatively treated as
potentially belonging to either channel. It can overflag; use separate masks
for an exact audit. Raw image means remain in their original acquisition frame.

Each new cell row includes `sat_flag_legacy` and `qc_pass_legacy_saturation`.
The latter holds all *other current* QC settings fixed and changes only the
saturation rule. To evaluate actual historical images, rerun the same manifest,
parameters and segmentation configuration into a new output, then run:

```bash
python scripts/run_pipeline.py --manifest images.csv --controls controls.csv \
  --green-extraction G --red-extraction R --out results/reprocessed.csv
python scripts/audit_saturation_registration.py \
  --cells results/reprocessed.csv --out results/saturation-impact.json
```

The audit reports changed clipping flags, QC cell counts and median ratios per
sample/timepoint. It is a saturation-only counterfactual under current processing,
not a claim of byte-for-byte reproduction of a historical environment. It does
not recalculate old posteriors. Original spatial masks cannot be recovered from
old per-cell CSVs; raw images are necessary. No historical experimental TIFFs
were available in the implementation checkout, so their numerical impact is
**not yet determined**. Synthetic tests cover integer, fractional, zero and
opposite-direction shifts, including both newly and no-longer flagged cells.

### 2. Safe output and provenance

`process_experiment` refuses existing CSVs or sidecar manifests. It writes in a
sibling temporary directory and publishes a CSV plus `cells.csv.manifest.json`
with no-clobber links. Ordinary failures clean up unpublished work. The manifest
is the completion marker: a process kill between the two publications can leave
an incomplete CSV without a manifest, which must not be accepted as a full run.
Concurrent cooperative writers are excluded by an exclusive `.lock` file; a
crash can leave a lock that must be inspected before manual removal.

The manifest records input file and array hashes, effective settings, QC
thresholds, git commit/status, source-file hashes and package versions. Hashes
detect changes; they do not authenticate biological provenance. API-supplied
opaque objects (e.g. a loaded Cellpose model) are marked nonserializable; retain
external model weights/configuration separately for exact reproduction.

The integrated workflow stages its whole output directory and publishes it only
when complete. Existing output directories are rejected. Input files are checked
for changes during processing. Publication assumes a local filesystem with link
and rename support; output bundles are not multi-file database transactions.

### 3. Focus and contrast-to-noise

The per-cell focus implementation is integrated from `feat/per-cell-focus-qc`
through commit `5cdaa83`. It scores green and registered red **before** background
clipping/unmixing, preserving the noise floor. Scores, CNR, channel agreement,
status and mode are exported. Low CNR, unknown CNR, unscorable cells, low focus
scores and channel disagreement have distinct statuses. Low CNR is not labeled
as proven defocus.

Legacy/API default `FocusPolicy(mode="report")` measures but does not reject.
Explicit `mode="enforce"` combines these checks with saturation, border, field
focus, registration and lamp QC in `qc_pass`; it requires `min_score` and rejects
indeterminate measurements. Scores are not corrections to intensities. Thresholds
are acquisition-dependent and need local validation. A noiseless/quantized ring
can have unknown CNR and is rejected in enforce mode, rather than silently passed.

```bash
python scripts/run_pipeline.py --manifest images.csv --controls controls.csv \
  --green-extraction G --red-extraction R --focus-mode enforce \
  --focus-min-score YOUR_VALIDATED_SCORE --focus-min-cnr YOUR_VALIDATED_CNR \
  --out results/cells.csv
```

### 4. Paired summaries and reference contract

`summary_method="paired_cell_medians_v1"` means the median of individual cell
ratios and the median corrected mean green intensity, computed from the **same**
QC-passing, finite, above-floor cells. It is not the ratio of population medians.
Quartiles describe the cell distribution; cells and fields are not independent
biological replicates. Samples below `min_cells` have missing summaries and an
explicit `insufficient_cells` status.

The adapter requires explicit specimen/culture IDs and consistent acquisition
metadata. Image manifests support `specimen_id`, `biological_replicate_id` and
`gain_setting` in addition to the existing fields. Acquisition JSON retains
objective, exposure, ND settings, extraction methods, gain and green units.
Paths in image/control manifests resolve relative to the manifest file.

Build snapshot calibration with `summary_method="paired_cell_medians_v1"` from
references measured with that same estimator, image processing, QC and acquisition
settings. Every noise-training row must also declare this method. Older artifacts
without a declared estimator remain valid for the old direct APIs, but the new
adapter rejects them: rebuild them from appropriate references instead of simply
relabeling incompatible measurements. In particular, an ideal model-generated
pixel intensity is not automatically the post-segmentation mean image intensity.

### 5. Single CLI

Copy `examples/workflow.json` and `examples/images.csv` into your experiment
folder and provide real input paths/artifacts. The example is a **template**:
`min_score=null` deliberately requires a measured threshold. Its extraction,
exposure, thresholds and sample count are illustrative, not recommendations for
the real microscope. Paths in the JSON resolve relative to that JSON.

```bash
python scripts/run_workflow.py --config experiment/workflow.json \
  --out-dir results/experiment-01
```

For posterior mode supply `table`, `calibration`, `noise`, `prior` and an explicit
`max_mahalanobis_squared`. The population green floor must match noise fitting.
For compatibility mode omit `noise`/`prior` and supply `ratio_tolerance`,
`green_tolerance` and `green_floor`. One configuration uses one session-bound
calibration/noise artifact; run different calibrated sessions separately.

Experimental mode requires controls and explicit enforced focus thresholds.
Synthetic mode may explicitly use `bleed_default` and report-only focus.
Calibration independence and matching acquisition settings remain checked by the
snapshot/posterior APIs. The CLI consumes fitted calibration/noise artifacts; it
does not fit references from unknown targets or invent missing acquisition data.

Outputs:

- `cells.csv` and its processing manifest; `samples.csv` with paired medians.
- `saturation-impact.json`, comparing old/new saturation eligibility.
- One `decoded-XXXX.json` per summary, containing the posterior or rejection.
- `qc-report.json` and `report.html`, with counts, MAP ties, credible-set sizes,
  links to full posterior matrices and rejection reasons.
- Copies of the forward/calibration/noise/prior artifacts and `run-manifest.json`
  with hashes of all outputs and external inputs.

A below-QC sample is reported, not silently dropped; a field with no segmented
cells gets `no_segmented_cells`. A poor-fit/prior-only result remains explicitly
marked and must not be presented as a reliable estimate. Successful CLI exit
means the report was produced, not that every sample was accepted. Structural
configuration/input failures abort publication rather than leave a partial run.

### 6. Robustness evaluation

```bash
python scripts/evaluate_robustness.py --out-dir results/robustness \
  --variants ox er --trials 12 --noise-scales 0.5 1 2 --refine
```

Default decoding grids are 9 doses (0–800 µM) by 17 times (0–8 h). `--refine`
adds nested 17-by-33 grids with the same truths and noisy target draws. Each
variant uses on-grid targets plus off-grid targets simulated independently from
the ODE, never interpolated from the decoding table. Within-condition covariance
is fitted from separate synthetic culture labels at each noise level.

Scenarios apply ±10% ratio gain, ±0.1 normalized-green offset or ±10% green span
changes to targets while keeping decoder calibration and training fixed. The
same latent target noise is used across scenarios. These are controlled
miscalibration checks, not propagation of calibration uncertainty.

`trials.csv`, `evaluation.json`, all observations, forward truth tables and noise
artifacts record error ranges over MAP ties, fit failures, credible area and
inclusion. Off-grid truths are not discrete nodes: their reported metric is
inclusion of the nearest-node midpoint cell. It must **not** be called nominal
continuous 95% coverage. Rejections remain in the trial denominator. Small
sample counts are diagnostic; increase trials for stable estimates.

The earlier dense clearance-misspecification benchmark is also available in
`scripts/verify_dense_posterior.py` (integrated from `d7c2fdd`); see
[dense posterior verification](dense_posterior_verification.md).

Experimental holdout validation, empirical QC thresholds, model/calibration
uncertainty propagation and hierarchical longitudinal noise remain open.
