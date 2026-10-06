# Synthetic end-to-end verification

## Run the complete workflow

From the repository root, with the project environment active:

```bash
python scripts/verify_end_to_end.py --out-dir outputs/e2e-verification
```

Use a new output directory on every run. Existing output is never overwritten.
Exit code 0 means all eight checks passed; a failed check produces exit code 1
and a report to inspect. Errors encountered after directory creation leave the
partial artifacts and the failed execution check for diagnosis.

Read `outputs/e2e-verification/report.md` for a compact expected/actual table,
and `report.json` for structured checks, assumptions and file fingerprints.

## What runs, how and why

This executable verification joins the production model/protocol runner, image
processing functions, cell exporter, population-summary runner and calibrated
sample adapter. It uses idealized, noise-free synthetic images because their
underlying dose, timing and session scaling are known. The objective is to
verify that the stages compose correctly under controlled assumptions. This is
not a fitted experimental dataset, instrument calibration or biological
validation.

1. The protocol CLI runs the oxidative model at reference dose 25 uM, using the
   complete default profile, finite preincubation of 2 hours, and six observation
   times [0.25, 0.5, 0.75, 1, 1.5, 2] hours. It exports model, fluorescence and
   protocol metadata.
2. At each time, synthetic green/red/bright-field arrays contain three separated
   circular cells. Green intensity is 1000 inside cells and zero outside;
   red = 1.5 * model red/green ratio * green. Bright field is 0.2 inside and 0.5
   outside. The images therefore encode model ratios, not a physical camera or
   absolute concentration-to-intensity law. They have no noise, blur, nonzero
   background, clipping or bleed-through. Existing segmentation, correction,
   extraction and CSV export run on these arrays. Arrays are retained as NPZ.
3. Session calibration fits the extracted median image ratios against model
   ratios at alternating reference times (three fit points). The other three
   times are held out. The fitted coefficient should recover the known 1.5
   scale, and calibrated holdout ratios should recover model ratios. Acquisition
   settings are identical across reference and target samples.
4. A fresh model simulation generates target ratios for 75 uM at 0.5 and 1 hour.
   These create target images; cell extraction and the summary runner generate
   calibrated medians with three usable cells per time. The sample adapter then
   checks provenance and performs known-time discrete dose decoding.
5. Joint decoding uses relative spacing [0, 0.5], candidate elapsed times
   [0.5, 1], and dose grid [0, 25, 75]. The table covers 0.5, 1 and 1.5 hours.
   It should recover 75 uM and elapsed time 0.5 hour. A deliberately loose
   absolute tolerance of 10 produces multiple compatible pairs and no point
   estimate; this is an explicit ambiguity-handling test, not an empirically
   measured noise level.
6. Continuous known-time decoding uses a separate dose grid
   [25, 30, 35, 40, 45, 50] uM. Interpolation is checked by fresh simulations at
   quarter, half and three-quarter points per segment, with allowance 0.001 in
   model ratio units. Independently simulated 37.5 uM ratios generate another
   image sample. The recovered dose must be within 0.5 uM of 37.5. This allowance
   and target are fixed computational verification settings, not experimental
   error estimates or certified uniform interpolation bounds.
7. A low-green sample uses intensity 100 with the same ratios. The summary's
   explicitly synthetic image green floor is 500; no cells qualify. Decoding
   must reject the sample as having insufficient usable cells. This verifies
   propagation of filtering into abstention, not a calibrated detector floor.

The summary runner uses `min_cells=3`, original per-cell red/green ratios and
separate calibrated summaries. Session calibration and factor f are not applied
again by the decoders. Cells are not independent biological replicates.

## Retained artifacts

The output directory contains reference model/protocol exports, calibration
JSON and manifest, synthetic NPZ arrays, reference/target cell CSVs, sample
summary CSVs and metadata, dose tables, interpolation checks, each decoding
result, low-green rejection details, and reports. The JSON report includes the
verification script hash and SHA256 hashes of the stage artifacts. Dose tables
and calibration references additionally retain their existing model/profile,
solver, initialization and exposure provenance. Reports do not hash themselves.
No results need to be added to git; the script and documentation are reproducible
project sources.

## Recorded example run

An initial run on 2026-10-06 produced all eight checks as passing:

| Check | Expected | Actual |
| --- | --- | --- |
| Reference extraction | 18 QC-passing cells | 18 |
| Session coefficient | 1.5 | 1.4999999999999996 |
| Holdout maximum ratio error | < 1e-10 | 3.33e-16 |
| Known-time dose | 75 uM | 75 uM |
| Joint pair | 75 uM / 0.5 hour | 75 uM / 0.5 hour |
| Loose-tolerance ambiguity | No point estimate | Ambiguous |
| Continuous dose | 37.5 uM within 0.5 uM | 37.5442912200363 uM |
| Low-green sample | Zero usable cells, reject | Zero usable cells, reject |

Small numerical differences across solver/runtime versions are possible; the
checks use explicit tolerances. The test suite invokes the executable workflow,
checks artifact fingerprints and verifies refusal to overwrite existing output.

## Scope and remaining validation

This run verifies the oxidative, default-profile, finite-preincubation case
under fixed, noise-free assumptions. It does not exhaustively validate ER,
other profiles, equilibrium initialization, arbitrary clearance, model mismatch,
real image artifacts or the entire dose–time domain. It does not establish
biological accuracy, camera linearity, statistical uncertainty or joint
continuous identifiability. Experimental reference images and measurements are
needed for those claims. The broader test suite continues to cover additional
components and input contracts separately.

## Visual report for users

Render an existing verification run without repeating model simulations:

```bash
python scripts/render_e2e_report.py \
  --input-dir outputs/e2e-verification \
  --out-dir outputs/e2e-visual
```

On macOS, open the report with:

```bash
open outputs/e2e-visual/index.html
```

The local HTML report includes the expected/actual verification table and four
PNG panels per image: GFP, mCherry, bright field and segmentation boundaries.
Each image is identified by sample and time; extracted/QC-passing/usable cell
counts and original/calibrated median ratios are shown where available. The
reference images indicate the fit/holdout split. Low-green samples clearly show
that calibrated summaries are unavailable when no cells pass the green floor.

A fixed display range is used per channel across the whole run, with numerical
ranges printed in every caption. Low-green intensity therefore remains visibly
lower than target intensity. Green and red have separate ranges, so channel
brightness must not be used to estimate their numerical ratio; read the printed
ratio summaries. Units are arbitrary synthetic intensity units. These previews
are not calibrated concentrations or original microscope acquisitions. The NPZ
arrays remain the numerical source of truth.

Segmentation is recomputed with the same deterministic default Otsu function
used by this synthetic pipeline. Its cell count must match the cell CSV; the
viewer does not claim to load archived segmentation masks. Input artifacts are
checked against the verification report fingerprints before output creation.
Existing visual directories are not overwritten. `preview_metadata.json`
records source report identity, display scales, image/sample/time mapping and
segmentation convention. This renderer is scoped to the generated synthetic
verification runs, not arbitrary experimental images.

### Model curves

The visual report also exports `time-course.png` and `dose-response.png`.
The time plot shows the reference model's observed green, mature red reporter
and red/green ratio at the six saved times. Calibrated reference image ratios
are overlaid on the ratio plot, distinguishing fit and holdout times. Absolute
model observable values are not calibrated camera intensities.

The dose plots show the stored response tables at each observation time,
alongside calibrated synthetic target ratios at their known generation doses.
These marker dose coordinates are ground truth for this verification, not
estimated dose outputs. Lines connect saved grid points; no new simulations or
continuous time fits are implied. Decoder estimates remain in the expected/actual
results table. Dense model trajectories, uncertain parameters and experimental
observations require separate inputs and validation.

The renderer requires Matplotlib, included in the project's `image` extra.
Existing users can install missing plotting dependencies with
`python -m pip install -e '.[model,image,dev]'`.
