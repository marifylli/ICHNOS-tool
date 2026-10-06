# ICHNOS tool

Research software for fluorescence image processing, ox/er model simulation
and conditional decoding of ratio and green fluorescence observations.

Ratio-only decoders support discrete/checked continuous dose at known times
and discrete dose–elapsed-time compatibility. A separate two-observable
snapshot decoder now uses ratio plus reference-normalized green to retain
compatible dose/time grid pairs. It requires explicit reference calibration
and observation metadata. A conditional two-observable grid posterior is available; experimental coverage and calibration uncertainty remain unvalidated.
A unique grid match is not proof of identifiability.
Synthetic verification checks software behavior under stated assumptions;
experimental recovery accuracy has not been established.

The research repositories (`Ichnos_PULSE`, `ICHNOS-ablation`,
`ichnos-fisher`, `DryLabTool`) retain the supporting studies.
Migrated code and model sources are recorded in
[docs/PROVENANCE.md](docs/PROVENANCE.md).

## Status

| Component | State |
| --- | --- |
| Image pipeline (`ichnos_image`) | Implemented; computational tests pass |
| RGB extraction | Explicit CLI/API selection; experimental mapping unresolved |
| Saturation QC | Checked before RGB extraction and forwarded through the pipeline |
| Shared schema (`ichnos.schema`) | Migrated; schema versioning pending |
| SBML models and builder (`ichnos.build`) | Migrated for ox and er |
| Shared simulator (`ichnos.simulate`) | Implemented with solver settings and observable selections |
| Parameter profiles (`ichnos.params`) | Packaged; value, unit and provenance checks implemented |
| Model exports (`ichnos.io`) | Unique archives with manifests and overwrite protection |
| Experimental exposure protocol | Constant-stress baseline and optional first-order clearance with an explicit rate in h^-1; observation times supplied in hours; computational equilibration or explicit finite preincubation with initialization assumptions recorded |
| Measurement mapping and session calibration | Model fluorescence exports, relative ratio scaling and two-reference green normalization implemented; independent experimental calibration remains open |
| Ratio-only decoders | Discrete dose/time and checked continuous dose implemented |
| Two-observable snapshot decoder | Ratio + reference-normalized green compatibility on a discrete grid; explicit calibration/observation artifacts, no posterior |
| End-to-end and visual reports | Implemented for a fixed noise-free synthetic oxidative scenario |
| Estimator uncertainty | Conditional grid posterior with fitted biological-replicate covariance; calibration/model/covariance uncertainty not propagated |

Passing tests verify software behaviour under controlled assumptions.
They do not establish experimental accuracy for dose or time estimation.


See [the decoder scientific review](docs/decoder_scientific_review.md) for
reproduced counterexamples, calibration diagnostics and the required two-observable design.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[image,model,dev]"
python -m pytest -q
```

The optional Cellpose backend requires the `cellpose` extra.

## Image pipeline CLI

```bash
python scripts/run_pipeline.py --help
```

For RGB acquisitions, explicitly select the extraction method for each
channel using `--green-extraction` and `--red-extraction`.
The same choices are applied to samples and crosstalk controls.

The CLI saturation threshold defaults to the instrument profile and can
be overridden with `--saturation-value`. It is applied to stored image
components before fluorescence extraction.

Explicit extraction choices are configuration options, not evidence
that those choices have been experimentally calibrated.

### Experimental timing in image manifests

`timepoint` is an integer point identifier, not elapsed time.

Image manifests optionally accept:

- `sampling_time_hours`: actual sampling time from stress onset.
- `measurement_time_hours`: actual measurement time from stress onset.

Both fields use hours and accept fractional values, including 0.5 and 0.75.
Missing columns or blank entries remain unknown. Supplied values must be
finite and non-negative.

The two times are preserved separately in the per-cell CSV. No automatic
conversion between sampling and measurement time is performed.

The time used for comparison with the simulator must reflect the sample
handling protocol, including whether biological activity continues between
sampling and imaging.

## Run a stress protocol

Install the model dependencies in the active virtual environment:

```bash
python -m pip install -e ".[model]"
```

Run a finite zero-stress preincubation scenario followed by stress:

```bash
python scripts/run_protocol.py \
  --variant ox \
  --profile default \
  --dose 75 \
  --dose-units uM \
  --times-hours 0.5 1 2 3 \
  --initialization finite-preincubation \
  --preincubation-hours 2 \
  --out-dir outputs/ox_75uM_preincubation_2h
```

The preincubation duration is an explicit model assumption. The model
does not represent glucose-to-galactose switching or a validated
experimental initial state.

For computational equilibration instead, use:

```bash
python scripts/run_protocol.py \
  --variant er \
  --profile default \
  --dose 250 \
  --dose-units uM \
  --times-hours 0.75 2 4 \
  --initialization equilibrium \
  --out-dir outputs/er_250uM_equilibrium
```

Equilibrium mode requires convergence of the checked readouts over the
zero-stress preparation interval. This does not establish experimental
equilibrium or convergence of every model state.

Observation times are supplied in hours from stress onset. Supply the
actual elapsed times appropriate to the measurement, rather than assuming
that nominal sampling times and image acquisition times are identical.

Stress remains constant unless `--clearance-rate-per-hour` is supplied.
That option requires an explicit positive finite rate in h^-1 and adds
assumed first-order decay. No experimentally calibrated default is provided.

Each run creates a new output directory containing:

- `results.csv`: model outputs and `time_hours`.
- `fluorescence.csv`: `time_hours`, `green` (`Observed_Green`), `red`
  (`Reporter_red`, mature mCherry) and `ratio_red_green`
  (`Measured_Ratio_RG`). The channels are model concentrations in nM,
  not camera intensity units. `Total_red_pool` includes immature forms
  and is not the red fluorescence signal.
- `metadata.json`: parameter profile and provenance, exposure,
  initialization, solver settings, runtime versions and hashes, plus
  the fluorescence mapping and the model's `f` and `eps` values.
- `model.sbml`: the loaded model before preparation and dose application.

The recorded profile includes its original parameter values; the exposure
metadata records the dose applied by the protocol.

Existing output paths are rejected. File writes are not a single atomic
transaction; a failed write can leave a partial output directory.

These exports contain model observables, not experimentally calibrated
fluorescence or decoder estimates.

The ratio direction is red/green (mCherry/GFP). The exported model ratio
is `f * Reporter_red / (Observed_Green + eps)`; FRET and `f` are already
included and are not applied again. Session calibration is not yet applied.

## Implemented workflows and remaining scientific work

Checked continuous dose interpolation at known times is documented in
[docs/continuous_dose_decoder.md](docs/continuous_dose_decoder.md). It retains
all compatible dose regions, checks intermediate simulations and supports
the existing sample-summary adapter. These are compatibility regions,
not experimentally validated confidence intervals.

The checked adapter from population summaries to the discrete dose decoder
is documented in [docs/sample_decoder_linkage.md](docs/sample_decoder_linkage.md).
It selects one sample, requires an explicit elapsed-time field and checks
session-calibration and model/protocol provenance before decoding.

The first model-generated discrete dose decoder at known observation times
is documented in [docs/dose_decoder.md](docs/dose_decoder.md). It reports
all compatible grid doses and abstains on ambiguity or incompatibility.
This is computational verification, not experimental dose validation.

Relative session calibration and its synthetic verification procedure are
documented in [docs/session_calibration.md](docs/session_calibration.md).
The synthetic artifacts verify software scale recovery; they are not an
experimental microscope calibration. Image-pipeline sample identifiers,
QC-filtered population summaries and optional session calibration are
documented in [docs/population_summary.md](docs/population_summary.md).

1. Align experimental metadata with the simulator: variant, initial dose,
   actual elapsed times, exposure history and initial model state.
2. Evaluate the exposure assumption and constrain clearance rates using
   experimental data; no experimentally calibrated clearance default exists.
3. Extend interpolation verification and grid refinement to the intended
   operating domain and sharp-response regions.
4. Validate two-reference green calibration independently, test snapshot
   identifiability across the intended domain and add replicate-based likelihoods.
5. Validate the workflow against experimental references for each variant.

Copper/CuSO4 modelling is outside the current implementation scope.
Its images can be processed without a copper-specific model or decoder.

See [docs/limitations.md](docs/limitations.md) for limitations and
unresolved measurement decisions.

### Joint discrete dose and elapsed time

`python scripts/run_joint_decoder.py` checks every dose/elapsed-time pair in a
model-generated grid using already calibrated red/green ratios and known spacing
between measurements. Elapsed time is measured from exposure onset to the first
measurement. It retains ambiguous pairs and does not interpolate or claim
statistical confidence intervals or experimental validation. See
[the procedure and synthetic verification](docs/joint_dose_time_decoder.md).

The sample-summary adapter also supports `--mode joint` with an explicit
`--elapsed-time-grid-hours` grid. It uses recorded times as measurement spacing
and retains the calibration/model provenance checks. See
[sample joint decoder linkage](docs/sample_joint_decoder_linkage.md).

### Complete synthetic verification

Run `python scripts/verify_end_to_end.py --out-dir outputs/e2e-verification`
with a new output directory to exercise reference calibration, synthetic image
processing, population summaries, discrete/continuous dose and joint dose/time
decoding, ambiguity and low-green rejection. Read the generated `report.md` and
`report.json` for expected/actual checks. See
[the complete procedure and limitations](docs/end_to_end_verification.md).

To view channel previews, cell outlines and image-linked summaries, run
`python scripts/render_e2e_report.py --input-dir outputs/e2e-verification --out-dir outputs/e2e-visual`
and open `outputs/e2e-visual/index.html` in a browser. The report uses fixed
channel display scales and explicitly identifies the images as synthetic.
The visual report includes reference model time courses and dose–response grid
plots, with calibrated synthetic image ratios overlaid for comparison.

### Ratio and green snapshot decoding

Run `python scripts/verify_snapshot_decoder.py --out-dir outputs/two-observable-verification`
to verify a late oxidative snapshot with model-dependent synthetic green.
See [the two-observable decoder](docs/two_observable_decoder.md) for reference
contracts, CLI usage and the remaining uncertainty/experimental validation work.

### Conditional snapshot posterior

See [snapshot posterior](docs/snapshot_posterior.md) for correlated replicate-noise fitting, explicit priors, discrete credible sets and a synthetic held-out check. This is conditional inference under a fixed model and calibration, not experimental validation.
