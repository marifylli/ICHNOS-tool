# Session calibration: method and verification record

Implementation date: 2026-10-05.

## Purpose and scope

This step implements a relative red/green scale per imaging session. It
does not calibrate absolute GFP/mCherry intensities or estimate stress dose.
The convention is:

`image_ratio_red_green = c_session * Measured_Ratio_RG`

To compare an original image ratio with the model, divide it by `c_session`.
`Measured_Ratio_RG = f * Reporter_red / (Observed_Green + eps)` already
includes the reporter's FRET effect and `f`. Neither is applied again.
`c_session` is relative to that specific forward model, not an independent
estimate of `f`. Changes to the model, profile, initialization, protocol,
solver or observation times change the reference identifier.

## Why synthetic data were used

Experimental session reference measurements were not supplied for this
implementation. Synthetic data provide known ground truth, allowing a
check that the software recovers an injected scale and correctly transforms
previously unseen ratios. A successful check demonstrates computational
behaviour under the assumed scale model. It does not demonstrate biological
validity or quantitative accuracy of the actual microscope.

These are synthetic numeric ratios generated from model outputs, not
synthetic microscopy images. They do not test segmentation, optical
crosstalk, illumination or camera linearity.

## Reproducible procedure

1. Run the existing protocol simulator with explicit variant, parameter
   profile, dose, initialization and observation times. The example below
   uses ox, the default profile, constant stress of 75 uM, assumed two-hour
   zero-stress preincubation and nine observation times. These are software
   test settings, not a recommended experimental protocol.
2. Read `Measured_Ratio_RG` from the mapped `fluorescence.csv`.
3. Generate `image_ratio = 1.5 * model_ratio * (1 + noise)`. The known
   scale is 1.5. Seed is 20261005. The exact case uses zero noise; the noisy
   case uses independent uniform relative perturbations from -0.02 to +0.02.
   This bounded noise is a software test scenario, not a measured camera
   noise distribution.
4. Split alternating observation times into five fitting pairs and four
   holdout pairs. Holdout values are not used to estimate the scale.
5. Fit a single multiplicative coefficient by unweighted least squares
   through the origin: `c_session = sum(model_ratio * image_ratio) /
   sum(model_ratio ** 2)`. Store all fitted pairs and the fitting RMSE.
6. Divide the original synthetic image ratios by the fitted coefficient.
   Record the relative coefficient error and holdout RMSE against known
   model ratios. Keep original and corrected ratios in separate columns.
7. Save `calibration.json`, `verification.json` and `synthetic_ratios.csv`.
   Artifacts explicitly record `source_kind: synthetic` and
   `experimental_instrument_validated: false`, plus seed, generation
   equation, noise, data split and the source protocol metadata.

Run from the repository root in the configured environment:

```bash
python scripts/run_protocol.py \
  --variant ox --profile default --dose 75 --dose-units uM \
  --times-hours 0 0.25 0.5 0.75 1 1.5 2 3 4 \
  --initialization finite-preincubation --preincubation-hours 2 \
  --out-dir outputs/session-calibration-protocol

python scripts/verify_session_calibration.py \
  --protocol-dir outputs/session-calibration-protocol \
  --out-dir outputs/session-calibration-exact \
  --scale 1.5 --relative-noise 0 --seed 20261005

python scripts/verify_session_calibration.py \
  --protocol-dir outputs/session-calibration-protocol \
  --out-dir outputs/session-calibration-noisy \
  --scale 1.5 --relative-noise 0.02 --seed 20261005
```

Output directories must be new. This avoids replacing a prior verification.

## Tests and limits

Verification run on 2026-10-05 with the settings above:
The full test suite passed: 288 passed, 1 skipped, with FutureWarning treated
as an error. The whitespace check was clean.

| Case | Known scale | Recovered scale | Relative scale error | Holdout RMSE |
| --- | --- | --- | --- | --- |
| No noise | 1.5 | 1.5000000000000002 | 2.22e-16 | 9.61e-17 |
| Uniform relative noise +/-2% | 1.5 | 1.4926906799394253 | 0.004873 (0.4873%) | 0.004529 |

RMSE is in model ratio units. These figures describe one seeded synthetic
example, not experimental accuracy or an uncertainty confidence interval.
The exact case recovers the coefficient to floating-point precision.

Tests cover exact scale recovery, correction on unseen ratios, seeded noise,
nonunit model `f`, invalid/no-signal references, session/model mismatches,
JSON persistence and rejected artifact relabelling. The CLI test exercises
real protocol simulation and both synthetic verification scenarios.

At least three positive matched reference ratios are required by the fitting
API as a software guard, not as evidence that three measurements establish
experimental calibration. For real data, references must match biological
condition and time, pass QC and use acquisition settings represented by the
calibration. Cells must not be treated as independent biological replicates.
Population aggregation and image-pipeline integration are implemented as a
separate summary step; see `population_summary.md`. Real acquisition
settings and calibration references still require experimental verification.

The method assumes a constant multiplicative scale and zero intercept.
It does not establish detector linearity, detection limits, model validity,
experimental uncertainty or a valid scale across changed acquisition settings.
A synthetic calibration must never be used as the calibration of real images.
Even an artifact fitted from experimental references is not automatically
labelled as experimentally validated.
