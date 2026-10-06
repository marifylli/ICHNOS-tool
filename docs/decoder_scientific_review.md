# Decoder scientific review — 2026-10-06

## Finding and scope correction

The implemented decoders use `Measured_Ratio_RG` alone. Green enters image QC
and optional detection-floor checks; it is not an inferential observable.
This does not implement the intended two-observable (ratio, green) snapshot
inverse problem. Passing synthetic workflow tests establishes computational
consistency for selected cases, not general dose/time identifiability, estimator
uncertainty, or experimental accuracy.

`decode_joint_dose_time` accepts a single observation at relative time zero.
For multiple observations it requires known relative spacing. A single accepted
grid pair may reflect the chosen grid/tolerance; it is not evidence that the
ratio uniquely identifies both dose and time. Time courses can also be ambiguous.
The claim supported today is **ratio-based compatibility under a fixed model,
protocol, initialization and explicit tolerance**.

## Reproduced model check

Using `scripts/run_protocol.py`, default profiles, equilibrium initialization,
constant stress and observation times 0.5, 1 and 8 hours:

| Variant | Dose (uM) | Ratio 0.5 h | Ratio 1 h | Ratio 8 h | Green 8 h |
| --- | --- | --- | --- | --- | --- |
| ox | 10 | 0.821136 | 0.811409 | 0.822489 | 7.055040 |
| ox | 200 | 0.527005 | 0.457822 | 0.818912 | 27.115954 |
| ox | 400 | 0.502045 | 0.435539 | 0.818398 | 32.186359 |
| ox | 800 | 0.497467 | 0.431264 | 0.818270 | 33.376963 |
| er | 10 | 0.822799 | 0.822799 | 0.822799 | 2.356958 |
| er | 50 | 0.822799 | 0.822798 | 0.822786 | 2.357170 |
| er | 100 | 0.822797 | 0.822787 | 0.822611 | 2.360255 |
| er | 800 | 0.817848 | 0.777549 | 0.844047 | 3.137814 |

For example:

```bash
python scripts/run_protocol.py --variant ox --dose 800 --dose-units uM \
  --times-hours 0.5 1 8 --initialization equilibrium --out-dir outputs/probe-ox-800
```

The oxidative late ratio carries much less dose contrast than green. The ER
10–100 uM ratio is nearly constant, not exactly constant at full precision.
Green is also nearly flat in that ER range, so adding green cannot by itself
repair an insensitive forward model. Dose range, model profile and the
experimental regime must agree. The ER 800 uM ratio exceeds baseline by 8 hours.
These comparisons are model outputs, not measured biosensor behavior.

## Two-reference green normalization: required contract

For independently acquired reference samples in the same session, a candidate
mapping is `z_image = (G_image - G_low) / (G_high - G_low)`. Apply the analogous
mapping to `Observed_Green` using the matching reference protocol, dose and
reference time. Reference times and acquisition settings must be fixed and
recorded. If reference values depend on time, that dependence must be modeled
explicitly; do not normalize using the unknown target onset time.

The low/high samples must be separate from decoded targets. "High" must be a
biologically characterized reference, not automatically the largest grid dose
or a saturated camera measurement. A near-zero reference span must be rejected.
Normalization places the reference anchors at 0 and 1; noisy or extrapolating
targets can lie outside that interval and must not be silently clipped.

Required provenance includes session, biological culture/replicate, condition,
reference/target role, reference dose and elapsed time, exposure/gain/filter,
background controls, and independent calibration/validation split. Image
intensity must be generated from model green in future synthetic checks: the
existing end-to-end images hold green at 1000 and encode only the model ratio.
They cannot test a green-sensitive inverse model.

## Noise and posterior: required contract

Estimate measurement variation from independent biological culture replicates
at matched conditions/time, not from the number of segmented cells. Repeated
images/cells are nested observations. Estimate covariance as well as marginal
variances: ratio and green share a denominator and are generally correlated.
Do not pool different expected responses as if their biological variation were
noise. Estimation of sample-mean versus individual-culture uncertainty must be
consistent with the quantity being decoded.

A joint likelihood on a declared transform can then be combined with an
explicit prior over the dose/time grid. Priors must account for grid cell size
if they represent a continuous uniform density on an irregular grid. A
posterior region is conditional on that noise model, calibration, parameter
profile, protocol and prior; it is not an automatically validated error bound.
Use independent holdout cultures and recovery/coverage checks, and propagate
calibration/model uncertainty where relevant.

`log(green_norm)` is undefined at the low reference (zero) and for negative
normalized observations. Do not add an arbitrary epsilon to hide that problem.
Choose a likelihood appropriate for normalized green (possibly on its linear
scale) or a separately calibrated strictly positive green measurement model.
A posterior must retain disconnected modes and report prior/domain dependence.

## Ratio calibration and background

The present calibration fits `image_ratio = c_session * model_ratio` through
zero. This is conditional on model reference predictions. Independent targets
and an independent measurement calibration/validation dataset are necessary to
avoid using synthetic self-consistency as evidence that the biology is correct.
Three fitted reference pairs provide very little diagnostic information.

`ichnos.calibration_diagnostics.compare_ratio_calibration` now compares the
origin-constrained fit with an affine fit, reports reference range, both RMSEs
and residual degrees of freedom, and optionally evaluates held-out paired
references. It does not automatically replace the calibration or establish
reference independence. The CLI examines the existing artifact's fit pairs:

```bash
python scripts/audit_ratio_calibration.py \
  --calibration outputs/e2e-verification/calibration.json \
  --out outputs/ratio-calibration-audit.json
```

Improvement in training RMSE from an extra parameter is not validation. A ratio
intercept is only a diagnostic: offsets in the two channels generally produce
`(a_R R + b_R)/(a_G G + b_G)`, not a constant offset to the ratio. Dark/background
controls and cellular autofluorescence controls are distinct requirements.
The pipeline does call per-channel background subtraction; it does not establish
that cellular autofluorescence has been removed.

## Photobleaching

`correct.correct_photobleaching` exists but is not called by the pipeline.
Do not fit and divide out the biological response trajectory to activate it:
that risks removing the response being decoded. A correction needs independent
controls and illumination/acquisition history; destructive samples from the
same culture are not the same cells repeatedly illuminated. This remains open
and exported data must continue to be described as not bleaching-corrected.

## Repository publication and next gate

At review start `origin/main` was nine commits behind
`origin/feat/end-to-end-verification`, which pointed to `80b92a0` and contained
no end-to-end scripts. The end-to-end and visualization patches existed locally
and are included in the review integration branch. README and limitations have
been corrected to describe the actual implementation.

The next scientific implementation gate is a traceable two-observable mapping
and snapshot compatibility test suite over the intended ox/er domain, followed
by a biological-replicate noise artifact and posterior inversion. These are
outstanding work, not capabilities delivered by this documentation/diagnostic
repair. Do not describe the original snapshot timer design as complete.

### Publish the tested local work from the user's checkout

The assistant's environment could fetch public branches but could not push:
Git reported no available GitHub username/credentials. No remote branch was
changed by the review. After applying the review patch in the checkout that
already has the end-to-end/visualization work, commit the explicit source files:

```bash
git add README.md pyproject.toml docs/limitations.md \
  docs/end_to_end_verification.md docs/decoder_scientific_review.md \
  ichnos/calibration_diagnostics.py scripts/audit_ratio_calibration.py \
  scripts/verify_end_to_end.py scripts/render_e2e_report.py \
  tests/model/test_calibration_diagnostics.py \
  tests/model/test_end_to_end_verification.py

git diff --cached --check &&
git commit -m "Publish verification workflow and correct decoder scope" &&
git push origin feat/end-to-end-verification
```

On success, fast-forward main; stop if any command fails:

```bash
git switch main &&
git pull --ff-only origin main &&
git merge --ff-only feat/end-to-end-verification &&
git push origin main
```

Do not force-push or reset away local work if branches have diverged. These
commands integrate the stacked history and documentation; they do not implement
the outstanding two-observable posterior decoder.
