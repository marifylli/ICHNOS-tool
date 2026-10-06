# Conditional snapshot posterior

## Purpose and scope

`ichnos.posterior` adds a probabilistic alternative to tolerance-based snapshot
compatibility. It combines ratio and reference-normalized green to assign mass
to each stored (dose, elapsed time) pair. Ambiguous modes remain visible. This
is a **discrete conditional posterior**, not continuous inversion or proof of
experimental accuracy. The older ratio-only sample adapters are unchanged.

## Observation and replicate contract

Coordinates are `(log(calibrated red/green ratio), green_norm)`. Green is linear
because two-reference normalization legitimately gives zero or negative values;
logging it would exclude valid observations. Ratios must be strictly positive.
The same session, acquisition, reference independence, source-kind and detection
floor checks used by the snapshot decoder apply. Green is never clipped.

Each noise-training row follows the snapshot observation JSON contract and adds:

```json
{
  "observation_unit": "biological_replicate_summary",
  "condition_id": "ox-400-8h",
  "dose_uM": 400,
  "time_hours": 8
}
```

Retain all required snapshot fields, including session, specimen and biological
replicate IDs, acquisition settings, image-scale ratio, corrected green and
source. Each row represents **one independent biological culture summary**, not
one cell, field or repeated measurement. Use the same summary estimator and
sampling procedure for training and the target. Unequal precision from different
cell counts is not modeled. The software checks declared IDs; it cannot verify
actual biological independence.

Training IDs cannot repeat, including across conditions. Calibration references,
noise-fitting cultures and the target must be distinct. Paired longitudinal
cultures require a different hierarchical model and are currently rejected.
Conditions must name exact table nodes. At least three cultures per condition
are required numerically; this is not a recommendation that three are sufficient.

## Noise fitting and likelihood

For each condition, subtract its empirical mean in transformed coordinates.
Pool residual outer products and divide by `N - number_of_conditions` to obtain
a full 2 x 2 sample covariance. Between-condition signal therefore does not
inflate noise. The off-diagonal term retains ratio/green correlation. This is
particularly relevant because the ratio already contains green.

The covariance is assumed constant across the grid. Singular or ill-conditioned
fits are rejected without automatic regularization. Recorded empirical means
and their differences from model predictions support bias diagnostics; inference
does not silently subtract those discrepancies.

At each grid node, evaluate a bivariate Gaussian centered on the model's
transformed response, multiply by the declared prior mass and normalize using
log-sum-exp. The fitted covariance, model and calibration are held fixed.
Uncertainty in these quantities is **not propagated**. No biological noise is
inferred from synthetic cells or assigned from their count.

## Priors and output interpretation

A prior JSON is mandatory:

- `{"kind":"uniform_grid"}` assigns equal mass to each stored node.
- `{"kind":"uniform_density","dose_bounds_uM":[0,800],"time_bounds_hours":[0,8]}`
  assigns masses proportional to midpoint-cell areas within explicit bounds.
  Bounds must enclose the grid. This approximates uniform density on physical
  dose/time coordinates; it differs from uniform node mass on an irregular grid.
- `{"kind":"explicit_mass","masses":[...]}` supplies a nonnegative matrix in
  dose-row/time-column order, with positive total mass. Zeros exclude nodes.

Outputs include the full posterior matrix, both marginal distributions, all
MAP ties and a discrete credible set. The set accumulates highest **node masses**
until the requested mass is reached, retaining all cutoff ties. It may contain
disconnected alternatives and exceed 95% mass. It is not a continuous HPD region
or a frequentist confidence interval. Results depend on the grid and prior.

`prior_only` flags a likelihood that is flat on prior support.
`poor_model_fit` flags a minimum squared Mahalanobis distance exceeding the
explicit user threshold. The normalized conditional posterior is still returned
but must not be interpreted as evidence that the model fits. This guard is not a
calibrated p-value. Below-floor observations produce no posterior.

## Commands

Starting with the artifacts from `verify_snapshot_decoder.py`:

```bash
python scripts/verify_posterior.py \
  --snapshot-dir outputs/two-observable-verification \
  --out-dir outputs/posterior-verification
```

The output directory must be new. It contains training summaries, fitted noise,
prior, a separate target, held-out observations, verification JSON and
`posterior.png`. The heatmap shows categorical grid-node masses, not a continuous
physical density surface.

The general CLI can reuse these example files:

```bash
python scripts/run_posterior.py fit-noise \
  --table outputs/posterior-verification/table.json \
  --calibration outputs/posterior-verification/calibration.json \
  --replicates outputs/posterior-verification/replicates.json \
  --green-floor 100 --out outputs/refitted-noise.json

python scripts/run_posterior.py decode \
  --table outputs/posterior-verification/table.json \
  --calibration outputs/posterior-verification/calibration.json \
  --noise outputs/refitted-noise.json \
  --observation outputs/posterior-verification/observation.json \
  --prior outputs/posterior-verification/prior.json \
  --max-mahalanobis-squared 25 --credible-mass 0.95 \
  --out outputs/refitted-posterior.json
```

Files are not overwritten. Noise artifacts retain training rows, fitted
covariance, condition diagnostics and table/calibration fingerprints. Validation
recomputes the fit from those rows. Fingerprints detect changes, not authenticity.

## Why synthetic verification, and what it establishes

No independent experimental replicate dataset is supplied. The verifier therefore
explicitly generates summary-level Gaussian perturbations with known standard
deviations 0.02 in log-ratio and 0.035 in normalized green, correlation -0.35.
These are illustrative assumptions, not measured instrument or biological noise.
It uses 12 independent synthetic culture labels at each of three conditions,
then separate target and held-out draws. It does not simulate noisy microscopy
or validate segmentation under experimental noise.

With the default seed, the illustrative 400 uM / 8 h target gives approximately
53.5% mass at 400 uM / 8 h and 46.5% at 800 uM / 8 h. Both are retained in the
95% credible set. This illustrates loss of certainty when noise is introduced,
not successful unique dose recovery.

A small held-out experiment samples truths from the declared grid prior and
records inclusion in credible sets. The default run includes 24/24 truths;
this small, same-model synthetic result does not establish nominal 95% coverage.
`computational_checks_passed` checks normalization and accumulated posterior
mass only. Experimental validation, condition-dependent noise, calibration and
covariance uncertainty, model discrepancy, off-grid truth and domain-wide
identifiability remain open.
