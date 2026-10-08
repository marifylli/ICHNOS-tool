# Dense posterior verification

## What this benchmark tests

The original 15-node example checks mechanics but cannot establish uncertainty
coverage. A discrete 95% set can accumulate almost 100% of posterior mass, and
24 included truths are not evidence of experimentally calibrated uncertainty.

`verify_dense_posterior.py` evaluates the production joint decoder at every
node of a denser oxidative grid, with separate targets for each repeat. Defaults
are 17 doses from 0 to 800 uM and 33 times from 0 to 8 h (561 nodes), three draws
per node and two generating scenarios. This is 3,366 independent synthetic
summary draws, each compared using ratio alone and ratio plus green.

- **Matched:** default model generates and decodes the data.
- **Clearance misspecification:** generation uses exponential stress clearance
  at 0.2 per hour; decoding still assumes constant exposure. Dose refers to the
  initial dose. Calibration and fitted noise remain those of the default model;
  they are not refitted to conceal the mismatch.

The covariance is fitted from 36 separate synthetic culture labels. Its assumed
Gaussian generating values remain illustrative. Targets are generated directly
in `(log ratio, green_norm)` space; this benchmark does not simulate microscopy
noise or exercise segmentation. The earlier noiseless snapshot image workflow
is a separate check and does not fill that gap.

The prior approximates uniform density on physical dose/time coordinates, with
explicit bounds [0,800] and [0,8]. Midpoint-cell widths determine node masses.
The ratio-only comparison uses the *marginal* one-dimensional Gaussian variance
of log-ratio. It does not condition on the unobserved green channel. Both methods
receive exactly the same noisy observation and prior within each comparison.

## Run

After producing `outputs/two-observable-verification`:

```bash
python scripts/verify_dense_posterior.py \
  --snapshot-dir outputs/two-observable-verification \
  --out-dir outputs/dense-posterior-verification
```

Progress is printed per dose. The output directory must not already exist.
Optional `--dose-nodes`, `--time-nodes`, `--repeats`, `--seed` and
`--clearance-rate` make the scenario explicit. An odd number of dose nodes
retains the illustrative 400 uM point. More repeats improve the precision of
per-node diagnostics but increase runtime. Three repeats per node are only a
small diagnostic sample, not sufficient local coverage estimation.

## Outputs and interpretation

- `trials.csv`: individual MAP dose/time absolute errors, number of tied MAP
  nodes, credible-set node count/fraction, actual accumulated mass, truth inclusion
  and the production joint decoder's poor-fit flag.
- `per-node.csv`: means of those numeric diagnostics by truth node, scenario
  and method. Actual posterior mass and set size must be read alongside inclusion.
- `verification.json`: aggregate diagnostics, grid, prior, seed and assumptions.
  Equal trial counts per node mean aggregate inclusion is **not** prior-predictive
  coverage under the cell-area prior. There is no pass/fail assertion of 95% coverage.
- `posterior-comparison.png` and `.pdf`: ratio-only versus ratio plus green for
  the same 400 uM / 8 h target within each scenario. Both panels use the same
  color scale; each scenario has its own independently drawn target. The red
  cross is the generating truth. Colors are node masses, not density per uM/hour.
- JSON artifacts retain both forward tables, calibration, noise-training rows,
  held-out observations, illustrative production posterior results and prior.

MAP ties are retained. Error minima and maxima are calculated separately over
the tied nodes for each coordinate; no arbitrary first node is presented as a
unique estimate. Zero-dose and initial-time ambiguities remain in the benchmark.
A credible set can be disconnected; its node count is not an interval width.

The Mahalanobis threshold is the existing explicit value 25, not a calibrated
p-value. A misspecified model can produce a plausible but wrong posterior without
triggering this threshold. Evaluate bias, set size and inclusion as well as flags.

## Remaining scientific limits

All truths are on-grid. Grid refinement, off-grid targets, heteroscedastic noise,
model-parameter and calibration uncertainty, ER coverage and real biological
replicates remain necessary for stronger claims. A finer grid alone does not
resolve biological non-identifiability. Wiki figures from this workflow must
be labeled **in silico demonstration under assumed noise and calibration**.

The demo uses green references at 0 and 800 uM, both at 8 h. Those values are not
hard-coded requirements of the calibration API: reference dose/time coordinates
must be known, represented in the forward table and measured consistently.
High-dose saturation is not established. The numerical minimum of three
replicates per condition is not a validated wet-lab sample-size recommendation.
