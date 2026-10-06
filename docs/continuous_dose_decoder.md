# Continuous dose decoder at known observation times

Implemented on 2026-10-06. Extends the discrete decoder with checked
piecewise-linear interpolation in dose. Onset time is still known, not
estimated. All outputs remain model-based and not experimentally validated.

## What was implemented, how and why

`ichnos/continuous_decoder.py` interpolates the vector of predicted ratios
between adjacent grid doses. For a segment `[d0, d1]`, write
`d = d0 + alpha * (d1 - d0)`, with `alpha` in `[0, 1]`. At each observation
time, the interpolated ratio is `R0 + alpha * (R1 - R0)`.

Piecewise-linear interpolation is used because it is explicit, reproducible
and introduces no spline overshoot. It does not assume a monotonic global
response. It still approximates the nonlinear forward model, so intermediate
doses require computational checks against fresh simulations.

For each segment, solve every timepoint's ratio compatibility inequality
for alpha and intersect the resulting intervals. Merge touching dose
intervals across segments. This preserves all compatible regions rather
than selecting the first local optimum or nearest grid dose.

Within one compatible region, report the dose that minimizes the weighted
squared ratio residual on the linear segments. The weights use the explicit
effective compatibility tolerances. If several equal minima exist, return
ambiguity instead of an arbitrary point. An exactly flat compatible segment
also causes abstention from a point estimate.

The reported regions are compatibility regions, not confidence intervals
or evidence of biological identifiability. A broad connected region can
indicate weak discrimination even if a numerical best fit is reported.
Inspect region width, not only the point estimate.

## Interpolation verification

Before continuous decoding, run a direct forward simulation at 25%, 50%
and 75% of every grid segment. Re-simulate the original grid anchors too,
checking the complete profile, prepared-model identity, solver and anchor
predictions against the original artifact.

Record the tested doses, direct ratio predictions, interpolated predictions,
maximum absolute error per segment/time, permitted error and pass/fail
status in a separate interpolation-check artifact. The artifact is bound to
the exact dose table and protected by its own content hash.

The caller must supply a positive interpolation error allowance in ratio
units. If any checked point exceeds it, the decoder returns
`interpolation_check_failed` with no dose estimate. Refine the grid or
reconsider the permitted approximation error deliberately; do not simply
increase the allowance to make an inadequate grid appear accurate.

For a passing check, the effective ratio compatibility margin is:

`measurement_ratio_tolerance + interpolation_error_allowance`

The second term accounts for the explicitly accepted approximation in the
compatibility calculation. Neither term is automatically a measured noise
level, and their sum is not a statistical confidence interval.

Three points per segment provide a sampled computational check, not a
certified uniform error bound at every possible intermediate dose. Local
features between checked doses can be missed. The output explicitly records
this limitation. Grid refinement and additional independent checks remain
necessary when the response varies sharply or the desired accuracy is tighter.

## Status and domain behaviour

| Status | Meaning |
| --- | --- |
| `single_compatible_region` | One connected dose region, with an interpolated residual-minimizing estimate |
| `ambiguous` | Separate compatible regions or several equally good minima; no point estimate |
| `flat_response` | Compatible constant response cannot distinguish doses; no point estimate |
| `interpolation_check_failed` | Sampled interpolation error exceeds the specified allowance |
| `no_continuous_match` | No interpolated trajectory jointly meets the ratio margins |
| `out_of_response_domain` | Signal exceeds the interpolated response envelope plus margins |
| `below_detection_floor` | Supplied green values fail the explicitly supplied floor |

No extrapolation beyond the grid is performed. Boundary contact is recorded.
No match does not prove a particular cause or a true dose beyond the domain.
Calibration error, wrong protocol, model error or measurement noise can also
produce incompatibility. Existing detection-floor checks remain available.

## Commands

First build a dose table using the existing builder. Illustrative software
test settings for ox, not a recommended experimental design:

```bash
python scripts/run_decoder.py build \
  --variant ox --profile default \
  --doses-uM 25 30 35 40 45 50 --times-hours 0.5 1 \
  --initialization finite-preincubation --preincubation-hours 2 \
  --out continuous-dose-table.json

python scripts/run_continuous_decoder.py check \
  --table continuous-dose-table.json --allowed-error 0.001 \
  --out interpolation-check.json

python scripts/run_continuous_decoder.py decode \
  --table continuous-dose-table.json --validation interpolation-check.json \
  --observations observations.csv --ratio-tolerance 0.00001 \
  --data-kind synthetic --out continuous-dose.json
```

The observations CSV uses `time_hours` and
`calibrated_ratio_red_green` in the model ratio scale. The example error
allowances are software-test settings, not experimentally calibrated errors.
All output filenames must be new.

For the checked sample-summary workflow, add these options to the existing
`scripts/decode_sample.py` command:

```bash
--mode continuous --interpolation-validation interpolation-check.json
```

All sample, time, calibration and protocol checks are retained. The calibrated
sample ratio is not divided by the session coefficient again. Discrete mode
remains the default. A validation artifact supplied with discrete mode, or
continuous mode without one, is rejected.

## Why and how synthetic data were used

Verification on 2026-10-06: 382 passed, 1 skipped in the full suite,
with FutureWarning treated as an error. The whitespace check was clean.

Controlled linear tables check exact off-grid recovery. Nonmonotonic and
constant tables check ambiguity and flat-response abstention. These cases
test the mathematics without relying on accidental properties of the current
biological model.

Separately, an ox grid at 25, 30, 35, 40, 45 and 50 uM, observed at 0.5 and
1 hour after assumed two-hour preincubation, is checked with direct intermediate
simulations. Independently simulated 37.5 and 42.5 uM trajectories, absent
from the original grid, are decoded with error below 0.5 uM in this test.
The test also checks that the known dose falls in the compatible region.

CLI tests exercise check-artifact generation, continuous decoding and output
overwrite rejection. Sample-linkage tests confirm that the same provenance
checks and calibration handling work in continuous mode.

Known synthetic doses provide a reference for software verification. They
do not establish real experimental accuracy, noise coverage or sensitivity.
The tested accuracy applies to this grid and these model/protocol settings,
not arbitrary doses, variants or instruments. Joint dose/time estimation
and experimentally justified uncertainty remain separate subsequent work.
