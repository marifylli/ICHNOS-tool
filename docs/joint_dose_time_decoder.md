# Joint dose–time decoder

## What it estimates and why

A fluorescence ratio can be compatible with a small dose measured late or a
larger dose measured early. A dose decoder at a known time does not address this
confounding. This stage checks all pairs of dose and elapsed time on explicitly
specified discrete grids. It reports every compatible pair and abstains from a
point estimate if there is more than one.

Elapsed time means hours **from exposure onset to the first measurement**.
Observation input times are relative to the first measurement: zero, then known
positive increasing intervals. For a candidate elapsed time tau, predictions are
read at tau + relative_time. Thus the exposure onset would be tau hours before
the first measurement. This is not an estimate of an absolute calendar timestamp.
The model assumes the single onset/exposure protocol recorded in the table;
unknown repeated exposures, unknown clearance and uncertain initial conditions
are not fitted.

## How it works

`ichnos/joint_decoder.py` reuses `DoseTable` and the existing discrete dose
compatibility validation. For each candidate elapsed time it selects the stored
model predictions at shifted observation times. All required times must be
present in the response table; an incomplete table is rejected rather than
silently dropping candidate times. Time matching accepts floating point rounding
within 1e-12 hours, with exactly one table match required. There is no time or dose
interpolation and no extrapolation. The existing table validator requires at
least one positive observation time in every selected trajectory.

A pair is compatible if every absolute ratio residual is at most the supplied
positive tolerance. A scalar tolerance applies to all observations; the Python
API also accepts a vector. Ratios must already be calibrated to model scale,
red/green (mCherry/GFP). Session calibration and model factor f are not applied
again. Optional green measurements and a detector floor must share units; a
measurement at or below that floor produces an abstention.

Statuses:

| Status | Meaning |
| --- | --- |
| `unique_grid_pair` | One compatible grid pair; both point estimates returned |
| `ambiguous` | Multiple compatible pairs; all retained, no point estimates |
| `no_grid_pair_match` | No pair explains all observations within tolerance |
| `out_of_response_domain` | At least one observed ratio lies outside all candidate predictions plus tolerance |
| `below_detection_floor` | Green signal fails the explicitly provided detector floor |

A unique grid pair does not prove uniqueness between grid points or outside the
searched domain. Tolerance is a compatibility allowance, not a confidence
interval. An out-of-response result does not establish whether the cause is
noise, a dose outside the grid, an incorrect protocol, or model mismatch.

Outputs record all candidate pairs, grids, ratio tolerances, time convention,
table artifact ID, and decoder source hash. CLI output additionally records the
input CSV SHA256 and declared synthetic/experimental data kind. The table retains
model, profile, initialization, solver and exposure provenance.

## Synthetic verification: how and why

The checks deliberately use synthetic data because no complete experimental
validation of dose–time recovery is available. Synthetic inputs provide known
answers and controlled ambiguity, allowing us to verify the implementation.
They do not calibrate an instrument and do not validate the biological model.

Analytic synthetic tables test a unique pair, dose–time confounding from one
measurement, time-flat responses, mutually inconsistent observations, response
envelope rejection, malformed grids, missing shifted times, and detector-floor
abstention. Tests also verify CLI provenance and refusal to overwrite output.

A separate actual-model check builds an oxidative response table with doses
25 and 75 uM at 0.5, 1 and 1.5 hours, using the complete default profile and
finite preincubation of 2 hours. A fresh simulation at 75 uM supplies the ratios
at 1 and 1.5 hours. With relative times [0, 0.5], candidate elapsed times
[0.5, 1] and absolute ratio tolerance 1e-8, the decoder recovers the grid pair
75 uM and 1 hour. These exact noise-free, same-model measurements verify
computational consistency under fixed assumptions, not robustness to real noise
or independent biological validation. The small tolerance is a numerical test
setting, not a recommended experimental measurement error.

## Commands

Build the table with the existing runner, including every shifted time:

```bash
python scripts/run_decoder.py build \
  --variant ox --doses-uM 25 75 --times-hours 0.5 1 1.5 \
  --initialization finite-preincubation --preincubation-hours 2 \
  --out joint-table.json
```

Create an observation CSV with columns `relative_time_hours` and
`calibrated_ratio_red_green`. For two measurements half an hour apart, the first
column contains 0 and 0.5. Supply your already calibrated ratios in the second
column; choose a tolerance appropriate to your measurement assumptions.

```bash
python scripts/run_joint_decoder.py \
  --table joint-table.json --observations observations.csv \
  --elapsed-times-hours 0.5 1 --ratio-tolerance 0.01 \
  --data-kind synthetic --out joint-result.json
```

This command starts from calibrated ratios. The automatic sample-summary
adapter still supports known-time discrete/continuous dose decoding; it is not
yet extended to unknown onset. Joint continuous interpolation, statistically
calibrated uncertainty, experimental recovery tests and unknown protocol
parameter estimation remain future work.
