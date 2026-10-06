# First dose decoder: discrete compatibility at known times

Implemented on 2026-10-06. This first version determines which doses in an
explicit model-generated grid are compatible with ratios at known elapsed
times. It does not yet interpolate a continuous dose, estimate onset time,
fit a biological response curve, or establish experimental accuracy.

## What was implemented, how and why

`ichnos/decoder.py` builds a dose-response table by running each requested
dose with a fresh runner through the shared simulator and stress protocol.
All runs use the same complete parameter profile, initialization assumption,
observation times and exposure model. A fresh runner for every dose prevents
the previous simulation's state from contaminating the next one.

For each dose, store the trajectory of `Measured_Ratio_RG` at the requested
times. This observable already includes FRET and model `f`; they are not
applied again. Calibration here means a computational forward lookup table,
not an experimentally fitted microscope or biological calibration.

The JSON artifact contains dose and time grids, ratio predictions, complete
parameter profile, baseline SBML and its hash, exposure and initialization
records for every dose, actual solver settings, package versions, hashes of
the relevant code and an artifact identifier. The loader checks the artifact
hash and validates grid ordering, times, shape, values, ratio convention and
model-generated provenance. Existing files are never overwritten.

This record makes the forward map inspectable and reproducible. A decoder
must invert the same model and protocol used to generate its table.

## Input contract

Input ratios must already be on the model's `Measured_Ratio_RG` scale, using
the appropriate session calibration if they originate from images.
The observation CSV contains one sample trajectory, in increasing time:

```csv
time_hours,calibrated_ratio_red_green
0.5,MODEL_SCALE_VALUE_AT_0_5_H
1.0,MODEL_SCALE_VALUE_AT_1_H
```

The placeholders are not numeric data. Keep samples and sessions separate.
Select the actual elapsed times appropriate to the biological readout;
sampling and acquisition time are not automatically interchangeable.
The decoder requires exact agreement with the table's recorded time grid.
No automatic adapter from `samples.csv` is implemented in this step.
The caller must select one sample, the correct time field and the matching
variant, profile, protocol and session-calibration reference.

The caller must specify a positive absolute ratio tolerance, either one
scalar for all times or one value per time via the Python API. No tolerance
is silently inferred from cell quartiles or labelled as experimental noise.

## Decision rule

For dose `d`, prediction `R(d,t)`, observation `y(t)` and tolerance `delta(t)`:

`compatible(d) = all(abs(R(d,t) - y(t)) <= delta(t))`

Retain every compatible grid dose. There is no monotonicity assumption,
nearest-neighbour tie breaking or extrapolation.

| Status | Meaning | Dose estimate |
| --- | --- | --- |
| `unique_grid_match` | Exactly one grid dose meets all tolerances | That grid dose only |
| `ambiguous` | Several grid doses meet all tolerances | Missing; return all candidates |
| `out_of_response_domain` | At least one observation lies outside the grid prediction envelope plus tolerance | Missing |
| `no_grid_match` | Observations lie within the envelope, but no grid trajectory jointly matches them | Missing |
| `below_detection_floor` | A supplied green measurement is at or below the supplied floor | Missing |

A unique grid match does not establish uniqueness between grid doses or a
precise continuous concentration. Adjacent doses could still be compatible
on a finer grid. The tolerance and candidate set are not confidence intervals.
An out-of-response-domain result does not prove that the true dose is outside
the dose grid: wrong protocol, noise, calibration or model error could also
cause disagreement. No grid match likewise does not identify a cause.

## Detection floor

Optionally supply green measurements and a floor in the same corrected image
units. All supplied green values must exceed that floor. The software does
not infer the floor from model concentrations or camera dtype. If it is not
supplied, record `detection_floor_checked: false`; do not imply that the signal
passed experimental detectability checks. A synthetic floor only tests the
abstention path and is not a measured detector sensitivity.

## Commands

Example software-test grid and initialization, not an experimental design:

```bash
python scripts/run_decoder.py build \
  --variant ox --profile default \
  --doses-uM 0 25 75 150 300 \
  --times-hours 0.5 1 \
  --initialization finite-preincubation --preincubation-hours 2 \
  --out dose-table.json

python scripts/run_decoder.py decode \
  --table dose-table.json --observations observations.csv \
  --ratio-tolerance 0.01 --data-kind synthetic \
  --out decoded.json
```

The tolerance 0.01 is illustrative, not an experimentally calibrated error.
For `--green-floor`, the CSV additionally needs `corrected_mean_green`.
Equilibrium initialization and explicitly supplied first-order clearance
are supported through the same protocol functions. There is no fitted
default clearance rate. New output filenames are required.

## Synthetic verification and why it was used

Verification on 2026-10-06: the full suite passed with 342 passed,
1 skipped and FutureWarning treated as an error. The 30 new decoder tests
passed, including independent known-grid-dose recovery for ox and er.
The whitespace check was clean.

Known-dose model outputs provide ground truth for checking grid recovery.
Tests independently simulate a known grid dose through the shared protocol
and decode its ratio trajectory. Both ox and er are checked, with the first
CLI example using ox. CLI tests also check serialization and output identity.

Controlled tables exercise identical responses at different doses,
off-grid observations, inconsistent multi-timepoint signals, out-of-envelope
signals, explicit detection-floor failure, invalid values and mismatched
times. These controlled cases check abstention and ambiguity without relying
on an accidental property of the current biological model.

This is computational verification, not independent biological validation.
Model-generated data necessarily share the model assumptions. Experimental
response fitting, calibrated noise/uncertainty coverage, grid refinement,
continuous interpolation and joint dose/time estimation remain future work.
