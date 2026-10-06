# Calibrated sample summaries to joint dose–time decoding

## What changed and why

`decode_sample` and `scripts/decode_sample.py` now accept `mode="joint"`.
The existing sample selection, CSV fingerprint, cell count, calibration identity,
model/profile, solver, exposure and initialization compatibility checks run before
decoding. This joins the image-derived population summaries to the joint grid
search without applying session calibration or the model factor f again.

The caller chooses `sampling_time_hours` or `measurement_time_hours` explicitly.
In joint mode these recorded times supply **only measurement spacing**: the first
selected time is subtracted from all subsequent times. For example, recorded
sampling times [0.5, 1] and measurement times [0.75, 1.25] both imply spacing
[0, 0.5]. Candidate elapsed times are provided separately; they mean time since
exposure onset at the first selected measurement. Recorded timestamps are not
assumed to identify onset. The selected time field must contain finite,
nonnegative, strictly increasing values and at least one positive value, as in
the existing time validator.

Every elapsed time plus relative measurement time must be present in the
response table. Missing shifted times are errors, not silently skipped
candidates. Output retains every compatible dose–time pair. A unique pair is a
result on the specified discrete grid, not continuous identifiability or a
statistical confidence interval. Interpolation validation is accepted only for
continuous known-time mode, not for joint mode.

The sample linkage records the original selected timestamps and explicitly marks
`recorded_times_used_as_spacing_only=true`, alongside session, sample, condition,
calibration identity, coefficient, CSV fingerprint and usable cell counts.
The summary's image green floor remains recorded as a filtering setting; this
stage does not invent a calibrated instrument detection floor.

## Use

```bash
python scripts/decode_sample.py \
  --samples-csv summary/samples.csv \
  --summary-metadata summary/metadata.json \
  --dose-table joint-table.json \
  --calibration-reference reference/metadata.json \
  --session-id s1 --sample-id sample-a \
  --time-field sampling_time_hours \
  --mode joint --elapsed-time-grid-hours 0.5 1 \
  --ratio-tolerance 0.01 --out sample-joint-result.json
```

Here the response table must cover both shifted trajectories. Output files are
created exclusively; existing results are not overwritten. Default `discrete`
and `continuous` modes retain their known-time behavior. The elapsed grid is
required for joint mode and rejected in other modes.

## Synthetic verification: how and why

Tests reuse the image-processing integration fixture: a model-generated 25 uM
reference defines a synthetic session factor of 1.5. Three simplified synthetic
cell images per time are generated using the 75 uM model ratio at 0.5 and 1 hour,
then segmented, corrected, extracted, exported and summarized with that session
calibration. The adapter recovers the pair 75 uM / 0.5 hour from relative spacing
[0, 0.5] and candidate elapsed grid [0.5]. Selecting measurement timestamps
[0.75, 1.25] instead gives the same result because spacing, not the original time
origin, is used. This intentionally small grid verifies the adapter; ambiguity
across candidate times is tested in the joint decoder tests, not claimed absent
throughout the model domain.

Contract tests reject absent/irrelevant elapsed grids, joint interpolation
artifacts and missing response times. CLI tests verify correct output and
refusal to overwrite existing results. Existing sample provenance and calibration
checks remain active for joint decoding.

These controlled synthetic inputs are used to check software composition with a
known answer. The reference is not a real instrument calibration; shared-model,
noise-free recovery is not independent biological validation. Experimental
accuracy, uncertain clearance/initial states, joint continuous interpolation and
statistical uncertainty remain outside this stage.
