# Two-observable snapshot decoder

## Implemented scope

The new `ichnos.snapshot_decoder` module uses a **single snapshot**, consisting
of an original corrected image red/green ratio and a corrected green summary.
It searches every dose/time pair in a forward-model table containing both
`Measured_Ratio_RG` and `Observed_Green`. No target elapsed time or inter-sample
spacing is supplied. The original ratio-only decoders retain their existing
interfaces.

This first two-observable implementation uses two explicit compatibility
allowances. It is not a likelihood or posterior. It returns every compatible
grid pair, and abstains from a point estimate on ambiguity, missing joint
agreement or insufficient green signal. A unique pair is conditional on the
chosen discrete grid, model, protocol, calibration and tolerances. It does not
establish domain-wide structural or practical identifiability.

## Forward model and artifacts

`build_observable_table` shares the production forward simulation with
`build_dose_table`. Each dose uses a fresh protocol runner, the complete selected
profile and the declared initialization/clearance. The new table stores both
observable arrays at the same dose/time grid and carries the existing
model/profile/solver/protocol provenance. FRET and the model factor f are
already included; the decoder does not apply them again.

A new table schema is used. Existing ratio-only artifacts are not silently
upgraded: they do not contain green. The new table and calibration have content
hashes and are checked on loading/use. Calibration binds to the complete table
identity; rebuilding/changing the table requires re-binding and verifying the
reference configuration.

## Two-reference green calibration

The calibration configuration declares an unstressed reference (dose zero) and
a high-dose reference. Each has a specimen ID, biological replicate ID, source,
known reference dose, known reference elapsed time and corrected image green.
Both reference dose/time coordinates must be present in the model table.
The corresponding model green values are retrieved from that table.

For image and model values respectively:

```
green_norm_image = (image_green - image_green_low) / (image_green_high - image_green_low)
green_norm_model = (Observed_Green - model_green_low) / (model_green_high - model_green_low)
```

Anchors use fixed, explicitly recorded reference times. Their selection does
not depend on the unknown target time. The caller supplies strictly positive
minimum acceptable reference spans in image and model units. Reversed or
insufficient spans are rejected. This threshold is an explicit experimental
choice, not an estimated calibration uncertainty.

The anchor values map to 0 and 1. Target observations and model predictions may
fall outside that interval; no clipping is performed. Negative normalized
values are supported and are not logged. A separate explicit green detection
floor uses the original corrected image units.

The high anchor is not automatically a proven biological saturation condition,
and camera saturation must already have been excluded by QC. The artifact
records `high_reference_saturation_validated=false`. Two anchors determine an
affine mapping exactly; they do not test linearity or estimate uncertainty.

## Ratio calibration, acquisition and independence

The configuration embeds a validated existing ratio-calibration artifact and
its protocol reference metadata. Model, profile, initialization, solver,
exposure and f/eps compatibility checks run before decoding. The existing ratio
calibration is still a model-conditioned multiplicative scale. This patch does
not independently validate that scale or implement channel autofluorescence
correction.

Acquisition settings must be recorded in a dictionary containing objective,
exposure_ms_green, exposure_ms_red, green_units, extraction_green,
extraction_red and gain_setting. Additional relevant settings can be included.
The same full dictionary must accompany the target. Values describe the actual
processing/acquisition convention; matching labels do not establish instrument
linearity or physically verify settings.

The target specimen and biological replicate IDs must differ from both green
references and all declared ratio-reference IDs. Both lists of ratio-reference
specimens and biological replicates are required. Reusing references as targets
is rejected even when different specimen names share the same replicate ID.
These checks enforce declared separation, not external proof of independent
cultures. Synthetic reference labels are explicitly synthetic.

## Observation JSON

Each observation represents one QC-filtered specimen summary. Use original
per-cell-ratio and corrected-green summaries from the same accepted population,
with a consistent aggregation method. The caller is responsible for extracting
that summary and satisfying QC; the existing population-summary adapter does
not yet automate this new schema.

Required fields are:

| Field | Meaning |
| --- | --- |
| `session_id` | Session used for reference calibration |
| `specimen_id` | Target specimen ID, distinct from references |
| `biological_replicate_id` | Target culture/replicate ID, distinct from references |
| `data_kind` | `synthetic` or `experimental`, consistent with calibration source |
| `acquisition` | Exact acquisition/processing settings dictionary |
| `ratio_scale` | Must be `image`, preventing an explicitly pre-calibrated input |
| `image_ratio_red_green` | Original corrected image ratio, before session scale |
| `corrected_green` | Green summary in the same units as the reference anchors |

The ratio is divided by c_session once; green is normalized once. There is no
input field for a known target dose/time. The result includes both individual
observable candidate counts and all joint candidates, normalized observations,
reference/table identities, tolerances, floor, and any out-of-anchor-range flag.

## Reproducible synthetic image verification

```bash
python scripts/verify_snapshot_decoder.py --out-dir outputs/two-observable-verification
```

This creates a usable `table.json`, `calibration-config.json`, `calibration.json`,
`observation.json`, `decoded.json`, and `verification.json`, plus reference
protocol exports, image NPZs and extracted cell CSVs. Existing directories are
rejected. A failed recovery check returns a nonzero exit status.

The verification deliberately addresses late oxidative ratio degeneracy:

- Default ox profile, equilibrium initialization, constant stress.
- Table doses [0, 10, 200, 400, 800] uM and times [0.5, 1, 8] hours.
- Fixed low/high green references: 0 and 800 uM, both at 8 hours.
- Separate ratio-reference images at 25 uM and 0.5, 1, 8 hours.
- Every cell's generated green is `20 + 100 * Observed_Green`.
- Generated red is `1.5 * Measured_Ratio_RG * generated_green`.
- Production image segmentation/correction/extraction yields three cells per
  idealized scene. Reference and target values are medians of these cells.
- A fresh simulation supplies the target at 400 uM / 8 hours. Its original
  dose/time is not supplied to the inverse decoder.
- Ratio tolerance 0.005, normalized-green tolerance 0.002, green floor 100;
  minimum image/model anchor spans 100 and 1 respectively. These are declared
  synthetic check settings, not recommended experimental uncertainties.

The observed check retained eight ratio-only grid pairs and one joint pair,
correctly recovering 400 uM / 8 hours. Green varies with the model in these
images, unlike the earlier ratio-only end-to-end scenes with constant green.
No noise is added. Separate IDs and a fresh simulation do not make this an
independent biological validation; all data are generated by the same model.
The three cells are not three biological replicates.

## CLI with retained inputs

The verification's files can be decoded again to a new output:

```bash
python scripts/run_snapshot_decoder.py decode \
  --table outputs/two-observable-verification/table.json \
  --calibration outputs/two-observable-verification/calibration.json \
  --observation outputs/two-observable-verification/observation.json \
  --ratio-tolerance 0.005 --green-tolerance 0.002 --green-floor 100 \
  --out outputs/two-observable-verification/decoded-again.json
```

For another model table use `run_snapshot_decoder.py build --help`; for binding
reference configuration use `calibrate --help`. The emitted configuration is a
complete example. Replace its references, acquisition settings and data-kind
metadata with the actual independently measured references for real work;
synthetic coefficients must not be relabeled as experimental calibration.

## Verification and remaining work

Tests include cases in which each observable separately admits several pairs
but the pair resolves a snapshot, and cases in which both observables remain
ambiguous. They cover reference/target overlap, source/session/settings mismatch,
flat reference spans, missing provenance, fingerprinted artifact changes, floor
abstention, missing joint matches, no clipping, fresh-model recovery, CLI
roundtrips, output protection and legacy forward-result agreement.

Posterior inference, correlated biological-replicate noise, calibration
uncertainty, independent experimental references and full ox/er domain coverage
remain open. The green mapping does not rescue an ER range where both model
observables are insensitive. The existing sample-summary adapter and visual
report still serve the prior ratio-only workflow; this snapshot CLI uses its
own explicit observation/calibration artifacts.
