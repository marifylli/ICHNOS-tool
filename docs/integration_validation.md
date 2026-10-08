# Integration verification — 2026-10-07

Baseline reviewed: `main` at `491fab6`. No experimental image dataset was present
in the checkout. Tests and benchmarks below establish computational behavior,
not biological accuracy.

## Verification performed

- Saturation tests: integer/fractional/zero registration, interpolation support,
  newly clipped and no-longer-clipped cells, and isolated QC impact accounting.
- Storage tests: refusal to overwrite, failure during processing, failure during
  manifest publication, preservation of earlier CSVs, and output hash matching.
- Focus tests: score behavior with defocus/gain, CNR guards, indeterminate states,
  extraction/export propagation and enforced QC rejection.
- Adapter tests: paired medians, identity/acquisition/estimator contracts and
  equality with direct snapshot/posterior API results.
- Real CLI test: on-disk TIFFs and manifests through segmentation, summaries,
  posterior and HTML/JSON reporting; rejected samples retained; overwrite refused.
- ox/er benchmark tests: on/off-grid targets, noise, calibration shifts and nested
  refinement. Off-grid membership is explicitly separated from discrete coverage.

The final complete suite passed **494 tests with one optional Cellpose skip**,
with FutureWarnings treated as errors (`python -m pytest -q -W error::FutureWarning`).
`git diff --check` also passed.

## Robustness run

Command:

```bash
python scripts/evaluate_robustness.py --out-dir outputs/robustness-integration \
  --trials 12 --refine
```

Seed: 20261007. Variants: ox and er. Grids: 9×17 and nested 17×33. Domain:
0–800 µM, 0–8 h. Three noise scales (0.5, 1, 2); seven calibration scenarios.
Twelve on-grid and twelve off-grid truths per variant, reused across scenarios
and grids. There were **2,016 posterior evaluations**, not 2,016 independent
biological observations. The scenarios share latent noise to support paired
comparisons. Thirty-six independent synthetic training labels per variant/noise
level fit the covariance separately from the targets.

Selected results below use the refined grid, baseline noise scale and 12 off-grid
targets per row. Dose/time errors are mean absolute MAP errors. These selected
rows had no MAP ties. "Cell inclusion" refers to the midpoint cell containing
the continuous truth, **not continuous 95% coverage**.

| Variant | Target calibration scenario | Dose error (µM) | Time error (h) | Truth-cell inclusion | Poor-fit flag |
| --- | --- | ---: | ---: | ---: | ---: |
| ox | matched | 166.4 | 1.48 | 11/12 | 0/12 |
| ox | ratio gain −10% | 179.0 | 2.70 | 0/12 | 0/12 |
| ox | green offset +0.1 reference span | 188.5 | 0.99 | 5/12 | 0/12 |
| er | matched | 98.0 | 2.75 | 12/12 | 0/12 |
| er | ratio gain +10% | 97.3 | 3.91 | 8/12 | 5/12 |
| er | green offset +0.1 reference span | 103.5 | 2.66 | 5/12 | 0/12 |

These small same-model synthetic checks illustrate that calibration errors can
exclude the truth without triggering the explicit Mahalanobis threshold. They
do not calibrate that threshold or establish error rates for the microscope.
Matched cell inclusion also does not imply accurate point estimates or narrow
credible regions. Inspect the full `evaluation.json`/`trials.csv`, including
credible area, all error bounds and all scenarios, before interpreting recovery.

## Still requiring real data

The historical effect of the saturation fix is unknown until the original
images/manifests are reprocessed. Acquisition-specific focus/CNR thresholds,
calibration references, noise-training cultures and held-out experimental targets
must be supplied before judging experimental performance. No numerical focus
threshold or calibration artifact has been invented for the real instrument.
