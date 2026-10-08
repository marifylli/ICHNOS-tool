# ICHNOS PI Bridge — D3 Scientific Report

**Scope:** deterministic synthetic OX/ER forward modeling and cross-scenario Red/Green nearest-grid analysis. **Validation status:** not experimentally validated.

## Research question

How sensitive is a modeled Red/Green-to-dose/time nearest-grid assignment to physics-informed scenario assumptions, and how does that sensitivity differ between OX and ER?

## Methods

Five scenarios were evaluated: OX `frozen`, `ox_pi_reference`, `ox_pi_systematic`; ER `frozen`, `er_m2_n4`. All used six doses (25–50 µM, step 5 µM), finite 2 h zero-stress preincubation, and times measured after stress onset. D3a/D3b used 0.5 and 1 h; D3c expanded to 0.25, 0.5, 0.75, 1, 1.5, 2, 3 and 4 h.

Each truth-scenario dose/time point was matched to the smallest absolute Red/Green difference among all dose/time candidates under an assumed scenario. The procedure was deterministic, without measurement noise, residual threshold, or experimental calibration. Self-scenario exact matches are expected by construction and are **not** evidence of predictive validity.

## Results

### D3b: two observation times

| Metric | OX cross-scenario | ER cross-scenario |
|---|---:|---:|
| Comparisons | 72 | 24 |
| Nearest-grid time changed | 35/72 (48.6%) | 0/24 (0%) |
| Nearest-grid dose changed | 56/72 (77.8%) | 20/24 (83.3%) |

All cross-scenario assignments were `nearest_only` rather than exact matches or ties.

### D3c: eight observation times

| Metric | OX cross-scenario | ER cross-scenario |
|---|---:|---:|
| Comparisons | 288 | 96 |
| Nearest-grid time changed | 173/288 (60.1%) | 0/96 (0%) |
| Nearest-grid dose changed | 207/288 (71.9%) | 80/96 (83.3%) |

There were 624 comparisons in total including same-scenario assignments (432 OX, 192 ER). All cross-scenario matches were `nearest_only`; same-scenario exact matches were 48/48 per scenario by construction. D3b and D3c percentages are descriptive, **not** a controlled like-for-like improvement comparison: both the truth grid and candidate set changed.

### D3c cross-scenario nearest Red/Green residual ranges

| Variant | Truth → Assumed | Minimum | Maximum |
|---|---|---:|---:|
| OX | frozen → ox_pi_reference | 9.7722063247e-06 | 0.0516443578668 |
| OX | frozen → ox_pi_systematic | 5.34661862582e-05 | 0.071472651941 |
| OX | ox_pi_reference → frozen | 9.7722063247e-06 | 0.0415253930291 |
| OX | ox_pi_reference → ox_pi_systematic | 0.000207253697947 | 0.0198282940742 |
| OX | ox_pi_systematic → frozen | 5.34661862582e-05 | 0.0332217269117 |
| OX | ox_pi_systematic → ox_pi_reference | 0.000207253697947 | 0.0113286407026 |
| ER | frozen → er_m2_n4 | 2.97950109029e-10 | 1.61090249018e-06 |
| ER | er_m2_n4 → frozen | 2.97950109029e-10 | 2.02764950711e-05 |

These are model-predicted absolute Red/Green differences, **not** measured errors or evidence of observational indistinguishability.

### Time identifiability diagnostics on D3c grid

For **known dose**, the smallest same-dose difference between two distinct observation times was:

| Scenario | Minimum absolute Red/Green time separation |
|---|---:|
| OX frozen | 0.000228104942394 |
| OX ox_pi_reference | 0.00142648794373 |
| OX ox_pi_systematic | 0.0015975086525 |
| ER frozen | 0.0130061442382 |
| ER er_m2_n4 | 0.0130085833718 |

When **dose was also unknown**, the smallest cross-time separation across any dose pair was OX frozen `0.000228104942394`, OX reference `6.43216944529e-06`, OX systematic `4.35651387468e-05`, ER frozen `0.0130000364917`, ER M2 n4 `0.0129926759102`.

Strict, zero-tolerance monotonicity of Red/Green with time across eight points:

- OX frozen: 0/6 doses strictly increasing; 6/6 nonmonotonic or flat.
- OX reference: 2/6 strictly increasing; 4/6 nonmonotonic or flat.
- OX systematic: 4/6 strictly increasing; 2/6 nonmonotonic or flat.
- ER frozen and ER M2 n4: 6/6 strictly increasing.

“Nonmonotonic or flat” is a numerical classification with no tolerance for solver fluctuations.

## Interpretation

**OX:** Cross-scenario assumptions materially affect nearest-grid time assignments, and the modeled time trajectories are not uniformly monotonic. Some alternate matches have tiny residuals, while others differ substantially; a changed nearest match is not automatically a plausible explanation of the observation.

**ER:** Cross-scenario nearest-grid time assignments remained unchanged on the tested grid, while dose assignments changed frequently. Cross-scenario residuals were numerically very small. Without an experimental noise model, neither time identifiability nor dose ambiguity can be translated into validated measurement accuracy.

The two variants show different patterns of sensitivity to the assumed physics-informed scenario. These patterns are conditional on the specified models and grids.

## Reproducibility and checks

OX/ER scenario contracts have fixed canonical hashes and validated manifests. D3a and D3c shared-time forward predictions agreed within `rtol=1e-8`, `atol=1e-10`. The final artifact-integrity check passed for all five scenarios, four forward artifacts and two robustness artifacts (156 D3b plus 624 D3c comparisons). This verifies basic structural consistency, **not** full end-to-end integration or experimental validation.

## Limitations and completion statement

This work does not establish physical clearance rates, experimental dose/time decoder accuracy, sensor noise tolerance, camera calibration, or readiness for production use. D0 saturation-mask registration and output atomicity issues remain outside this research prototype's completed scope. No canonical SBML, PULSE science, or other Tool components were changed as part of this PI bridge workflow.

**Completion status:** D1–D3c deterministic synthetic PI bridge analyses and their basic integrity gate are complete. The deliverable is a **reproducible synthetic research prototype with documented limitations**, not an experimentally validated or production-integrated decoder.
