# ICHNOS Physics-Informed Bridge

**Status:** completed deterministic synthetic research prototype; **not** experimentally validated or production-integrated.

## Scope and boundaries

This isolated `pi_bridge/` workspace connects documented OX/ER physics-informed scenario assumptions to forward-model predictions and deterministic cross-scenario nearest-grid analyses. It does not modify canonical SBML, the underlying PULSE science, existing ICHNOS Tool code, or colleagues' sensitivity/Fisher/ablation analyses.

The modeled observable is Red/Green; model-generated Green is also retained. These values are **not** calibrated camera intensities. The effective `k_clear` parameter represents decay of sensor drive, not independently established physical H2O2 or DTT clearance.

## Frozen scenario contracts

| Variant | Scenario | Effective `k_clear` (h^-1) | Kinetic overrides |
|---|---|---:|---|
| OX | `frozen` | 0 | None |
| OX | `ox_pi_reference` | 1.95 | None |
| OX | `ox_pi_systematic` | 3.30 | None |
| ER | `frozen` | 0 | None |
| ER | `er_m2_n4` | 0.5032 | `K_act_er=929.5`, `n_er=4.0`, `k_on_er=3.474`, `k_off_er=13.674`, `d_x_er=1.755` |

Contract SHA-256 digests (canonical JSON):

- OX: `0e8ec503b0bdc84e4addab4f26e93a266d4e9b746a3fbb47701e7eb649c5b261`
- ER: `6e54921440af799f178275fc0d856f52d6fa7425e5bdded15fb24c0cc4e1f163`

Source scenario definition: `Ichnos_PULSE/python/pi_scenarios.py` at commit `f4f07176ddf8fce8f0af3c8374b4811b09c9e297` (as recorded during bridge development).

## Reproduction from repository root

Use the established Python environment with the ICHNOS Tool dependencies installed. Run in this order:

```powershell
python pi_bridge/scenarios/prepare_pi_bundles.py
python pi_bridge/forward/build_ox_pi_forward_tables.py
python pi_bridge/forward/build_er_pi_forward_tables.py
python pi_bridge/robustness/evaluate_cross_scenario.py
python pi_bridge/forward/build_d3c_extended_tables.py
python pi_bridge/robustness/evaluate_d3c_cross_scenario.py
```

**Important:** Some generators deliberately refuse to overwrite existing output files. On a clean checkout, create the missing artifacts once. For repeat runs, use a clean, separately backed-up workspace; do not indiscriminately delete existing scientific outputs. The scenario contract generator also refuses to overwrite divergent existing contracts.

### Grid and artifact inventory

- D3a: six doses `[25,30,35,40,45,50]` µM × times `[0.5,1.0]` h, for three OX and two ER scenarios.
- D3c: the same six doses × `[0.25,0.5,0.75,1.0,1.5,2.0,3.0,4.0]` h.
- Both use finite 2 h zero-stress preincubation; times are measured after stress onset. No convergence to steady state is asserted.
- D3a output: `pi_bridge/results/{ox,er}_pi_forward_tables.json`.
- D3b output: `pi_bridge/results/d3b_cross_scenario_v1.json`.
- D3c forward output: `pi_bridge/results/{ox,er}_pi_forward_tables_d3c.json`.
- D3c robustness output: `pi_bridge/results/d3c_cross_scenario_v1.json`.

OX and ER output schemas differ: OX uses `ratio_red_green` and `observed_green`, while ER uses `ratios` and `greens`. Each scenario has a six-by-time-point matrix for each observable.

The cross-scenario evaluator treats each truth-scenario dose/time prediction as an observation and chooses the candidate with minimum absolute Red/Green residual from an assumed scenario's grid. It uses exact floating-point equality for ties and **does not** introduce a noise model, acceptance tolerance, or validated decoder.

## Verification completed

- Scenario contract generation, manifest construction, and validation succeeded with expected hashes.
- Forward tables were generated successfully for all five scenarios.
- D3a/D3c values at shared times (0.5 and 1 h) agreed within `rtol=1e-8`, `atol=1e-10` for both modeled observables.
- Final read-only artifact gate passed: files exist, expected scenario IDs are present, all Red/Green and Green matrices have expected dimensions and finite values, robustness comparison counts are correct, residuals are finite/nonnegative, explicit `experimentally_validated=false` flags are present, and no tracked-file changes were outstanding.
- The final gate is **structural/basic integrity**, not a full integration, scientific validation, or independent provenance audit.

## Known limitations and follow-up

1. **No experimental validation:** all D3a–D3c results are deterministic model predictions, not measured decoder performance.
2. **No measurement noise or uncertainty model:** very small residuals must not be labeled indistinguishable without a justified experimental error model.
3. **Grid dependence:** nearest matches, dose/time changes, and apparent identifiability depend on sampled doses, times, and scenarios.
4. **Effective-parameter semantics:** `k_clear` is a sensor-drive decay assumption, not a direct biochemical clearance measurement.
5. **Known D0 issues deferred:** saturation-mask/red-image registration and output atomicity require separate engineering work before broader production use.
6. **No verified production integration:** the PI bridge is an isolated research prototype, not a replacement for Tool decoding or an end-to-end validated pipeline.

See `D3_SCIENTIFIC_REPORT.md` for quantitative results and interpretation.
