# Limitations and open items

Things this repository does **not** do, and claims it must not make.

## Implemented scope and unresolved inference

- Ratio-only dose decoders and a discrete joint dose/time decoder exist.
  A single snapshot is accepted, but its ratio alone generally does not identify
  dose and time. Multiple observations require known relative spacing.
- A separate two-observable snapshot decoder now uses ratio and two-reference
  normalized green for grid compatibility. It requires explicit observation
  and reference artifacts; the older sample adapter remains ratio-only; the new integrated workflow uses paired summaries and explicit acquisition/replicate metadata.
  A separate conditional posterior uses fitted biological-replicate covariance.
  Independent experimental calibration and empirical coverage validation remain open.
- Relative ratio calibration assumes a zero-intercept multiplicative mapping
  to model reference values. This is model-conditioned normalization, not
  independent validation of the model or camera. Two-reference green
  normalization is implemented separately; its experimental validation remains open.
- The noise-free synthetic end-to-end workflow and visual reports exist; they
  test a limited oxidative scenario and do not establish experimental accuracy.
- Background subtraction is called by the image pipeline, but this does not
  establish subtraction of cellular autofluorescence or validate zero intercept
  in the ratio calibration. An affine diagnostic is available separately.
- Optional first-order clearance is implemented, but its rate has not been experimentally calibrated. Constant stress remains an available baseline.
  The current constant-stress baseline is a model assumption, not a
  measurement of the experimental exposure history.
- Initial model states have not yet been validated against the experimental
  galactose pre-incubation procedure.
- Conditional discrete posterior credible sets are computed with a fixed Gaussian
  covariance estimated from independent biological summaries. Calibration, covariance
  estimation and model-parameter uncertainty are not propagated. CRLB results
  from `ichnos-fisher` remain theoretical bounds, not empirical coverage validation.
- No copper/CuSO4-specific model or decoder is included.

## Known weaknesses carried over from the migrated code

- The CLI uses the instrument-profile saturation threshold unless
  --saturation-value overrides it. Saturation is checked before RGB
  extraction and passed to per-cell QC as a raw saturation mask.
  The Python API still defaults ImageSet.saturation_value to 65535 for
  compatibility with existing 16-bit workflows. API callers must supply
  the appropriate threshold or a raw saturation mask for their images.
- `correct.py` implements a photobleaching correction that `pipeline.py` never
  calls. Exported CSVs are **not** bleaching-corrected.
- Flat fields are supported in the Python API but not exposed through the CLI.
- `timepoint` is an acquisition index, not elapsed biological time.
- `integrated_green` is a sum of pixel intensities, not a reporter
  concentration. Its experimental mapping to model concentration remains unresolved.

## Exposure model and clearance decision

The simulator supports two explicit exposure assumptions:

- Constant stress: clears=False, with no clearance rate supplied.
- First-order decay: clears=True, with an explicitly supplied positive,
  finite clearance_rate_per_hour in h^-1.

For first-order decay, dS/dt = -k_clear * S, so
S(t) = S0 * exp(-k_clear * t). Dose is expressed in uM and numerical
simulation time in hours.

Clearance is optional because the experiment measures fluorescence,
not stress concentration over time. A single initial addition does not
determine the subsequent concentration trajectory or its decay rate.
Fluorescence dynamics alone do not establish a unique clearance rate.

No numerical clearance default is supplied: an arbitrary value would
change predictions without experimental support. The constant-stress
baseline is retained for comparison; it is also an assumption, not
evidence that stress remains constant.

Existing parameter profiles are not automatically recalibrated when
clearance is enabled. Predictions under clearance remain conditional on
the selected rate, parameter profile and observation model.

load_protocol_model prepares a fresh model for the requested exposure.
Protocol execution rejects a loaded stress law that disagrees with the
requested configuration. Simulation results record exposure metadata.

Observation times are supplied explicitly in hours after stress onset.
The simulation includes the interval from onset to the first observation,
even when the first requested observation occurs later than zero.

Synthetic clearance rates used in tests verify implementation behaviour;
they are not recommended experimental values.


## Experimental data and model alignment

- The experiment specifies an initial stress addition followed by repeated
  sampling from the same culture tubes. It does not directly measure
  stress concentration over time.
- A fluorescence decline alone does not establish a clearance rate.
  Clearance must remain an explicit modelling assumption until its
  parameter is constrained under a stated observation model.
- Record actual stress-addition, sampling and image-acquisition times.
  The current timepoint field is an acquisition index, not a replacement
  for elapsed biological time.
- Record initial dose and model variant explicitly when linking image
  results to simulations.
- Use a consistent ratio convention: the current image code computes
  corrected red/green, corresponding to mCherry/GFP.
- Repeated sampling provides culture-level time-course observations;
  it does not track the same individual cells over time.
- Additional observation times do not by themselves guarantee unique
  parameter estimates. Identifiability depends on the selected parameters,
  model observables and measurement uncertainty.

## Found in Step 4

**Assignment-rule observables are absent from a default simulation.**
RoadRunner's default selections are floating species and rate-rule states.
`Observed_Green`, `Measured_Ratio_RG`, `Ratio_RG_FRET`, `Total_red_pool` and
`b_fret` are assignment-rule parameters -- between them the entire measured
readout of this circuit -- and do not appear in the output of a plain
`simulate()`. Code looking for them finds nothing. `ichnos.simulate.load_model`
now widens the selections when given the model.

**The default `S_er` is effectively zero stress.** The ER sensing module ships
with `S_er = 100 µM` against `K_act_er = 2345.3 µM`, `n_er = 3.97`. Measured at
12 h after onset:

| `S_er` (µM) | `Observed_Green` |
| --- | --- |
| 0 | 2.3570 |
| 100 (default) | 2.3626 |
| 500 | 2.6000 |
| 1000 | 3.2200 |
| 2345 | 5.5234 |
| 5000 | 6.6161 |
| 10000 | 6.7040 |

The default dose moves the readout by 0.24% against no stress at all. Any ER
run at the shipped default is a near-baseline run. Calibration must span a
dose range around `K_act_er`, not around the default.

**The oxidative transient is fast.** `A_ox` peaks at 0.8627 about 4 minutes
after onset and settles to 0.3052. A 200 h run over 500 points samples every
0.4 h and misses it entirely. `simulate.peak_summary()` returns an
`undersampled` flag rather than printing a warning.

## Corrected in Step 2

The migrated pipeline described an Olympus **SC30** camera. That is not this
lab's camera. The confirmed setup is a **QImaging MicroPublisher 3.3 RTV**
(3.45 µm sensor pixel, 2048×1536, 10-bit ADC, cooled, colour) behind an
Olympus **U-TV0.5XC-3** 0.5X adapter on an Olympus **CKX41**.

Two concrete consequences, both now fixed:

- `pixel_size_at_sample_um()` defaulted to a 3.2 µm sensor pixel and a 1.0X
  adapter. The real figures are 3.45 µm and 0.5X, so every default-path
  spatial scale was wrong by a factor of 2.16.
- The objective pixel-size table held `60X: 0.1076` and `100X: 0.1300`,
  measured on an unrelated public dataset's microscope. The real values are
  4X 1.725, 10X 0.690, 40X 0.1725, 60X 0.115 µm/px. This turret has no 100X,
  and the 40X the dry lab proposes for fluorescence was missing entirely, so
  `rolling_ball_radius_for_objective("40X")` raised while `"100X"` returned a
  plausible-looking radius for an objective that does not exist.

SC30 constants were removed rather than kept as an inactive profile: this lab
does not own that camera. They remain in git history and in DryLabTool.

## Still unresolved after Step 2

Hardware is identified; measurement calibration is not. Tracked in
`ichnos_image.instrument.UNRESOLVED`, and
`instrument.is_quantitatively_calibrated()` returns False while any entry is
open. The ones that block specific claims:

- **RGB extraction.** The camera is colour; each acquisition is a 24-bit RGB
  file, not a single-channel measurement. Which component or calibrated
  combination represents the signal per cube is undecided, so
  `image_io.extract_fluorescence_plane()` requires an explicit method and has
  no default. A generic grayscale conversion must never become the silent
  fallback.
- **Post-snap macro.** `CAM_APPLY_LUT` is confirmed on preset 3 (phase
  contrast) and unphotographed on preset 6 (fluorescence). One intensity
  value in five is absent from the sample fluorescence file, which is a
  rescaling signature. `gamma = 1` is therefore not sufficient to call these
  files raw linear measurements.
- **Clipped background.** `Offset = 0` puts 13% of red and 43% of blue pixels
  at exactly zero. Background estimates are biased upward and read noise is
  not symmetric about the background level.
- **Exposure.** "Adjust Exp for Binning" makes capture exposure 4× the
  preview exposure (100/400 ms and 300/1200 ms both observed in one preset).
  The capture value is written into the TIFF; the pipeline reads it and must
  not apply the factor again.
- **Channel name.** Image-Pro Plus writes `Alexa Fluor 560` into every file
  regardless of cube. It must never be used to identify a channel.
- **Fluorescence objective, lamp intensity, burner hours, per-channel
  exposure and gain, D460/50M placement, dark frames, PTC, flat fields,
  single-colour controls.** All open. Recorded as `None`, never as a default.

## Resolved morphology compatibility

Deprecated binary morphology calls have been replaced with erosion,
dilation and closing using mode="ignore". This retains border handling
for the symmetric, odd-sized disk footprints used by the pipeline.
The image extra requires scikit-image>=0.23.

The minimum-object-size rule is now explicit: an object containing at
least min_size pixels is retained.

The compatibility helper uses max_size=min_size-1 when the new API is
available, and min_size on older supported versions.

This restores the intended minimum-size rule. Compared with the observed
deprecated-call behaviour in scikit-image 0.26.0, objects exactly at the
threshold may now be retained.

Regression tests cover objects below, at and above the threshold.
The complete test suite passes with FutureWarnings treated as errors.

## Real-data analysis, October 2026

What the analysis of the seven imaged sessions (H2O2 transfected 5/10, 6/10,
7/10; DTT transfected 6/10, 7/10; non-transfected H2O2 8/10 and DTT 9/10)
can and cannot support. Numbers and commands: [run_log.md](run_log.md).
No further experiments will be run, so these limits are final for this data set.

### Measurement

- **Bleed-through** of green into red is set to 0.05, an indicative value:
  no single-colour controls were imaged.
- **Exposure linearity.** The green background rises linearly with exposure
  (pooled R² 0.70), which supports per-second normalisation of green. Red
  was never imaged at several exposures; its normalisation is assumed.
- **No flat-field correction.** Uneven illumination across the frame adds
  field-to-field noise. It does not bias within-day comparisons unless doses
  were imaged systematically at different positions.
- **Post-snap processing** (a possible look-up table on fluorescence frames)
  and the constant red background of ~50–60 grey levels remain unexplained
  (`instrument.UNRESOLVED`).
- **Signal size.** Median background-corrected cell signal is a few grey
  levels on an 8-bit scale, so single-cell ratios are noisy.

### Segmentation

- Cellpose on bright-field, default settings: recall 0.91, precision 0.90
  against 212 cells marked by one person on 7 crops. Points score detection
  and separation, not outline accuracy. Small dark cells are the main misses.
- 6 of 294 fields have no matching bright-field frame and are not analysed.
  6 dense fields were processed with the foreground cap raised from 20 % to
  50 % of the frame.

### Design

- Each day used a different transformation, and the non-transfected cultures
  were imaged on other days. Every effect is therefore measured within one
  day against that day's 0 µM. The autofluorescence baseline is one constant
  per stressor, not day-matched.
- Fields of one sample are not independent cultures, so the number of
  independent units is small: 3 days for H2O2 and 2 for DTT, one culture
  per dose and day. On 7/10 H2O2 only 2 stressed time points lie inside the
  0 µM time range.

### Findings that limit claims

- **Induction.** 12 of 66 within-day comparisons are significant (q < 0.05).
  The only dose-ordered response is H2O2 on 7/10: more reporter-positive cells
  (+14 percentage points at 600 µM) and higher population red. On that day the
  response appeared as an increase in the share of reporter-positive cells
  rather than in their brightness. This rests on one day, on a positive
  threshold taken from non-transfected cells imaged on another day at another
  exposure, and on the assumed exposure normalisation of red; it is a
  hypothesis for the next cycle. On 5/10 green fell with dose. Other days show
  no dose-ordered response; their detection limits are in `effects.csv`.
- **Peak time** cannot be given. On 7/10 H2O2 the 600 µM red signal is still
  rising at the last sample (7.6 h).
- **Model.** The default, uncalibrated model overpredicts the H2O2 response of
  reporter-positive cells 1.7–5× at ≥ 75 µM. Assumed fast H2O2 clearance
  narrows the gap but does not close it. At the DTT doses used (≤ 200 µM) the
  model predicts no response (K_act ≈ 2.3 mM). The absence of a DTT response
  is consistent with that range, but a null result cannot confirm it: the
  reporter did not respond to H2O2 on two of three days either, and K_act_er
  comes from a four-point fit, two points of which were read from the text of
  the source paper.
  Clearance rates are assumptions, not measurements.
- **Decoding.** Leave-one-day-out, no decoder recovers the dose of a single
  field better than chance, with or without a same-day 0 µM reference. The
  empirical decoders mostly stay near the prior, but they become overconfident
  when the held-out day lies outside the spread of the training days (both
  7/10 days). With one training day (DTT) the between-day spread cannot be
  estimated at all.
- **Time since the stress.** The tool decodes it into windows from the tandem
  timer, but in these data the timer ratio of reporter-positive cells does not
  change with time since the stress (slopes −0.08 to +0.02 per hour, p ≥ 0.12
  on every day). Leave-one-day-out, the time decoder scores accuracy 0.19
  against chance 0.31 for H2O2 and 0.20 against 0.28 for DTT. These data do
  not test the timer: no stress-induced burst of synthesis was detected. Basal
  promoter activity and day-dependent transformation are alternative
  explanations we cannot separate.
- **The model-based decoder must not be offered to users** in its current
  state. For H2O2 it puts every field at 0–75 µM with confidence (−22 bits
  against the prior). For DTT it returns the prior.

### What would remove these limits

Single-colour controls (bleed-through); one field imaged at several exposures
(red linearity); flat-field frames; at least 3–4 independent days per
stressor, each with its own 0 µM; DTT doses around the model's activation
range (≥ 1 mM) or a recalibrated ER module; a measurement of how fast H2O2
disappears from the medium; reference samples for a fluorescence-to-model
calibration.

## Licensing

`ICHNOS-ablation` is MIT. `DryLabTool`, `Ichnos_PULSE` and `ichnos-fisher`
carry no licence file at the frozen commits. Permission or an upstream licence
is required before this repository is published.
