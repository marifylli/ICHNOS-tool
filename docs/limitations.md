# Limitations and open items

Things this repository does **not** do, and claims it must not make.

## Not implemented

- No decoder exists. Experimental images cannot yet be converted into
  dose or elapsed-time estimates.
- Measurement mapping, complete session calibration and a decoder
  calibration artifact remain unfinished.
- Optional clearance has not yet been implemented in the shared simulator.
  The current constant-stress baseline is a model assumption, not a
  measurement of the experimental exposure history.
- Initial model states have not yet been validated against the experimental
  galactose pre-incubation procedure.
- No estimator uncertainty is computed. CRLB results from `ichnos-fisher`
  are theoretical bounds under assumed measurement conditions, not
  confidence intervals for an implemented estimator.
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
  concentration. The mapping to model observables is unresolved.

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

## Licensing

`ICHNOS-ablation` is MIT. `DryLabTool`, `Ichnos_PULSE` and `ichnos-fisher`
carry no licence file at the frozen commits. Permission or an upstream licence
is required before this repository is published.
