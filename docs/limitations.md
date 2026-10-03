# Limitations and open items

Things this repository does **not** do, and claims it must not make.

## Not implemented

- There is no decoder. No dose and no elapsed-time estimate can be produced.
- There is no model, simulator or calibration artifact here yet.
- No uncertainty is computed. CRLB figures from `ichnos-fisher` are a local
  theoretical bound under an assumed noise model, not a confidence interval of
  any estimator, and are not reproduced here.

## Known weaknesses carried over from the migrated code

- `ImageSet.saturation_value` still defaults to 65535 and the CLI does not
  pass a value from the manifest, so 8-bit saturated cells can go unflagged.
  `image_io.saturation_value_for_file()` now returns the right value (255);
  wiring it through `ImageSet` and the CLI is outstanding.
- `correct.py` implements a photobleaching correction that `pipeline.py` never
  calls. Exported CSVs are **not** bleaching-corrected.
- Flat fields are supported in the Python API but not exposed through the CLI.
- `timepoint` is an acquisition index, not elapsed biological time.
- `integrated_green` is a sum of pixel intensities, not a reporter
  concentration. The mapping to model observables is unresolved.

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

## Open decision found during Step 1

`segment.segment_cells()` calls
`morphology.remove_small_objects(binary, min_size=min_size)`. In
scikit-image 0.26 `min_size` is deprecated in favour of `max_size`, and the
boundary semantics changed with it:

| version | `min_size=N` removes |
| --- | --- |
| < 0.26 | objects strictly smaller than N (an N-pixel object is **kept**) |
| 0.26 | objects of N pixels or fewer (an N-pixel object is **dropped**) |

Verified on 0.26: a single 5-pixel object is removed by both `min_size=5` and
`max_size=5`, and kept by `max_size=4`.

So the segmentation already changed by one pixel of threshold when
scikit-image was upgraded, silently. The call has deliberately **not** been
migrated yet, because the two possible replacements are not equivalent:

- `max_size=min_size` preserves today's (0.26) behaviour;
- `max_size=min_size - 1` preserves the original intent of the `min_size`
  argument ("keep objects of at least this size").

This needs a decision, a pinned scikit-image version, and a regression test
on a fixed scene before it is changed. Until then the deprecation warning is
left in place as the visible reminder.

## Licensing

`ICHNOS-ablation` is MIT. `DryLabTool`, `Ichnos_PULSE` and `ichnos-fisher`
carry no licence file at the frozen commits. Permission or an upstream licence
is required before this repository is published.
