# Limitations and open items

Things this repository does **not** do, and claims it must not make.

## Not implemented

- There is no decoder. No dose and no elapsed-time estimate can be produced.
- There is no model, simulator or calibration artifact here yet.
- No uncertainty is computed. CRLB figures from `ichnos-fisher` are a local
  theoretical bound under an assumed noise model, not a confidence interval of
  any estimator, and are not reproduced here.

## Known weaknesses carried over from the migrated code

- `ichnos/config.py` still describes an Olympus **SC30** camera. The team's
  confirmed camera is a **QImaging MicroPublisher 3.3 RTV** on an Olympus
  CKX41. Fixed in Step 2; until then the active defaults are wrong for the
  real setup.
- `ImageSet.saturation_value` defaults to 65535 and the CLI does not pass a
  value from the manifest, so 8-bit saturated cells can go unflagged.
- `correct.py` implements a photobleaching correction that `pipeline.py` never
  calls. Exported CSVs are **not** bleaching-corrected.
- Flat fields are supported in the Python API but not exposed through the CLI.
- `timepoint` is an acquisition index, not elapsed biological time.
- `integrated_green` is a sum of pixel intensities, not a reporter
  concentration. The mapping to model observables is unresolved.

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
