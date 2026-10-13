# Architecture

Two Python packages and a set of command-line scripts. Every script refuses
to overwrite an existing output and writes what it was run with (`run.json`
or a manifest next to its outputs). Every step run on real data is recorded,
PASS or FAIL, in [run_log.md](run_log.md).

| Package | What it holds |
| --- | --- |
| `ichnos` | The model: SBML sensing modules for ox (H2O2) and er (DTT), the merge into one circuit (`build`), parameter profiles (`params`), stress protocols (constant or first-order clearance, `protocol`), simulation (`simulate`), model fluorescence (`fluorescence`), and the model-based decoders (dose tables, snapshot compatibility, conditional grid posterior, `workflow`). |
| `ichnos_image` | Images to numbers: acquisition log and manifest, segmentation, bright-field alignment, per-cell extraction and QC, exposure normalisation, and the analysis of real data (time courses, induction test, model comparison, single-field decoder). |

## Route used for the October 2026 data

Each stage reads the previous stage's file. Commands are as run (see
run_log.md for dates, commits and results).

| # | Stage | Script (module) | Output |
| --- | --- | --- | --- |
| 1 | Wet-lab log → one row per imaged field | `build_acquisition_log.py` (`acquisition_log`) | `data/acquisition/acquisition_log.csv` |
| 2 | Match TIFFs to log rows by time | `build_image_manifest.py` (`acquisition_match`) | `data/acquisition/manifest/images.csv` |
| 3 | Find each field's bright-field frame and its shift onto the fluorescence | `align_brightfield.py --search-up 2 --write-manifest …` (`align`) | `outputs/align_20261011_v2/images_aligned.csv` |
| 4 | Segment on bright-field (Cellpose), extract per-cell green/red, QC | `run_pipeline_batched.py … -- --segmentation-method cellpose --segmentation-source brightfield --cellpose-gpu --bleed-default 0.05 --saturation-value 255` (`segment`, `pipeline`, `extract`, `focus`) | `outputs/full_20261012_cellpose/cells_all.csv` |
| 5 | Grey levels → grey levels per second | `normalize_exposure.py` (`exposure`) | `cells_all_per_s.csv` |
| 6 | Time courses, autofluorescence baseline, peak times | `fluorescence_timecourse.py` (`timecourse`) | `timecourse/` |
| 7 | Each dose against the same day's 0 µM, with detection limits | `detect_induction.py` (`induction`) | `induction/` |
| 8 | Model prediction against the measured fold over 0 µM | `compare_model.py` (`model_compare`, `ichnos.protocol`) | `model_compare/` |
| 9 | Dose from a single field, leave-one-day-out | `evaluate_decoder_lodo.py` (`lodo_decoder`) | `decoder_lodo/` |

Stages 1–5 take hours (stage 4: ~2 min per field on the Apple GPU, ~11 h for
294 fields) and are run one by one. Stages 6–9 run from one configuration:

```bash
python3 scripts/run_analysis.py --config examples/analysis_202610.json --out-dir outputs/analysis_<date>
```

which writes the four sub-folders and `run-manifest.json` (configuration,
input sha256, git commit, sha256 of every output). About 2–3 minutes.

Decoding a new culture: `decode_field.py` (`field_decoder`) trains the
paired or absolute dose decoder and the timer-based time-window decoder on
every reference day, decodes the new fields, and attaches to each answer a
verdict and that decoder's leave-one-day-out record
([user_guide.md](user_guide.md)). `make_worked_example.py` turns a reference
day into a "new culture" with hidden doses to try it.

The ichnos_vail application: `app/ichnos_vail.py` (Streamlit) over the `ichnos_vail` package. It
builds the scripts' own command lines (`commands`), runs them as detached
background jobs with status and log on disk (`jobs`, `chain`), checks a
manifest before a run (`manifest_check`) and displays the files the scripts
wrote (`results`). It computes nothing itself, so the scripts' tests cover
every result it shows.

Supporting tools: `score_segmentation.py` scores any segmentation method
against the hand-marked cells in `data/segmentation_reference/`;
`check_exposure_linearity.py` tests the camera background against exposure.

## Choices that shape the route

- **Masks from bright-field, not fluorescence.** Fluorescence-based masks
  scored recall 0.44 / precision 0.31 against hand-marked cells; Cellpose on
  bright-field 0.91 / 0.90. A mask is then independent of how bright the
  reporter is.
- **Unit of analysis.** A field is (session, condition, sample); `sample_id`
  repeats across time points. Fields of one sample and time are not
  independent cultures; they are pooled, never counted as replicates.
- **Comparisons within a day.** Each day used a different transformation and
  the non-transfected cultures were imaged on other days, so every effect is
  measured against the same day's 0 µM (stages 7–9).
- **Validation by day.** Anything fitted from data is fitted on other days
  and scored on the held-out day (decision of 2026-10-12, run_log.md).

## The model-based workflow (`run_workflow.py`)

`scripts/run_workflow.py` (`ichnos.workflow`) is the route the repository was
built around: image manifest → cells → paired population summaries →
snapshot compatibility or a conditional grid posterior over dose × time from
the model. Configuration format: [integrated_workflow.md](integrated_workflow.md),
example: `examples/workflow.json`.

It is **not** applied to the October data, for reasons that are properties of
the data, not of the code:

- In experimental mode it requires green→red crosstalk controls and forbids a
  default bleed-through value. No controls were imaged; 0.05 is indicative.
- It requires a session calibration from reference samples and a replicate
  covariance from independent cultures (≥ 3 per condition). Neither exists.
- Stages 8–9 show the uncalibrated model overpredicts the H2O2 response 1.7–5×
  and predicts none at the DTT doses used; a decoder built on it was
  confidently wrong on held-out days (run_log.md, 2026-10-13 08:07 and 08:27).

It stays in the repository, with its synthetic verifications, for data that
meet those requirements.
