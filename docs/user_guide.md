# User guide

From microscope images of a yeast culture carrying a red/green fluorescent
timer reporter to two answers:

1. **Does my reporter respond to the stress, and how large a response could
   my experiment have detected?** (analysis of a reference experiment)
2. **Which dose was this new culture exposed to, and how long ago?**
   (decoding new images against that reference, with a verdict for each
   answer saying whether it can be trusted)

The tool runs from the command line (Python 3.10 or later, macOS or Linux).
Steps 1–3 below take hours for a full experiment; everything after takes
minutes.

## What it will and will not tell you

- It **will** segment cells on bright-field, measure each cell's green and
  red, remove autofluorescence, and test every dose against the same day's
  unstressed culture, with confidence intervals and detection limits.
- It **will** give a dose for a new image only when the reference experiment
  supports one. Every answer comes with a verdict and the decoder's own
  score on held-out reference days.
- It **will not** give a dose from the ICHNOS model alone. That model is
  uncalibrated; on our data it was confidently wrong (docs/limitations.md).
- With the team's October 2026 data as reference, no new image gets the
  verdict "dose estimate": those data contain no dose signal that carries over
  to a new day. That is a property of the data, and the tool says so.

## The ichnos_vail application

Everything below can also be done from **ichnos_vail**, a local app in the browser:

```bash
pip install -e ".[ichnos_vail]"
streamlit run app/ichnos_vail.py        # macOS, for long runs: caffeinate -i streamlit run app/ichnos_vail.py
```

Pages follow the steps of this guide (check images, images → cells, analyse,
decode a new culture, jobs). Each page shows the exact command it is about to
run; long runs continue in the background if the browser is closed and are
found again under **Jobs** (kept in `outputs/ichnos_vail_jobs/`). The app is served to
this computer only (`.streamlit/config.toml`). It computes nothing itself:
results are those of the scripts below.

## 0. Install

```bash
git clone https://github.com/marifylli/ICHNOS-tool && cd ICHNOS-tool
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[image,model,dev]"
pip install cellpose          # bright-field segmentation; downloads a 1.15 GB model on first use
python3 -m pytest -q          # optional check; the Cellpose tests need that download
```

A GPU helps: Cellpose takes ~2 minutes per frame on an Apple GPU
(`--cellpose-gpu`), much longer on a CPU.

## 1. Plan the imaging

| You need | Why |
| --- | --- |
| The same reporter, objective, camera and settings for reference and new images | The baseline, thresholds and dose classes come from the reference |
| Per field: green, red and bright-field frames of the same position | Cells are found on bright-field, measured in green and red |
| An unstressed (0 µM) culture imaged on **every** day | All comparisons are made within a day; the paired decoder needs it |
| A non-transfected culture of each stressor, at least once | Autofluorescence baseline and the "reporter-positive" threshold |
| The time from adding the stress to each image | Time axis of every analysis |
| Reference: several doses on **at least 3 independent days** | Validation holds each day out; with 1–2 days the between-day spread cannot be estimated |
| Bleed-through controls (green-only cells imaged in the red channel) | Otherwise an indicative 0.05 is assumed |

## 2. Describe the images: the manifest

One CSV row per field (full column list: [data_contract.md](data_contract.md)):

```csv
session_id,sample_id,condition_id,timepoint,acquisition_order,green_path,red_path,bright_field_path,exposure_ms_green,exposure_ms_red,nd_filter_green,nd_filter_red,objective,burner_hours,lamp_warmup_minutes,measurement_time_hours
20261101_h2o2_tr,ctrl_2h_R1,h2o2_0uM_2h,0,0,img/g001.tif,img/r001.tif,img/bf001.tif,400,400,,,40X,,,2.0
20261101_h2o2_tr,d300_2h_R1,h2o2_300uM_2h,0,1,img/g002.tif,img/r002.tif,img/bf002.tif,400,400,,,40X,,,2.1
20261103_h2o2_nt,nt_2h_R1,h2o2_0uM_2h,0,0,img/g101.tif,img/r101.tif,img/bf101.tif,800,800,,,40X,,,2.0
```

- `session_id` = `<yyyymmdd>_<stressor>_<tr|nt>`; `condition_id` =
  `<stressor>_<dose>uM_<label>` for known doses.
- For a **new culture of unknown dose**, use any other condition name
  (`mystery_A`); keep `h2o2_0uM_<label>` for its same-day control.
- Paths are relative to the manifest's folder. Blank cells mean "not
  recorded", never zero.

If bright-field frames sit in separate folders or are shifted relative to
the fluorescence, let the tool find and align them (writes a new manifest):

```bash
python3 scripts/align_brightfield.py --manifest images.csv --out align/align.csv \
    --write-manifest align/images_aligned.csv --search-up 2
```

## 3. Images to cells

```bash
python3 scripts/run_pipeline_batched.py --manifest align/images_aligned.csv --out-dir run1 --batch-size 20 -- \
    --segmentation-method cellpose --segmentation-source brightfield --cellpose-gpu \
    --bleed-default 0.05 --saturation-value 255
python3 scripts/normalize_exposure.py --cells run1/cells.csv --out run1/cells_per_s.csv
```

- Use `--controls controls.csv` instead of `--bleed-default` if you imaged
  bleed-through controls.
- For colour (RGB) files add `--green-extraction G --red-extraction R`.
- `--saturation-value` is the camera's maximum (255 for 8-bit, 4095 for
  12-bit, …).
- Very dense slides may be refused by the foreground check; rerun those rows
  with `--max-foreground-fraction 0.5`.
- Check segmentation on your own images once: mark cells on a few crops
  and run `scripts/score_segmentation.py` (see `data/segmentation_reference/README.md`).

## 4. Analyse a reference experiment

Copy `examples/analysis_202610.json`, point `cells` at your
`cells_per_s.csv`, and run:

```bash
python3 scripts/run_analysis.py --config my_analysis.json --out-dir analysis1
```

| Look at | It tells you |
| --- | --- |
| `induction/effects.csv` | Each dose against the same day's 0 µM: effect, 95 % interval, q-value; `mde` = smallest effect this experiment could have detected |
| `induction/trend.csv` | Whether the effect grows with dose |
| `timecourse/timecourse_<stressor>.png`, `peaks.csv` | Curves and peak times (a peak at the first or last sample is only a bound) |
| `model_compare/compare.png` | Measured fold over 0 µM against the ICHNOS model |
| `decoder_lodo/metrics.csv`, `posterior.png` | Could a single image's dose be recovered on a day the decoder never saw? |
| `run-manifest.json` | Settings, input checksum, code version, output checksums |

## 5. Decode a new culture

```bash
python3 scripts/decode_field.py --reference-cells analysis_ref/cells_per_s.csv \
    --cells new_run/cells_per_s.csv --stressor h2o2 --out-dir new_run/decoded
```

`--decoder paired` (default) compares each image with the same-day 0 µM
culture; `--decoder absolute` uses the image alone and is sensitive to
anything that differs between days (transformation efficiency, lamp,
slide).

The printout starts with the decoder's **track record** on held-out
reference days, then one row per field and per sample:

| Verdict | Meaning | What to do |
| --- | --- | --- |
| `dose estimate` | Informative posterior from a decoder that beat chance on held-out days | Report `map_dose_uM` with the 80 % credible set `credible80_uM` |
| `not validated` | The posterior points somewhere, but on held-out days this decoder did no better than chance | Do not report a dose |
| `no dose information` | Posterior ≈ prior: the image cannot tell the doses apart | Do not report a dose |
| `unlike reference` | Farther from every reference dose than reference images ever are | Check reporter, settings and stressor; extend the reference |
| `needs a 0 µM image` | Paired decoder without a same-day control | Image a 0 µM culture on the same day, or use `--decoder absolute` |
| `too few cells` | Fewer usable cells than `--min-cells` (default 20) | Image more cells or a denser field |

**Time since the stress** is decoded next to the dose, into windows
(`--time-windows 0,2,4,6,8` by default: 0–2, 2–4, 4–6, 6–8 h), from the
tandem timer: red/green of the reporter-positive cells, which does not depend
on how much reporter a cell carries. Leave `measurement_time_hours` blank in
the manifest when you do not know when the stress came. The time decoder has
its own track record and its own verdict:

| Time verdict | Meaning | What to do |
| --- | --- | --- |
| `time estimate` | Informative posterior from a time decoder that beat chance on held-out days | Report `time_map` with `time_credible80` |
| `not validated` / `no time information` / `unlike reference` | As for the dose | Do not report a time |
| `too few reporter cells` | Fewer than 5 reporter-positive cells: the timer cannot be read | Image more cells |
| `unstressed control` | A 0 µM culture has no time since stress | — |
| `time not available` | The reference cannot train a time decoder (e.g. its stressed images cover fewer than 2 windows) | Extend the reference |

A reference field whose window no other reference day covers cannot be
scored on held-out days; it is reported as unscorable, never as right or
wrong. Windows can only be told apart if the reference has stressed images
in each of them, on more than one day.

Full numbers: `decoded.csv` (posteriors `p_<dose>uM` and `pt_<window>`),
`validation.csv` and `validation_time.csv` (track records), `summary.json`.

## Worked example (team data, October 2026)

The 7/10 H2O2 day plays the "new culture": its 16 stressed samples are
renamed `sample01…sample16` so that their doses are hidden (`key.csv` keeps
the answers), its 0 µM samples keep their names as the same-day control, and
the day is left out of the reference.

```bash
CELLS=outputs/full_20261012_cellpose/cells_all_per_s.csv
python3 scripts/make_worked_example.py --cells $CELLS --session 20261007_h2o2_tr --out-dir outputs/example_user
python3 scripts/decode_field.py --reference-cells $CELLS --cells outputs/example_user/cells_per_s.csv \
    --stressor h2o2 --exclude-reference-session 20261007_h2o2_tr --out-dir outputs/example_user/decoded
```

Result: the track record on the two remaining H2O2 days is at chance
(accuracy 0.19 against 0.20, −1.78 bits per field). Of 48 fields (39 hidden,
9 controls): 23 `unlike reference` (7/10 has far more reporter-positive cells
than the reference days), 15 `not validated`, 10 `no dose information`, none
`dose estimate`. The most probable dose would have been right for 15 % of the
hidden fields, below chance, so declining to name a dose is the correct
answer on these data ([run_log.md](run_log.md), 2026-10-13).

With `--hide-time` the example also blanks the hours since stress
(`key.csv` keeps them). The time decoder then has two reference days, 5/10
(images at 1–4 h) and 6/10 (6–7 h), which never cover each other's windows,
so it has no track record; of the 39 hidden fields, 20 get `no time
information`, 7 `not validated`, 7 `unlike reference`, 5 `too few reporter
cells`, none `time estimate`. With all three H2O2 days as reference, its
leave-one-day-out record is accuracy 0.19 against chance 0.31 (−0.95 bits per
field); for DTT 0.20 against 0.28 (−0.21 bits): the timer ratio does not carry
the time since stress in these data, because the stress did not set off new
reporter synthesis in the cells that carry it.

## When something goes wrong

| Message | Cause |
| --- | --- |
| `output directory already exists` | Every run writes a new folder; choose another `--out-dir` |
| `lacks [...]; run normalize_exposure.py first` | The cell table is not per-second yet |
| `session_id not of the form ...` / `condition_id does not start with ...` | Reference identifiers break the naming rules above |
| `reference has no <stressor> sessions with non-transfected cells` | Add a non-transfected session for the baseline |
| Fields refused for foreground fraction | Dense slide; see step 3 |

Further reading: [architecture.md](architecture.md) (how the stages fit
together), [data_contract.md](data_contract.md) (every file and column),
[limitations.md](limitations.md) (what the results can and cannot support).
