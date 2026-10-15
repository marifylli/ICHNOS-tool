"""ichnos_vail — the ICHNOS application: a local graphical front end.

    pip install -e ".[ichnos_vail]"
    streamlit run app/ichnos_vail.py

Runs on this computer only (see .streamlit/config.toml). Every button runs
the same script a user could run by hand; the exact command is shown before
it starts, and long runs continue in the background if the browser is closed.
"""
from __future__ import annotations

import importlib.util
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

from ichnos_vail import commands, jobs, manifest_check, results  # noqa: E402

st.set_page_config(page_title="ichnos_vail", page_icon="🔬", layout="wide")

PAGES = ["Start", "1 · Check images", "2 · Images → cells", "3 · Analyse an experiment",
         "4 · Decode a new culture", "Jobs"]
STATE_LABEL = {"created": "⏳ created", "running": "🔄 running", "stopping": "⏹ stopping", "finished": "✅ finished",
               "failed": "❌ failed", "stopped": "⏹ stopped", "interrupted": "⚠️ interrupted"}
DEFAULT_REFERENCE = ROOT / "data" / "reference" / "cells_all_per_s.csv.gz"


# ---------------------------------------------------------------- helpers

def stamp() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M")


def as_path(text: str) -> Path | None:
    text = (text or "").strip()
    if not text:
        return None
    path = Path(text).expanduser()
    return path if path.is_absolute() else ROOT / path


def path_input(label: str, key: str, default: str = "", help: str | None = None) -> Path | None:
    """A path field whose value survives switching pages (widget state alone does not)."""
    saved = f"saved_{key}"
    if saved not in st.session_state:
        st.session_state[saved] = default
    value = st.text_input(label, value=st.session_state[saved], key=f"w_{key}", help=help)
    st.session_state[saved] = value
    return as_path(value)


def set_path(key: str, value) -> None:
    """Fill a path field from code (e.g. the output of the previous step)."""
    st.session_state[f"saved_{key}"] = str(value)
    st.session_state.pop(f"w_{key}", None)


def show_commands(cmds: list[list[str]]) -> None:
    st.caption("This runs exactly the following (copy it to run without ichnos_vail):")
    st.code("\n".join(commands.as_shell(c) for c in cmds), language="bash")


def launch(name: str, cmds: list[list[str]], errors: list[str], key: str, *, outputs=(), progress=None) -> None:
    if errors:
        for error in errors:
            st.error(error)
        return
    job = jobs.create_job(name, cmds, outputs=[str(o) for o in outputs], progress=progress)
    jobs.start_job(job)
    st.session_state[key] = str(job)
    st.success(f"Started. The job keeps running if you close this tab — find it again under **Jobs**. ({job.name})")


@st.fragment(run_every="3s")
def job_panel(key: str) -> None:
    """Live status of the job stored under session key ``key``."""
    if not st.session_state.get(key):
        return
    job = Path(st.session_state[key])
    status, spec = jobs.read_status(job), jobs.read_job(job)
    state = status.get("state", "?")
    steps = len(spec.get("commands", []))
    started = status.get("started")
    elapsed = ((status.get("ended") or time.time()) - started) / 60 if started else 0
    cols = st.columns([2, 1, 1, 1])
    cols[0].markdown(f"**{spec.get('name', job.name)}** — {STATE_LABEL.get(state, state)}")
    cols[1].markdown(f"step {status.get('step', 0)} / {steps}")
    cols[2].markdown(f"{elapsed:.0f} min")
    if state == "running" and cols[3].button("Stop", key=f"stop-{key}"):
        jobs.stop_job(job)
    fraction = jobs.progress(job)
    if fraction is not None:
        st.progress(fraction, text=f"{fraction:.0%} of batches done")
    with st.expander("Log", expanded=state in ("running", "failed", "interrupted")):
        st.code(jobs.tail(job, 40) or "(nothing printed yet)", language="text")
    if state == "failed":
        st.error("The job failed — the end of the log says why. Fix the input and start again.")
    if state == "interrupted":
        st.warning("The job's process is gone (computer slept or restarted?). Segmentation runs can be "
                   "resumed: start the same run with the same output folder.")


def job_state(key: str) -> str | None:
    if not st.session_state.get(key):
        return None
    return jobs.read_status(Path(st.session_state[key])).get("state")


def verdict_style(value: str) -> str:
    colour = results.VERDICTS.get(value, ("#8a8a8a", ""))[0]
    return f"background-color: {colour}; color: white; font-weight: 600"


# ---------------------------------------------------------------- pages

def page_start() -> None:
    st.title("ichnos_vail")
    st.markdown(
        "From microscope images of a yeast culture carrying the red/green timer reporter to two answers: "
        "**does the reporter respond to the stress** (analysis of an experiment), and **which dose was a new "
        "culture exposed to** — with a verdict saying whether that answer can be trusted.\n\n"
        "Work through the steps in the sidebar. Each step shows the command it runs, so every result can be "
        "reproduced without ichnos_vail.")
    st.subheader("This computer")

    @st.cache_resource
    def environment() -> dict:
        found = {name: importlib.util.find_spec(module) is not None for name, module in
                 (("Cellpose (segmentation)", "cellpose"), ("PyTorch", "torch"),
                  ("libRoadRunner (model)", "roadrunner"), ("scikit-image", "skimage"))}
        gpu = "none found"
        if found["PyTorch"]:
            import torch
            if torch.cuda.is_available():
                gpu = "CUDA"
            elif getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
                gpu = "Apple GPU (MPS)"
        return {"python": platform.python_version(), "found": found, "gpu": gpu}

    env = environment()
    cols = st.columns(len(env["found"]) + 1)
    for col, (name, ok) in zip(cols, env["found"].items()):
        col.metric(name, "installed" if ok else "missing")
    cols[-1].metric("GPU", env["gpu"])
    if not env["found"]["Cellpose (segmentation)"]:
        st.warning("Cellpose is missing: `pip install cellpose` (bright-field segmentation needs it).")
    if env["gpu"] == "none found":
        st.info("No GPU: Cellpose will run on the CPU, several times slower than ~2 min per frame.")
    st.subheader("Read before trusting a result")
    st.markdown(
        "- In the repository's `docs` folder: `user_guide.md` (how to use the tool), `data_contract.md` "
        "(what each file contains), `limitations.md` (what results can and cannot support).\n"
        "- A dose is only reported when the verdict is **dose estimate**. With the team's October 2026 data as "
        "reference no image reaches it: those data contain no dose signal that carries over to a new day.")


def page_manifest() -> None:
    st.title("1 · Check images")
    st.markdown("Point to the manifest: one CSV row per imaged field (columns: `docs/data_contract.md`). "
                "Problems are found here, not hours into a run.")
    manifest = path_input("Manifest (CSV)", "manifest_path", "",
                          help="Absolute, or relative to the repository folder. Image paths inside it are "
                               "relative to the manifest's own folder.")
    if manifest and st.button("Check", type="primary"):
        st.session_state["manifest_checked"] = str(manifest)
    if st.session_state.get("manifest_checked") and manifest and str(manifest) == st.session_state["manifest_checked"]:
        table, issues, summary = manifest_check.check_manifest(manifest)
        if table is not None:
            c = st.columns(4)
            c[0].metric("Fields", summary.get("rows", 0))
            c[1].metric("Sessions", len(table["session_id"].unique()) if "session_id" in table else 0)
            c[2].metric("Errors", summary.get("errors", 0))
            c[3].metric("Warnings", summary.get("warnings", 0))
        if issues:
            frame = pd.DataFrame(issues)
            st.dataframe(frame.style.map(lambda v: "color: #c23b3b; font-weight: 600" if v == "error" else
                                         ("color: #b26b00" if v == "warning" else ""), subset=["level"]),
                         hide_index=True)
        else:
            st.success("No problems found.")
        if summary.get("sessions"):
            st.markdown("**Sessions**")
            st.dataframe(pd.DataFrame(summary["sessions"]).T)
        if table is not None:
            with st.expander("First rows"):
                st.dataframe(table.head(20))
            if summary.get("errors", 1) == 0:
                set_path("pipeline_manifest", manifest)

    st.divider()
    st.subheader("Find and align bright-field frames (optional)")
    st.markdown("Use this when bright-field frames sit in other folders or are shifted relative to the "
                "fluorescence. It writes a new manifest; the original is never changed.")
    if manifest:
        out = path_input("Output folder", "align_out", f"outputs/align_{stamp()}")
        search_up = st.number_input("Folders to search upwards", 0, 5, 2)
        cmds, errors = commands.align(manifest, out, int(search_up))
        show_commands(cmds)
        if st.button("Start alignment", disabled=job_state("align_job") == "running"):
            launch("align bright-field", cmds, errors, "align_job", outputs=[out / "images_aligned.csv"])
            set_path("pipeline_manifest", out / "images_aligned.csv")
        job_panel("align_job")
        if job_state("align_job") == "finished":
            st.success(f"Aligned manifest: `{out / 'images_aligned.csv'}` — step 2 will use it.")


def page_pipeline() -> None:
    st.title("2 · Images → cells")
    st.markdown("Segments every field on bright-field, measures green and red per cell, applies QC, then "
                "converts to grey levels per second. Long: about 2 minutes per field with Cellpose on a GPU. "
                "Runs in batches; an interrupted run resumes when started again with the same output folder.")
    manifest = path_input("Manifest", "pipeline_manifest", "")
    out = path_input("Output folder", "pipeline_out", f"outputs/cells_{stamp()}")
    c = st.columns(3)
    segmentation = c[0].selectbox("Segmentation", commands.SEGMENTATION, index=0)
    source = c[1].selectbox("Segment on", ["brightfield", "sum", "green", "red"], index=0)
    gpu = c[2].checkbox("Use the GPU (Cellpose)", value=True)
    c = st.columns(3)
    bleed_mode = c[0].radio("Green → red bleed-through", ["indicative value", "controls file"], horizontal=False)
    bleed = controls = None
    if bleed_mode == "indicative value":
        bleed = c[1].number_input("Bleed-through", 0.0, 1.0, 0.05, 0.01,
                                  help="Used when no single-colour controls were imaged (stated as a limitation).")
    else:
        with c[1]:
            controls = path_input("Controls CSV", "controls_path", "")
    saturation = c[2].number_input("Camera maximum (saturation)", 1, 65535, 255,
                                   help="255 for 8-bit, 4095 for 12-bit, 65535 for 16-bit images")
    with st.expander("Advanced"):
        c = st.columns(4)
        green = c[0].selectbox("Green from RGB", commands.EXTRACTION, help="Only for colour (RGB) files")
        red = c[1].selectbox("Red from RGB", commands.EXTRACTION)
        batch = c[2].number_input("Fields per batch", 1, 200, 20)
        dense = c[3].checkbox("Dense slides (allow masks over 50 % of the frame)")
    if not (manifest and out):
        return
    cmds, errors = commands.images_to_cells(manifest, out, batch_size=int(batch), segmentation=segmentation,
                                            source=source, gpu=gpu, bleed_default=bleed, controls=controls,
                                            saturation=float(saturation), green_extraction=green,
                                            red_extraction=red, max_foreground=0.5 if dense else None)
    rows = 0
    if manifest.is_file():
        rows = sum(1 for _ in manifest.open(encoding="utf-8", errors="replace")) - 1
        minutes = rows * (2 if (segmentation == "cellpose" and gpu) else 0.5)
        st.caption(f"{rows} fields; rough duration {minutes / 60:.1f} h. Keep the computer awake "
                   f"(macOS: start ichnos_vail with `caffeinate -i streamlit run app/ichnos_vail.py`).")
    show_commands(cmds)
    if st.button("Start", type="primary", disabled=job_state("pipeline_job") == "running"):
        launch("images to cells", cmds, errors, "pipeline_job", outputs=[out / "cells_per_s.csv"],
               progress={"kind": "batches", "total": max(1, -(-rows // int(batch)))})
    job_panel("pipeline_job")
    if job_state("pipeline_job") == "finished":
        st.success(f"Cell table: `{out / 'cells_per_s.csv'}`")
        set_path("analysis_cells", out / "cells_per_s.csv")
        set_path("decode_cells", out / "cells_per_s.csv")


def show_analysis(folder: Path) -> None:
    figures, tables = results.analysis_figures(folder), results.analysis_tables(folder)
    if not figures and not tables:
        st.info("No analysis results in this folder yet.")
        return
    manifest = results.read_json(folder / "run-manifest.json")
    if manifest:
        st.caption(f"Run status: {manifest.get('status')} · code {str(manifest.get('code', {}).get('commit'))[:7]} · "
                   f"input sha256 {manifest.get('input_sha256', '')[:12]}…")
    for title, frame in tables.items():
        with st.expander(title, expanded=title.startswith("Detected")):
            st.dataframe(frame, hide_index=True)
    for title, image in figures:
        st.markdown(f"**{title}**")
        st.image(str(image), width="stretch")


def page_analysis() -> None:
    st.title("3 · Analyse an experiment")
    st.markdown("Time courses, each dose against the same day's 0 µM, the model comparison and the "
                "leave-one-day-out decoder test. A few minutes.")
    cells = path_input("Cell table (cells_per_s.csv)", "analysis_cells",
                       str(DEFAULT_REFERENCE.relative_to(ROOT)) if DEFAULT_REFERENCE.exists() else "")
    out = path_input("Output folder", "analysis_out", f"outputs/analysis_{stamp()}")
    steps = st.multiselect("Steps", ["timecourse", "induction", "model_compare", "decoder_lodo"],
                           default=["timecourse", "induction", "model_compare", "decoder_lodo"])
    with st.expander("Advanced"):
        c = st.columns(4)
        n_boot = c[0].number_input("Resamples", 50, 10000, 1000, 50)
        seed = c[1].number_input("Seed", 0, 10_000, 0)
        quantile = c[2].number_input("Reporter-positive quantile", 0.5, 0.99, 0.9, 0.01)
        min_cells = c[3].number_input("Smallest field (cells)", 1, 1000, 20)
        clear_text = st.text_input("Model clearance scenarios (per hour; 'none' = constant stress)", "none, 1, 3")
        decoders = st.multiselect("Decoders", ["absolute", "paired", "model"], default=["absolute", "paired", "model"])
    if cells and out and steps:
        try:
            clearance = [None if t.strip().lower() == "none" else float(t) for t in clear_text.split(",") if t.strip()]
        except ValueError:
            st.error("Clearance scenarios: numbers or 'none', separated by commas.")
            return
        config = out.parent / f"{out.name}_config.json"
        cmds, errors = commands.analysis(config, cells, out)
        show_commands(cmds)
        if st.button("Start", type="primary", disabled=job_state("analysis_job") == "running"):
            if not errors:
                commands.write_analysis_config(config, cells, steps=steps, n_boot=int(n_boot), seed=int(seed),
                                               positive_quantile=float(quantile), clearance_per_h=clearance,
                                               decoders=decoders, min_cells=int(min_cells))
            launch("analysis", cmds, errors, "analysis_job", outputs=[out / "run-manifest.json"])
            set_path("analysis_view", out)
        job_panel("analysis_job")
    st.divider()
    view = path_input("Show results of an analysis folder", "analysis_view", "")
    if view and view.exists():
        show_analysis(view)


def show_decoded(folder: Path) -> None:
    summary, decoded = results.decode_result(folder)
    if decoded.empty:
        st.info("No decoded results in this folder yet.")
        return
    fields = decoded[decoded["unit"] == "field"]

    def chips(column: str) -> None:
        counts = fields[column].value_counts()
        spans = "".join(f"<span style='{verdict_style(v)};padding:6px 10px;border-radius:6px;"
                        f"white-space:nowrap'>{v}: {n}</span>" for v, n in counts.items())
        st.markdown(f"<div style='display:flex;flex-wrap:wrap;gap:6px'>{spans}</div>", unsafe_allow_html=True)

    left, right = st.columns(2)
    with left:
        st.markdown("#### Dose")
        gain = summary.get("validated_info_gain_bits")
        (st.success if gain is not None and gain > 0 else st.warning)(
            f"**Track record.** {results.track_record_text(summary)}")
        chips("verdict")
    with right:
        st.markdown("#### Time since the stress")
        gain_t = summary.get("time_validated_info_gain_bits")
        (st.success if gain_t is not None and gain_t > 0 else st.warning)(
            f"**Track record.** {results.time_record_text(summary)}")
        if "time_verdict" in fields:
            chips("time_verdict")
    with st.expander("What the verdicts mean", expanded=False):
        st.dataframe(pd.DataFrame([{"verdict": v, "what to do": t} for v, (_, t) in results.VERDICTS.items()]),
                     hide_index=True)
    unit = st.radio("Show", ["field", "sample"], horizontal=True,
                    help="sample = the fields of one condition averaged")
    shown = decoded[decoded["unit"] == unit]
    front = ["condition_id", "sample_id", "n_cells", "verdict", "map_dose_uM", "credible80_uM",
             "time_verdict", "time_map", "time_credible80", "known_dose_uM", "known_time_h"]
    front = [c for c in front if c in shown.columns]
    shown = shown[front + [c for c in shown.columns if c.startswith(("p_", "pt_"))]].copy()
    for column in ("credible80_uM", "time_map", "time_credible80"):
        if column in shown:
            shown[column] = shown[column].fillna("")
    verdicts = [c for c in ("verdict", "time_verdict") if c in shown.columns]
    st.dataframe(shown.style.map(verdict_style, subset=verdicts).format(precision=2), hide_index=True)
    st.caption("Read each verdict before its answer. map_dose_uM / time_map: the most probable dose / time window; "
               "credible80: the doses / windows holding 80 % of the posterior; p_… and pt_…: the full posteriors. "
               "known_…: what the table itself says (blank when unknown).")


def page_decode() -> None:
    st.title("4 · Decode a new culture")
    st.markdown("Which dose was my culture exposed to, and how long ago? Both decoders are trained on a reference "
                "experiment, and their scores on reference days they never saw are shown with every answer. Name "
                "the same-day unstressed culture `<stressor>_0uM_<label>` in the manifest; give unknown samples any "
                "other name. Leave `measurement_time_hours` blank if you do not know when the stress came.")
    reference = path_input("Reference cell table", "decode_reference",
                           str(DEFAULT_REFERENCE.relative_to(ROOT)) if DEFAULT_REFERENCE.exists() else "")
    cells = path_input("New culture's cell table", "decode_cells", "")
    stressor, exclude = "", []
    if reference and reference.is_file():
        try:
            sessions = results.sessions(reference)
        except (ValueError, KeyError) as exc:
            st.error(f"Cannot read session ids from the reference: {exc}")
            return
        stressors = sorted(sessions["stressor"].dropna().unique())
        c = st.columns(3)
        stressor = c[0].selectbox("Stressor", stressors) if stressors else ""
        tr = sessions[(sessions["stressor"] == stressor) & (sessions["strain"] == "tr")]["session_id"].tolist()
        exclude = c[1].multiselect("Leave reference days out", tr,
                                   help="e.g. the day you are decoding, when it is also in the reference")
        decoder = c[2].radio("Decoder", ["paired", "absolute"], help="paired needs the same-day 0 µM culture")
        c = st.columns(2)
        min_cells = c[0].number_input("Smallest field (cells)", 1, 1000, 20)
        windows = c[1].text_input("Time windows (hours since stress, edges)", "0,2,4,6,8")
        out = path_input("Output folder", "decode_out", f"outputs/decoded_{stamp()}")
        if cells and out:
            cmds, errors = commands.decode(reference, cells, stressor, out, decoder=decoder, exclude=exclude,
                                           min_cells=int(min_cells), time_windows=windows)
            show_commands(cmds)
            if st.button("Decode", type="primary", disabled=job_state("decode_job") == "running"):
                launch("decode", cmds, errors, "decode_job", outputs=[out / "decoded.csv"])
                set_path("decode_view", out)
            job_panel("decode_job")
    st.divider()
    view = path_input("Show a decoded folder", "decode_view", "")
    if view and view.exists():
        show_decoded(view)


def page_jobs() -> None:
    st.title("Jobs")
    rows = jobs.list_jobs()
    if not rows:
        st.info(f"No jobs yet. They are kept in `{jobs.default_jobs_dir()}`.")
        return
    table = pd.DataFrame([{**r, "state": STATE_LABEL.get(r["state"], r["state"]), "path": r["path"].name}
                          for r in rows])
    st.dataframe(table[["created", "name", "state", "step", "steps", "path"]], hide_index=True)
    choice = st.selectbox("Open", [r["path"].name for r in rows])
    st.session_state["jobs_selected"] = str(jobs.default_jobs_dir() / choice)
    spec = jobs.read_job(Path(st.session_state["jobs_selected"]))
    show_commands(spec.get("commands", []))
    job_panel("jobs_selected")


# ---------------------------------------------------------------- layout

page = st.sidebar.radio("Steps", PAGES, key="page")
with st.sidebar:
    st.divider()
    st.caption("Recent jobs")
    for row in jobs.list_jobs(limit=5):
        st.caption(f"{STATE_LABEL.get(row['state'], row['state'])} · {row['name']} · {row['created'][5:16]}")

{"Start": page_start, "1 · Check images": page_manifest, "2 · Images → cells": page_pipeline,
 "3 · Analyse an experiment": page_analysis, "4 · Decode a new culture": page_decode,
 "Jobs": page_jobs}[page]()
