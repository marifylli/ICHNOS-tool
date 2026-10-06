"""Reproducible noise-free synthetic image-to-decoder verification."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
from ichnos.calibration import fit_session_calibration, model_reference_id, save_session_calibration
from ichnos.continuous_decoder import validate_interpolation, save_interpolation_validation
from ichnos.decoder import build_dose_table, save_dose_table
from ichnos.sample_decoder import decode_sample
from ichnos_image.pipeline import ImageSet, process_image_set
from ichnos_image.export import export_csv

ROOT = Path(__file__).resolve().parents[1]
SCALE = 1.5
MODEL = dict(variant="ox", initialization="finite-preincubation", preincubation_hours=2)


def write_json(path, data):
    encoded = json.dumps(data, indent=2, allow_nan=False) + "\n"
    with path.open("x", encoding="utf-8") as handle:
        handle.write(encoded)


def image_cells(directory, name, times, ratios, *, intensity=1000.):
    """Persist idealized arrays and run the production image processing path."""
    yy, xx = np.mgrid[:128, :128]
    mask = np.zeros((128, 128), dtype=bool)
    for cy, cx in ((32, 32), (64, 96), (96, 32)):
        mask |= (yy-cy)**2 + (xx-cx)**2 < 10**2
    records = []
    for index, (time, ratio) in enumerate(zip(times, ratios)):
        green = np.where(mask, intensity, 0.)
        red = SCALE * ratio * green
        bright = np.where(mask, .2, .5)
        np.savez_compressed(directory / f"{name}-{index}.npz", green=green, red=red, bright_field=bright)
        cells = process_image_set(ImageSet(
            green=green, red=red, bright_field=bright, session_id="synthetic-e2e",
            sample_id=name, condition_id=name, timepoint=index, acquisition_order=index,
            exposure_ms_green=100, exposure_ms_red=200, nd_filter_green=0,
            nd_filter_red=0, objective="40X", burner_hours=1, lamp_warmup_minutes=30,
            sampling_time_hours=float(time), measurement_time_hours=float(time+.25),
        ), bleed_green_to_red=0)
        records.extend(cells)
    export_csv(records, directory / f"{name}-cells.csv")
    return records


def summarize(directory, name):
    runpy.run_path(str(ROOT / "scripts/summarize_cells.py"))["run"](SimpleNamespace(
        cells=directory / f"{name}-cells.csv", out_dir=directory / f"{name}-summary",
        data_kind="synthetic", calibrations=directory / "calibrations.csv",
        min_cells=3, green_floor=500.,
    ))


def run(out_dir):
    if out_dir.exists():
        raise FileExistsError(f"output directory already exists: {out_dir}")
    out_dir.mkdir(parents=True)
    checks = []
    def check(name, passed, expected, actual):
        checks.append(dict(name=name, passed=bool(passed), expected=expected, actual=actual))
    try:
        reference_dir = out_dir / "reference"
        command = [sys.executable, str(ROOT / "scripts/run_protocol.py"),
                   "--variant", "ox", "--dose", "25", "--dose-units", "uM",
                   "--times-hours", ".25", ".5", ".75", "1", "1.5", "2",
                   "--initialization", "finite-preincubation", "--preincubation-hours", "2",
                   "--out-dir", str(reference_dir)]
        completed = subprocess.run(command, capture_output=True, text=True)
        if completed.returncode:
            raise ValueError(completed.stderr)
        reference = json.loads((reference_dir / "metadata.json").read_text())
        fluorescence = pd.read_csv(reference_dir / "fluorescence.csv")
        reference_ratios = fluorescence.ratio_red_green.to_numpy()
        reference_times = fluorescence.time_hours.to_numpy()
        records = image_cells(out_dir, "calibration-reference", reference_times, reference_ratios)
        frame = pd.read_csv(out_dir / "calibration-reference-cells.csv")
        medians = frame.groupby("timepoint").ratio_red_green.median().to_numpy()
        check("reference_cell_extraction", len(records)==18 and frame.qc_pass.all(), 18, len(records))
        calibration = fit_session_calibration(reference_ratios[::2], medians[::2],
            session_id="synthetic-e2e", reference_id=model_reference_id(reference),
            source_kind="synthetic", source="idealized synthetic images; alternating fit/holdout reference times")
        save_session_calibration(calibration, out_dir / "calibration.json")
        pd.DataFrame([dict(session_id=calibration.session_id, reference_id=calibration.reference_id,
                           calibration_path="calibration.json")]).to_csv(out_dir / "calibrations.csv", index=False)
        holdout_error = float(np.max(np.abs(medians[1::2]/calibration.c_session-reference_ratios[1::2])))
        check("session_scale", abs(calibration.c_session-SCALE)<1e-10, SCALE, calibration.c_session)
        check("calibration_holdout", holdout_error<1e-10, "max absolute error < 1e-10", holdout_error)
        table = build_dose_table(**MODEL, doses_uM=[0,25,75], times_hours=[.5,1,1.5])
        save_dose_table(table, out_dir / "table.json")
        # Fresh forward simulation supplies target ground truth, not table lookup.
        truth = build_dose_table(**MODEL, doses_uM=[75,100], times_hours=[.5,1])
        image_cells(out_dir, "dose75", truth.times_hours, truth.ratios[0])
        summarize(out_dir, "dose75")
        def options(name, table_name="table.json"):
            return dict(samples_csv=out_dir/f"{name}-summary/samples.csv",
                        summary_metadata=out_dir/f"{name}-summary/metadata.json",
                        calibration_reference=reference_dir/"metadata.json", dose_table=out_dir/table_name,
                        session_id="synthetic-e2e", sample_id=name,
                        time_field="sampling_time_hours", ratio_tolerance=1e-8)
        # Known-time dose decoding needs precisely the selected observation times.
        known = build_dose_table(**MODEL, doses_uM=[0,25,75], times_hours=[.5,1])
        save_dose_table(known, out_dir/"known-table.json")
        discrete = decode_sample(**options("dose75", "known-table.json"))
        write_json(out_dir/"discrete.json", discrete)
        check("known_time_dose", discrete["dose_estimate_uM"]==75, 75, discrete["dose_estimate_uM"])
        joint = decode_sample(**options("dose75"), mode="joint", elapsed_time_grid_hours=[.5,1])
        write_json(out_dir/"joint.json", joint)
        check("joint_dose_time", joint["status"]=="unique_grid_pair" and joint["dose_estimate_uM"]==75 and joint["elapsed_time_estimate_hours"]==.5,
              dict(dose_uM=75,elapsed_hours=.5), dict(status=joint["status"],dose_uM=joint["dose_estimate_uM"],elapsed_hours=joint["elapsed_time_estimate_hours"]))
        ambiguous_options = options("dose75")
        ambiguous_options["ratio_tolerance"] = 10.
        ambiguous = decode_sample(**ambiguous_options, mode="joint", elapsed_time_grid_hours=[.5,1])
        write_json(out_dir/"ambiguous.json", ambiguous)
        check("loose_tolerance_ambiguity", ambiguous["status"]=="ambiguous" and ambiguous["dose_estimate_uM"] is None,
              "ambiguous, no point estimate", ambiguous["status"])
        continuous_table = build_dose_table(**MODEL, doses_uM=[25,30,35,40,45,50], times_hours=[.5,1])
        save_dose_table(continuous_table,out_dir/"continuous-table.json")
        validation = validate_interpolation(continuous_table, allowed_error=.001)
        save_interpolation_validation(validation,out_dir/"interpolation.json")
        offgrid = build_dose_table(**MODEL, doses_uM=[37.5,42.5], times_hours=[.5,1])
        image_cells(out_dir,"dose37.5",offgrid.times_hours,offgrid.ratios[0])
        summarize(out_dir,"dose37.5")
        continuous = decode_sample(**options("dose37.5","continuous-table.json"),mode="continuous",
                                   interpolation_validation=out_dir/"interpolation.json")
        write_json(out_dir/"continuous.json",continuous)
        estimate = continuous["dose_estimate_uM"]
        check("continuous_dose", estimate is not None and abs(estimate-37.5)<.5,
              "37.5 uM within 0.5 uM",estimate)
        image_cells(out_dir,"low-green",truth.times_hours,truth.ratios[0],intensity=100.)
        summarize(out_dir,"low-green")
        low = pd.read_csv(out_dir/"low-green-summary/samples.csv")
        try:
            decode_sample(**options("low-green","known-table.json"))
            rejected = False
            reason = "unexpected accepted sample"
        except ValueError as exc:
            reason = str(exc)
            rejected = "insufficient or invalid usable cell counts" in reason
        write_json(out_dir/"low-green.json",dict(rejected=rejected,reason=reason))
        check("low_green_rejection", rejected and (low.n_cells_used==0).all(),
              "zero usable cells, decoding rejected",dict(usable_cells=low.n_cells_used.tolist(),reason=reason))
    except (ValueError,KeyError,OSError) as exc:
        check("execution_error",False,"complete execution",str(exc))
    report = dict(schema_version=1,passed=all(row["passed"] for row in checks),checks=checks,
                  data_kind="synthetic",experimentally_validated=False,
                  assumptions=dict(variant="ox",profile="default",preincubation_hours=2,
                                   known_session_scale=SCALE,noise="none",bleed_green_to_red=0,
                                   summary_green_floor=500.,min_cells=3),
                  limitations=["idealized images, fixed model and protocol; not experimental accuracy",
                               "three cells are not independent biological replicates",
                               "loose-tolerance ambiguity is deliberately constructed",
                               "interpolation checks are sampled, not a uniform error certificate"],
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    report["files_sha256"] = {str(path.relative_to(out_dir)):hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in sorted(out_dir.rglob("*")) if path.is_file()}
    write_json(out_dir/"report.json",report)
    lines = ["# Synthetic end-to-end verification", "", f"Passed: {report['passed']}", "",
             "| Check | Passed | Expected | Actual |", "| --- | --- | --- | --- |"]
    for row in checks:
        lines.append(f"| {row['name']} | {row['passed']} | {json.dumps(row['expected'])} | {json.dumps(row['actual'])} |")
    lines.extend(["", "Noise-free synthetic software verification; not biological or instrument validation."])
    (out_dir/"report.md").write_text("\n".join(lines)+"\n")
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir",type=Path,required=True)
    args=parser.parse_args()
    try:
        report=run(args.out_dir)
    except (ValueError,KeyError,OSError) as exc:
        parser.exit(1,f"Error: {exc}\n")
    print(f"Verification passed={report['passed']} -> {args.out_dir/'report.md'}")
    if not report["passed"]:
        parser.exit(1,"Verification failed; inspect report.json\n")


if __name__=="__main__":
    main()
