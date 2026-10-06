import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from ichnos.calibration import load_session_calibration


ROOT = Path(__file__).resolve().parents[2]


def test_synthetic_verification_from_real_protocol_outputs(tmp_path):
    protocol_dir = tmp_path / "protocol"
    command = [
        sys.executable, str(ROOT / "scripts/run_protocol.py"),
        "--variant", "ox", "--dose", "75", "--dose-units", "uM",
        "--times-hours", "0", "0.25", "0.5", "0.75", "1", "1.5", "2", "3", "4",
        "--initialization", "finite-preincubation", "--preincubation-hours", "2",
        "--out-dir", str(protocol_dir),
    ]
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    for noise in ("0", "0.02"):
        output = tmp_path / f"verification-{noise}"
        verification_command = [
            sys.executable, str(ROOT / "scripts/verify_session_calibration.py"),
            "--protocol-dir", str(protocol_dir), "--out-dir", str(output),
            "--relative-noise", noise,
        ]
        completed = subprocess.run(
            verification_command, cwd=ROOT, capture_output=True, text=True
        )
        assert completed.returncode == 0, completed.stderr
        report = json.loads((output / "verification.json").read_text())
        calibration = load_session_calibration(output / "calibration.json")
        assert calibration.source_kind == report["source_kind"] == "synthetic"
        assert report["experimental_instrument_validated"] is False
        assert report["n_fit"] == 5
        assert report["n_holdout"] == 4
        assert report["scale_relative_error"] < .02
        with (output / "synthetic_ratios.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        fit_rows = [row for row in rows if row["split"] == "fit"]
        np.testing.assert_array_equal(
            calibration.model_ratios,
            [float(row["model_ratio_red_green"]) for row in fit_rows],
        )
        if noise == "0":
            assert report["scale_relative_error"] < 1e-12
            assert report["holdout_rmse_model_ratio_units"] < 1e-12
        previous = (output / "calibration.json").read_bytes()
        completed = subprocess.run(
            verification_command, cwd=ROOT, capture_output=True, text=True
        )
        assert completed.returncode != 0
        assert "already exists" in completed.stderr
        assert (output / "calibration.json").read_bytes() == previous
