import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "run_protocol.py"


def _command(variant, out_dir):
    return [
        sys.executable,
        str(SCRIPT),
        "--variant", variant,
        "--dose", "75",
        "--dose-units", "uM",
        "--times-hours", "0", "0.5", "1",
        "--initialization", "finite-preincubation",
        "--preincubation-hours", "2",
        "--clearance-rate-per-hour", "0.5",
        "--out-dir", str(out_dir),
    ]


@pytest.mark.parametrize("variant", ["ox", "er"])
def test_cli_exports_results_and_assumptions(tmp_path, variant):
    out_dir = tmp_path / variant
    completed = subprocess.run(
        _command(variant, out_dir),
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, (
        completed.stdout + completed.stderr
    )

    metadata = json.loads(
        (out_dir / "metadata.json").read_text(encoding="utf-8")
    )
    assert metadata["parameter_profile"]["variant"] == variant
    assert metadata["parameter_profile"]["parameters"]
    assert metadata["exposure"]["initial_dose"] == 75.0
    assert metadata["exposure"]["model"] == "first_order_decay"
    assert metadata["exposure"]["clearance_rate_per_hour"] == 0.5
    assert metadata["initialization"]["method"] == (
        "finite_preincubation_assumption"
    )
    assert metadata["initialization"]["zero_stress_duration_hours"] == 2.0
    assert metadata["initialization"][
        "experimental_initial_state_validated"
    ] is False

    sbml = (out_dir / "model.sbml").read_text(encoding="utf-8")
    assert metadata["model_sbml_sha256"] == hashlib.sha256(
        sbml.encode("utf-8")
    ).hexdigest()

    with (out_dir / "results.csv").open(
        encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == metadata["csv_columns"]
        rows = list(reader)

    times = np.array([float(row["time_hours"]) for row in rows])
    stress = np.array([float(row[f"S_{variant}"]) for row in rows])
    np.testing.assert_allclose(times, [0.0, 0.5, 1.0])
    np.testing.assert_allclose(
        stress,
        75.0 * np.exp(-0.5 * times),
        rtol=1e-6,
    )
    assert "Observed_Green" in rows[0]


def test_cli_rejects_existing_output_directory(tmp_path):
    out_dir = tmp_path / "existing"
    out_dir.mkdir()
    sentinel = out_dir / "keep.txt"
    sentinel.write_text("keep", encoding="utf-8")

    completed = subprocess.run(
        _command("ox", out_dir),
        cwd=ROOT,
        capture_output=True,
        text=True,
    )

    assert completed.returncode != 0
    assert "already exists" in completed.stderr
    assert sentinel.read_text(encoding="utf-8") == "keep"
    assert sorted(path.name for path in out_dir.iterdir()) == [
        "keep.txt"
    ]