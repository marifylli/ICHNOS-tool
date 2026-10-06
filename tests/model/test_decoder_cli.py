import json
from pathlib import Path
import subprocess
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/run_decoder.py"


def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], cwd=ROOT, capture_output=True, text=True)


def test_cli_build_decode_and_reject_overwrite(tmp_path):
    artifact = tmp_path / "table.json"
    completed = run(
        "build", "--variant", "ox", "--doses-uM", "0", "25", "75",
        "--times-hours", ".5", "1", "--initialization", "finite-preincubation",
        "--preincubation-hours", "2", "--out", artifact,
    )
    assert completed.returncode == 0, completed.stderr
    data = json.loads(artifact.read_text())
    observations = tmp_path / "observations.csv"
    pd.DataFrame({
        "time_hours": data["times_hours"],
        "calibrated_ratio_red_green": data["ratios"][1],
    }).to_csv(observations, index=False)
    output = tmp_path / "decoded.json"
    completed = run(
        "decode", "--table", artifact, "--observations", observations,
        "--ratio-tolerance", "1e-9", "--data-kind", "synthetic", "--out", output,
    )
    assert completed.returncode == 0, completed.stderr
    result = json.loads(output.read_text())
    assert result["status"] == "unique_grid_match"
    assert result["dose_estimate_uM"] == 25
    assert result["artifact_id"] == data["artifact_id"]
    previous = output.read_bytes()
    completed = run(
        "decode", "--table", artifact, "--observations", observations,
        "--ratio-tolerance", "1e-9", "--data-kind", "synthetic", "--out", output,
    )
    assert completed.returncode != 0
    assert output.read_bytes() == previous
