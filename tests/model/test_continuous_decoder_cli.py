import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

from ichnos.decoder import build_dose_table, save_dose_table


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/run_continuous_decoder.py"


def command(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], cwd=ROOT, capture_output=True, text=True)


def test_cli_checks_interpolation_and_decodes_offgrid_data(tmp_path):
    table = build_dose_table(
        variant="ox", doses_uM=[25, 30, 35, 40, 45, 50], times_hours=[.5, 1],
        initialization="finite-preincubation", preincubation_hours=2,
    )
    path = tmp_path / "table.json"
    save_dose_table(table, path)
    validation = tmp_path / "validation.json"
    completed = command("check", "--table", path, "--allowed-error", ".001", "--out", validation)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(validation.read_text())["passed"] is True
    predicted = build_dose_table(
        variant="ox", doses_uM=[37.5, 42.5], times_hours=[.5, 1],
        initialization="finite-preincubation", preincubation_hours=2,
    )
    observations = tmp_path / "observations.csv"
    pd.DataFrame({"time_hours": [.5, 1], "calibrated_ratio_red_green": predicted.ratios[0]}).to_csv(observations, index=False)
    output = tmp_path / "result.json"
    args = ("decode", "--table", path, "--validation", validation,
            "--observations", observations, "--ratio-tolerance", "1e-5",
            "--data-kind", "synthetic", "--out", output)
    completed = command(*args)
    assert completed.returncode == 0, completed.stderr
    result = json.loads(output.read_text())
    assert result["status"] == "single_compatible_region"
    assert abs(result["dose_estimate_uM"] - 37.5) < .5
    assert result["experimentally_validated"] is False
    previous = output.read_bytes()
    assert command(*args).returncode != 0
    assert output.read_bytes() == previous
