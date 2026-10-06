"""Check benchmark metrics and an independent small CLI scenario."""
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from ichnos.decoder import DoseTable
from ichnos.snapshot_decoder import ObservableTable
from scripts.verify_dense_posterior import ratio_posterior, summarize_mass

ROOT = Path(__file__).resolve().parents[2]


def test_ratio_marginal_ignores_green_and_matches_gaussian_odds():
    table = ObservableTable(DoseTable([0, 1], [0, 1], [[1, 2], [1, 2]], {"source_kind":"model_generated", "ratio_direction":"red/green", "experimentally_validated":False}), [[1, 100], [3, 9]])
    mass = ratio_posterior(1., table, .5, {'kind':'uniform_grid'})
    assert np.isclose(mass.sum(), 1)
    assert np.isclose(mass[0,0]/mass[0,1], np.exp(np.log(2)**2))
    assert np.allclose(mass[0], mass[1])


def test_credible_and_map_ties_are_not_broken_arbitrarily():
    result = summarize_mass(np.array([[.5, 0], [0, .5]]), np.array([0,800]), np.array([0,8]), (0,0))
    assert result['map_tie_count'] == 2
    assert result['map_dose_abs_error_min_uM'] == 0
    assert result['map_dose_abs_error_max_uM'] == 800
    assert result['map_time_abs_error_max_hours'] == 8
    assert result['credible_node_count'] == 2
    assert result['credible_mass_achieved'] == 1


def test_dense_cli_records_both_scenarios_and_independent_targets(tmp_path):
    snapshot = tmp_path/'snapshot'
    command = [sys.executable, str(ROOT/'scripts/verify_snapshot_decoder.py'), '--out-dir', str(snapshot)]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    output = tmp_path/'dense'
    command = [sys.executable, str(ROOT/'scripts/verify_dense_posterior.py'), '--snapshot-dir', str(snapshot),
        '--out-dir', str(output), '--dose-nodes','3','--time-nodes','3','--repeats','1']
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    trials = pd.read_csv(output/'trials.csv')
    assert len(trials) == 2*2*3*3
    assert set(trials.scenario) == {'matched','clearance_misspecified'}
    assert (trials.credible_mass_achieved >= .95-1e-12).all()
    assert trials.groupby(['scenario','method']).size().eq(9).all()
    noise = json.loads((output/'noise.json').read_text())
    targets = json.loads((output/'observations.json').read_text())
    assert set(noise['biological_replicate_ids']).isdisjoint(x['observation']['biological_replicate_id'] for x in targets)
    default = json.loads((output/'table.json').read_text())
    clearance = json.loads((output/'clearance-truth-table.json').read_text())
    assert not np.allclose(default['green'], clearance['green'])
    assert (output/'posterior-comparison.png').is_file()
    assert (output/'posterior-comparison.pdf').is_file()
    assert subprocess.run(command, capture_output=True).returncode != 0
