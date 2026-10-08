import json
import numpy as np
import pytest
from ichnos.evaluation import posterior_metrics, run_evaluation


def test_offgrid_metric_is_explicit_and_tied_maps_retained():
    result=dict(posterior_computed=True,status='posterior_available',credible_mass_achieved=1.,
        map_grid_pairs=[dict(dose_uM=0.,elapsed_time_hours=0.),dict(dose_uM=10.,elapsed_time_hours=2.)],
        credible_set=[dict(dose_uM=0.,elapsed_time_hours=0.)])
    metrics=posterior_metrics(result,np.array([0.,10.]),np.array([0.,2.]),1.,.2,False)
    assert metrics['truth_in_discrete_set'] is None
    assert metrics['truth_cell_in_set']
    assert metrics['map_dose_error_min']==1 and metrics['map_dose_error_max']==9
    assert metrics['credible_area_fraction']==.25


def test_benchmark_ox_er_offgrid_noise_calibration_and_refinement(tmp_path):
    report=run_evaluation(tmp_path/'evaluation',dose_nodes=3,time_nodes=3,trials=2,
        noise_scales=(.5,2.),scenarios=('matched','ratio_gain_plus_10pct'),refine=True)
    assert report['n_trials']==2*2*2*2*2*2
    assert {s['variant'] for s in report['summaries']}=={'ox','er'}
    assert {s['truth_kind'] for s in report['summaries']}=={'on_grid','off_grid'}
    assert all(s['discrete_set_inclusion_fraction'] is None for s in report['summaries'] if s['truth_kind']=='off_grid')
    assert (tmp_path/'evaluation'/'manifest.json').exists()
    with pytest.raises(FileExistsError):
        run_evaluation(tmp_path/'evaluation',trials=2)
