import numpy as np
import pytest
from ichnos.calibration_diagnostics import compare_ratio_calibration


def test_nonzero_offset_is_detected_and_checked_on_holdout():
    result=compare_ratio_calibration([1,2,4],[2.5,4.5,8.5],holdout_model=[1.5,3],holdout_image=[3.5,6.5])
    assert result['affine']['intercept']==pytest.approx(.5)
    assert result['affine']['slope']==pytest.approx(2)
    assert result['affine']['holdout_rmse']<1e-12
    assert result['through_origin']['holdout_rmse']>.05
    assert result['affine']['residual_degrees_of_freedom']==1
    assert result['independent_reference_provenance_verified'] is False
    assert result['corrections_applied'] is False


def test_proportional_references_do_not_invent_offset():
    result=compare_ratio_calibration([1,2,4],[1.5,3,6])
    assert abs(result['affine']['intercept'])<1e-12
    assert result['through_origin']['rmse']<1e-12
    assert result['heldout_comparison_available'] is False


@pytest.mark.parametrize('x,y,kwargs',[
    ([1,1,1],[2,2,2],{}),([1,2],[2,4],{}),([1,2,3],[2,4],{}),
    ([1,2,3],[2,4,6],dict(holdout_model=[4])),
    ([1,2,3],[2,4,6],dict(holdout_model=[4],holdout_image=[8,9])),
    ([1,2,np.nan],[2,4,6],{}),
])
def test_invalid_or_unidentifiable_diagnostic_rejected(x,y,kwargs):
    with pytest.raises(ValueError):
        compare_ratio_calibration(x,y,**kwargs)
