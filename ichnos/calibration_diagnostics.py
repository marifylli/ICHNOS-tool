"""Compare ratio calibration assumptions; diagnostics do not apply corrections."""
import numpy as np
from .calibration import _ratios


def compare_ratio_calibration(model_ratios, image_ratios, *, holdout_model=None, holdout_image=None):
    """Compare origin-constrained and affine fits on identical reference pairs.

    The caller must supply independent held-out references when available.
    In-sample improvement with one additional parameter is not validation.
    A ratio intercept is not an estimate of per-channel autofluorescence.
    """
    x = _ratios(model_ratios, "model_ratios", positive=True)
    y = _ratios(image_ratios, "image_ratios")
    if x.shape != y.shape or len(x) < 3:
        raise ValueError("need at least three paired references")
    if np.ptp(x) <= np.finfo(float).eps * max(1., float(np.max(x))):
        raise ValueError("references do not span a resolvable ratio range")
    scale = float(np.dot(x,y)/np.dot(x,x))
    centered = x-x.mean()
    slope = float(np.dot(centered,y-y.mean())/np.dot(centered,centered))
    intercept = float(y.mean()-slope*x.mean())
    result = dict(
        n_reference_pairs=len(x), ratio_range=[float(x.min()),float(x.max())],
        through_origin=dict(slope=scale,intercept=0.,rmse=float(np.sqrt(np.mean((y-scale*x)**2))),residual_degrees_of_freedom=len(x)-1),
        affine=dict(slope=slope,intercept=intercept,rmse=float(np.sqrt(np.mean((y-slope*x-intercept)**2))),residual_degrees_of_freedom=len(x)-2),
        corrections_applied=False, independent_reference_provenance_verified=False,
        inference="model-conditioned diagnostic; in-sample RMSE is not validation",
        heldout_comparison_available=False,
    )
    if (holdout_model is None) != (holdout_image is None):
        raise ValueError("both held-out arrays must be supplied")
    if holdout_model is not None:
        hx = _ratios(holdout_model,"holdout_model",positive=True)
        hy = _ratios(holdout_image,"holdout_image")
        if hx.shape != hy.shape:
            raise ValueError("held-out arrays must be paired")
        result.update(heldout_comparison_available=True,n_holdout_pairs=len(hx))
        for key,a,b in [('through_origin',scale,0.),('affine',slope,intercept)]:
            result[key]['holdout_rmse']=float(np.sqrt(np.mean((hy-a*hx-b)**2)))
    return result
