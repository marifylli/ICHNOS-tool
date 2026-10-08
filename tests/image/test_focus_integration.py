import numpy as np
from dataclasses import asdict
from ichnos_image.focus import CellFocus, FocusPolicy, evaluate_cell
from ichnos_image.pipeline import ImageSet, process_image_set
import pytest


def test_low_cnr_and_unknown_are_not_called_defocus():
    good = CellFocus(1, 5., 10., True)
    policy = FocusPolicy('enforce', 1.5)
    for other, reason in ((CellFocus(1, 5., .2, True), 'focus_low_cnr'),
                          (CellFocus(1, 5., float('nan'), True), 'focus_cnr_unknown')):
        status, agreement, passed = evaluate_cell(good, other, policy)
        assert status == reason and agreement is None and not passed
    with pytest.raises(ValueError):
        FocusPolicy('enforce')


def test_pipeline_exports_scores_and_rejects_low_cnr(monkeypatch):
    from ichnos_image import pipeline
    labels = np.zeros((30,30), int); labels[10:20,10:20] = 1
    green = np.where(labels, 100., 10.)
    monkeypatch.setattr(pipeline.segment, 'segment_cells', lambda *a, **k: labels)
    monkeypatch.setattr(pipeline.segment, 'focus_score', lambda *a: 1.)
    monkeypatch.setattr(pipeline.correct, 'estimate_registration_shift', lambda *a: (0,0))
    monkeypatch.setattr(pipeline, 'score_cells', lambda *a: [CellFocus(1,3.,.1,True)])
    image = ImageSet(green,green,'s',0,0,100,100,0,0,'40X',1,30)
    record, = process_image_set(image, bleed_green_to_red=0, focus_policy=FocusPolicy('enforce',1.5))
    assert not record.qc_pass
    assert record.qc_reasons == 'focus_low_cnr'
    assert record.contrast_to_noise_green == .1
    assert record.focus_score_red == 3.
    report, = process_image_set(image, bleed_green_to_red=0)
    assert report.qc_pass and report.focus_status == 'focus_low_cnr'
