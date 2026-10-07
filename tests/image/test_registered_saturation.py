import numpy as np
import pandas as pd
import pytest
from ichnos_image.correct import register_saturation_mask
from ichnos_image.pipeline import ImageSet, process_image_set
from ichnos_image.saturation_audit import saturation_impact
from dataclasses import asdict


def test_fractional_shift_flags_all_interpolation_contributors():
    mask = np.zeros((10, 10), bool)
    mask[4, 4] = True
    moved = register_saturation_mask(mask, (.5, -.5))
    assert set(map(tuple, np.argwhere(moved))) == {(4,3),(4,4),(5,3),(5,4)}
    np.testing.assert_equal(register_saturation_mask(mask, (0,0)), mask)


@pytest.mark.parametrize('shift,red_pixel,expected,legacy', [
    ((0,-2), (8,10), True, False),
    ((0,2), (8,9), False, True),
    ((0,-.5), (8,10), True, False),
])
def test_pipeline_moves_red_saturation_and_reports_impact(monkeypatch, shift, red_pixel, expected, legacy):
    labels = np.zeros((20,20), int); labels[6:10,6:10] = 1
    green = np.where(labels,100.,10.); red = green.copy()
    red_mask = np.zeros_like(labels, bool); red_mask[red_pixel] = True
    from ichnos_image import pipeline
    monkeypatch.setattr(pipeline.segment, 'segment_cells', lambda *a, **k: labels)
    monkeypatch.setattr(pipeline.segment, 'focus_score', lambda *a: 1.)
    monkeypatch.setattr(pipeline.correct, 'estimate_registration_shift', lambda *a: shift)
    image = ImageSet(green,red,'s',0,0,100,100,0,0,'40X',1,30,
                     sample_id='target', green_saturation_mask=np.zeros_like(labels,bool),
                     red_saturation_mask=red_mask)
    records = process_image_set(image, bleed_green_to_red=0)
    assert records[0].sat_flag is expected
    assert records[0].sat_flag_legacy is legacy
    report = saturation_impact(pd.DataFrame([asdict(r) for r in records]))
    assert report['samples'][0]['n_qc_after'] == int(not expected)
    assert report['samples'][0]['n_qc_before'] == int(not legacy)
