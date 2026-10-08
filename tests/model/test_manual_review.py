import pandas as pd
import pytest
from test_population import cell_frame
from ichnos.manual_review import IDENTITY, review_summaries


def flags(cells, ids):
    frame = cells.loc[cells.cell_id.isin(ids), IDENTITY].copy()
    frame['review_reason'] = 'uncertain boundary'
    return frame


def test_identity_isolates_field_and_preserves_baseline():
    cells = pd.concat([cell_frame(), cell_frame(order=2)], ignore_index=True)
    original = cells.copy(deep=True)
    review = flags(cells[cells.acquisition_order.eq(1)], [3])
    baseline, sensitivity, comparison, marked = review_summaries(
        cells, review, data_kind='experimental', min_cells=1)
    assert baseline.iloc[0].n_cells_used == 6
    assert sensitivity.iloc[0].n_cells_used == 5
    assert marked.review_reason.notna().sum() == 1
    pd.testing.assert_frame_equal(cells, original)


@pytest.mark.parametrize('kind', ['unknown', 'duplicate'])
def test_invalid_flags_are_rejected(kind):
    cells = cell_frame()
    review = flags(cells, [1])
    if kind == 'unknown':
        review['cell_id'] = 999
    else:
        review = pd.concat([review, review])
    with pytest.raises(ValueError):
        review_summaries(cells, review, data_kind='experimental')


def test_all_flagged_preserves_group_with_insufficient_status():
    cells = cell_frame()
    baseline, sensitivity, comparison, _ = review_summaries(
        cells, flags(cells, [1, 2, 3, 4]), data_kind='experimental')
    assert sensitivity.iloc[0].n_cells_used == 0
    assert sensitivity.iloc[0].summary_status == 'insufficient_cells'
    assert pd.isna(comparison.iloc[0].ratio_red_green_median_change_percent)
