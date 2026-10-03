"""End-to-end smoke test on synthetic data: segment -> correct -> extract -> export.

Uses method="otsu" so this runs without cellpose/torch installed.
"""
import numpy as np
import pandas as pd

from ichnos_image import segment, correct, extract, export


def _make_synthetic_field(n_cells=3, size=128, seed=0, bleed_green_to_red=0.05):
    rng = np.random.default_rng(seed)
    bf = np.full((size, size), 0.5)
    green = np.zeros((size, size), dtype=np.float64)
    red = np.zeros((size, size), dtype=np.float64)
    yy, xx = np.mgrid[0:size, 0:size]

    centers = rng.integers(25, size - 25, size=(n_cells, 2))
    for cy, cx in centers:
        blob = ((yy - cy) ** 2 + (xx - cx) ** 2) < 10**2
        bf[blob] -= 0.3
        green[blob] = 1000.0
        red[blob] = 100.0

    red_with_bleed = red + bleed_green_to_red * green
    return bf, green, red_with_bleed


def test_end_to_end_smoke(tmp_path):
    bleed = 0.05
    bf, green_raw, red_raw = _make_synthetic_field(bleed_green_to_red=bleed)

    labels = segment.segment_cells(bf, method="otsu", min_size=30)
    assert labels.max() >= 1

    green_bg, _ = correct.subtract_background(green_raw)
    red_bg, _ = correct.subtract_background(red_raw)
    green_corr, red_corr = correct.unmix_crosstalk(green_bg, red_bg, bleed_green_to_red=bleed)

    features = extract.extract_per_cell(labels, green_raw, red_raw, green_corr, red_corr, saturation_value=65535.0)
    assert len(features) >= 1
    assert all(f.corrected_mean_red < f.raw_mean_red for f in features)

    edge_ids = segment.border_touching_labels(labels)
    records = export.build_records(
        features,
        session_id="test-session",
        timepoint=0,
        edge_flagged_ids=edge_ids,
        focus_score=8000.0,
        registration_shift_px=0.1,
        exposure_ms_green=800.0,
        exposure_ms_red=2000.0,
        nd_filter_green=0.0,
        nd_filter_red=0.0,
        objective="60x",
        burner_hours=10.0,
        lamp_warmup_minutes=30.0,
        acquisition_order=1,
    )
    assert len(records) == len(features)
    assert all(r.qc_pass for r in records)

    out_csv = tmp_path / "out.csv"
    export.export_csv(records, out_csv)
    assert out_csv.exists()

    df = pd.read_csv(out_csv)
    assert list(df.columns) == export.CSV_COLUMNS
    assert (df["ratio_red_green"] > 0).all()
