from dataclasses import asdict
import json
import runpy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from ichnos.calibration import fit_session_calibration, save_session_calibration
from ichnos_image.pipeline import ImageSet, process_image_set
from ichnos_image.export import export_csv


ROOT = Path(__file__).resolve().parents[2]


def image_records():
    yy, xx = np.mgrid[:128, :128]
    mask = np.zeros((128, 128), dtype=bool)
    for cy, cx in ((32, 32), (64, 96), (96, 32)):
        mask |= (yy - cy) ** 2 + (xx - cx) ** 2 < 10 ** 2
    green = np.where(mask, 1000., 0.)
    red = .6 * green
    bright_field = np.where(mask, .2, .5)
    images = ImageSet(
        green=green, red=red, bright_field=bright_field,
        session_id="001", sample_id="002", condition_id="003",
        timepoint=1, acquisition_order=1,
        exposure_ms_green=100, exposure_ms_red=200,
        nd_filter_green=0, nd_filter_red=0, objective="40X",
        burner_hours=1, lamp_warmup_minutes=30,
        sampling_time_hours=.5, measurement_time_hours=.75,
    )
    return process_image_set(images, bleed_green_to_red=0)


def test_images_to_cells_to_calibrated_sample_summary(tmp_path):
    records = image_records()
    assert len(records) == 3
    assert all(record.qc_pass for record in records)
    assert all(record.sample_id == "002" for record in records)
    original_csv = tmp_path / "cells.csv"
    export_csv(records, original_csv)
    original_bytes = original_csv.read_bytes()
    calibration = fit_session_calibration(
        [.1, .2, .3], [.15, .3, .45], session_id="001",
        reference_id="synthetic-model", source_kind="synthetic",
        source="known scale test",
    )
    save_session_calibration(calibration, tmp_path / "calibration.json")
    manifest = tmp_path / "calibrations.csv"
    pd.DataFrame([{
        "session_id": "001", "reference_id": "synthetic-model",
        "calibration_path": "calibration.json",
    }]).to_csv(manifest, index=False)
    run = runpy.run_path(str(ROOT / "scripts/summarize_cells.py"))["run"]
    args = SimpleNamespace(
        cells=original_csv, out_dir=tmp_path / "summary", data_kind="synthetic",
        calibrations=manifest, min_cells=3, green_floor=0,
    )
    run(args)
    output = pd.read_csv(args.out_dir / "samples.csv", dtype={"sample_id": str})
    assert len(output) == 1
    assert output.iloc[0].sample_id == "002"
    assert output.iloc[0].n_cells_used == 3
    assert output.iloc[0].ratio_red_green_median == pytest.approx(.6)
    assert output.iloc[0].calibrated_ratio_red_green_median == pytest.approx(.4)
    metadata = json.loads((args.out_dir / "metadata.json").read_text())
    assert metadata["calibrations"]["001"]["source_kind"] == "synthetic"
    assert metadata["experimental_instrument_validated"] is False
    assert metadata["f_applied_again"] is False
    assert original_csv.read_bytes() == original_bytes
    with pytest.raises(FileExistsError):
        run(args)
    args.out_dir = tmp_path / "rejected"
    args.data_kind = "experimental"
    with pytest.raises(ValueError, match="synthetic calibration"):
        run(args)
    assert not args.out_dir.exists()


@pytest.mark.parametrize("name", ["sample_id", "condition_id"])
def test_cell_records_reject_blank_identifiers(name):
    fields = asdict(image_records()[0])
    # ImageSet and CellRecord share these identifier checks.
    from ichnos.schema import CellRecord, validate
    fields[name] = " "
    with pytest.raises(ValueError, match=name):
        validate(CellRecord(**fields))
