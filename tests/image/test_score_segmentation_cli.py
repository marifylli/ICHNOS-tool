"""scripts/score_segmentation.py end to end on a tiny reference set."""
import runpy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from PIL import Image

from ichnos_image import segment

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "score_segmentation.py"


def reference(tmp_path):
    frame = np.full((60, 80), 150, dtype=np.uint8)
    Image.fromarray(frame).save(tmp_path / "bf.tif")
    Image.fromarray(frame[10:50, 20:70]).save(tmp_path / "a.png")
    pd.DataFrame([{"crop": "a.png", "session_id": "s", "sample_id": "x",
                   "bright_field_path": str(tmp_path / "bf.tif"), "y0": 10, "x0": 20,
                   "height": 40, "width": 50, "pixel_size_um": 0.1725}]).to_csv(tmp_path / "crops.csv", index=False)
    pd.DataFrame({" ": [1, 2], "X": [10.0, 40.0], "Y": [10.0, 30.0]}).to_csv(tmp_path / "a.points.csv", index=False)


class FakeCellpose:
    """Stands in for the 1.1 GB model: one mask around the first point only."""

    calls = 0

    def eval(self, image, diameter=None, channels=None):
        FakeCellpose.calls += 1
        masks = np.zeros(image.shape, dtype=np.int32)
        masks[15:25, 25:35] = 1   # frame coordinates of crop point (10, 10)
        return masks, None, None


def test_cellpose_is_loaded_once_and_scored(tmp_path, monkeypatch):
    reference(tmp_path)
    loads = []
    monkeypatch.setattr(segment, "load_cellpose_model", lambda gpu=False: loads.append(gpu) or FakeCellpose())
    main = runpy.run_path(str(SCRIPT))["main"]
    out = tmp_path / "score.csv"
    assert main(["--reference", str(tmp_path), "--method", "cellpose", "--out", str(out)]) == 0
    assert loads == [False] and FakeCellpose.calls == 1
    row = pd.read_csv(out).set_index("crop").loc["a"]
    assert (row.marked, row.separated, row.false_masks) == (2, 1, 0)
    assert "cellpose_model" not in Path(str(out) + ".json").read_text()


def test_existing_output_is_refused(tmp_path):
    reference(tmp_path)
    out = tmp_path / "score.csv"
    out.write_text("keep")
    main = runpy.run_path(str(SCRIPT))["main"]
    with pytest.raises(SystemExit, match="already exists"):
        main(["--reference", str(tmp_path), "--method", "transmitted", "--out", str(out)])
