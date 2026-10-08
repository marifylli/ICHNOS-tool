"""The bright-field frames have to be found before they can be used.

A manifest built before anyone noticed them has no bright_field_path, and a
folder can hold several bright-field files with nothing in their names to
say which field each belongs to. The script discovers the candidates, keeps
the one that fits, and leaves the rows where none fits without one -- which
is the real shape of the team's session: one bright-field per folder, so the
folders imaged at two fields have a field it does not cover.
"""
import runpy
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from PIL import Image
from scipy import ndimage as ndi
from skimage import morphology

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "align_brightfield.py"


@pytest.fixture(scope="module")
def script():
    return runpy.run_path(str(SCRIPT))


def cells_at(centres, shape=(192, 240), radius=8):
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    mask = np.zeros(shape, dtype=bool)
    for cy, cx in centres:
        mask |= (yy - cy) ** 2 + (xx - cx) ** 2 <= radius**2
    return mask


def brightfield_of(cells, seed=0):
    rng = np.random.default_rng(seed)
    rim = (ndi.binary_dilation(cells, morphology.disk(2))
           & ~ndi.binary_erosion(cells, morphology.disk(1)))
    halo = ndi.binary_dilation(cells, morphology.disk(4)) & ~cells
    image = np.full(cells.shape, 150.0)
    image[halo] += 25
    image[rim] -= 55
    return image + rng.normal(0, 1.5, cells.shape)


def fluorescence_of(cells, seed=0, amplitude=80.0):
    rng = np.random.default_rng(seed)
    image = 40 + rng.normal(0, 2.0, cells.shape)
    image[cells] += amplitude
    return image


def scattered(seed, shape=(192, 240), n=22, radius=8):
    """A different arrangement of cells for every field.

    Not a shared grid with an offset: fields that look alike make the
    'unrelated' comparison score as high as the matching one, and the test
    would then be measuring the fixture rather than the code.
    """
    rng = np.random.default_rng(seed)
    centres = [(int(rng.integers(radius * 2, shape[0] - radius * 2)),
                int(rng.integers(radius * 2, shape[1] - radius * 2))) for _ in range(n)]
    return cells_at(centres, shape=shape, radius=radius)


def build_session(root: Path):
    """Two folders. The second holds two fields and one bright-field, which
    belongs to the second of them."""
    rows = []
    for folder, fields, shift in (("d0", ["01a"], (7, -5)), ("d1", ["01a", "01b"], (0, 0))):
        (root / folder).mkdir(parents=True)
        for n, field in enumerate(fields):
            # A fixed seed per field, not hash(): string hashing is salted
            # per process, which would make this fixture -- and so the test
            # -- quietly different on every run.
            cells = scattered(seed=sum(ord(c) for c in folder + field))
            lit = np.roll(np.roll(cells, shift[0], axis=0), shift[1], axis=1)
            def save(array, name):
                Image.fromarray(np.clip(array, 0, 255).astype(np.uint8)).save(
                    root / folder / name)
                return f"{folder}/{name}"
            green = save(fluorescence_of(lit, seed=n), f"[{field}][Alexa 488].tif")
            red = save(fluorescence_of(lit, seed=n + 9, amplitude=95),
                       f"[{field}][Alexa 594].tif")
            if field == fields[-1]:
                save(brightfield_of(cells, seed=n), "Image000.tif")
            rows.append(dict(
                session_id="s", timepoint=len(rows), acquisition_order=len(rows),
                sample_id=f"{folder}-{field}", green_path=green, red_path=red,
                exposure_ms_green=100.0, exposure_ms_red=100.0, nd_filter_green=0.0,
                nd_filter_red=0.0, objective="60X", burner_hours=5.0,
                lamp_warmup_minutes=30.0))
    path = root / "images.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_candidates_exclude_the_fluorescence_pair(tmp_path, script):
    """By path, not by name: a glob wide enough to match the fluorescence
    files must not let an image become its own mask source."""
    (tmp_path / "f").mkdir()
    paths = []
    for name in ("[01a][Alexa 488].tif", "[01a][Alexa 594].tif", "Image000.tif"):
        path = tmp_path / "f" / name
        Image.fromarray(np.zeros((8, 8), dtype=np.uint8)).save(path)
        paths.append(path)
    found = script["candidates_for"](paths[:2], "*.tif")
    assert [p.name for p in found] == ["Image000.tif"]


def test_it_finds_the_frames_a_manifest_never_mentioned(tmp_path):
    manifest = build_session(tmp_path)
    assert "bright_field_path" not in pd.read_csv(manifest).columns

    done = subprocess.run(
        [sys.executable, str(SCRIPT), "--manifest", str(manifest),
         "--out", str(tmp_path / "aligned.csv"), "--min-size", "80", "--write-back"],
        capture_output=True, text=True,
    )
    assert done.returncode == 0, done.stderr
    written = pd.read_csv(manifest)
    assert "bright_field_path" in written.columns

    aligned = pd.read_csv(tmp_path / "aligned.csv")
    assert len(aligned) == 3
    assert (aligned.n_candidates == 1).all()


def test_the_field_without_its_own_brightfield_is_left_without_one(tmp_path):
    """Giving it the folder's bright-field anyway would put masks from
    another field on it, and nothing downstream would say so."""
    manifest = build_session(tmp_path)
    subprocess.run(
        [sys.executable, str(SCRIPT), "--manifest", str(manifest),
         "--out", str(tmp_path / "aligned.csv"), "--min-size", "80", "--write-back"],
        check=True, capture_output=True,
    )
    written = pd.read_csv(manifest).set_index("sample_id")
    kept = written.bright_field_path.notna() & (written.bright_field_path != "")
    assert kept["d0-01a"] and kept["d1-01b"]
    assert not kept["d1-01a"]
    assert pd.isna(written.loc["d1-01a", "brightfield_shift_dy"])


def test_a_glob_that_matches_nothing_says_so_instead_of_failing_obscurely(tmp_path):
    manifest = build_session(tmp_path)
    done = subprocess.run(
        [sys.executable, str(SCRIPT), "--manifest", str(manifest),
         "--out", str(tmp_path / "aligned.csv"), "--brightfield-glob", "nothing*.tif"],
        capture_output=True, text=True,
    )
    assert done.returncode != 0
    assert "--brightfield-glob" in (done.stdout + done.stderr)
