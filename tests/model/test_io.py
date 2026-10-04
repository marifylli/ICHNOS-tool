import json
from datetime import datetime
from pathlib import Path
from uuid import UUID

import pytest

from ichnos import io


class FixedDatetime(datetime):
    @classmethod
    def now(cls, tz=None):
        return cls(2026, 10, 4, 12, 0, 0, tzinfo=tz)


@pytest.fixture
def frozen_export(monkeypatch):
    monkeypatch.setattr(io, "datetime", FixedDatetime)
    monkeypatch.setattr(
        io, "_collect_shared_param_snapshot", lambda model: {}
    )


@pytest.mark.parametrize("variant", ["ox", "er"])
def test_exports_at_identical_time_preserve_both_archives(
    tmp_path, frozen_export, variant
):
    first = io.save_merged_sbml(
        "first model", variant, None, {}, tmp_path
    )
    second = io.save_merged_sbml(
        "second model", variant, None, {}, tmp_path
    )

    assert first[0] == second[0]
    assert first[1] != second[1]
    assert first[2] != second[2]

    assert Path(first[1]).read_text(encoding="utf-8") == "first model"
    assert Path(second[1]).read_text(encoding="utf-8") == "second model"
    assert Path(second[0]).read_text(encoding="utf-8") == "second model"

    for paths in (first, second):
        manifest = json.loads(
            Path(paths[2]).read_text(encoding="utf-8")
        )
        assert manifest["variant"] == variant
        assert manifest["built_at"] == "2026-10-04T12:00:00.000000"


@pytest.mark.parametrize("variant", ["ox", "er"])
def test_archive_collision_preserves_existing_files(
    tmp_path, frozen_export, monkeypatch, variant
):
    monkeypatch.setattr(io, "uuid4", lambda: UUID(int=1))

    paths = io.save_merged_sbml(
        "original model", variant, None, {}, tmp_path
    )
    original_contents = [
        Path(path).read_bytes() for path in paths
    ]

    with pytest.raises(FileExistsError):
        io.save_merged_sbml(
            "replacement model", variant, None, {}, tmp_path
        )

    assert [
        Path(path).read_bytes() for path in paths
    ] == original_contents
