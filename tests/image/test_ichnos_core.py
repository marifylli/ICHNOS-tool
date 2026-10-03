"""Core ichnos package: shared schema + config are importable, and
ichnos_image actually uses them (not duplicated copies)."""
from ichnos import config
from ichnos.schema import CellRecord, CSV_COLUMNS
from ichnos_image import export
from ichnos_image.schema import CellRecord as ImageCellRecord


def test_ichnos_image_reexports_the_same_schema_class():
    assert ImageCellRecord is CellRecord


def test_export_defaults_match_config():
    import inspect

    sig = inspect.signature(export.build_records)
    assert sig.parameters["focus_score_threshold"].default == config.FOCUS_SCORE_THRESHOLD
    assert (
        sig.parameters["registration_shift_threshold_px"].default
        == config.REGISTRATION_SHIFT_THRESHOLD_PX
    )
    assert (
        sig.parameters["lamp_warmup_threshold_minutes"].default
        == config.LAMP_WARMUP_THRESHOLD_MINUTES
    )


def test_csv_columns_match_dataclass_fields():
    from dataclasses import fields

    assert CSV_COLUMNS == [f.name for f in fields(CellRecord)]
