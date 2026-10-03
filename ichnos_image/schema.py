"""Re-exports the shared schema from the core `ichnos` package.

CellRecord/CSV_COLUMNS/validate used to be defined here as a placeholder,
before `ichnos.schema` existed. They now live there (the real, shared
definition); this module just re-exports them so existing
`from ichnos_image.schema import ...` / `from ichnos_image import CellRecord`
call sites keep working unchanged.
"""
from ichnos.schema import CellRecord, CSV_COLUMNS, validate

__all__ = ["CellRecord", "CSV_COLUMNS", "validate"]
