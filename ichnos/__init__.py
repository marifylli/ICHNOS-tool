"""ichnos — core package shared across the ICHNOS pipeline's per-domain
packages (ichnos_image today; a future ratio/decoder package for Stage 6
would live alongside it here).

    schema.py   CellRecord / CSV_COLUMNS: the one definition of what a row
                of the final per-cell dataset looks like -- ichnos_image
                (and anything downstream that reads its CSV output) imports
                this rather than defining its own copy.
    config.py   Calibrated pipeline constants (QC thresholds, protocol
                numbers) with their provenance, so a number isn't just a
                hardcoded literal somewhere inside a function.

Deliberately excludes Stage 6 (ratio/FRET decoder) for now -- out of scope
until that stage is tackled.
"""
