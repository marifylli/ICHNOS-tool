"""ichnos -- the core package: the shared data contract, the model, and
(once they exist) calibration and decoding.

    schema.py        CellRecord / CSV_COLUMNS: the one definition of a row of
                     the per-cell dataset. ichnos_image writes these; anything
                     reading that CSV imports this rather than re-deriving it
                     from the header.
    config.py        Model configuration -- which SBML files make up each
                     variant, shared-parameter rules, unsupported variants.
    build.py         Merges TIP-TetR, the reporter and one sensing module into
                     one SBML model.
    merge_checks.py  Unit, collision and shared-parameter checks run during
                     the merge.
    io.py            Saving a merged model with its manifest.
    models/          The four packaged SBML sources plus manifest.json.

Nothing is imported eagerly here. ichnos_image needs only ichnos.schema, and
importing it must not pull in libsbml or a solver.

Not implemented yet: calibrate.py, decode.py, uncertainty.py. There is no
decoder, so no dose or elapsed-time estimate can be produced.
"""
