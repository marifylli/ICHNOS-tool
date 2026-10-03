"""ichnos_image — image processing package for the ICHNOS pipeline.

bright-field/DIC + GFP + mCherry images -> per-cell, per-timepoint CSV rows,
covering pipeline Stages 2-5 and 8:

    segment.py   Stage 2  bright-field/DIC -> per-cell label masks
    correct.py   Stages 1,3,4  background, flat-field, co-registration, crosstalk, photobleaching
    extract.py   Stage 5  masks + corrected channels -> per-cell features (+ raw ratio)
    export.py    Stage 8  per-cell records -> CSV, validated against schema.CellRecord
    synthesize.py          real single-GFP images -> synthetic dual-channel data with
                            known crosstalk/ratio ground truth, for training/validating
                            correct.unmix_crosstalk() before real dual-channel wet-lab
                            data exists
    pipeline.py             orchestrates segment/correct/extract/export over a batch of
                            image sets into one combined experiment CSV
    metadata.py             reads acquisition metadata (exposure, OME-XML) from a real
                            microscope TIFF via tifffile (optional dependency)

Stage 6 (FRET/maturation-kinetics ratio correction, c_session calibration)
and Stage 7 (dataset-level QC policy) live partly here (raw ratio, per-image
QC flags) and partly in the core `ichnos` package once it exists.
"""
from .segment import segment_cells, load_cellpose_model, border_touching_labels, focus_score
from .correct import (
    subtract_background,
    pixel_size_at_sample_um,
    suggest_rolling_ball_radius,
    rolling_ball_radius_for_objective,
    estimate_flat_field,
    flat_field_correct,
    estimate_registration_shift,
    apply_shift,
    unmix_crosstalk,
    estimate_crosstalk_coefficient,
    calibrate_crosstalk_from_control,
    correct_photobleaching,
)
from .extract import extract_per_cell, compute_ratio, CellFeatures
from .export import build_records, export_csv
from .schema import CellRecord, CSV_COLUMNS
from .image_io import (
    load_image,
    load_png16,
    save_png16,
    extract_fluorescence_plane,
    plane_for_channel,
    saturation_value_for_file,
    UncalibratedExtractionError,
)
from .synthesize import make_synthetic_pair, build_dataset
from . import instrument
from .pipeline import (
    ImageSet,
    process_image_set,
    process_experiment,
    CrosstalkCalibrationWarning,
)
from .metadata import read_tiff_metadata

__all__ = [
    "segment_cells",
    "load_cellpose_model",
    "border_touching_labels",
    "focus_score",
    "subtract_background",
    "pixel_size_at_sample_um",
    "suggest_rolling_ball_radius",
    "rolling_ball_radius_for_objective",
    "estimate_flat_field",
    "flat_field_correct",
    "estimate_registration_shift",
    "apply_shift",
    "unmix_crosstalk",
    "estimate_crosstalk_coefficient",
    "calibrate_crosstalk_from_control",
    "correct_photobleaching",
    "extract_per_cell",
    "compute_ratio",
    "CellFeatures",
    "build_records",
    "export_csv",
    "CellRecord",
    "CSV_COLUMNS",
    "load_png16",
    "save_png16",
    "make_synthetic_pair",
    "build_dataset",
    "ImageSet",
    "process_image_set",
    "process_experiment",
    "read_tiff_metadata",
]
