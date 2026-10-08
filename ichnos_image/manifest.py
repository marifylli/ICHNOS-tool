"""Image and control manifest readers shared by both CLIs."""
import json
from pathlib import Path
from skimage.color import rgb2gray

import numpy as np
import pandas as pd
from ichnos_image.instrument import SATURATION_VALUE
from ichnos.schema import validate_elapsed_hours, validate_optional_identifier

from ichnos_image import ImageSet
from ichnos_image.correct import calibrate_crosstalk_from_control
from ichnos_image.image_io import (
    load_image,
    plane_for_channel,
    saturation_mask_for_image,
)

def _load_bright_field(path: str | Path):
    """Convert bright-field images to a 2D segmentation plane."""
    image = load_image(path)

    if image.ndim == 2:
        return image

    if image.ndim == 3 and image.shape[-1] in (3, 4):
        return rgb2gray(image[..., :3])

    raise ValueError(f"unsupported bright-field shape: {image.shape}")

def _optional_elapsed_hours(row, name):
    value = getattr(row, name, None)
    if value is None or pd.isna(value):
        return None

    return validate_elapsed_hours(value, field_name=name)


def _optional_identifier(row, name):
    value = getattr(row, name, None)
    if value is None or pd.isna(value):
        return None
    return validate_optional_identifier(str(value), field_name=name)

_FALSE_WORDS = {"false", "0", "no", "n", "blank", "none"}
_TRUE_WORDS = {"true", "1", "yes", "y"}


def _optional_expect_cells(row) -> bool:
    """Read the optional expect_cells column; default True.

    Spelled out in words rather than left to pandas' truthiness because a
    manifest is hand-edited: "no" and "false" and an empty cell all turn up,
    and the quiet failure -- a blank control read as expecting cells --
    looks exactly like a segmentation failure in the QC output.
    """
    value = getattr(row, "expect_cells", None)
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return True
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip().lower()
    if text == "":
        return True
    if text in _FALSE_WORDS:
        return False
    if text in _TRUE_WORDS:
        return True
    raise ValueError(
        f"expect_cells must be one of {sorted(_TRUE_WORDS | _FALSE_WORDS)} or empty, got {value!r}"
    )


def _build_image_sets(
    manifest_path: Path,
    *,
    green_extraction: str | None = None,
    red_extraction: str | None = None,
    saturation_value: float = SATURATION_VALUE,
    image_loader=load_image,
) -> list[ImageSet]:
    manifest = pd.read_csv(
        manifest_path,
        dtype={"session_id": str, "sample_id": str, "condition_id": str, "specimen_id": str, "biological_replicate_id": str},
    )
    image_sets = []

    for row in manifest.itertuples():
        expect_cells = _optional_expect_cells(row)
        sample_id = _optional_identifier(row, "sample_id")
        condition_id = _optional_identifier(row, "condition_id")
        sampling_time_hours = _optional_elapsed_hours(
            row, "sampling_time_hours"
        )
        measurement_time_hours = _optional_elapsed_hours(
            row, "measurement_time_hours"
        )
        bright_field = (
            _load_bright_field(Path(manifest_path).parent / row.bright_field_path)
            if getattr(row, "bright_field_path", "")
            and pd.notna(row.bright_field_path)
            else None
        )
        def resolve(value):
            path = Path(value)
            return path if path.is_absolute() else Path(manifest_path).parent / path
        raw_green = image_loader(resolve(row.green_path))
        raw_red = image_loader(resolve(row.red_path))

        green_saturation = saturation_mask_for_image(
            raw_green, saturation_value
        )
        red_saturation = saturation_mask_for_image(
            raw_red, saturation_value
        )

        if green_saturation.shape != red_saturation.shape:
            raise ValueError("green and red images must have matching shapes")

        raw_saturation = green_saturation | red_saturation

        image_sets.append(
            ImageSet(
                source_paths=tuple(str(resolve(p)) for p in (row.green_path, row.red_path,
                    getattr(row, "bright_field_path", None)) if isinstance(p, str) and p),
                green=plane_for_channel(
                    raw_green,
                    "green",
                    method=green_extraction,
                ),
                red=plane_for_channel(
                    raw_red,
                    "red",
                    method=red_extraction,
                ),
                saturation_value=saturation_value,
                raw_saturation_mask=raw_saturation,
                green_saturation_mask=green_saturation,
                red_saturation_mask=red_saturation,
                bright_field=bright_field,
                session_id=str(row.session_id),
                timepoint=int(row.timepoint),
                acquisition_order=int(row.acquisition_order),
                exposure_ms_green=float(row.exposure_ms_green),
                exposure_ms_red=float(row.exposure_ms_red),
                nd_filter_green=float(row.nd_filter_green),
                nd_filter_red=float(row.nd_filter_red),
                objective=str(row.objective),
                burner_hours=float(row.burner_hours),
                lamp_warmup_minutes=float(row.lamp_warmup_minutes),
                sampling_time_hours=sampling_time_hours,
                measurement_time_hours=measurement_time_hours,
                sample_id=sample_id,
                condition_id=condition_id,
                expect_cells=expect_cells,
                specimen_id=_optional_identifier(row, 'specimen_id'),
                biological_replicate_id=_optional_identifier(row, 'biological_replicate_id'),
                acquisition_json=json.dumps(dict(
                    objective=str(row.objective), green_units='corrected_image_units',
                    extraction_green=green_extraction or ('scalar' if raw_green.ndim == 2 else None),
                    extraction_red=red_extraction or ('scalar' if raw_red.ndim == 2 else None),
                    gain_setting=_optional_identifier(row, 'gain_setting'),
                    exposure_ms_green=float(row.exposure_ms_green),
                    exposure_ms_red=float(row.exposure_ms_red),
                    nd_filter_green=float(row.nd_filter_green), nd_filter_red=float(row.nd_filter_red),
                ), sort_keys=True, allow_nan=False),
            )
        )

    return image_sets


def _calibrate_bleed_per_session(
    controls_path: Path | None,
    sessions: set[str],
    bleed_default: float | None,
    *,
    green_extraction: str | None = None,
    red_extraction: str | None = None,
    image_loader=load_image, calibrator=calibrate_crosstalk_from_control,
) -> dict:
    bleed_by_session = {}
    if controls_path is not None:
        controls = pd.read_csv(controls_path, dtype={"session_id": str})
        if controls.session_id.duplicated().any():
            raise ValueError("duplicate control session_id")
        for row in controls.itertuples():
            control_green = plane_for_channel(
                image_loader(Path(controls_path).parent / row.control_green_path),
                "green",
                method=green_extraction,
            )
            control_red = plane_for_channel(
                image_loader(Path(controls_path).parent / row.control_red_path),
                "red",
                method=red_extraction,
            )
            calibration = calibrator(control_green, control_red)
            bleed_by_session[str(row.session_id)] = calibration["bleed_green_to_red"]
            print(
                f"calibrated session={row.session_id}: "
                f"bleed_green_to_red={calibration['bleed_green_to_red']:.4f} "
                f"residual_check={calibration['residual_check']:.2f}"
            )

    missing = sessions - set(bleed_by_session)
    if missing:
        if bleed_default is None:
            raise SystemExit(
                f"no control and no --bleed-default for session(s): {sorted(missing)}. "
                "Add them to --controls or pass --bleed-default."
            )
        for session_id in missing:
            bleed_by_session[session_id] = bleed_default

    return bleed_by_session
