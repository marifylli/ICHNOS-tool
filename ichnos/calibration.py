"""Relative red/green session calibration, separate from image corrections.

image_ratio = c_session * Measured_Ratio_RG.
Measured_Ratio_RG already includes the reporter's FRET and f.
This module neither estimates f nor calibrates absolute channel intensities.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from .fluorescence import RATIO_DIRECTION


def _nonempty(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value


def _ratios(values, name, *, positive=False):
    raw = np.asarray(values, dtype=object)
    if raw.ndim != 1 or not raw.size:
        raise ValueError(f"{name} must be a non-empty one-dimensional sequence")
    if any(isinstance(value, (bool, np.bool_)) for value in raw):
        raise ValueError(f"{name} must not contain booleans")
    try:
        array = raw.astype(float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must contain numbers") from exc
    if not np.isfinite(array).all() or (array < 0).any():
        raise ValueError(f"{name} must contain finite non-negative ratios")
    if positive and (array <= 0).any():
        raise ValueError(f"{name} must contain strictly positive reference ratios")
    return array


def model_reference_id(metadata: dict) -> str:
    """Bind calibration to the forward model, f, solver and protocol.

    Ignore timestamps, file paths and software host details. Keep observation
    times, initialization, exposure and the complete parameter profile.
    """
    reference = {
        key: metadata[key]
        for key in (
            "model_sbml_sha256", "parameter_profile", "protocol",
            "observation_times_hours", "initialization", "solver_used",
            "fluorescence",
        )
    }
    encoded = json.dumps(reference, sort_keys=True, allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True)
class SessionCalibration:
    session_id: str
    reference_id: str
    c_session: float
    source_kind: str
    source: str
    model_ratios: tuple[float, ...]
    image_ratios: tuple[float, ...]
    rmse: float

    def __post_init__(self):
        for name in ("session_id", "reference_id", "source"):
            _nonempty(getattr(self, name), name)
        if self.source_kind not in {"synthetic", "experimental_reference"}:
            raise ValueError("source_kind must be synthetic or experimental_reference")
        for name in ("c_session", "rmse"):
            raw = getattr(self, name)
            if isinstance(raw, (bool, np.bool_)):
                raise ValueError(f"{name} must be numeric")
            value = float(raw)
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and non-negative")
            object.__setattr__(self, name, value)
        if self.c_session <= 0:
            raise ValueError("c_session must be strictly positive")
        x = _ratios(self.model_ratios, "model_ratios", positive=True)
        y = _ratios(self.image_ratios, "image_ratios")
        if x.shape != y.shape or x.size < 3:
            raise ValueError("need at least three matched reference pairs")
        object.__setattr__(self, "model_ratios", tuple(map(float, x)))
        object.__setattr__(self, "image_ratios", tuple(map(float, y)))

    def to_dict(self) -> dict:
        return {
            "schema_version": 1,
            **asdict(self),
            "ratio_direction": RATIO_DIRECTION,
            "model_observable": "Measured_Ratio_RG",
            "equation": "image_ratio = c_session * Measured_Ratio_RG",
            "method": "unweighted least squares through origin",
            "n_reference_pairs": len(self.model_ratios),
            "f_already_in_model_ratio": True,
            "experimental_instrument_validated": False,
            "purpose": (
                "computational verification with synthetic reference ratios; "
                "not microscope calibration"
                if self.source_kind == "synthetic" else
                "relative session scale from experimental references; "
                "instrument linearity and model validity require separate validation"
            ),
        }


def fit_session_calibration(
    model_ratios, image_ratios, *, session_id: str,
    reference_id: str, source_kind: str, source: str,
) -> SessionCalibration:
    """Fit one multiplicative scale from matched reference ratios.

    Reference pairs must have the same biological condition and elapsed time.
    Use QC-passing sample summaries; do not treat cells as independent replicates.
    Positive model references avoid fitting an undefined/no-signal ratio.
    """
    x = _ratios(model_ratios, "model_ratios", positive=True)
    y = _ratios(image_ratios, "image_ratios")
    if x.shape != y.shape or x.size < 3:
        raise ValueError("need at least three matched reference pairs")
    # Scale before dot products to avoid squaring large reference values.
    scale = float(x.max())
    normalized = x / scale
    coefficient = float(np.dot(normalized, y) / np.dot(normalized, normalized) / scale)
    residual = y - coefficient * x
    rmse = float(np.linalg.norm(residual / math.sqrt(len(x))))
    return SessionCalibration(
        session_id, reference_id, coefficient, source_kind, source,
        tuple(x), tuple(y), rmse,
    )


def calibrate_image_ratios(
    image_ratios, calibration: SessionCalibration, *,
    session_id: str, reference_id: str,
) -> np.ndarray:
    """Convert uncalibrated image ratios to the recorded model ratio scale.

    Pass original image ratios, never a previously calibrated result.
    The input is preserved and the caller retains both columns.
    """
    if session_id != calibration.session_id:
        raise ValueError("session_id differs from the calibration session")
    if reference_id != calibration.reference_id:
        raise ValueError("reference_id differs from the calibration forward model")
    ratios = _ratios(image_ratios, "image_ratios")
    corrected = ratios / calibration.c_session
    if not np.isfinite(corrected).all():
        raise ValueError("calibrated ratios are not finite")
    return corrected


def save_session_calibration(calibration: SessionCalibration, path) -> Path:
    """Write an exclusive JSON artifact with all fitted reference pairs."""
    path = Path(path)
    encoded = json.dumps(calibration.to_dict(), indent=2, allow_nan=False) + "\n"
    with path.open("x", encoding="utf-8") as handle:
        handle.write(encoded)
    return path


def load_session_calibration(path) -> SessionCalibration:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return session_calibration_from_dict(data)


def session_calibration_from_dict(data: dict) -> SessionCalibration:
    """Validate an embedded artifact using the same checks as file loading."""
    if data.get("schema_version") != 1 or data.get("ratio_direction") != RATIO_DIRECTION:
        raise ValueError("unsupported calibration schema or ratio direction")
    calibration = SessionCalibration(**{
        name: data[name] for name in SessionCalibration.__dataclass_fields__
    })
    expected = calibration.to_dict()
    # Reject inconsistent method/provenance rather than silently relabelling it.
    for key, value in expected.items():
        if key not in {"model_ratios", "image_ratios"} and data.get(key) != value:
            raise ValueError(f"inconsistent calibration metadata: {key}")
    fitted = fit_session_calibration(
        calibration.model_ratios, calibration.image_ratios,
        session_id=calibration.session_id, reference_id=calibration.reference_id,
        source_kind=calibration.source_kind, source=calibration.source,
    )
    if not math.isclose(fitted.c_session, calibration.c_session, rel_tol=1e-12):
        raise ValueError("c_session differs from the recorded reference fit")
    if not math.isclose(fitted.rmse, calibration.rmse, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("rmse differs from the recorded reference fit")
    return calibration
