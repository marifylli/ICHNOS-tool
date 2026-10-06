"""Discrete dose compatibility at known times; no extrapolation or time fit."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import platform
from importlib.metadata import version

import numpy as np

from .simulate import validate_observation_times


def _vector(values, name, *, positive=False):
    raw = np.asarray(values, dtype=object)
    if raw.ndim != 1 or not raw.size:
        raise ValueError(f"{name} must be a non-empty one-dimensional sequence")
    if any(isinstance(value, (bool, np.bool_)) for value in raw):
        raise ValueError(f"{name} must not contain booleans")
    try:
        array = raw.astype(float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must contain numeric values") from exc
    if not np.isfinite(array).all() or (array < 0).any():
        raise ValueError(f"{name} must be finite and non-negative")
    if positive and (array <= 0).any():
        raise ValueError(f"{name} must be strictly positive")
    return array


def _digest(data):
    encoded = json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode()).hexdigest()


@dataclass
class DoseTable:
    doses_uM: np.ndarray
    times_hours: np.ndarray
    ratios: np.ndarray
    provenance: dict

    def __post_init__(self):
        self.doses_uM = _vector(self.doses_uM, "doses_uM")
        if self.doses_uM.size < 2 or (np.diff(self.doses_uM) <= 0).any():
            raise ValueError("dose grid must contain at least two strictly increasing doses")
        self.times_hours = validate_observation_times(self.times_hours)
        raw = np.asarray(self.ratios, dtype=object)
        if raw.shape != (len(self.doses_uM), len(self.times_hours)):
            raise ValueError("ratios must have shape (n_doses, n_times)")
        self.ratios = _vector(raw.ravel(), "ratios").reshape(raw.shape)
        if self.provenance.get("source_kind") != "model_generated":
            raise ValueError("table source_kind must be model_generated")
        if self.provenance.get("ratio_direction") != "red/green":
            raise ValueError("table ratio_direction must be red/green")
        if self.provenance.get("experimentally_validated") is not False:
            raise ValueError("this model-generated table is not experimentally validated")
        # Validate serializability before a run can be exported.
        json.dumps(self.provenance, allow_nan=False)

    def to_dict(self):
        payload = {
            "schema_version": 1, "doses_uM": self.doses_uM.tolist(),
            "times_hours": self.times_hours.tolist(), "ratios": self.ratios.tolist(),
            "provenance": self.provenance,
        }
        return {**payload, "artifact_id": _digest(payload)}


def save_dose_table(table: DoseTable, path) -> Path:
    path = Path(path)
    encoded = json.dumps(table.to_dict(), indent=2, allow_nan=False) + "\n"
    with path.open("x", encoding="utf-8") as handle:
        handle.write(encoded)
    return path


def load_dose_table(path) -> DoseTable:
    data = json.loads(Path(path).read_text())
    artifact_id = data.pop("artifact_id")
    if data.get("schema_version") != 1 or _digest(data) != artifact_id:
        raise ValueError("unsupported or modified dose artifact")
    return DoseTable(
        data["doses_uM"], data["times_hours"], data["ratios"], data["provenance"]
    )


def build_dose_table(
    *, variant: str, doses_uM, times_hours, initialization: str,
    profile_name: str = "default", preincubation_hours=None,
    clearance_rate_per_hour=None, equilibration_horizon_hours=50.,
    equilibration_tolerance=1e-6,
) -> DoseTable:
    """Generate each grid dose with a fresh runner and the shared protocol."""
    import libsbml
    from . import build, naming, params, protocol
    from .fluorescence import model_fluorescence

    doses = _vector(doses_uM, "doses_uM")
    if doses.size < 2 or (np.diff(doses) <= 0).any():
        raise ValueError("dose grid must contain at least two strictly increasing doses")
    times = validate_observation_times(times_hours)
    if initialization not in {"equilibrium", "finite-preincubation"}:
        raise ValueError("initialization must be equilibrium or finite-preincubation")
    if initialization == "equilibrium" and preincubation_hours is not None:
        raise ValueError("preincubation_hours requires finite-preincubation")
    if initialization == "finite-preincubation":
        duration = _vector([preincubation_hours], "preincubation_hours", positive=True)[0]
    else:
        duration = None
    profile = params.load_profile(variant, profile_name)
    sbml = build.build_variant_sbml_string(variant, save_sbml=False)
    document = libsbml.readSBMLFromString(sbml)
    model = document.getModel()
    params.apply_profile(model, profile)
    baseline = libsbml.writeSBMLToString(document)
    ratios = []
    records = []
    for dose in doses:
        exposure = protocol.StressProtocol(
            variant=variant, dose=float(dose),
            clears=clearance_rate_per_hour is not None,
            clearance_rate_per_hour=clearance_rate_per_hour,
            equilibration_horizon_hours=equilibration_horizon_hours,
            equilibration_tolerance=equilibration_tolerance,
        )
        runner, loaded_model = protocol.load_protocol_model(baseline, exposure)
        prepared_model_sha256 = hashlib.sha256(
            runner.getCurrentSBML().encode()
        ).hexdigest()
        names = naming.id_to_name_map(loaded_model)
        if initialization == "equilibrium":
            result, _ = protocol.run_protocol_at_times(
                runner, exposure, times_hours=times, id_to_name=names
            )
        else:
            result, _ = protocol.run_protocol_after_preincubation(
                runner, exposure, times_hours=times, id_to_name=names,
                preincubation_hours=float(duration),
            )
        ratios.append(model_fluorescence(result)["ratio_red_green"])
        records.append({
            "dose_uM": float(dose), "exposure": result.exposure,
            "initialization": result.initialization,
            "solver_used": asdict(result.solver),
            "prepared_model_sbml_sha256": prepared_model_sha256,
        })
    provenance = {
        "source_kind": "model_generated", "experimentally_validated": False,
        "variant": variant, "parameter_profile": asdict(profile),
        "model_sbml_sha256": hashlib.sha256(baseline.encode()).hexdigest(),
        "baseline_sbml": baseline,
        "ratio_direction": "red/green", "observable": "Measured_Ratio_RG",
        "f_already_in_model_ratio": True,
        "f": profile.parameters["f"].value,
        "eps": profile.parameters["eps"].value,
        "runs": records,
        "runtime": {
            "python": platform.python_version(),
            "packages": {name: version(name) for name in (
                "numpy", "python-libsbml", "libroadrunner", "tellurium", "PyYAML"
            )},
        },
        "code_sha256": {
            name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ("decoder.py", "protocol.py", "simulate.py", "fluorescence.py")
        },
        "scope": "discrete dose grid at known observation times; no fitted experimental response",
    }
    return DoseTable(doses, times, np.asarray(ratios), provenance)


def decode_dose(
    table: DoseTable, *, times_hours, calibrated_ratios, ratio_tolerance,
    green_measurements=None, green_floor=None,
) -> dict:
    """Keep every grid dose within explicit absolute ratio tolerances.

    A unique candidate is a discrete grid result, not proof of continuous
    identifiability. Tolerances are compatibility bounds, not confidence levels.
    Optional green values/floor must share the caller's measurement units.
    """
    times = validate_observation_times(times_hours)
    if not np.array_equal(times, table.times_hours):
        raise ValueError("observation times differ from the dose table")
    observations = _vector(calibrated_ratios, "calibrated_ratios")
    if observations.shape != times.shape:
        raise ValueError("one calibrated ratio is required per observation time")
    raw_tolerance = np.asarray(ratio_tolerance, dtype=object)
    if raw_tolerance.ndim == 0:
        tolerances = np.repeat(_vector([ratio_tolerance], "ratio_tolerance", positive=True), len(times))
    else:
        tolerances = _vector(ratio_tolerance, "ratio_tolerance", positive=True)
    if tolerances.shape != times.shape:
        raise ValueError("ratio_tolerance must be scalar or one value per time")
    if (green_measurements is None) != (green_floor is None):
        raise ValueError("green_measurements and green_floor must be supplied together")
    result = {
        "artifact_id": table.to_dict()["artifact_id"],
        "times_hours": times.tolist(), "calibrated_ratios": observations.tolist(),
        "ratio_tolerance": tolerances.tolist(), "dose_units": "uM",
        "grid_range_uM": [float(table.doses_uM[0]), float(table.doses_uM[-1])],
        "candidate_doses_uM": [], "dose_estimate_uM": None,
        "discrete_grid_only": True, "time_estimated": False,
        "experimentally_validated": False,
        "tolerance_is_confidence_interval": False,
        "detection_floor_checked": green_floor is not None,
    }
    if green_floor is not None:
        green = _vector(green_measurements, "green_measurements")
        floor = _vector([green_floor], "green_floor")[0]
        if green.shape != times.shape:
            raise ValueError("one green measurement is required per observation time")
        result.update(green_measurements=green.tolist(), green_floor=float(floor))
        if (green <= floor).any():
            return {**result, "status": "below_detection_floor"}
    residuals = np.abs(table.ratios - observations)
    compatible = (residuals <= tolerances).all(axis=1)
    candidates = table.doses_uM[compatible].tolist()
    result["candidate_doses_uM"] = candidates
    result["max_absolute_ratio_residual_per_grid_dose"] = residuals.max(axis=1).tolist()
    if len(candidates) == 1:
        return {**result, "status": "unique_grid_match", "dose_estimate_uM": candidates[0]}
    if len(candidates) > 1:
        return {**result, "status": "ambiguous"}
    outside = ((observations < table.ratios.min(axis=0) - tolerances)
               | (observations > table.ratios.max(axis=0) + tolerances)).any()
    return {**result, "status": "out_of_response_domain" if outside else "no_grid_match"}
