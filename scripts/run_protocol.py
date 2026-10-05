"""Run one stress protocol and export model outputs and assumptions."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
from dataclasses import asdict
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import libsbml

from ichnos import build, naming, params, protocol, simulate


def _package_version(name):
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def _sha256(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def make_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", required=True, choices=["ox", "er"])
    parser.add_argument("--profile", default="default")
    parser.add_argument("--dose", required=True, type=float)
    parser.add_argument("--dose-units", required=True, choices=["uM"])
    parser.add_argument(
        "--times-hours",
        required=True,
        nargs="+",
        type=float,
        help="Strictly increasing observation times after stress onset",
    )
    parser.add_argument(
        "--initialization",
        required=True,
        choices=["equilibrium", "finite-preincubation"],
    )
    parser.add_argument("--preincubation-hours", type=float)
    parser.add_argument(
        "--equilibration-horizon-hours",
        type=float,
        default=50.0,
        help="Computational convergence horizon; not a culture duration",
    )
    parser.add_argument(
        "--equilibration-tolerance",
        type=float,
        default=1e-6,
    )
    parser.add_argument(
        "--clearance-rate-per-hour",
        type=float,
        help="Explicit assumed first-order clearance rate in h^-1",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        type=Path,
        help="New output directory; existing paths are rejected",
    )
    return parser


def run(args):
    if args.out_dir.exists():
        raise FileExistsError(
            f"output path already exists: {args.out_dir}"
        )

    if args.initialization == "finite-preincubation":
        if args.preincubation_hours is None:
            raise ValueError(
                "finite-preincubation requires --preincubation-hours"
            )
        duration = args.preincubation_hours
        if not (0 < duration < float("inf")):
            raise ValueError(
                "--preincubation-hours must be positive and finite"
            )
    elif args.preincubation_hours is not None:
        raise ValueError(
            "--preincubation-hours requires finite-preincubation"
        )

    times = simulate.validate_observation_times(args.times_hours)

    stress_protocol = protocol.StressProtocol(
        variant=args.variant,
        dose=args.dose,
        dose_units=args.dose_units,
        clears=args.clearance_rate_per_hour is not None,
        clearance_rate_per_hour=args.clearance_rate_per_hour,
        equilibration_horizon_hours=args.equilibration_horizon_hours,
        equilibration_tolerance=args.equilibration_tolerance,
    )

    profile = params.load_profile(args.variant, args.profile)
    baseline = build.build_variant_sbml_string(
        args.variant,
        save_sbml=False,
    )
    document = libsbml.readSBMLFromString(baseline)
    model = document.getModel()
    if model is None:
        raise ValueError("builder returned SBML without a model")

    # Apply the complete profile before adding an optional stress rate rule.
    params.apply_profile(model, profile)
    baseline = libsbml.writeSBMLToString(document)

    runner, loaded_model = protocol.load_protocol_model(
        baseline,
        stress_protocol,
    )
    names = naming.id_to_name_map(loaded_model)

    # Export stress even when it is a constant parameter.
    stress_id = naming.resolve_id(
        loaded_model,
        stress_protocol.stress_parameter(),
        kinds=("parameter",),
    )
    selections = list(runner.timeCourseSelections)
    if stress_id not in selections:
        selections.append(stress_id)
    runner.timeCourseSelections = selections

    # Snapshot before preparation or dose application.
    prepared_sbml = runner.getCurrentSBML()

    if args.initialization == "finite-preincubation":
        result, _ = protocol.run_protocol_after_preincubation(
            runner,
            stress_protocol,
            preincubation_hours=args.preincubation_hours,
            times_hours=times,
            id_to_name=names,
        )
    else:
        result, _ = protocol.run_protocol_at_times(
            runner,
            stress_protocol,
            times_hours=times,
            id_to_name=names,
            require_equilibrium=True,
        )

    metadata = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "parameter_profile": asdict(profile),
        "protocol": asdict(stress_protocol),
        "observation_times_hours": result.time.tolist(),
        "exposure": result.exposure,
        "initialization": result.initialization,
        "solver_used": asdict(result.solver),
        "integrator": runner.getIntegrator().getName(),
        "csv_columns": ["time_hours", *result.columns],
        "model_sbml_sha256": _sha256(prepared_sbml),
        "runner_script_sha256": hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        "runtime": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "packages": {
                name: _package_version(name)
                for name in (
                    "ichnos-tool",
                    "numpy",
                    "python-libsbml",
                    "libroadrunner",
                    "tellurium",
                    "PyYAML",
                )
            },
        },
        "interpretation": {
            "outputs": "model observables, not calibrated fluorescence",
            "initial_state": "source SBML initial conditions",
            "profile_initial_state": (
                "documented profile values; no separate initial-state "
                "override is applied by this script"
            ),
            "galactose_switching_explicitly_modelled": False,
            "experimentally_validated": False,
        },
    }

    # Reject non-standard JSON numbers before creating output files.
    metadata_text = json.dumps(
        metadata,
        indent=2,
        sort_keys=True,
        allow_nan=False,
    ) + "\n"

    args.out_dir.parent.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir()

    with (args.out_dir / "model.sbml").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        handle.write(prepared_sbml)

    with (args.out_dir / "results.csv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_hours", *result.columns])
        for time, values in zip(result.time, result.values):
            writer.writerow([float(time), *map(float, values)])

    with (args.out_dir / "metadata.json").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        handle.write(metadata_text)

    print(f"Saved protocol run to {args.out_dir}")


def main():
    parser = make_parser()
    args = parser.parse_args()
    try:
        run(args)
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()