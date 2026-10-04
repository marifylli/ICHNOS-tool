# ICHNOS tool

Software under development for converting fluorescence microscopy images
of the ICHNOS yeast stress biosensor into estimates of stress dose and
time since stress onset, with quality-control flags and uncertainty.

The current implementation processes images and simulates the ox and er
models. It does not yet estimate dose or time from experimental images.

The research repositories (`Ichnos_PULSE`, `ICHNOS-ablation`,
`ichnos-fisher`, `DryLabTool`) retain the supporting studies.
Migrated code and model sources are recorded in
[docs/PROVENANCE.md](docs/PROVENANCE.md).

## Status

| Component | State |
| --- | --- |
| Image pipeline (`ichnos_image`) | Implemented; computational tests pass |
| RGB extraction | Explicit CLI/API selection; experimental mapping unresolved |
| Saturation QC | Checked before RGB extraction and forwarded through the pipeline |
| Shared schema (`ichnos.schema`) | Migrated; schema versioning pending |
| SBML models and builder (`ichnos.build`) | Migrated for ox and er |
| Shared simulator (`ichnos.simulate`) | Implemented with solver settings and observable selections |
| Parameter profiles (`ichnos.params`) | Packaged; value, unit and provenance checks implemented |
| Model exports (`ichnos.io`) | Unique archives with manifests and overwrite protection |
| Experimental exposure protocol | Constant-stress baseline; clearance support pending |
| Measurement mapping and session calibration | Not complete |
| Calibration artifact and decoder | Not implemented |
| Estimator uncertainty | Not implemented |

Passing tests verify software behaviour under controlled assumptions.
They do not establish experimental accuracy for dose or time estimation.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[image,model,dev]"
python -m pytest -q
```

The optional Cellpose backend requires the `cellpose` extra.

## Image pipeline CLI

```bash
python scripts/run_pipeline.py --help
```

For RGB acquisitions, explicitly select the extraction method for each
channel using `--green-extraction` and `--red-extraction`.
The same choices are applied to samples and crosstalk controls.

The CLI saturation threshold defaults to the instrument profile and can
be overridden with `--saturation-value`. It is applied to stored image
components before fluorescence extraction.

Explicit extraction choices are configuration options, not evidence
that those choices have been experimentally calibrated.

## Next implementation steps

1. Align experimental metadata with the simulator: variant, initial dose,
   actual elapsed times, exposure history and initial model state.
2. Add optional clearance with explicit assumptions and parameter provenance.
3. Implement the mapping from model observables to measured fluorescence,
   session calibration and population grouping.
4. Build a calibration artifact and decoder for one known variant.
5. Evaluate dose/time ambiguity and estimator uncertainty using an
   appropriate measurement-noise model.
6. Extend the validated workflow to the second variant.

Copper/CuSO4 modelling is outside the current implementation scope.
Its images can be processed without a copper-specific model or decoder.

See [docs/limitations.md](docs/limitations.md) for limitations and
unresolved measurement decisions.
