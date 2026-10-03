# ICHNOS tool

Turns fluorescence microscopy images of the ICHNOS yeast stress biosensor into
an estimate of **stress dose** and **time since stress onset**, with QC flags
and stated uncertainty.

This repository is the *product*. The research repositories
(`Ichnos_PULSE`, `ICHNOS-ablation`, `ichnos-fisher`, `DryLabTool`) keep the
studies and stay unchanged; code is migrated here with its provenance
recorded in [`docs/PROVENANCE.md`](docs/PROVENANCE.md).

## Status

| Stage | State |
| --- | --- |
| Image pipeline (`ichnos_image`) | migrated, tests green |
| Shared schema (`ichnos.schema`) | migrated as-is, versioning not yet added |
| SBML model + builder (`ichnos.build`) | not yet migrated — Step 3 |
| Calibration / decoder (`ichnos.calibrate`, `ichnos.decode`) | **not implemented** |
| Uncertainty | **not implemented** |

No decoder exists yet. Nothing in this repository estimates dose or time
today; the name describes the target, not the current capability.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[image,model,dev]"
python -m pytest -q
```

## Implementation order

1. **Step 1 — baseline migration** *(done)*: image pipeline and shared schema
   moved here verbatim, test-to-script dependency removed, suite green.
2. Step 2 — split image config into `instrument.py`, extract the image loader
   into `image_io.py`, apply the confirmed QImaging/CKX41 instrument profile.
3. Step 3 — migrate the four SBML files and the builder from `Ichnos_PULSE`,
   with no parameter changes.
4. Step 4 — one shared simulator, naming, parameter profiles, protocol.
5. Step 5 — measurement mapping and session calibration, population grouping.
6. Step 6 — calibration artifact and decoder for one known variant.
7. Step 7 — uncertainty with a validated noise model; second variant.

Model or protocol changes always go in a separate commit from file moves.
