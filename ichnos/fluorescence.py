"""Model-to-fluorescence mapping.

Model channels are concentration-based signals, not camera intensity units.
The reporter already applies FRET and the ratio scale f.
"""

from __future__ import annotations

import numpy as np

from .simulate import SimulationResult


RATIO_DIRECTION = "red/green"

FLUORESCENCE_OBSERVABLES = (
    "Observed_Green",
    "Reporter_red",
    "Measured_Ratio_RG",
)


def model_fluorescence(result: SimulationResult) -> dict[str, np.ndarray]:
    """Return model signals without applying FRET or f a second time."""
    missing = [
        name for name in FLUORESCENCE_OBSERVABLES
        if not result.has(name)
    ]
    if missing:
        raise ValueError(
            f"missing fluorescence observables: {missing}"
        )

    mapping = {
        "green": "Observed_Green",
        "red": "Reporter_red",
        "ratio_red_green": "Measured_Ratio_RG",
    }

    signals = {}
    for channel, observable in mapping.items():
        values = np.asarray(result[observable], dtype=float)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(
                f"{observable} must contain finite non-negative values"
            )
        signals[channel] = values.copy()

    return signals