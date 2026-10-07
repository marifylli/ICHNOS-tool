"""Diagnostic reproduction of saturation-mask / registration alignment.

This script does not modify the ICHNOS image pipeline.
It demonstrates whether a red-channel saturation mask remains spatially
consistent after red-channel registration.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import shift as ndi_shift

from ichnos_image import correct


def main() -> None:
    shape = (32, 32)

    # Reference/green coordinate system.
    green = np.zeros(shape, dtype=float)
    green[10:18, 10:18] = 100.0

    # Moving/red channel is displaced four columns to the right.
    red = np.zeros(shape, dtype=float)
    red[10:18, 14:22] = 100.0

    # One saturated pixel belonging to the displaced red object.
    red[13, 20] = 255.0

    red_saturation = np.zeros(shape, dtype=bool)
    red_saturation[13, 20] = True

    # Known shift needed to align red to green.
    shift_rc = (0.0, -4.0)

    registered_red = correct.apply_shift(red, shift_rc)

    # Current problematic behaviour:
    # red is registered, but its saturation mask remains unregistered.
    current_mask = red_saturation

    # Expected geometry:
    # apply the same spatial transform to the discrete mask using
    # nearest-neighbour interpolation rather than image interpolation.
    registered_mask = ndi_shift(
        red_saturation.astype(np.uint8),
        shift=shift_rc,
        order=0,
        mode="constant",
        cval=0,
    ).astype(bool)

    cell_mask = np.zeros(shape, dtype=bool)
    cell_mask[10:18, 10:18] = True

    current_sat_flag = bool(current_mask[cell_mask].any())
    expected_sat_flag = bool(registered_mask[cell_mask].any())

    print("ICHNOS saturation / registration diagnostic")
    print("--------------------------------------------")
    print(f"Applied red shift:             {shift_rc}")
    print(f"Original saturated coordinate: {np.argwhere(red_saturation).tolist()}")
    print(f"Registered saturated coord.:   {np.argwhere(registered_mask).tolist()}")
    print()
    print(f"Current unregistered mask says saturated: {current_sat_flag}")
    print(f"Registered mask says saturated:           {expected_sat_flag}")
    print()

    if current_sat_flag == expected_sat_flag:
        print("RESULT: diagnostic did NOT reproduce the alignment bug.")
    else:
        print("RESULT: BUG REPRODUCED.")
        print(
            "The registered red image and the unregistered red saturation "
            "mask disagree in coordinate space."
        )

    # Keep this variable alive deliberately: this confirms that the same
    # existing ICHNOS image transform executes successfully.
    assert registered_red.shape == shape


if __name__ == "__main__":
    main()