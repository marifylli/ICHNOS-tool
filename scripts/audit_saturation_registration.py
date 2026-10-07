"""Audit a reprocessed CSV; original images are required to recover spatial masks."""
import argparse
import json
from pathlib import Path
import pandas as pd
from ichnos_image.saturation_audit import saturation_impact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cells', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = saturation_impact(pd.read_csv(args.cells, dtype={'session_id': str, 'sample_id': str}))
    with args.out.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)


if __name__ == '__main__':
    main()
