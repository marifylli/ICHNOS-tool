"""Inspect origin and affine reference fits without changing calibration."""
import argparse
import json
from pathlib import Path
from ichnos.calibration import load_session_calibration
from ichnos.calibration_diagnostics import compare_ratio_calibration


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calibration',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    try:
        calibration=load_session_calibration(args.calibration)
        result=compare_ratio_calibration(calibration.model_ratios,calibration.image_ratios)
        result.update(session_id=calibration.session_id,reference_id=calibration.reference_id,
                      source_kind=calibration.source_kind,holdout_note='CLI examines stored fit pairs only; no independent holdout supplied')
        encoded=json.dumps(result,indent=2,allow_nan=False)+'\n'
        with args.out.open('x',encoding='utf-8') as handle:
            handle.write(encoded)
        print(args.out)
    except (ValueError,KeyError,OSError) as exc:
        parser.exit(1,f'Error: {exc}\n')


if __name__=='__main__':
    main()
