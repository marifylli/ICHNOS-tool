"""Build a two-observable table, bind reference calibration, or decode one snapshot."""
import argparse
import hashlib
import json
from pathlib import Path
from ichnos.snapshot_decoder import (
    build_observable_table, build_snapshot_calibration, decode_snapshot,
    load_observable_table, save_artifact,
)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    build=commands.add_parser('build')
    build.add_argument('--variant',choices=['ox','er'],required=True)
    build.add_argument('--profile-name',default='default')
    build.add_argument('--doses-uM',type=float,nargs='+',required=True)
    build.add_argument('--times-hours',type=float,nargs='+',required=True)
    build.add_argument('--initialization',choices=['equilibrium','finite-preincubation'],required=True)
    build.add_argument('--preincubation-hours',type=float)
    build.add_argument('--clearance-rate-per-hour',type=float)
    build.add_argument('--out',type=Path,required=True)
    calibration=commands.add_parser('calibrate')
    calibration.add_argument('--table',type=Path,required=True)
    calibration.add_argument('--config',type=Path,required=True)
    calibration.add_argument('--out',type=Path,required=True)
    decode=commands.add_parser('decode')
    decode.add_argument('--table',type=Path,required=True)
    decode.add_argument('--calibration',type=Path,required=True)
    decode.add_argument('--observation',type=Path,required=True)
    decode.add_argument('--ratio-tolerance',type=float,required=True)
    decode.add_argument('--green-tolerance',type=float,required=True)
    decode.add_argument('--green-floor',type=float,required=True)
    decode.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    try:
        if args.out.exists():
            raise FileExistsError(f'output path already exists: {args.out}')
        if args.command=='build':
            options=vars(args).copy()
            options.pop('command')
            options.pop('out')
            result=build_observable_table(**options).to_dict()
        elif args.command=='calibrate':
            config=json.loads(args.config.read_text())
            result=build_snapshot_calibration(load_observable_table(args.table),**config)
        else:
            result=decode_snapshot(
                load_observable_table(args.table),json.loads(args.calibration.read_text()),
                json.loads(args.observation.read_text()),ratio_tolerance=args.ratio_tolerance,
                green_tolerance=args.green_tolerance,green_floor=args.green_floor,
            )
            result['observation_file_sha256']=hashlib.sha256(args.observation.read_bytes()).hexdigest()
        save_artifact(result,args.out)
        print(f"{result.get('status',args.command)} -> {args.out}")
    except (ValueError,KeyError,TypeError,OSError) as exc:
        parser.exit(1,f'Error: {exc}\n')


if __name__=='__main__':
    main()
