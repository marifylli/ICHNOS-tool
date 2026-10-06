"""Fit biological-replicate covariance or compute a conditional snapshot posterior."""
import argparse
import hashlib
import json
from pathlib import Path
from ichnos.posterior import fit_replicate_noise, decode_posterior
from ichnos.snapshot_decoder import load_observable_table, save_artifact


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    for name in ('fit-noise','decode'):
        command=commands.add_parser(name)
        command.add_argument('--table',type=Path,required=True)
        command.add_argument('--calibration',type=Path,required=True)
        command.add_argument('--out',type=Path,required=True)
        if name=='fit-noise':
            command.add_argument('--replicates',type=Path,required=True)
            command.add_argument('--green-floor',type=float,required=True)
            command.add_argument('--min-replicates-per-condition',type=int,default=3)
        else:
            command.add_argument('--noise',type=Path,required=True)
            command.add_argument('--observation',type=Path,required=True)
            command.add_argument('--prior',type=Path,required=True)
            command.add_argument('--max-mahalanobis-squared',type=float,required=True)
            command.add_argument('--credible-mass',type=float,default=.95)
    args=parser.parse_args()
    try:
        if args.out.exists():
            raise FileExistsError(f'output path already exists: {args.out}')
        table=load_observable_table(args.table)
        calibration=json.loads(args.calibration.read_text())
        if args.command=='fit-noise':
            result=fit_replicate_noise(table,calibration,json.loads(args.replicates.read_text()),
                green_floor=args.green_floor,min_replicates_per_condition=args.min_replicates_per_condition)
        else:
            result=decode_posterior(table,calibration,json.loads(args.noise.read_text()),
                json.loads(args.observation.read_text()),prior=json.loads(args.prior.read_text()),
                max_mahalanobis_squared=args.max_mahalanobis_squared,credible_mass=args.credible_mass)
            result['observation_file_sha256']=hashlib.sha256(args.observation.read_bytes()).hexdigest()
        save_artifact(result,args.out)
        print(f"{result.get('status',args.command)} -> {args.out}")
    except (ValueError,KeyError,TypeError,OSError) as exc:
        parser.exit(1,f'Error: {exc}\n')


if __name__=='__main__':
    main()
