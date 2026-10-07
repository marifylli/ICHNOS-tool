"""Evaluate ox/er on/off-grid recovery under noise and calibration perturbations."""
import argparse
from pathlib import Path
from ichnos.evaluation import run_evaluation, SCENARIOS


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir',type=Path,required=True)
    parser.add_argument('--variants',nargs='+',choices=['ox','er'],default=['ox','er'])
    parser.add_argument('--dose-nodes',type=int,default=9)
    parser.add_argument('--time-nodes',type=int,default=17)
    parser.add_argument('--trials',type=int,default=12)
    parser.add_argument('--noise-scales',type=float,nargs='+',default=[.5,1.,2.])
    parser.add_argument('--scenarios',nargs='+',choices=list(SCENARIOS),default=list(SCENARIOS))
    parser.add_argument('--seed',type=int,default=20261007)
    parser.add_argument('--refine',action='store_true')
    args=parser.parse_args()
    try:
        report=run_evaluation(**vars(args))
    except (ValueError,KeyError,OSError) as exc:
        parser.exit(1,f'Error: {exc}\n')
    print(f"{report['n_trials']} synthetic trials -> {args.out_dir}")


if __name__=='__main__':
    main()
