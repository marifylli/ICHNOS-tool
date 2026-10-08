"""Process an image manifest through QC, paired summaries and snapshot/posterior inference."""
import argparse
import json
from pathlib import Path
from ichnos.workflow import run_workflow


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    try:
        report = run_workflow(args.config, args.out_dir)
    except (ValueError, KeyError, TypeError, OSError) as exc:
        parser.exit(1, f'Error: {exc}\n')
    print(json.dumps({'output':str(args.out_dir),'status_counts':report['status_counts']},indent=2))


if __name__ == '__main__':
    main()
