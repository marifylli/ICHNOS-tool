"""Publish baseline and manual-review sensitivity summaries in a new directory."""
import argparse
import json
from pathlib import Path
import pandas as pd
from ichnos.artifacts import new_output_directory, code_provenance, sha256_file
from ichnos.manual_review import review_summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cells', type=Path, required=True)
    parser.add_argument('--review', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--data-kind', choices=['synthetic', 'experimental'], required=True)
    parser.add_argument('--min-cells', type=int, default=3)
    parser.add_argument('--green-floor', type=float, default=0.0)
    args = parser.parse_args()
    try:
        hashes = {'cells': sha256_file(args.cells), 'review': sha256_file(args.review)}
        dtypes = {'session_id': str, 'sample_id': str}
        cells = pd.read_csv(args.cells, dtype=dtypes)
        review = pd.read_csv(args.review, dtype=dtypes)
        outputs = review_summaries(cells, review, data_kind=args.data_kind,
                                  min_cells=args.min_cells, green_floor=args.green_floor)
        with new_output_directory(args.out_dir) as staged:
            for name, frame in zip(['baseline', 'without-flagged', 'comparison', 'annotated-cells'], outputs):
                frame.to_csv(staged / (name + '.csv'), index=False)
            if hashes != {'cells': sha256_file(args.cells), 'review': sha256_file(args.review)}:
                raise ValueError('input files changed during processing')
            metadata = dict(schema_version=1, code=code_provenance(), input_sha256=hashes,
                data_kind=args.data_kind, min_cells=args.min_cells, green_floor=args.green_floor,
                flagged_objects=len(review), baseline_qc_unchanged=True,
                interpretation='Manual flags are uncertainty annotations, not confirmed segmentation errors.',
                experimentally_validated=False)
            (staged / 'metadata.json').write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n')
        print(f'{len(review)} flagged objects -> {args.out_dir}')
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f'Error: {error}\n')


if __name__ == '__main__':
    main()
