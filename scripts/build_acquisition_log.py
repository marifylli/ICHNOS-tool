"""Build acquisition_log.csv from the wet-lab microscopy log.

    python scripts/build_acquisition_log.py --pdf log.pdf --out data/acquisition/acquisition_log.csv

The PDF is converted with ``pdftotext -layout`` (poppler) when it is
installed, otherwise with pypdf (``python3 -m pip install pypdf``).
``--text`` accepts an already-converted file instead.
Warnings (a heading not found, an addition time not stated in its section)
are printed and make the script exit with status 1 unless --allow-warnings
is given, so a changed document cannot pass unnoticed.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ichnos_image.acquisition_log import parse_log_text, write_rows


def pdf_to_text(pdf: Path) -> str:
    """Layout-preserving text: poppler's pdftotext if installed, otherwise pypdf."""
    exe = shutil.which("pdftotext")
    if exe is not None:
        done = subprocess.run([exe, "-layout", "--", str(pdf), "-"], check=True, capture_output=True)
        return done.stdout.decode("utf-8")
    try:
        from pypdf import PdfReader
    except ImportError:
        raise SystemExit("need pdftotext (brew install poppler) or pypdf (python3 -m pip install pypdf)") from None
    return "\n\f".join(page.extract_text(extraction_mode="layout") for page in PdfReader(str(pdf)).pages)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--pdf", type=Path)
    source.add_argument("--text", type=Path, help="output of pdftotext -layout")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--allow-warnings", action="store_true")
    args = ap.parse_args(argv)

    text = pdf_to_text(args.pdf) if args.pdf else args.text.read_text(encoding="utf-8")
    result = parse_log_text(text)
    for warning in result.warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    if result.warnings and not args.allow_warnings:
        print("not written: resolve the warnings or pass --allow-warnings", file=sys.stderr)
        return 1
    write_rows(result.rows, args.out)
    print(f"{len(result.rows)} rows -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
