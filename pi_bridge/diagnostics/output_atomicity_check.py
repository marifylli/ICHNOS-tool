"""Diagnostic reproduction of premature output deletion.

This script does not modify the ICHNOS pipeline.
It checks whether an existing valid output survives a processing failure.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from ichnos_image import pipeline


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        out_csv = Path(tmp) / "existing_results.csv"

        # Simulate an already-valid result that must not be lost
        # if a replacement run fails.
        original_content = "valid,existing,result\n1,2,3\n"
        out_csv.write_text(original_content, encoding="utf-8")

        print("ICHNOS output atomicity diagnostic")
        print("----------------------------------")
        print(f"Existing output before run: {out_csv.exists()}")

        # We only need one object in image_sets. process_image_set is mocked
        # to fail before it needs to inspect the object.
        from types import SimpleNamespace

        dummy_image_set = SimpleNamespace(session_id="diagnostic_session")
        try:
            with patch(
                "ichnos_image.pipeline.process_image_set",
                side_effect=RuntimeError("intentional diagnostic failure"),
            ):
                pipeline.process_experiment(
                    [dummy_image_set],
                    out_csv,
                    bleed_green_to_red=0.0,
                )
        except Exception as exc:
            print(f"Controlled processing failure: {type(exc).__name__}: {exc}")

        exists_after = out_csv.exists()

        print()
        print(f"Existing output after failed run: {exists_after}")

        if exists_after:
            preserved = out_csv.read_text(encoding="utf-8") == original_content
            print(f"Original content preserved:       {preserved}")
        else:
            preserved = False
            print("Original content preserved:       False")

        print()

        if not exists_after:
            print("RESULT: BUG REPRODUCED.")
            print(
                "The previous valid output was deleted before the new "
                "processing completed successfully."
            )
        elif not preserved:
            print("RESULT: BUG REPRODUCED.")
            print(
                "The previous output survived as a path but its original "
                "content was not preserved."
            )
        else:
            print("RESULT: diagnostic did NOT reproduce the atomicity bug.")


if __name__ == "__main__":
    main()
