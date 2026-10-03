"""Reading acquisition metadata straight from a real microscope TIFF file
(exposure time, objective, pixel size, instrument state), instead of typing
it into a manifest CSV by hand.

Deliberately minimal and NOT validated against a real Olympus TIFF (none
available yet) -- extracts whatever tifffile can find (OME-XML if present,
the ImageDescription tag, raw TIFF tags) and hands it back as-is, rather
than guessing at specific field names/paths that may not match what the
real acquisition software actually writes. Once a real Olympus TIFF exists,
inspect read_tiff_metadata()'s output on it first, then extend this with a
proper field-extraction function (e.g. exposure_ms_from_metadata()) built
against what's actually there.
"""
from __future__ import annotations

from pathlib import Path


def read_tiff_metadata(path: str | Path) -> dict:
    """Best-effort dump of a TIFF's embedded metadata via `tifffile`.

    Returns a dict with:
    - shape, dtype: of the first page/frame
    - ome_xml: the embedded OME-XML metadata string, if the file has any
      (most scientific TIFF writers -- Bio-Formats, MetaMorph, etc. -- embed
      this; a raw Olympus .oib/.vsi export may or may not)
    - image_description: the TIFF ImageDescription tag's raw value, if set
      (some acquisition software puts instrument settings here as free text
      instead of/alongside OME-XML)
    - tags: every other TIFF tag on the first page, name -> value, for
      manual inspection when neither of the above has what you need
    """
    try:
        import tifffile
    except ImportError as exc:
        raise ImportError("tifffile is required for read_tiff_metadata(); pip install tifffile") from exc

    with tifffile.TiffFile(str(path)) as tif:
        page = tif.pages[0]
        result: dict = {"shape": page.shape, "dtype": str(page.dtype)}

        if tif.ome_metadata:
            result["ome_xml"] = tif.ome_metadata

        description_tag = page.tags.get("ImageDescription")
        if description_tag is not None:
            result["image_description"] = description_tag.value

        result["tags"] = {
            tag.name: tag.value
            for tag in page.tags
            if tag.name not in ("ImageDescription",)
        }

        return result
