#!/usr/bin/env python3
"""Measure whether a morphology channel marks nuclei or washes across the section.

The per-cell flat-file metadata cannot answer this. `Mean.<channel>` is computed
only INSIDE segmented cells, so a wash of signal sitting everywhere -- including
the extracellular space -- is never sampled, and if it lifts every cell equally it
cancels out of any within-cell ratio. That is precisely the "haze of watercolor"
failure: signal above baseline across the whole slide that does not mark nuclei
specifically.

Answering it needs the images. `CompartmentLabels_F*.tif` labels every pixel as
background / nuclear / cytoplasm / membrane, so the morphology TIFF can be scored
against it directly:

  localization_index      nuclear mean / background mean. THE headline. A stain
                          that marks nuclei is high; a wash approaches 1.
  nuclear_over_cytoplasm  does the signal sit in NUCLEI, or merely inside cells?
                          A nuclear stain is >1; DAPI bound to cytoplasmic RNA
                          drags this toward (or under) 1.
  background_share        background mean / nuclear mean, i.e. how much of the
                          nuclear reading is just the wash underneath it.
  saturated_frac          fraction of pixels at the dtype ceiling. A clipped stain
                          also loses contrast, and looks hazy, for an unrelated
                          reason -- worth excluding before blaming chemistry.

Compartment codes are read from the label image itself: every distinct value gets
a row, with the conventional names applied when they match. An unexpected encoding
still yields data instead of a wrong answer.

Usage:
    uv run python pipeline/python/morphology_localization.py \
        --fov-dir staged_fovs --slide-id 7583G27583G7 \
        --output localization_by_fov.csv
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

# Page order within a CosMx morphology TIFF, matching napari-cosmx's convention.
CHANNEL_ORDER = ["B", "G", "Y", "R", "U"]
DEFAULT_MARKERS = ["Histone", "Empty", "rRNA", "GFAP", "DAPI"]

# Conventional CompartmentLabels encoding. Applied only to values that appear.
COMPARTMENT_NAMES = {0: "background", 1: "nuclear", 2: "cytoplasm", 3: "membrane"}
BACKGROUND = "background"
NUCLEAR = "nuclear"
CYTOPLASM = "cytoplasm"

MORPHOLOGY_PATTERN = re.compile(r".*C902_P99_N99_F(?P<fov>\d+)\.TIF$", re.IGNORECASE)
COMPARTMENT_PATTERN = re.compile(r"CompartmentLabels_F(?P<fov>\d+)\.tif$", re.IGNORECASE)

# Below this many pixels a compartment is segmentation noise, not a measurement.
MIN_COMPARTMENT_PIXELS = 1000


def _pillow_page(path: Path, index: int) -> np.ndarray:
    from PIL import Image

    with Image.open(path) as image:
        image.seek(index)
        return np.asarray(image)


def read_page(path: Path, index: int = 0,
              tif: tifffile.TiffFile | None = None) -> np.ndarray:
    """Decode one TIFF page, falling back to Pillow when tifffile lacks a codec.

    tifffile delegates LZW and friends to imagecodecs, which the pipeline
    container does not ship, and CompartmentLabels images are LZW-compressed.
    Pillow decodes them natively. Reading tags never needs a codec, so the
    channel metadata still comes from tifffile either way.
    """
    try:
        if tif is not None:
            return tif.pages[index].asarray()
        return tifffile.imread(path, key=index)
    except ValueError as exc:
        message = str(exc)
        if "imagecodecs" not in message and "COMPRESSION" not in message:
            raise
        return _pillow_page(path, index)


class MissingCompartmentError(Exception):
    """A FOV's morphology image has no matching CompartmentLabels image."""


class ChannelCountError(Exception):
    """The morphology TIFF's page count does not match the channel convention."""


def channel_markers(tif: tifffile.TiffFile, overrides: dict[str, str]) -> list[str]:
    """Marker name per page, from MorphologyKit metadata where it is present.

    The kit metadata is known to mislabel targets and to disagree between studies
    for the same physical panel, so --channel-name always wins.
    """
    markers = list(DEFAULT_MARKERS)
    try:
        description = tif.pages[0].tags["ImageDescription"].value
        reagents = json.loads(description)["MorphologyKit"]["MorphologyReagents"]
        by_channel = {r["Fluorophore"]["ChannelId"]: r["BiologicalTarget"].replace("/", "_")
                      for r in reagents}
        markers = [by_channel.get(c, d) for c, d in zip(CHANNEL_ORDER, DEFAULT_MARKERS)]
    except (KeyError, ValueError, TypeError):
        pass  # fall back to the conventional order
    return [overrides.get(c, m) for c, m in zip(CHANNEL_ORDER, markers)]


def compartment_masks(labels: np.ndarray) -> dict[str, np.ndarray]:
    """Boolean mask per compartment actually present in the label image."""
    masks: dict[str, np.ndarray] = {}
    for value in np.unique(labels):
        mask = labels == value
        if int(mask.sum()) < MIN_COMPARTMENT_PIXELS:
            continue
        masks[COMPARTMENT_NAMES.get(int(value), f"code_{int(value)}")] = mask
    return masks


def channel_rows(plane: np.ndarray, marker: str, channel: str,
                 masks: dict[str, np.ndarray]) -> dict:
    """Localization metrics for one channel of one FOV."""
    ceiling = np.iinfo(plane.dtype).max if np.issubdtype(plane.dtype, np.integer) else None
    row: dict = {
        "channel": channel,
        "marker": marker,
        "saturated_frac": (float((plane == ceiling).mean()) if ceiling is not None
                           else float("nan")),
        "whole_image_mean": float(plane.mean()),
    }
    for name, mask in masks.items():
        row[f"{name}_mean"] = float(plane[mask].mean())
        row[f"{name}_frac_px"] = float(mask.mean())

    background = row.get(f"{BACKGROUND}_mean")
    nuclear = row.get(f"{NUCLEAR}_mean")
    cytoplasm = row.get(f"{CYTOPLASM}_mean")

    # A zero background would make the ratio infinite rather than merely large;
    # NaN keeps an unmeasurable FOV out of the summaries instead of dominating them.
    row["localization_index"] = (nuclear / background
                                 if background and nuclear is not None and background > 0
                                 else float("nan"))
    row["nuclear_over_cytoplasm"] = (nuclear / cytoplasm
                                     if cytoplasm and nuclear is not None and cytoplasm > 0
                                     else float("nan"))
    row["background_share"] = (background / nuclear
                               if nuclear and background is not None and nuclear > 0
                               else float("nan"))
    return row


def fov_pairs(fov_dir: Path) -> list[tuple[str, Path, Path]]:
    """(fov, morphology tif, compartment tif) for every FOV with both images."""
    compartments = {}
    for path in fov_dir.rglob("CompartmentLabels_F*.tif"):
        match = COMPARTMENT_PATTERN.search(path.name)
        if match:
            compartments[match.group("fov").lstrip("0") or "0"] = path

    pairs = []
    for path in sorted(fov_dir.rglob("*.TIF")) + sorted(fov_dir.rglob("*.tif")):
        match = MORPHOLOGY_PATTERN.match(path.name)
        if not match:
            continue
        fov = match.group("fov").lstrip("0") or "0"
        if fov not in compartments:
            print(f"WARN: no CompartmentLabels for FOV {fov}, skipping {path.name}",
                  file=sys.stderr)
            continue
        pairs.append((fov, path, compartments[fov]))
    return pairs


def analyse_fov(slide_id: str, fov: str, morphology: Path, compartment: Path,
                overrides: dict[str, str]) -> list[dict]:
    labels = read_page(compartment)
    masks = compartment_masks(labels)
    if NUCLEAR not in masks or BACKGROUND not in masks:
        print(f"WARN: FOV {fov} has compartments {sorted(masks)}; "
              f"needs both '{NUCLEAR}' and '{BACKGROUND}' for the ratios",
              file=sys.stderr)

    rows = []
    with tifffile.TiffFile(morphology) as tif:
        n_pages = len(tif.pages)
        if n_pages != len(CHANNEL_ORDER):
            raise ChannelCountError(
                f"{morphology.name} has {n_pages} pages, expected {len(CHANNEL_ORDER)} "
                f"({'/'.join(CHANNEL_ORDER)})")
        markers = channel_markers(tif, overrides)
        for index, (channel, marker) in enumerate(zip(CHANNEL_ORDER, markers)):
            plane = read_page(morphology, index, tif)
            if plane.shape != labels.shape:
                raise ChannelCountError(
                    f"FOV {fov}: morphology plane {plane.shape} does not match "
                    f"CompartmentLabels {labels.shape}; binning mismatch")
            row = {"slide_id": slide_id, "fov": int(fov)}
            row.update(channel_rows(plane, marker, channel, masks))
            rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fov-dir", type=Path, required=True,
                        help="Directory holding staged morphology and CompartmentLabels TIFFs")
    parser.add_argument("--slide-id", required=True)
    parser.add_argument("--channel-name", action="append", default=[],
                        help="Override the kit's channel-to-marker map, CH=MARKER "
                             "with CH in B/G/Y/R/U. Repeatable.")
    parser.add_argument("--max-fovs", type=int,
                        help="Analyse at most this many FOVs (they are ~150MB each)")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    overrides = {}
    for pair in args.channel_name:
        channel, separator, name = pair.partition("=")
        if not separator or channel not in CHANNEL_ORDER:
            parser.error(f"invalid --channel-name {pair!r}: expected CH=MARKER "
                         f"with CH in {CHANNEL_ORDER}")
        overrides[channel] = name

    pairs = fov_pairs(args.fov_dir)
    if not pairs:
        print(f"ERROR: no morphology/CompartmentLabels pairs under {args.fov_dir}",
              file=sys.stderr)
        sys.exit(1)
    if args.max_fovs:
        pairs = pairs[:args.max_fovs]
    print(f"Analysing {len(pairs)} FOV(s) for {args.slide_id}")

    rows = []
    for fov, morphology, compartment in pairs:
        print(f"  FOV {fov}: {morphology.name}")
        rows.extend(analyse_fov(args.slide_id, fov, morphology, compartment, overrides))

    frame = pd.DataFrame(rows)
    frame.to_csv(args.output, index=False)
    print(f"Wrote {args.output}")

    columns = [c for c in ("localization_index", "nuclear_over_cytoplasm",
                           "background_share", "saturated_frac")
               if c in frame.columns]
    print(f"\nMedians across {frame['fov'].nunique()} FOV(s):")
    print(frame.groupby("marker")[columns].median().round(3).to_string())


if __name__ == "__main__":
    main()
