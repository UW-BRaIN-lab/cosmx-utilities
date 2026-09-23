#!/usr/bin/env python3
"""Per-FOV donor / block / region annotations, joined to the cell ids the pipeline uses.

A CosMx slide carries two tissue pieces from two DIFFERENT cases, and this cohort's FOV
numbering runs continuously 1-200 across both, so a slide is not a donor and the pieces cannot
be separated from the cell ids alone. Only the AtoMx annotation reference can split them, which
makes this join a prerequisite for any per-donor claim.

Two files are needed because they key differently: the annotation reference uses canonical slide
names (`7495 G3 7302 G3`) while the flat files — and therefore our stage-1 cell ids — use AtoMx
export folder names (`7495G37302G3`). The crosswalk maps between them.

Committed copies live in pipeline/reference/; see that README for provenance and the curation
traps (orientation half-swaps on three slides, irregular FOV ranges, known export gaps).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

_REFERENCE = Path(__file__).resolve().parents[1] / "reference"
DEFAULT_ANNOTATIONS = _REFERENCE / "gbm_fov_annotations.csv"
DEFAULT_CROSSWALK = _REFERENCE / "gbm_slide_name_crosswalk.csv"

ANNOTATION_COLUMNS = {"Slide", "Case", "Block", "Region", "FOVs"}
CROSSWALK_FOLDER = "AtoMx_flatfile_folder"
CROSSWALK_CANONICAL = "Canonical_slide_name"
# Stage-1 cell ids are "<slide folder>_F<fov>_C<cell>", e.g. 7104A297104A23_F1_C2.
CELL_ID_RE = re.compile(r"^(?P<slide>.+)_F(?P<fov>\d+)_C\d+$")


def load_fov_annotations(annotations: Path = DEFAULT_ANNOTATIONS,
                         crosswalk: Path = DEFAULT_CROSSWALK,
                         extra_columns: tuple[str, ...] = ()) -> pd.DataFrame:
    """Per-FOV Case/Block/Region indexed by unit key "<folder>:F<fov>".

    `extra_columns` names further crosswalk columns to carry through — `Export_source` is the
    useful one, since it identifies which of the five flat-file exports a slide came from and so
    turns any grouping into a batch test.
    """
    # The committed copy is plain UTF-8; the upstream AtoMx file carries a BOM that would
    # otherwise corrupt the first header. utf-8-sig reads both.
    ann = pd.read_csv(annotations, encoding="utf-8-sig")
    missing = ANNOTATION_COLUMNS - set(ann.columns)
    if missing:
        sys.exit(f"ERROR: {annotations} is missing {sorted(missing)}.")

    xw = pd.read_csv(crosswalk)
    for col in (CROSSWALK_FOLDER, CROSSWALK_CANONICAL, *extra_columns):
        if col not in xw.columns:
            sys.exit(f"ERROR: {crosswalk} is missing {col}.")

    keep = [CROSSWALK_CANONICAL, CROSSWALK_FOLDER, *extra_columns]
    ann = ann.merge(xw[keep], left_on="Slide", right_on=CROSSWALK_CANONICAL, how="left")
    unmapped = ann[CROSSWALK_FOLDER].isna()
    if unmapped.any():
        names = sorted(ann.loc[unmapped, "Slide"].unique())
        print(f"WARNING: {len(names)} annotated slide(s) absent from the crosswalk, dropped: "
              f"{', '.join(names[:5])}{' ...' if len(names) > 5 else ''}")
        ann = ann[~unmapped]

    ann["unit"] = ann[CROSSWALK_FOLDER] + ":F" + ann["FOVs"].astype(int).astype(str)
    if ann["unit"].duplicated().any():
        sys.exit("ERROR: the annotation reference has more than one row per (slide, FOV). "
                 "A FOV must map to exactly one case.")
    out = ann.set_index("unit")[["Case", "Block", "Region", *extra_columns]]
    out.insert(0, "slide", ann.set_index("unit")[CROSSWALK_FOLDER])
    return out


def units_from_cell_ids(cell_ids: pd.Index) -> pd.Series:
    """Map "<folder>_F<fov>_C<cell>" cell ids to the "<folder>:F<fov>" unit key."""
    parsed = pd.Series(cell_ids, index=cell_ids).str.extract(CELL_ID_RE)
    bad = parsed["slide"].isna()
    if bad.all():
        sys.exit("ERROR: no cell id matched '<slide>_F<fov>_C<cell>'. "
                 f"First id: {cell_ids[0]!r}")
    if bad.any():
        print(f"WARNING: {int(bad.sum()):,} cell ids did not parse; dropped.")
    parsed = parsed[~bad]
    return parsed["slide"] + ":F" + parsed["fov"].astype(int).astype(str)


def annotate_cells(cell_ids: pd.Index, **kwargs) -> pd.DataFrame:
    """Per-CELL Case/Block/Region (+extras), indexed by cell id. Unannotated cells are dropped."""
    units = units_from_cell_ids(cell_ids)
    ann = load_fov_annotations(**kwargs)
    out = pd.DataFrame({"unit": units}).join(ann, on="unit")
    unannotated = int(out["Case"].isna().sum())
    if unannotated:
        print(f"NOTE: {unannotated:,}/{len(out):,} cells have no FOV annotation and are dropped "
              f"from the grouped view.")
        out = out[out["Case"].notna()]
    return out
