#!/usr/bin/env python3
"""Drop cells in FOVs the SME ruled out of analysis, by (slide, FOV).

Shared by the stages that turn cells into something downstream trusts -- the InSituType /
InSituTree input (prep_insitutype_inputs.py) and the de novo component tables
(diagnose_denovo_components.py) -- because an exclusion applied in one and not the other
produces results that disagree for reasons nobody can see.

Why not obs['Region_excluded']: flatfiles_to_anndata.py writes that column, but only when
stage 1 was run with --fov-annotations. A cohort assembled before the SME's exclusions
existed carries the stale upstream Region and no exclusion flag at all, so the annotation
CSV has to be read directly if the exclusion is to apply to data already on disk.

See pipeline/reference/FOV_ANNOTATIONS.md for the annotation file itself.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

SLIDE_COL = "slide_id"
FOV_COL = "fov"
# Columns of fov_annotations.csv.
FOV_SLIDE_COL, FOV_FOV_COL, FOV_EXCLUDE_COL = "slide_id", "fov", "exclude"


def excluded_fov_keys(annotations_csv: Path) -> set[tuple[str, int]]:
    """The (slide, FOV) pairs marked exclude=1, validated."""
    annotations = pd.read_csv(annotations_csv, dtype={FOV_SLIDE_COL: str})
    for col in (FOV_SLIDE_COL, FOV_FOV_COL, FOV_EXCLUDE_COL):
        if col not in annotations.columns:
            sys.exit(f"ERROR: {annotations_csv} has no '{col}' column")
    excluded = annotations.loc[annotations[FOV_EXCLUDE_COL].astype(int) == 1]
    return set(zip(excluded[FOV_SLIDE_COL].astype(str).str.strip(),
                   excluded[FOV_FOV_COL].astype(int)))


def excluded_mask(obs: pd.DataFrame, annotations_csv: Path) -> pd.Series:
    """Boolean over obs: True where the cell sits in an excluded FOV.

    Exits if the keys match no cell at all. Matching nothing means the slide ids or FOV
    numbering disagree, and the caller would then carry on with contaminated data believing
    it had been filtered -- the precise failure this exists to prevent.
    """
    missing = [c for c in (SLIDE_COL, FOV_COL) if c not in obs.columns]
    if missing:
        sys.exit(f"ERROR: obs is missing {missing}, needed to apply FOV exclusions")

    keys = excluded_fov_keys(annotations_csv)
    if not keys:
        print(f"NOTE: {annotations_csv} marks no FOV excluded; nothing dropped")
        return pd.Series(False, index=obs.index)

    mask = pd.Series(
        list(zip(obs[SLIDE_COL].astype(str).str.strip(), obs[FOV_COL].astype(int))),
        index=obs.index).isin(keys)
    if not mask.any():
        sys.exit(
            f"ERROR: none of the {len(keys)} excluded (slide, FOV) pairs matched any cell.\n"
            f"  obs slides:        {sorted(obs[SLIDE_COL].astype(str).unique())[:6]}\n"
            f"  annotation slides: {sorted({slide for slide, _ in keys})[:6]}\n"
            f"  If this input was ALREADY filtered upstream -- e.g. an InSituTree result "
            f"built from an insitutype_input.h5 that stage 70 filtered -- then there is "
            f"nothing left to exclude and you should omit the FOV annotations. Otherwise "
            f"the slide ids or FOV numbering disagree.")
    return mask


def drop_excluded_fovs(df: pd.DataFrame, annotations_csv: Path) -> pd.DataFrame:
    """df without the cells in excluded FOVs, reporting how many went."""
    mask = excluded_mask(df, annotations_csv)
    if not mask.any():
        return df
    print(f"Excluding {int(mask.sum()):,} cells in "
          f"{len(excluded_fov_keys(annotations_csv))} SME-excluded FOVs "
          f"({mask.mean() * 100:.1f}% of the cohort)")
    return df[~mask]
