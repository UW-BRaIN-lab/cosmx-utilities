#!/usr/bin/env python3
"""The per-cell batch label stage 3a hands the patient-batch PCA, and the identifier hygiene it needs.

The PCA corrects at the patient level (`Case`). The GBM surgical specimens are the same patients
as some of the postmortem tissue, but were resected, fixed and stored differently, so they carry a
batch effect of their own and must not share a batch with that patient's postmortem tissue. Cells
whose `Region` is one of the split regions therefore get `<Case>_<Region>` (e.g. `7488_Surgical`);
everyone else keeps their plain `Case`.

`Case` and `Block` are identifiers, but the AtoMx flat files make them numbers on most slides and
text on any slide carrying a non-numeric value (the nerve FOVs are `Nerve`). Concatenated, that is a
column of mixed ints and strings, which anndata refuses to write, so they are coerced to text per
slide before the concat.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

IDENTIFIER_COLUMNS = ("Case", "Block")
REGION_COL = "Region"
PCA_BATCH_COL = "pca_batch"
SPLIT_SEPARATOR = "_"


def _identifier_text(value):
    if pd.isna(value):
        return value
    # A slide with a missing Case reads the column as float, so 7706 arrives as 7706.0.
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return str(int(value))
    return str(value).strip()


def identifiers_as_text(obs: pd.DataFrame) -> pd.DataFrame:
    """obs with Case/Block (where present) as text, missing values left missing."""
    obs = obs.copy()
    for col in IDENTIFIER_COLUMNS:
        if col in obs.columns:
            obs[col] = obs[col].astype(object).map(_identifier_text)
    return obs


def pca_batch_labels(obs: pd.DataFrame, batch_col: str, split_regions: list[str]) -> pd.Series:
    """Per-cell PCA batch: obs[batch_col] as text, split off as `<batch>_<Region>` in split_regions.

    Region values are matched exactly. A split region absent from the data is reported, not an
    error: the default (`Surgical`) applies to every study and most have no surgical tissue.
    """
    batch = obs[batch_col].astype(str)
    if not split_regions:
        return batch
    if REGION_COL not in obs.columns:
        print(f"NOTE: no '{REGION_COL}' column; PCA batch is plain '{batch_col}' "
              f"(split regions {split_regions} not applied)")
        return batch

    region = obs[REGION_COL].astype(str)
    split = region.isin(split_regions)
    absent = [r for r in split_regions if r not in set(region)]
    if absent:
        print(f"NOTE: split region(s) {absent} not present in obs['{REGION_COL}']; "
              f"those batches are not split")
    if not split.any():
        return batch

    suffix = region.str.replace(" ", SPLIT_SEPARATOR, regex=False)
    batch = batch.where(~split, batch + SPLIT_SEPARATOR + suffix)
    split_counts = batch[split].value_counts().sort_index()
    print(f"PCA batch: {int(split.sum()):,} cells in {split_regions} split into "
          f"{len(split_counts)} batch(es) of their own:")
    for label, count in split_counts.items():
        print(f"  {label:28s} {count:9,d}")
    return batch
