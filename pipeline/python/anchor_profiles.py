#!/usr/bin/env python3
"""Shared reader for the `$profiles` matrix inside an InSituType result h5.

`R/insitutype_typing.R` writes every run's converged reference profiles alongside the
per-cell calls (/profiles, /profile_genes, /profile_types). Several tools downstream want
that matrix — prep_insitutree_profiles.py (build the InSituTree reference),
prep_supervised_profiles.py (named-only profiles for a fully supervised re-score) — so the
read + orientation + named/de-novo split live here rather than being copied into each.

Everything in this module is pure: no I/O beyond the one h5 read, no argparse.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

# De-novo clusters = InSituType cluster_name_pool = 1-2 lowercase letters (a..z, aa..; K>26
# overflows into two letters). Named types always carry an uppercase/digit/separator.
DENOVO_LABEL_RE = re.compile(r"^[a-z]{1,2}$")


def decode(arr: np.ndarray) -> list[str]:
    """hdf5r variable-length UTF-8 comes back as bytes-or-str; normalize to str."""
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in arr]


def is_denovo(label: str) -> bool:
    return bool(DENOVO_LABEL_RE.match(str(label)))


def split_named_denovo(types: list[str]) -> tuple[list[str], list[str]]:
    """Partition a run's cell-type labels into (named reference types, de-novo letters)."""
    named = [t for t in types if not is_denovo(t)]
    denovo = [t for t in types if is_denovo(t)]
    return named, denovo


def read_profiles(h5_path: Path) -> pd.DataFrame:
    """Read /profiles into a genes x types DataFrame, orienting robustly.

    R/insitutype_typing.R writes profiles as an R matrix (genes x types); crossing the
    R->HDF5->numpy boundary can transpose it, so we key orientation off the known
    gene/type vector lengths rather than trusting the stored axis order.
    """
    with h5py.File(h5_path, "r") as f:
        for key in ("profiles", "profile_genes", "profile_types"):
            if key not in f:
                print(f"ERROR: /{key} missing in {h5_path}. This run's h5 predates the "
                      f"profiles writeback — regenerate it, or dump $profiles from the "
                      f"companion insitutype_result.rds instead.", file=sys.stderr)
                sys.exit(1)
        mat = f["profiles"][()]
        genes = decode(f["profile_genes"][()])
        types = decode(f["profile_types"][()])

    n_genes, n_types = len(genes), len(types)
    if mat.shape == (n_genes, n_types):
        pass
    elif mat.shape == (n_types, n_genes):
        mat = mat.T
    else:
        print(f"ERROR: /profiles shape {mat.shape} matches neither "
              f"({n_genes} genes, {n_types} types) nor its transpose.", file=sys.stderr)
        sys.exit(1)

    return pd.DataFrame(mat, index=pd.Index(genes, name="gene"), columns=types)


def read_cell_calls(h5_path: Path) -> pd.DataFrame:
    """Read the per-cell calls (/cell_id, /cell_type, /prob) as a DataFrame."""
    with h5py.File(h5_path, "r") as f:
        for key in ("cell_id", "cell_type"):
            if key not in f:
                sys.exit(f"ERROR: /{key} missing in {h5_path}.")
        out = pd.DataFrame({
            "cell_id": decode(f["cell_id"][()]),
            "cell_type": decode(f["cell_type"][()]),
        })
        if "prob" in f:
            out["prob"] = np.asarray(f["prob"][()], dtype=float)
    return out
