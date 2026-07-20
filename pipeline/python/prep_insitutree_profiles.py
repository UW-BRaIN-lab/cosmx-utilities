#!/usr/bin/env python3
"""Build the InSituTree reference profile matrix: Core-L4 GBmap + our de-novo tumor states.

InSituTree (stage 4, hierarchical typing) is SUPERVISED — unlike InSituType it does NOT
platform-rescale the profiles it is given (its runInSituTree -> collapseProfiles ->
supervisedSubcluster -> insitutypeML path does no scRNA->CosMx correction). So its
`full_profiles` matrix must already be in CosMx expression scale, AND every column must
share one consistent scaling or the per-branch likelihoods are incomparable.

We get both properties for free by sourcing the whole matrix from ONE InSituType run's
post-rescale `$profiles`: our original Core-L4 rescale run (Kopah `stage4`), whose
profiles hold the 54 Core GBmap level-4 named types + the de-novo clusters (a..l), all
produced by a single insitutype(rescale=TRUE) call. We keep:
  - every NAMED Core-L4 type (the well-anchored non-malignant + Neftel compartments), and
  - the de-novo TUMOR-state clusters we validated from that run's markers, renamed to
    stable, collision-free names (see DEFAULT_DENOVO_RENAME; identities in
    reference/denovo_annotations/stage4.csv).
De-novo clusters that merely duplicate a conserved GBmap type (a=Myeloid, c=Vascular,
e=Astrocyte, g/i=Neuronal, l=Oligo) are dropped — GBmap's versions are better anchored.

The kept de-novo names are the LEAVES of the Malignant branch (+ a low-signal sink) in
insitutree_hierarchy.json, so the rename map here and that JSON MUST stay in lockstep.

Reads (stage the h5 down from Kopah first, e.g.:
    s5cmd cp s3://$KOPAH_BUCKET/$KOPAH_PREFIX/stage4/insitutype_result.h5 .):
  --profiles-h5   an InSituType run's result h5 (written by R/insitutype_typing.R) with
                  datasets /profiles (genes x types), /profile_genes, /profile_types.
  --annotations   that run's de-novo annotation CSV (reference/denovo_annotations/*.csv):
                  columns denovo_label, annotation, ... — used only to print/verify which
                  identity each kept letter carries (the actual rename is DEFAULT_DENOVO_RENAME).

Writes:
  --output        genes x cell-types CSV (gene = row index), linear scale, ready as the
                  InSituTree `full_profiles`. Columns = all named types + renamed de-novo.

Usage:
    uv run python pipeline/python/prep_insitutree_profiles.py \\
        --profiles-h5 insitutype_result.h5 \\
        --annotations pipeline/reference/denovo_annotations/stage4.csv \\
        --output pipeline/reference/insitutree_profiles.csv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

# De-novo letters (from the Core-L4 rescale run, reference/denovo_annotations/stage4.csv)
# to KEEP, mapped to stable leaf names for the InSituTree hierarchy. The `_denovo` suffix
# guarantees no collision with any of the 54 Core-L4 named types. Keep in lockstep with
# pipeline/reference/insitutree_hierarchy.json.
DEFAULT_DENOVO_RENAME = {
    "b": "Stress_denovo",       # Heat-shock/stress (cross-cutting)
    "d": "MES_AClike_denovo",   # MES/AC-like tumor
    "f": "Hypoxia_denovo",      # Hypoxia/angiogenic
    "h": "OPClike_denovo",      # OPC-like tumor
    "k": "MESlike_denovo",      # MES-like (mixed/weak)
    "j": "Low_signal_denovo",   # Low-signal/generic sink (housekeeping only)
}

DENOVO_LABEL_RE = re.compile(r"^[a-z]$")  # de-novo clusters are single lowercase letters


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profiles-h5", type=Path, required=True,
                   help="InSituType run result h5 with /profiles, /profile_genes, "
                        "/profile_types (from R/insitutype_typing.R).")
    p.add_argument("--annotations", type=Path, required=True,
                   help="De-novo annotation CSV for that run (denovo_label, annotation, ...).")
    p.add_argument("--output", type=Path, required=True,
                   help="Output InSituTree profile CSV (genes x cell types).")
    return p.parse_args()


def _decode(arr: np.ndarray) -> list[str]:
    """hdf5r variable-length UTF-8 comes back as bytes-or-str; normalize to str."""
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in arr]


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
        genes = _decode(f["profile_genes"][()])
        types = _decode(f["profile_types"][()])

    n_genes, n_types = len(genes), len(types)
    if mat.shape == (n_genes, n_types):
        pass
    elif mat.shape == (n_types, n_genes):
        mat = mat.T
    else:
        print(f"ERROR: /profiles shape {mat.shape} matches neither "
              f"({n_genes} genes, {n_types} types) nor its transpose.", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(mat, index=pd.Index(genes, name="gene"), columns=types)
    return df


def main() -> None:
    args = parse_args()

    profiles = read_profiles(args.profiles_h5)
    all_types = list(profiles.columns)
    named = [t for t in all_types if not DENOVO_LABEL_RE.match(t)]
    denovo = [t for t in all_types if DENOVO_LABEL_RE.match(t)]
    print(f"Profiles: {profiles.shape[0]:,} genes x {len(all_types)} types "
          f"({len(named)} named + {len(denovo)} de-novo {denovo})")

    # Confirm every de-novo letter we intend to keep is actually in the run.
    missing = sorted(set(DEFAULT_DENOVO_RENAME) - set(denovo))
    if missing:
        print(f"ERROR: de-novo letters {missing} requested but absent from this run's "
              f"profile_types {denovo}. Wrong run, or edit DEFAULT_DENOVO_RENAME.",
              file=sys.stderr)
        sys.exit(1)

    # Cross-check identities against the annotation table (informational).
    ann = pd.read_csv(args.annotations).set_index("denovo_label")
    for letter, new_name in DEFAULT_DENOVO_RENAME.items():
        identity = ann.loc[letter, "annotation"] if letter in ann.index else "(no annotation)"
        print(f"  keep de-novo '{letter}' -> '{new_name}'   [{identity}]")
    dropped = sorted(set(denovo) - set(DEFAULT_DENOVO_RENAME))
    print(f"  dropping de-novo {dropped} (conserved types better covered by GBmap)")

    kept = named + list(DEFAULT_DENOVO_RENAME)
    out = profiles[kept].rename(columns=DEFAULT_DENOVO_RENAME)

    # Guards: renamed de-novo must not collide with any named type, and no dup columns.
    collisions = set(DEFAULT_DENOVO_RENAME.values()) & set(named)
    assert not collisions, f"renamed de-novo collide with named types: {collisions}"
    assert not out.columns.duplicated().any(), "duplicate columns in output profiles"

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output)
    print(f"Wrote {args.output} ({out.shape[0]:,} genes x {out.shape[1]} cell types: "
          f"{len(named)} named + {len(DEFAULT_DENOVO_RENAME)} de-novo)")


if __name__ == "__main__":
    main()
