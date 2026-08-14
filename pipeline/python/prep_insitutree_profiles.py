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
  --annotations   that run's de-novo annotation CSV (reference/denovo_annotations/*.csv).
                  Two schemas are supported:
                    * DATA-DRIVEN (preferred): has `keep` + `new_name` columns — the kept
                      de-novo clusters (keep truthy) are renamed to their new_name. This is
                      the single source of truth; edit the CSV, not this file.
                    * LEGACY: only denovo_label/annotation — falls back to the hardcoded
                      DEFAULT_DENOVO_RENAME (the original Core-L4 pilot map).

Writes:
  --output        genes x cell-types CSV (gene = row index), linear scale, ready as the
                  InSituTree `full_profiles`. Columns = all named types + renamed de-novo.

The kept new_names are the LEAVES of the Malignant branch (+ Low_signal sink) in
insitutree_hierarchy.json, so whatever this run keeps MUST be reflected there in lockstep.

Usage:
    uv run python pipeline/python/prep_insitutree_profiles.py \\
        --profiles-h5 anchor_typing.h5 \\
        --annotations pipeline/reference/denovo_annotations/fullcohort_pruned_k27.csv \\
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

# FALLBACK rename map for LEGACY annotation CSVs that lack keep/new_name columns (the original
# Core-L4 pilot run). Preferred path is data-driven: a `keep`/`new_name` annotations CSV. The
# `_denovo` suffix guarantees no collision with any Core-L4 named type. Keep in lockstep with
# pipeline/reference/insitutree_hierarchy.json.
DEFAULT_DENOVO_RENAME = {
    "b": "Stress_denovo",       # Heat-shock/stress (cross-cutting)
    "d": "MES_AClike_denovo",   # MES/AC-like tumor
    "f": "Hypoxia_denovo",      # Hypoxia/angiogenic
    "h": "OPClike_denovo",      # OPC-like tumor
    "k": "MESlike_denovo",      # MES-like (mixed/weak)
    "j": "Low_signal_denovo",   # Low-signal/generic sink (housekeeping only)
}

# De-novo clusters = InSituType cluster_name_pool = 1-2 lowercase letters (a..z, aa..; K>26
# overflows into two letters). Named types always carry an uppercase/digit/separator.
DENOVO_LABEL_RE = re.compile(r"^[a-z]{1,2}$")

TRUTHY = {"true", "t", "1", "yes", "y"}


def load_rename_map(annotations_path: Path, denovo_labels: list[str]) -> dict[str, str]:
    """Return {denovo_label -> new_name} for the clusters to KEEP.

    Data-driven when the CSV has keep + new_name columns (kept = keep truthy, renamed to
    new_name); otherwise falls back to DEFAULT_DENOVO_RENAME. Validates that every kept label
    exists in this run, has a non-empty unique new_name, and (data-driven) that at least one
    row is kept.
    """
    ann = pd.read_csv(annotations_path)
    if "denovo_label" not in ann.columns:
        sys.exit(f"ERROR: {annotations_path} has no 'denovo_label' column.")
    ann["denovo_label"] = ann["denovo_label"].astype(str).str.strip()

    if {"keep", "new_name"}.issubset(ann.columns):
        kept = ann[ann["keep"].astype(str).str.strip().str.lower().isin(TRUTHY)].copy()
        kept["new_name"] = kept["new_name"].astype(str).str.strip()
        blank = kept[kept["new_name"].isin(("", "nan"))]["denovo_label"].tolist()
        if blank:
            sys.exit(f"ERROR: kept de-novo {blank} have an empty new_name in {annotations_path}.")
        rename = dict(zip(kept["denovo_label"], kept["new_name"]))
        if not rename:
            sys.exit(f"ERROR: no rows marked keep in {annotations_path}.")
        print(f"Rename map: DATA-DRIVEN from {annotations_path.name} "
              f"({len(rename)} keep, {len(ann) - len(rename)} drop)")
    else:
        rename = dict(DEFAULT_DENOVO_RENAME)
        print(f"Rename map: FALLBACK DEFAULT_DENOVO_RENAME (no keep/new_name in "
              f"{annotations_path.name}); {len(rename)} keep")

    missing = sorted(set(rename) - set(denovo_labels))
    if missing:
        sys.exit(f"ERROR: kept de-novo {missing} absent from this run's de-novo {denovo_labels}. "
                 f"Wrong run, or fix the annotations CSV.")
    dups = [n for n in set(rename.values()) if list(rename.values()).count(n) > 1]
    if dups:
        sys.exit(f"ERROR: duplicate new_name(s) in the rename map: {sorted(dups)}.")
    return rename


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

    # Which de-novo to keep and what to rename them (data-driven from the annotations CSV).
    rename = load_rename_map(args.annotations, denovo)

    # Print each keep with its identity (annotation column, if present) for the run log.
    ann = pd.read_csv(args.annotations)
    ann["denovo_label"] = ann["denovo_label"].astype(str).str.strip()
    ann_id = ann.set_index("denovo_label")["annotation"] if "annotation" in ann.columns else None
    for letter, new_name in rename.items():
        identity = ann_id.get(letter, "(no annotation)") if ann_id is not None else "(no annotation)"
        print(f"  keep de-novo '{letter}' -> '{new_name}'   [{identity}]")
    dropped = sorted(set(denovo) - set(rename))
    print(f"  dropping de-novo {dropped} (conserved types better covered by GBmap, or junk)")

    kept = named + list(rename)
    out = profiles[kept].rename(columns=rename)

    # Guards: renamed de-novo must not collide with any named type, and no dup columns.
    collisions = set(rename.values()) & set(named)
    assert not collisions, f"renamed de-novo collide with named types: {collisions}"
    assert not out.columns.duplicated().any(), "duplicate columns in output profiles"

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output)
    print(f"Wrote {args.output} ({out.shape[0]:,} genes x {out.shape[1]} cell types: "
          f"{len(named)} named + {len(rename)} de-novo)")
    print("REMINDER: add the kept new_names as leaves in insitutree_hierarchy.json (lockstep).")


if __name__ == "__main__":
    main()
