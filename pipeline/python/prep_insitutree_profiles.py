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
import sys
from pathlib import Path

import pandas as pd

from anchor_profiles import read_profiles, split_named_denovo

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
    # Several letters MAY share one new_name: that is a merge, and it is how a de-novo
    # split that turns out to be technical rather than biological is undone. Retina's
    # c/k/i partitioned one cortical population by depth -- their profiles differ only on
    # cytoskeletal and housekeeping genes, and every laminar marker sits below 0.16
    # expected counts with CUX2 and FEZF2 off-panel entirely.
    merged = {n: sorted(l for l, v in rename.items() if v == n)
              for n in sorted(set(rename.values()))}
    for new_name, letters in merged.items():
        if len(letters) > 1:
            print(f"  MERGE de-novo {letters} -> '{new_name}'")
    return rename


def select_columns(named: list[str], rename: dict[str, str], drop_named: list[str]) -> list[str]:
    """Columns to carry into the InSituTree matrix: the named reference types minus `drop_named`, then the kept de novo.

    A reference type is dropped when its profile is not a usable leaf (SORL1: the Brain_6k L6 neuron profile is 0.986
    cosine to L4 and absorbs shallow glia). Naming a type that is not in the run is an error, so a typo cannot
    silently keep it.
    """
    unknown = sorted(set(drop_named) - set(named))
    if unknown:
        sys.exit(f"ERROR: --drop-named {unknown} are not named types in this run: {named}")
    return [name for name in named if name not in drop_named] + list(rename)


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
    p.add_argument("--drop-named", action="append", default=[], metavar="TYPE",
                   help="A named reference type to leave out of the matrix (repeatable), e.g. 'L6 neuron'. "
                        "Its leaf must then also be absent from the hierarchy JSON.")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    profiles = read_profiles(args.profiles_h5)
    all_types = list(profiles.columns)
    named, denovo = split_named_denovo(all_types)
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
    print(f"  dropping de-novo {dropped} (not kept; see the annotations CSV for each one's disposition)")

    kept = select_columns(named, rename, args.drop_named)
    out = profiles[kept].rename(columns=rename)
    if out.columns.duplicated().any():
        # Unweighted mean of the source columns, as merge_reference_types.py does. These
        # come from ONE run at one scale, and the clusters being merged are of similar size
        # (retina c/k/i: 8,282 / 7,589 / 7,481), so weighting by cell count would move the
        # result by about 2% -- not worth a dependency on per-cell labels being present.
        collapsed = {}
        for new_name in dict.fromkeys(out.columns):
            block = out.loc[:, out.columns == new_name]
            if block.shape[1] > 1:
                letters = sorted(l for l, v in rename.items() if v == new_name)
                print(f"  '{new_name}': mean of {len(letters)} de-novo profiles {letters}")
            collapsed[new_name] = block.mean(axis=1)
        out = pd.DataFrame(collapsed, index=out.index)

    # Guards: renamed de-novo must not collide with any named type, and no dup columns.
    collisions = set(rename.values()) & set(named)
    assert not collisions, f"renamed de-novo collide with named types: {collisions}"
    assert not out.columns.duplicated().any(), "duplicate columns in output profiles"

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output)
    print(f"Wrote {args.output} ({out.shape[0]:,} genes x {out.shape[1]} cell types: "
          f"{len(kept) - len(rename)} named + {len(set(rename.values()))} de-novo)")
    print("Every kept new_name must be a leaf of the hierarchy this matrix will be "
          "used with, and vice versa. 83_prep_insitutree_profiles.sh checks that "
          "and fails if it does not hold.")


if __name__ == "__main__":
    main()
