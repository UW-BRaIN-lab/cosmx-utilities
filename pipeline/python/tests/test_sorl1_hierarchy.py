#!/usr/bin/env python3
"""The SORL1 InSituTree hierarchy must match the profile matrix columns exactly.

R/insitutree_typing.R stops if a hierarchy leaf is not a profile column or a leaf is listed twice,
and ignores a profile column that is in no leaf (so a type would silently never be called). The
profile matrix is built by prep_insitutree_profiles.py as: every named reference type plus the de novo
clusters marked keep in the annotation CSV, renamed to new_name. This rebuilds that expected set from the
two committed inputs (Brain_6k panel reference header, stage4_k21.csv) and compares it to the JSON in
BOTH directions, so a rename in one file cannot drift from the other.

Run:  uv run --with pandas python pipeline/python/tests/test_sorl1_hierarchy.py
"""
import json
from pathlib import Path

import pandas as pd

REFERENCE_DIR = Path(__file__).resolve().parents[2] / "reference"
TRUTHY = {"true", "t", "1", "yes", "y"}


def leaves(node):
    if isinstance(node, list):
        return list(node)
    return [leaf for child in node.values() for leaf in leaves(child)]


def render(node, indent=0):
    lines = []
    for name, child in node.items():
        if isinstance(child, list):
            lines.append("  " * indent + f"{name}  [{len(child)} leaves]")
            lines += ["  " * (indent + 1) + f"- {leaf}" for leaf in child]
        else:
            lines.append("  " * indent + name)
            lines += render(child, indent + 1)
    return lines


def expected_columns():
    named = list(pd.read_csv(REFERENCE_DIR / "brain6k_panel.csv", nrows=0).columns[1:])
    ann = pd.read_csv(REFERENCE_DIR / "denovo_annotations" / "stage4_k21.csv", dtype=str, keep_default_na=False)
    kept = ann[ann["keep"].str.strip().str.lower().isin(TRUTHY)]["new_name"].str.strip().tolist()
    return named, kept


def test_hierarchy_leaves_match_profile_columns_both_ways():
    hierarchy = json.load(open(REFERENCE_DIR / "sorl1_hierarchy.json"))
    named, kept = expected_columns()
    print("\n".join(render(hierarchy)))
    got = leaves(hierarchy)
    dup = sorted({x for x in got if got.count(x) > 1})
    print(f"{len(got)} leaves; expected {len(named)} named + {len(kept)} kept de novo")
    assert not dup, f"leaves listed more than once: {dup}"
    missing = sorted(set(named + kept) - set(got))
    extra = sorted(set(got) - set(named + kept))
    assert not missing, f"profile columns in no leaf (never called): {missing}"
    assert not extra, f"leaves with no profile column (R will stop): {extra}"
    assert len(named) == 13 and len(kept) == 15


def test_kept_names_are_unique_and_collision_free():
    named, kept = expected_columns()
    assert len(set(kept)) == len(kept), "duplicate new_name"
    assert not set(kept) & set(named), "a de novo new_name collides with a reference type"
    assert all(name.endswith("_denovo") for name in kept)


def test_no_empty_branches_and_sink_is_top_level():
    hierarchy = json.load(open(REFERENCE_DIR / "sorl1_hierarchy.json"))

    def walk(node, path):
        if isinstance(node, list):
            assert node, f"empty leaf list at {path}"
            return
        assert node, f"empty branch at {path}"
        for k, v in node.items():
            walk(v, path + "/" + k)
    walk(hierarchy, "")
    assert "Low_signal" in hierarchy and isinstance(hierarchy["Low_signal"], list)


if __name__ == "__main__":
    test_hierarchy_leaves_match_profile_columns_both_ways()
    test_kept_names_are_unique_and_collision_free()
    test_no_empty_branches_and_sink_is_top_level()
    print("OK")
