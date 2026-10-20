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
    assert len(named) == 13 and len(kept) == 9


def test_only_distinct_de_novo_leaves_are_kept():
    """A de novo leaf must differ from the leaves it competes with (>= 10 genes at 2-fold over its siblings).
    q, h, r, d (excitatory), n, t (glial) failed that and are not leaves; l is kept by judgment (real VIP markers).
    o ALSO failed it (5 neuronal genes) but was restored after the first InSituTree run: without it the reference
    Oligodendrocyte profile claimed only 57% of the Leiden-2 oligodendrocytes and 23% went to OPC.
    n was tested as OPC_low_count_denovo beside OPC (v3) and removed: the split ran the wrong way and gained nothing
    (it has 0 genes at 2-fold over reference OPC); OPC calls are handled by the opc_confident flag instead."""
    _, kept = expected_columns()
    assert set(kept) == {"Excitatory_PCP4_L5_denovo", "Neuron_NR4A2_SYNPR_denovo", "Excitatory_RORB_L4_denovo",
                         "Excitatory_BCL11B_deep_denovo", "Neuron_SEMA3E_HS3ST4_denovo", "VIP_interneuron_denovo",
                         "Oligodendrocyte_denovo", "Low_signal_neuronal_denovo", "Low_signal_generic_denovo"}


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
    # Top-level sink is the generic flat cluster; the neuron-leaning flat cluster sits in the Neuron branch,
    # where it competes with the neuron types it resembles (cosine 0.93 to Excitatory_CCK, 0.67 to the generic sink).
    assert hierarchy["Low_signal"] == ["Low_signal_generic_denovo"]
    assert hierarchy["Neuron"]["Neuron_low_signal"] == ["Low_signal_neuronal_denovo"]


def test_excitatory_branch_is_flat():
    """Cosine similarity does not support sub-groups among the excitatory leaves (within 0.941 vs between 0.924)."""
    hierarchy = json.load(open(REFERENCE_DIR / "sorl1_hierarchy.json"))
    excitatory = hierarchy["Neuron"]["Excitatory neuron"]
    assert isinstance(excitatory, list) and len(excitatory) == 8


# ---- v4: the revision after Maddie's spatial review (drops the Brain_6k L6 reference, adds h, r, t, and sink leaves) ----
V4_DROPPED_REFERENCE = ["L6 neuron"]
V4_SINK_LEAVES = {f"Neuron_unresolved_{x}_denovo" for x in "djqfb"} | {"Low_signal_neuronal_denovo"}


def v4_expected_columns():
    named, _ = expected_columns()
    named = [n for n in named if n not in V4_DROPPED_REFERENCE]
    ann = pd.read_csv(REFERENCE_DIR / "denovo_annotations" / "stage4_k21_v4.csv", dtype=str, keep_default_na=False)
    kept = ann[ann["keep"].str.strip().str.lower().isin(TRUTHY)]["new_name"].str.strip().tolist()
    return named, kept


def test_v4_leaves_match_profile_columns_and_the_built_matrix():
    hierarchy = json.load(open(REFERENCE_DIR / "sorl1_hierarchy_v4.json"))
    named, kept = v4_expected_columns()
    got = leaves(hierarchy)
    print("\n".join(render(hierarchy)))
    print(f"{len(got)} leaves; {len(named)} named + {len(kept)} kept de novo")
    assert sorted(got) == sorted(named + kept) and len(got) == len(set(got)) == 29
    matrix = list(pd.read_csv(REFERENCE_DIR / "sorl1_insitutree_profiles_v4.csv", nrows=0).columns[1:])
    assert sorted(matrix) == sorted(got), "the committed v4 profile matrix drifted from the hierarchy"


def test_v4_structure_follows_the_review():
    h = json.load(open(REFERENCE_DIR / "sorl1_hierarchy_v4.json"))
    got = set(leaves(h))
    assert "L6 neuron" not in got                                   # the L6 reference was a glia/shallow-cell attractor
    assert len(h["Neuron"]["Excitatory_deep"]) == 6                 # c, m, r, h, p, g
    assert {"Excitatory_CARTPT_GAP43_denovo", "Excitatory_SNCA_CD24_denovo"} <= set(h["Neuron"]["Excitatory_deep"])
    assert "Astrocyte_NEAT1_HSPA1A_denovo" in h["Glia"]["Astrocyte"]
    assert set(h["Neuron"]["Neuron_low_signal"]) == V4_SINK_LEAVES  # d, j, q, f, b + i
    assert not any("OPC_denovo" in leaf or "OPC_low" in leaf for leaf in got)       # n stays out (v3 rejected it)


def test_v4_collapse_map_covers_exactly_the_neuron_sink():
    collapse = pd.read_csv(REFERENCE_DIR / "sorl1_hierarchy_v4_collapse.csv", dtype=str)
    assert set(collapse["leaf"]) == V4_SINK_LEAVES and set(collapse["collapsed_label"]) == {"Neuron_low_signal"}


def test_select_columns_drops_named_and_rejects_typos():
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import prep_insitutree_profiles as prep
    assert prep.select_columns(["A", "L6 neuron", "B"], {"x": "X_denovo"}, ["L6 neuron"]) == ["A", "B", "x"]
    assert prep.select_columns(["A", "B"], {}, []) == ["A", "B"]
    try:
        prep.select_columns(["A", "B"], {}, ["L6 neurn"])
    except SystemExit as err:
        assert "L6 neurn" in str(err)
    else:
        raise AssertionError("a typo in --drop-named must stop the build, not silently keep the type")


if __name__ == "__main__":
    test_hierarchy_leaves_match_profile_columns_both_ways()
    test_only_distinct_de_novo_leaves_are_kept()
    test_kept_names_are_unique_and_collision_free()
    test_no_empty_branches_and_sink_is_top_level()
    test_excitatory_branch_is_flat()
    test_v4_leaves_match_profile_columns_and_the_built_matrix()
    test_v4_structure_follows_the_review()
    test_v4_collapse_map_covers_exactly_the_neuron_sink()
    test_select_columns_drops_named_and_rejects_typos()
    print("OK")
