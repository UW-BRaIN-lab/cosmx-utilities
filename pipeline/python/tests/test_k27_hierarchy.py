#!/usr/bin/env python3
"""The full-cohort InSituTree reference (k=27 anchor fit + the PI's annotations) must stay in lockstep.

Three committed files describe one reference: the annotation table (which k=27 letters are kept and their names),
the profile matrix built from it by prep_insitutree_profiles.py, and the hierarchy JSON. R/insitutree_typing.R stops
if a leaf has no profile column and silently never calls a profile column that is in no leaf, so the three are
compared in both directions here. The expected column set is rebuilt from the Core-L4 panel header and the
annotation table, so a rename in one file cannot drift from the others.

The hierarchy follows profile similarity (InSituTree's plotProfileSimilarity metric: cosine on the linear profiles),
not lineage, as the InSituTree authors advise. Two consequences are pinned below because they look like mistakes:
GBmap Astrocyte and the reactive astrocytes sit with the malignant AC/MES states, and normal OPC (GBmap's and `j`)
competes with the malignant OPC/NPC states.

Run:  uv run --no-project --with pandas --with numpy pytest -s pipeline/python/tests/test_k27_hierarchy.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

REFERENCE_DIR = Path(__file__).resolve().parents[2] / "reference"
HIERARCHY = REFERENCE_DIR / "insitutree_hierarchy_k27.json"
PROFILES = REFERENCE_DIR / "insitutree_profiles_k27.csv"
ANNOTATIONS = REFERENCE_DIR / "denovo_annotations" / "fullcohort_pruned_k27.csv"
NEURON_MAP = REFERENCE_DIR / "insitutree_hierarchy_k27_neuron_ei.csv"
TRUTHY = {"true", "t", "1", "yes", "y"}
NEAR_DUPLICATE_COSINE = 0.95
# Reference pairs already above the threshold; none of them involves a de novo leaf.
KNOWN_NEAR_DUPLICATES = {
    frozenset({"CD4_rest", "Mast"}),
    frozenset({"CD8_EM", "CD8_cytotoxic"}),
    frozenset({"TAM-MG_pro-infl_I", "TAM-MG_prolif"}),
    frozenset({"TAM-BDM_anti-infl", "TAM-MG_pro-infl_I"}),
    frozenset({"TAM-MG_aging_sig", "TAM-MG_pro-infl_I"}),
}


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


def load_hierarchy():
    return json.load(open(HIERARCHY))


def expected_columns():
    named = list(pd.read_csv(REFERENCE_DIR / "gbmap_level4_panel.csv", nrows=0).columns[1:])
    ann = pd.read_csv(ANNOTATIONS, dtype=str, keep_default_na=False)
    kept = ann[ann["keep"].str.strip().str.lower().isin(TRUTHY)]["new_name"].str.strip().tolist()
    return named, kept


def cosine_matrix():
    profiles = pd.read_csv(PROFILES, index_col=0)
    unit = profiles.values.T / np.linalg.norm(profiles.values.T, axis=1, keepdims=True)
    return pd.DataFrame(unit @ unit.T, index=profiles.columns, columns=profiles.columns)


def test_leaves_match_profile_columns_both_ways():
    hierarchy = load_hierarchy()
    named, kept = expected_columns()
    profile_columns = list(pd.read_csv(PROFILES, nrows=0).columns[1:])
    print("\n".join(render(hierarchy)))
    got = leaves(hierarchy)
    print(f"{len(got)} leaves; expected {len(named)} named + {len(kept)} kept de novo; matrix has {len(profile_columns)}")
    dup = sorted({x for x in got if got.count(x) > 1})
    assert not dup, f"leaves listed more than once: {dup}"
    assert sorted(profile_columns) == sorted(named + kept), "profile matrix is stale: re-run prep_insitutree_profiles.py"
    missing = sorted(set(profile_columns) - set(got))
    extra = sorted(set(got) - set(profile_columns))
    assert not missing, f"profile columns in no leaf (never called): {missing}"
    assert not extra, f"leaves with no profile column (R will stop): {extra}"
    assert len(named) == 54 and len(kept) == 22


def test_pi_decisions():
    """The PI's calls from the H&E review (2026-10-30): five drops, and every neuron letter kept for the E/I split."""
    ann = pd.read_csv(ANNOTATIONS, dtype=str, keep_default_na=False).set_index("denovo_label")
    dropped = set(ann.index[~ann["keep"].str.lower().isin(TRUTHY)])
    print(f"dropped: {sorted(dropped)}")
    assert dropped == {"aa", "h", "l", "r", "u"}
    assert ann.loc["j", "new_name"] == "OPC_denovo"
    assert ann.loc["z", "new_name"] == "Reactive_astrocyte_denovo"
    assert ann.loc["k", "new_name"] == "AClike_nonEGFRamp_denovo"
    assert (ann["pi_annotation"] != "").all(), "every letter carries the PI's call (b says it was not annotated)"


def test_kept_names_are_unique_collision_free_and_new():
    named, kept = expected_columns()
    assert len(set(kept)) == len(kept), "duplicate new_name"
    assert not set(kept) & set(named), "a de novo new_name collides with a reference type"
    assert all(name.endswith("_denovo") for name in kept)
    # The pilot reference (insitutree_hierarchy.json) used these names for DIFFERENT populations of the Wenyu run.
    # Low_signal_denovo is the one deliberate reuse: it names the same thing, the flat sink.
    pilot = {"Stress_denovo", "MES_AClike_denovo", "Hypoxia_denovo", "OPClike_denovo", "MESlike_denovo"}
    assert not set(kept) & pilot, f"reuses a pilot name for a different population: {set(kept) & pilot}"


def test_no_empty_branches_and_sink_is_top_level():
    hierarchy = load_hierarchy()

    def walk(node, path):
        if isinstance(node, list):
            assert node, f"empty leaf list at {path}"
            return
        assert node, f"empty branch at {path}"
        for name, child in node.items():
            walk(child, f"{path}/{name}")
    walk(hierarchy, "")
    assert hierarchy["Low_signal"] == ["Low_signal_denovo"]


def test_placements_follow_expression_not_lineage():
    hierarchy = load_hierarchy()
    astro_ac_mes = leaves(hierarchy["Astro_AC_MES"])
    assert {"Astrocyte", "Reactive_astrocyte_denovo"} <= set(astro_ac_mes)
    assert {"OPC", "OPC_denovo", "NPC-like_OPC"} <= set(hierarchy["OPC_NPC"])
    cos = cosine_matrix()
    print(f"OPC_denovo -> NPC-like_OPC {cos.loc['OPC_denovo', 'NPC-like_OPC']:.3f}; "
          f"Astrocyte -> Reactive_astrocyte_denovo {cos.loc['Astrocyte', 'Reactive_astrocyte_denovo']:.3f}")


def test_neuron_branch_is_flat_and_ei_map_covers_it():
    """Inhibitory s is as close to the excitatory letters (0.87-0.91) as they are to each other (0.91-0.93), so profile
    similarity does not support E/I sub-nodes. The split is made after typing with the label map instead."""
    hierarchy = load_hierarchy()
    neurons = hierarchy["Neuron"]
    assert isinstance(neurons, list) and len(neurons) == 5
    mapping = pd.read_csv(NEURON_MAP, dtype=str).set_index("leaf")["collapsed_label"]
    print(mapping.to_string())
    assert set(mapping.index) == set(neurons)
    assert set(mapping) == {"Excitatory_neuron", "Inhibitory_neuron"}


def test_no_de_novo_leaf_is_a_near_duplicate():
    cos = cosine_matrix()
    names = list(cos.columns)
    pairs = {frozenset({a, b}): cos.loc[a, b] for i, a in enumerate(names) for b in names[i + 1:]
             if cos.loc[a, b] >= NEAR_DUPLICATE_COSINE}
    for pair, value in sorted(pairs.items(), key=lambda kv: -kv[1]):
        print(f"  {' ~ '.join(sorted(pair))}  {value:.3f}")
    assert set(pairs) == KNOWN_NEAR_DUPLICATES
