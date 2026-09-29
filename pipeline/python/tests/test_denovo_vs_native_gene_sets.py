#!/usr/bin/env python3
"""Tests for --gene-set scoring in denovo_vs_native_pseudobulk.py.

Reuses the planted STRESS world from test_denovo_vs_native_siblings (every forced cell carries a
shared stress programme; B's natives carry some of it; A's none). A set of the stress genes must
score highest on the forced cells and lowest on A's natives, while a control set of genes nobody
expresses differently must score ~0 everywhere -- otherwise the score would just be measuring
depth. Also: keyword expansion, custom sets, bad input, and off-panel genes being reported.

Run: python pipeline/python/tests/test_denovo_vs_native_gene_sets.py   (or under pytest)
"""

import sys
from pathlib import Path

import pandas as pd

PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PY_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import denovo_vs_native_pseudobulk as dvn  # noqa: E402
import test_denovo_vs_native_siblings as world  # noqa: E402

STRESS_SET = "stress=" + ",".join(world.STRESS_NAMES)
CONTROL_SET = "control=" + ",".join(f"G{i}" for i in range(30, 36))


def test_parse_keywords_custom_and_errors():
    got = dvn.parse_gene_sets(["heat-shock", "Hypoxia", "mine=A, B"])
    assert got["heat-shock"] == dvn.HEAT_SHOCK_GENES
    assert got["hypoxia"] == dvn.HYPOXIA_GENES
    assert got["mine"] == ("A", "B")
    for bad in ("nonsense", "=A,B", "x="):
        try:
            dvn.parse_gene_sets([bad])
        except ValueError:
            continue
        raise AssertionError(f"{bad!r} must be rejected")


def test_stress_set_scores_forced_cells_high_and_control_flat():
    out = world._run(world._stress_world(), "--gene-set", STRESS_SET, "--gene-set", CONTROL_SET,
                     top_n=6)
    sc = pd.read_csv(out / "gene_set_scores.csv")
    print(sc[sc.destination == world.DEST_B].to_string(index=False))
    b = sc[(sc.destination == world.DEST_B) & (sc.gene_set == "stress")].set_index("group")["score"]
    forced_b = b[f"{world.LETTER}->{world.DEST_B}"]
    nat_b = b[f"{world.DEST_B} [native, unpinned]"]
    assert forced_b > nat_b > 1.0, (forced_b, nat_b)      # B's natives share some, forced have more
    a = sc[(sc.destination == world.DEST_A) & (sc.gene_set == "stress")].set_index("group")["score"]
    # A's baseline is B's natives (which share some stress), so A's forced cells sit only a little
    # above it, but far above A's own natives, which carry none.
    assert a[f"{world.LETTER}->{world.DEST_A}"] > 0.3, a.to_dict()
    assert a[f"{world.LETTER}->{world.DEST_A}"] > a[f"{world.DEST_A} [native, unpinned]"] + 2.0, a.to_dict()
    control = sc[sc.gene_set == "control"]["score"]
    # Not exactly zero: per-cell library-size normalisation pulls every unrelated gene down a little
    # in cells that carry a large extra programme. It must stay small next to the real effect.
    assert control.abs().max() < 0.4, control.describe()
    assert control.abs().max() < 0.15 * (forced_b - nat_b + 3.0), control.describe()
    assert (sc["n_genes_used"] == sc["n_genes_listed"]).all()


def test_off_panel_genes_are_reported_and_counted():
    spec = "ghost=G20,G21,NOTAGENE"
    out = world._run(world._stress_world(), "--gene-set", spec, top_n=6)
    sc = pd.read_csv(out / "gene_set_scores.csv")
    assert set(sc["n_genes_used"]) == {2} and set(sc["n_genes_listed"]) == {3}, sc.head()


def test_without_gene_set_no_file_is_written():
    out = world._run(world._stress_world(), top_n=6)
    assert not (out / "gene_set_scores.csv").exists()


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all tests passed")
