#!/usr/bin/env python3
"""Tests for reference_collinearity.py.

The fixture encodes the claim the script exists to support: the gene set decides the answer.
One pair of types differs ONLY in genes that gene selection drops, so on the full panel they
are separable and on the pruned panel they are not -- and the script must say so both times.

    uv run python pipeline/python/tests/test_reference_collinearity.py
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PY_DIR))
_spec = importlib.util.spec_from_file_location(
    "reference_collinearity", _PY_DIR / "reference_collinearity.py")
rc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rc)

N_SHARED = 60          # genes every type expresses the same way
N_PRIVATE = 12         # the genes that separate the "only via dropped genes" pair
KEPT = [f"KEPT{i:03d}" for i in range(N_SHARED)]
DROPPED = [f"DROP{i:03d}" for i in range(N_PRIVATE)]


def build_reference(seed: int = 3) -> pd.DataFrame:
    """Four types: a near-identical pair, a pair split only by dropped genes, one loner."""
    rng = np.random.default_rng(seed)
    base = rng.lognormal(mean=1.0, sigma=1.2, size=N_SHARED)
    loner = rng.lognormal(mean=1.0, sigma=1.2, size=N_SHARED)

    frame = pd.DataFrame(index=KEPT + DROPPED, dtype=float)
    # Twins: the same profile with 0.5% jitter. Inseparable on any gene set.
    frame["twin_a"] = np.concatenate([base, np.full(N_PRIVATE, 0.05)])
    frame["twin_b"] = np.concatenate([base * rng.normal(1.0, 0.005, N_SHARED),
                                      np.full(N_PRIVATE, 0.05)])
    # Split only by the dropped genes: identical on KEPT, opposite on DROP.
    frame["hidden_a"] = np.concatenate([base, np.concatenate([np.full(N_PRIVATE // 2, 40.0),
                                                              np.full(N_PRIVATE // 2, 0.05)])])
    frame["hidden_b"] = np.concatenate([base, np.concatenate([np.full(N_PRIVATE // 2, 0.05),
                                                              np.full(N_PRIVATE // 2, 40.0)])])
    frame["loner"] = np.concatenate([loner, np.full(N_PRIVATE, 0.05)])
    return frame


def write_fixture(tmp: Path) -> tuple[Path, Path]:
    reference = build_reference()
    ref_csv = tmp / "reference.csv"
    reference.to_csv(ref_csv, index_label="gene")
    keep = tmp / "kept_genes.txt"
    keep.write_text("\n".join(KEPT) + "\n")
    return ref_csv, keep


def test_pruning_decides_the_pair():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        ref_csv, keep = write_fixture(tmp)
        reference = rc.load_reference(ref_csv)

        full = rc.cosine_matrix(reference)
        pruned, unmatched = rc.restrict(reference, keep)
        pruned_cosines = rc.cosine_matrix(pruned)
        assert unmatched == 0

        hidden_full = full.loc["hidden_a", "hidden_b"]
        hidden_pruned = pruned_cosines.loc["hidden_a", "hidden_b"]
        print(f"hidden pair: full panel {hidden_full:.4f} -> pruned {hidden_pruned:.4f}")
        assert hidden_full < 0.99, "the dropped genes should separate them on the full panel"
        assert hidden_pruned > 0.9999, "on the kept genes they are the same profile"

        twins = pruned_cosines.loc["twin_a", "twin_b"]
        assert twins > 0.99, f"twins should be collinear on any gene set, got {twins}"
        for other in ("twin_a", "hidden_a"):
            assert pruned_cosines.loc["loner", other] < 0.99, "the loner must stay separable"


def test_components_group_the_inseparable_and_flag_chains():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        ref_csv, keep = write_fixture(tmp)
        pruned, _ = rc.restrict(rc.load_reference(ref_csv), keep)
        matrix = rc.cosine_matrix(pruned)

        groups = rc.components(matrix, 0.99)
        flat = {t for group in groups for t in group}
        assert "loner" not in flat, "a separable type must not be in a merge group"
        # On the kept genes all four of the others collapse into one chained component.
        assert any(set(group) >= {"twin_a", "twin_b"} for group in groups)
        assert any(set(group) >= {"hidden_a", "hidden_b"} for group in groups)

        # A bar above every pair proposes nothing.
        assert rc.components(matrix, 1.001) == []


def test_cli_writes_four_tables_with_the_delta_column():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        ref_csv, keep = write_fixture(tmp)
        out = tmp / "out"
        assert rc.main(["--reference", str(ref_csv), "--keep-genes", str(keep),
                        "--output-dir", str(out), "--max-component-size", "2"]) == 0

        for name in ("reference_collinearity_matrix.csv", "reference_collinearity_nearest.csv",
                     "reference_collinearity_pairs.csv", "reference_merge_candidates.csv"):
            assert (out / name).exists(), f"missing {name}"

        pairs = pd.read_csv(out / "reference_collinearity_pairs.csv")
        assert {"cosine_full_panel", "delta"} <= set(pairs.columns)
        hidden = pairs[(pairs.type_a.isin(["hidden_a", "hidden_b"])) &
                       (pairs.type_b.isin(["hidden_a", "hidden_b"]))]
        assert len(hidden) == 1
        assert hidden["delta"].iloc[0] > 0, "pruning made this pair LESS separable"
        print(f"delta for the hidden pair: {hidden['delta'].iloc[0]:+.4f}")

        nearest = pd.read_csv(out / "reference_collinearity_nearest.csv")
        assert nearest.iloc[0]["type"] == "loner", "most distinguishable type should sort first"

        candidates = pd.read_csv(out / "reference_merge_candidates.csv")
        # max-component-size 2 and a 4-way chain on the kept genes -> flagged, not proposed.
        assert candidates["chained"].any(), "an oversized group must be flagged as chained"


def test_rejects_a_reference_with_a_dead_profile():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        frame = build_reference()
        frame["dead"] = 0.0
        path = tmp / "dead.csv"
        frame.to_csv(path, index_label="gene")
        try:
            rc.cosine_matrix(rc.load_reference(path))
        except SystemExit as exc:
            assert "all-zero" in str(exc)
        else:
            raise AssertionError("an all-zero profile must not be scored")


if __name__ == "__main__":
    test_pruning_decides_the_pair()
    test_components_group_the_inseparable_and_flag_chains()
    test_cli_writes_four_tables_with_the_delta_column()
    test_rejects_a_reference_with_a_dead_profile()
    print("\nall reference_collinearity tests passed")
