#!/usr/bin/env python3
"""Tests for denovo_vs_reference_cosine.py on a synthetic InSituType result h5.

Plants one de novo cluster of each kind and requires the verdict to match:
  a = a noisy copy of named type Pvalb (known-like), b = the average of all named types (flat),
  c = expressed only on genes no named type uses (distinct), d = between two near-identical named
  types (known-like, ambiguous). Also pins the cosines, sizes and the reference-csv path.

Run:  uv run --with pandas --with numpy --with h5py python pipeline/python/tests/test_denovo_vs_reference_cosine.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "denovo_vs_reference_cosine.py"
N_GENES = 60
GENES = [f"G{i}" for i in range(N_GENES)]
CELLS_PER_TYPE = {"Pvalb": 50, "L4 IT": 40, "L4 IT B": 30, "Astro": 20, "a": 25, "b": 15, "c": 10, "d": 5}


def block_profile(start: int, stop: int, level: float = 10.0) -> np.ndarray:
    v = np.full(N_GENES, 0.05)
    v[start:stop] = level
    return v


def build_named() -> pd.DataFrame:
    rng = np.random.default_rng(0)
    l4 = block_profile(10, 20)
    return pd.DataFrame({
        "Pvalb": block_profile(0, 10),
        "L4 IT": l4,
        "L4 IT B": l4 + rng.normal(0, 0.05, N_GENES).clip(min=0),   # near-identical to L4 IT
        "Astro": block_profile(20, 30),
    }, index=GENES)


def build_denovo(named: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(1)
    return pd.DataFrame({
        "a": named["Pvalb"] + rng.normal(0, 0.3, N_GENES).clip(min=0),
        "b": named.mean(axis=1),
        "c": block_profile(40, 50),
        "d": (named["L4 IT"] + named["L4 IT B"]) / 2,
    }, index=GENES)


def write_h5(path: Path, profiles: pd.DataFrame) -> None:
    cell_types = [t for t, n in CELLS_PER_TYPE.items() for _ in range(n)]
    with h5py.File(path, "w") as f:
        f["profiles"] = profiles.to_numpy()                       # genes x types
        f["profile_genes"] = np.array(profiles.index, dtype="S")
        f["profile_types"] = np.array(profiles.columns, dtype="S")
        f["cell_id"] = np.array([f"c{i}" for i in range(len(cell_types))], dtype="S")
        f["cell_type"] = np.array(cell_types, dtype="S")


def run(*args: str) -> subprocess.CompletedProcess:
    cmd = [sys.executable, str(SCRIPT), *args]
    print("  $", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    print(res.stdout[-2500:], res.stderr[-1500:])
    return res


def main() -> None:
    named = build_named()
    denovo = build_denovo(named)
    # Fitted named columns are pulled toward the data under update_reference=TRUE; the reference
    # csv is the original. Make the fitted Pvalb deliberately different so the two modes differ.
    fitted = pd.concat([named, denovo], axis=1)
    fitted["Pvalb"] = fitted["Pvalb"] * 0 + block_profile(0, 10) + block_profile(40, 50, 6.0)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        h5, ref, out = tmp / "result.h5", tmp / "reference.csv", tmp / "out"
        write_h5(h5, fitted)
        named.to_csv(ref)

        print("test_verdicts_sizes_and_cosines_with_reference_csv")
        assert run("--result-h5", str(h5), "--reference-csv", str(ref),
                   "--output-dir", str(out)).returncode == 0
        t = pd.read_csv(out / "denovo_vs_reference.csv").set_index("cluster")
        assert list(t.index) == ["a", "b", "c", "d"], list(t.index)           # sorted by size
        assert t.loc["a", "n_cells"] == 25 and t.loc["d", "n_cells"] == 5
        assert t.loc["a", "best_named"] == "Pvalb" and t.loc["a", "best_cosine"] > 0.99
        assert t.loc["a", "verdict"] == "known-like", t.loc["a"]
        assert t.loc["b", "cosine_to_mean_named"] > 0.999 and t.loc["b", "verdict"] == "flat", t.loc["b"]
        assert t.loc["c", "verdict"] == "distinct" and t.loc["c", "best_cosine"] < 0.5, t.loc["c"]
        assert t.loc["d", "verdict"] == "known-like, ambiguous", t.loc["d"]
        assert t.loc["d", "margin"] < 0.02, t.loc["d"]
        assert abs(t["pct_cohort"].sum() - 100 * 55 / 195) < 0.05, t["pct_cohort"].sum()
        matrix = pd.read_csv(out / "denovo_vs_reference_matrix.csv", index_col=0)
        assert matrix.shape == (4, 4) and matrix.loc["a", "Pvalb"] == t.loc["a", "best_cosine"]

        print("test_fitted_columns_used_without_reference_csv")
        out2 = tmp / "out2"
        assert run("--result-h5", str(h5), "--output-dir", str(out2)).returncode == 0
        t2 = pd.read_csv(out2 / "denovo_vs_reference.csv").set_index("cluster")
        assert t2.loc["a", "best_cosine"] < t.loc["a", "best_cosine"], "fitted Pvalb was altered"

        print("test_reference_missing_a_named_type_is_refused")
        named.drop(columns=["Astro"]).to_csv(tmp / "short.csv")
        res = run("--result-h5", str(h5), "--reference-csv", str(tmp / "short.csv"),
                  "--output-dir", str(tmp / "out3"))
        assert res.returncode != 0 and "Astro" in (res.stdout + res.stderr)

        print("test_result_without_denovo_is_refused")
        write_h5(h5, named)
        res = run("--result-h5", str(h5), "--output-dir", str(tmp / "out4"))
        assert res.returncode != 0
    print("ALL PASSED")


if __name__ == "__main__":
    main()
