#!/usr/bin/env python3
"""Tests for marker_overlap_vs_reference.py on a synthetic InSituType result and reference.

Result types: a named `Oligodendrocyte` (markers = genes 0-9, the positive control), de novo `x`
(markers = genes 10-19), `y` (genes 50-59, which the reference has no cell type for) and `z` (only 3
genes above the expression floor, so it must be skipped, not scored). The reference gives types R_oligo,
R_astro, R_neuron and R_other genes 0-9, 10-19, 20-29 and 30-39 at a very different scale (x500) than the
fitted profiles, which is the platform gap this method must be immune to.

Run:  uv run --with pandas --with numpy --with h5py python pipeline/python/tests/test_marker_overlap_vs_reference.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from marker_overlap_vs_reference import hypergeom_sf  # noqa: E402

SCRIPT = Path(__file__).resolve().parents[1] / "marker_overlap_vs_reference.py"
N_GENES = 100
GENES = [f"G{i}" for i in range(N_GENES)]
CELLS = {"Oligodendrocyte": 40, "x": 30, "y": 20, "z": 10, "Neuron": 50}


def profile(high: range, level: float, background: float) -> np.ndarray:
    v = np.full(N_GENES, background)
    v[list(high)] = level
    return v


def build_result() -> pd.DataFrame:
    return pd.DataFrame({
        "Oligodendrocyte": profile(range(0, 10), 5.0, 0.05),
        "x": profile(range(10, 20), 5.0, 0.05),
        "y": profile(range(50, 60), 5.0, 0.05),
        "z": profile(range(70, 73), 5.0, 0.05),
        "Neuron": profile(range(80, 90), 5.0, 0.05),
    }, index=GENES)


def build_reference() -> pd.DataFrame:
    return pd.DataFrame({
        "R_oligo": profile(range(0, 10), 500.0, 5.0),
        "R_astro": profile(range(10, 20), 500.0, 5.0),
        "R_neuron": profile(range(20, 30), 500.0, 5.0),
        "R_other": profile(range(30, 40), 500.0, 5.0),
    }, index=GENES)


def write_h5(path: Path, profiles: pd.DataFrame) -> None:
    cell_types = [t for t, n in CELLS.items() for _ in range(n)]
    with h5py.File(path, "w") as f:
        f["profiles"] = profiles.to_numpy()
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
    print("test_hypergeom_matches_hand_computation")
    # universe 5, 2 successes, 2 draws: P(X>=1) = 1 - C(3,2)/C(5,2) = 0.7 ; P(X>=2) = 1/10
    assert abs(hypergeom_sf(1, 5, 2, 2) - 0.7) < 1e-12 and abs(hypergeom_sf(2, 5, 2, 2) - 0.1) < 1e-12
    assert hypergeom_sf(0, 5, 2, 2) == 1.0

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        h5, ref, out = tmp / "result.h5", tmp / "reference.csv", tmp / "out"
        write_h5(h5, build_result())
        build_reference().to_csv(ref)

        print("test_clusters_land_on_the_right_reference_type_despite_the_scale_gap")
        assert run("--result-h5", str(h5), "--reference-csv", str(ref), "--output-dir", str(out)).returncode == 0
        t = pd.read_csv(out / "marker_overlap.csv")
        best = t[t["rank"] == 1].set_index("cluster")
        assert best.loc["Oligodendrocyte", "reference_type"] == "R_oligo", best
        assert best.loc["x", "reference_type"] == "R_astro", best
        for c in ("Oligodendrocyte", "x"):
            assert best.loc[c, "concordance"] > 0.95 and best.loc[c, "overlap"] == 10, best.loc[c]
            assert best.loc[c, "overlap_p"] < 1e-6, best.loc[c]   # 10/10 hits; list padded to 20 -> p = 1.1e-8
        assert best.loc["Oligodendrocyte", "kind"] == "named" and best.loc["x", "kind"] == "de novo"
        assert best.loc["x", "n_cells"] == 30

        print("test_cluster_without_a_reference_counterpart_scores_near_chance")
        assert best.loc["y", "concordance"] < 0.7 and best.loc["y", "overlap"] == 0, best.loc["y"]
        assert best.loc["y", "overlap_p"] == 1.0

        print("test_cluster_with_too_few_markers_is_skipped")
        assert "z" not in set(t["cluster"])

        matrix = pd.read_csv(out / "marker_concordance_matrix.csv", index_col=0)
        assert matrix.loc["x", "R_astro"] == best.loc["x", "concordance"]
        assert matrix.shape == (4, 4), matrix.shape      # Oligodendrocyte, x, y, Neuron vs 4 reference types

        print("test_reference_with_other_gene_names_is_refused")
        build_reference().set_axis([f"ENSG{i}" for i in range(N_GENES)]).to_csv(tmp / "ens.csv")
        res = run("--result-h5", str(h5), "--reference-csv", str(tmp / "ens.csv"),
                  "--output-dir", str(tmp / "out2"))
        assert res.returncode != 0 and "gene naming" in (res.stdout + res.stderr)
    print("ALL PASSED")


if __name__ == "__main__":
    main()
