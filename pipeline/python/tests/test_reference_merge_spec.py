#!/usr/bin/env python3
"""The refit merge spec must be a COMPLETE, applicable ledger, and `keep` must be inert.

The spec is applied on another machine, against a reference regenerated from raw atlas
H5ADs, and the result is what every later stage types against. There is no opportunity to
discover there that a type was missing, spelled differently, or merged into a name something
else already uses -- so each of those is checked here against the committed reference.

    uv run python pipeline/python/tests/test_reference_merge_spec.py
"""
import csv
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

_PY_DIR = Path(__file__).resolve().parents[1]
_REFERENCE_DIR = _PY_DIR.parent / "reference"
SCRIPT = _PY_DIR / "merge_reference_types.py"
SPEC = _REFERENCE_DIR / "retina_type_merges_refit.csv"
REFERENCE = _REFERENCE_DIR / "retina_combined_panel.csv"

sys.path.insert(0, str(_PY_DIR))
_spec = importlib.util.spec_from_file_location("merge_reference_types", SCRIPT)
mrt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mrt)

EXPECTED_AFTER_REFIT = 53
RULED_OUT = ("Mona_T cell", "Mona_B cell", "Mona_natural killer cell",
             "Allen_Hippocampal CA1-3", "Allen_Medium spiny neuron",
             "Allen_Eccentric medium spiny neuron")


def reference_types() -> list[str]:
    with REFERENCE.open() as fh:
        return next(csv.reader(fh))[1:]


def run(reference: Path, spec: Path, output: Path, *extra: str):
    return subprocess.run([sys.executable, str(SCRIPT), "--reference", str(reference),
                           "--spec", str(spec), "--output", str(output), *extra],
                          capture_output=True, text=True)


def test_the_spec_covers_the_reference_exactly():
    spec = pd.read_csv(SPEC).fillna("")
    assert list(spec.columns)[:3] == ["action", "source_type", "merged_name"], list(spec.columns)
    types = reference_types()
    assert not spec["source_type"].duplicated().any()
    assert set(spec["source_type"]) == set(types), \
        sorted(set(types) ^ set(spec["source_type"]))
    assert (spec["decided_by"] != "").all() and (spec["note"] != "").all(), \
        "every row carries who decided it and why"
    assert set(spec["status"]) <= {"decided", "pending"}
    dropped = set(spec.loc[spec["action"] == "drop", "source_type"])
    for ruled_out in RULED_OUT:
        assert ruled_out in dropped, f"{ruled_out} was ruled out but is not dropped"
    print(f"  {len(types)} types, {len(spec)} rows, "
          f"{int((spec['status'] == 'decided').sum())} decided")


def test_no_retinal_type_moves_before_pepple_rules():
    spec = pd.read_csv(SPEC).fillna("")
    retina = spec[spec["source_type"].str.startswith("HRCA_")]
    moved = set(retina.loc[retina["action"] != "keep", "source_type"])
    # The one exception is the microglia duplicate, which 79 decided on tissue, not granularity.
    assert moved <= {"HRCA_microglial cell"}, moved
    pending_on_her = spec[spec["decided_by"] == "Kathryn Pepple"]
    assert len(pending_on_her) == 27, len(pending_on_her)
    assert (pending_on_her["action"] == "keep").all()
    print(f"  {len(pending_on_her)} retinal types held unchanged pending her answer")


def test_applying_it_gives_the_type_list_we_claim():
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "refit_panel.csv"
        result = run(REFERENCE, SPEC, out, "--require-complete")
        assert result.returncode == 0, result.stderr
        after = pd.read_csv(out, index_col=0)
        assert after.shape[1] == EXPECTED_AFTER_REFIT, after.shape
        assert "Excitatory neuron" in after.columns
        for gone in (*RULED_OUT, "HRCA_microglial cell", "Mona_fibroblast",
                     "Allen_Upper-layer intratelencephalic"):
            assert gone not in after.columns, gone
        for kept in ("HRCA_astrocyte", "Allen_Astrocyte", "Mona_astrocyte",
                     "HRCA_retinal pigment epithelial cell", "Allen_Microglia"):
            assert kept in after.columns, kept
        before = pd.read_csv(REFERENCE, index_col=0)
        assert after.index.equals(before.index), "merging must not touch the gene index"
        print(f"  64 -> {after.shape[1]} columns, gene index unchanged")


def _tiny_reference(path: Path):
    pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0], "c": [5.0, 6.0]},
                 index=["G1", "G2"]).to_csv(path, index_label="gene")


def _write_spec(path: Path, rows):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["action", "source_type", "merged_name"])
        w.writerows(rows)


def test_keep_is_inert_and_completeness_is_enforced():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        reference, out = tmp / "ref.csv", tmp / "out.csv"
        _tiny_reference(reference)

        # keep leaves the column exactly as it was.
        spec = tmp / "keep.csv"
        _write_spec(spec, [["keep", "a", ""], ["keep", "b", ""], ["drop", "c", ""]])
        assert run(reference, spec, out, "--require-complete").returncode == 0
        after = pd.read_csv(out, index_col=0)
        assert list(after.columns) == ["a", "b"]
        assert after["a"].tolist() == [1.0, 2.0]

        # --require-complete catches the type nobody ruled on.
        partial = tmp / "partial.csv"
        _write_spec(partial, [["drop", "c", ""]])
        assert run(reference, partial, out).returncode == 0, "silent without the flag"
        strict = run(reference, partial, out, "--require-complete")
        assert strict.returncode != 0 and "no spec row" in strict.stdout + strict.stderr

        # A keep row must not smuggle in a merged_name.
        bad = tmp / "bad.csv"
        _write_spec(bad, [["keep", "a", "Something"], ["keep", "b", ""], ["keep", "c", ""]])
        failed = run(reference, bad, out)
        assert failed.returncode != 0 and "must not set merged_name" in failed.stdout + failed.stderr

        # And a merge must not land on the name of a type that survives via keep.
        collide = tmp / "collide.csv"
        _write_spec(collide, [["merge", "a", "c"], ["merge", "b", "c"], ["keep", "c", ""]])
        clash = run(reference, collide, out)
        assert clash.returncode != 0 and "collide" in clash.stdout + clash.stderr
        print("  keep is inert; completeness, merged_name and collisions all enforced")


if __name__ == "__main__":
    for test in (test_the_spec_covers_the_reference_exactly,
                 test_no_retinal_type_moves_before_pepple_rules,
                 test_applying_it_gives_the_type_list_we_claim,
                 test_keep_is_inert_and_completeness_is_enforced):
        print(test.__name__)
        test()
    print("\nall reference merge spec tests passed")
