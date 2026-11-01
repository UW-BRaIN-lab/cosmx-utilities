"""The committed k=27 label-colour table must equal the Napari colour rule.

The UMAP job reads pipeline/reference/insitutree_k27_label_colors.csv because it cannot import
napari_colors; if the table drifted from the rule, the UMAP and Napari would colour a cell type
differently with no error anywhere.

Run:  uv run --with duckdb --with pytest python -m pytest scripts/tests/test_napari_label_colors.py -s
"""
import csv
import importlib.util
import json
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
REFERENCE = SCRIPTS.parent / "pipeline" / "reference"
HIERARCHY = REFERENCE / "insitutree_hierarchy_k27.json"
TABLE = REFERENCE / "insitutree_k27_label_colors.csv"


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


nc = _load("napari_colors", "napari_colors.py")
nlc = _load("napari_label_colors", "napari-label-colors.py")


def _table():
    with open(TABLE, newline="") as handle:
        return list(csv.DictReader(handle))


def test_hierarchy_names_are_depth_first_and_unique():
    names = nlc.hierarchy_names({"A": {"A1": ["x", "y"]}, "B": ["z", "x"]})
    print(names)
    assert names == ["A", "A1", "x", "y", "B", "z"]


def test_table_lists_every_hierarchy_name_once_in_order():
    rows = _table()
    names = nlc.hierarchy_names(json.loads(HIERARCHY.read_text()))
    print(f"{len(rows)} rows; first {[r['label'] for r in rows[:4]]}")
    assert [r["label"] for r in rows] == names
    assert {"Low_signal", "Low_signal_denovo", "Reactive_astrocyte_denovo", "Astrocyte"} <= set(names)


def test_table_colours_equal_the_napari_rule():
    mismatches = [(r["label"], r["hex_color"], nc.deterministic_color(r["label"]))
                  for r in _table() if r["hex_color"] != nc.deterministic_color(r["label"])]
    for mismatch in mismatches:
        print("  drift:", mismatch)
    assert not mismatches, "regenerate with scripts/napari-label-colors.py"
