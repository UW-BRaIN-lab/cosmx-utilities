#!/usr/bin/env python3
"""Write the Napari colour of every label in an InSituTree hierarchy, for figures that must match Napari.

A label's colour in every Napari `_metadata.csv` column is napari_colors.deterministic_color: a DuckDB
hash of its name, lightened where it would be too dark on Napari's black background. The cell-type UMAP
job (pipeline/slurm/95_celltype_umap.sh) runs in a container with no duckdb, so it reads colours from a
committed label,hex_color table (--color-table / COLOR_TABLE_BASENAME) instead of hashing.

pipeline/reference/napari_label_colors.csv was EXTRACTED from the live Napari columns, because 14 retina
leaves were written before the legibility rule and carry the raw hash there. A hierarchy whose Napari
columns were all written with the current rule (k=27: verified on the 3 H&E slides, 70/70 labels) can be
GENERATED instead, which also covers the leaves no slide shows yet. scripts/tests/test_napari_label_colors.py
keeps the committed k=27 table equal to the rule.

Every node and leaf name is written, in depth-first order: InSituTree can report a node's name (it does for
the single-child Low_signal node), and that label gets a Napari colour like any other.

Usage:
    uv run --with duckdb python scripts/napari-label-colors.py \\
        --hierarchy pipeline/reference/insitutree_hierarchy_k27.json \\
        --output pipeline/reference/insitutree_k27_label_colors.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from napari_colors import deterministic_color  # noqa: E402

COLUMNS = ("label", "hex_color")


def hierarchy_names(node) -> list[str]:
    """Every node and leaf name, depth-first, each once."""
    names: list[str] = []
    if isinstance(node, list):
        names.extend(node)
    else:
        for name, child in node.items():
            names.append(name)
            names.extend(hierarchy_names(child))
    return list(dict.fromkeys(names))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hierarchy", type=Path, required=True, help="InSituTree hierarchy JSON.")
    parser.add_argument("--output", type=Path, required=True, help="label,hex_color CSV to write.")
    args = parser.parse_args()

    names = hierarchy_names(json.loads(args.hierarchy.read_text()))
    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows((name, deterministic_color(name)) for name in names)
    print(f"Wrote {len(names)} label colours to {args.output}")


if __name__ == "__main__":
    main()
