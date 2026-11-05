#!/usr/bin/env python3
"""Generate Sankey diagrams comparing cell type assignments between two CosMx experiments.

Each cell (matched by cell_ID) is traced from its type in experiment A to its type
in experiment B.  The diagram shows how cells redistribute across types.

Usage:
    uv run --with plotly --with kaleido python scripts/cell-type-sankey.py \
        --bucket keene-cosmx-data \
        --prefix-a napari-stitched/CosMx-GBM/CosMx-GBM-segmentation-test-1.9.26/BsbneurosegmentationnoDAPIfixedGBmaptest4726_09_04_2026_13_07_19_475 \
        --prefix-b napari-stitched/CosMx-GBM/CosMx-GBM-segmentation-test-1.9.26/BsbneurosegmentationnoDAPItest33026_01_04_2026_15_05_56_903 \
        --output-dir histograms
"""

import argparse
import csv
import io
import os
from collections import Counter

import boto3
import duckdb
import plotly.graph_objects as go


def deterministic_color(value: str) -> str:
    """Generate a deterministic hex color from a string using DuckDB's hash()."""
    result = duckdb.sql(
        "SELECT printf('#%06X', abs(hash($1)) % 16777216)",
        params=[value],
    ).fetchall()
    return result[0][0]


MIN_FLOW_FRACTION = 0.005  # flows below this fraction of total cells are grouped into "Other"


def list_slides(s3, bucket: str, prefix: str) -> list[str]:
    resp = s3.list_objects_v2(
        Bucket=bucket, Prefix=prefix.rstrip("/") + "/", Delimiter="/",
    )
    return [
        p["Prefix"].rstrip("/").rsplit("/", 1)[-1]
        for p in resp.get("CommonPrefixes", [])
    ]


def load_metadata(s3, bucket: str, prefix: str, slide: str) -> dict[str, str]:
    """Return {cell_ID: cell_type} for a slide."""
    key = f"{prefix.rstrip('/')}/{slide}/_metadata.csv"
    obj = s3.get_object(Bucket=bucket, Key=key)
    body = obj["Body"].read().decode("utf-8")
    reader = csv.DictReader(io.StringIO(body))
    return {row["cell_ID"]: row["cell_type"] for row in reader}


def short_experiment_name(prefix: str) -> str:
    """Extract a readable short name from the S3 prefix."""
    run_dir = prefix.rstrip("/").rsplit("/", 1)[-1]
    # Take everything before the date/timestamp portion
    # Pattern: name_DD_MM_YYYY_HH_MM_SS_mmm
    parts = run_dir.split("_")
    # Walk backwards to find where the timestamp starts (6 numeric groups at the end)
    name_parts = []
    for part in parts:
        if part.isdigit() and len(name_parts) > 0:
            break
        name_parts.append(part)
    return "_".join(name_parts) if name_parts else run_dir[:30]


def build_sankey(
    slide: str,
    cells_a: dict[str, str],
    cells_b: dict[str, str],
    label_a: str,
    label_b: str,
    output_path: str,
) -> dict:
    """Build and save a Sankey diagram for one slide. Returns stats."""
    shared_ids = set(cells_a) & set(cells_b)

    # Count flows: (type_a, type_b) -> count
    flow_counts = Counter()
    for cid in shared_ids:
        flow_counts[(cells_a[cid], cells_b[cid])] += 1

    total_cells = len(shared_ids)
    min_flow = int(total_cells * MIN_FLOW_FRACTION)

    # Identify types that participate in at least one significant flow
    significant_types_a = set()
    significant_types_b = set()
    significant_flows = {}
    other_flow_total = 0

    for (ta, tb), count in flow_counts.items():
        if count >= min_flow:
            significant_types_a.add(ta)
            significant_types_b.add(tb)
            significant_flows[(ta, tb)] = count
        else:
            other_flow_total += count

    # Sort types by total cell count (descending) for better layout
    type_a_totals = Counter()
    type_b_totals = Counter()
    for (ta, tb), count in significant_flows.items():
        type_a_totals[ta] += count
        type_b_totals[tb] += count

    sorted_types_a = [t for t, _ in type_a_totals.most_common()]
    sorted_types_b = [t for t, _ in type_b_totals.most_common()]

    # Build node lists: left side (experiment A) then right side (experiment B)
    nodes_a = sorted_types_a
    nodes_b = sorted_types_b
    all_nodes = [f"{t} " for t in nodes_a] + list(nodes_b)  # trailing space distinguishes left from right

    node_index = {name: i for i, name in enumerate(all_nodes)}

    # Build links
    sources = []
    targets = []
    values = []
    for (ta, tb), count in significant_flows.items():
        sources.append(node_index[f"{ta} "])
        targets.append(node_index[tb])
        values.append(count)

    # Color nodes: deterministic color per cell type name so the same type
    # gets the same color on both sides of the diagram
    node_colors = [deterministic_color(name.strip()) for name in all_nodes]

    # Link colors: semi-transparent version of source node color
    link_colors = []
    for s in sources:
        base = node_colors[s].lstrip("#")
        r, g, b = int(base[:2], 16), int(base[2:4], 16), int(base[4:6], 16)
        link_colors.append(f"rgba({r},{g},{b},0.3)")

    # Clean display labels (remove trailing space from left nodes)
    display_labels = [n.rstrip() for n in all_nodes]

    changed = sum(
        count for (ta, tb), count in flow_counts.items() if ta != tb
    )

    fig = go.Figure(data=[go.Sankey(
        arrangement="snap",
        node=dict(
            pad=12,
            thickness=18,
            line=dict(color="white", width=0.5),
            label=display_labels,
            color=node_colors,
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=link_colors,
        ),
    )])

    fig.update_layout(
        title=dict(
            text=(
                f"Cell Type Reclassification — {slide}<br>"
                f"<sub>{label_a}  →  {label_b}  |  "
                f"{total_cells:,} cells  |  {changed:,} changed ({100*changed/total_cells:.1f}%)</sub>"
            ),
            font=dict(size=16),
        ),
        font=dict(size=10),
        width=1400,
        height=max(800, len(all_nodes) * 18),
    )

    fig.write_image(output_path, scale=2)

    return {
        "total_cells": total_cells,
        "changed": changed,
        "types_a": len(type_a_totals),
        "types_b": len(type_b_totals),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate Sankey diagrams comparing cell types between two CosMx experiments.",
    )
    parser.add_argument("--bucket", required=True)
    parser.add_argument("--prefix-a", required=True, help="S3 prefix for experiment A (left side)")
    parser.add_argument("--prefix-b", required=True, help="S3 prefix for experiment B (right side)")
    parser.add_argument("--label-a", default=None, help="Display label for experiment A (auto-detected if omitted)")
    parser.add_argument("--label-b", default=None, help="Display label for experiment B (auto-detected if omitted)")
    parser.add_argument("--output-dir", default="sankey")
    args = parser.parse_args()

    s3 = boto3.client("s3")
    os.makedirs(args.output_dir, exist_ok=True)

    label_a = args.label_a or short_experiment_name(args.prefix_a)
    label_b = args.label_b or short_experiment_name(args.prefix_b)
    print(f"Experiment A: {label_a}")
    print(f"Experiment B: {label_b}")

    slides_a = set(list_slides(s3, args.bucket, args.prefix_a))
    slides_b = set(list_slides(s3, args.bucket, args.prefix_b))
    shared_slides = sorted(slides_a & slides_b)

    if not shared_slides:
        print("ERROR: No shared slides between the two experiments")
        return

    only_a = slides_a - slides_b
    only_b = slides_b - slides_a
    if only_a:
        print(f"  Slides only in A: {only_a}")
    if only_b:
        print(f"  Slides only in B: {only_b}")

    print(f"Shared slides: {shared_slides}\n")

    for slide in shared_slides:
        print(f"  {slide} ...", end=" ", flush=True)
        cells_a = load_metadata(s3, args.bucket, args.prefix_a, slide)
        cells_b = load_metadata(s3, args.bucket, args.prefix_b, slide)

        output_path = os.path.join(args.output_dir, f"{slide}_sankey.png")
        stats = build_sankey(slide, cells_a, cells_b, label_a, label_b, output_path)
        print(
            f"{stats['total_cells']:,} cells, "
            f"{stats['changed']:,} changed ({100*stats['changed']/stats['total_cells']:.1f}%)"
        )

    print(f"\nSankey diagrams saved to {args.output_dir}/")


if __name__ == "__main__":
    main()
