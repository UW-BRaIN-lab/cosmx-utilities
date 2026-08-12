#!/usr/bin/env python3
"""Compare two cell-type columns in stitched _metadata.csv files.

Reads the multi-column _metadata.csv files produced by generate-slide-metadata.py
(each cell has several annotation columns, e.g. celltype_refit / celltype_norefit,
each paired with a <column>_color) and, for every run (study) under an S3
experiment prefix, compares two of those columns per cell:

  - a grouped abundance bar chart (proportion of cells per type, column A vs B)
  - a Sankey tracing each cell from its column-A type to its column-B type

Slides are grouped by their run folder (the directory just above each slide), so
e.g. an Initial-segmentation run and a Resegmentation run are summarised
separately rather than pooled (use --pool-all to pool everything into one group).

Usage:
    uv run --with plotly --with kaleido python scripts/compare-cell-typing.py \
        --bucket keene-cosmx-data \
        --prefix napari-stitched/CosMx-retina/CosMx-retina-brain-segmentation-test-4.1.26 \
        --output-dir compare_out

    # Different columns / labels, and draw every flow (no pooling of small ones):
    uv run --with plotly --with kaleido python scripts/compare-cell-typing.py \
        --bucket keene-cosmx-data --prefix napari-stitched/CosMx-Victoria/<run> \
        --column-a celltype_refit --label-a "Refit=TRUE" \
        --column-b celltype_norefit --label-b "Refit=FALSE" \
        --min-flow-fraction 0 --output-dir compare_out
"""

import argparse
import csv
import hashlib
import io
import os
import re
import sys
from collections import Counter, defaultdict

import boto3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go

COLUMN_A_SERIES = "#4C72B0"
COLUMN_B_SERIES = "#DD8452"


def deterministic_color(value: str) -> str:
    """Fallback color when a <column>_color column is absent: stable per value."""
    h = hashlib.md5(value.encode()).hexdigest()
    return f"#{h[:6]}"


def short_run_name(run: str) -> str:
    """Strip a trailing AtoMx date/timestamp for a readable run label."""
    stripped = re.sub(r"[_0-9]+$", "", run)
    return stripped or run


def discover_slides(s3, bucket: str, prefix: str) -> dict:
    """Find every slide _metadata.csv under the prefix, grouped by run folder.

    Returns {run_name: [(slide_name, key), ...]}. The run is the directory
    immediately above the slide (`.../<run>/<slide>/_metadata.csv`); slides
    sitting directly under the prefix are grouped under the prefix's last part.
    """
    groups = defaultdict(list)
    paginator = s3.get_paginator("list_objects_v2")
    base = prefix.rstrip("/")
    for page in paginator.paginate(Bucket=bucket, Prefix=base + "/"):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if not key.endswith("/_metadata.csv"):
                continue
            parts = key.split("/")
            slide = parts[-2]
            run = parts[-3] if len(parts) >= 3 else base.rsplit("/", 1)[-1]
            groups[run].append((slide, key))
    return groups


class Comparison:
    """Accumulates per-cell (column-A, column-B) pairs across a set of slides."""

    def __init__(self, col_a: str, col_b: str):
        self.col_a, self.col_b = col_a, col_b
        self.a_counts, self.b_counts = Counter(), Counter()
        self.pairs = Counter()
        self.a_color, self.b_color = {}, {}
        self.n_cells = 0

    def add_slide(self, s3, bucket: str, key: str) -> None:
        body = s3.get_object(Bucket=bucket, Key=key)["Body"].read().decode("utf-8")
        reader = csv.DictReader(io.StringIO(body))
        if self.col_a not in (reader.fieldnames or []) or \
           self.col_b not in (reader.fieldnames or []):
            print(f"  WARNING: {key} lacks {self.col_a}/{self.col_b} — skipped",
                  file=sys.stderr)
            return
        a_clr, b_clr = f"{self.col_a}_color", f"{self.col_b}_color"
        for row in reader:
            a, b = row.get(self.col_a, "").strip(), row.get(self.col_b, "").strip()
            if a:
                self.a_counts[a] += 1
                self.a_color.setdefault(a, row.get(a_clr, "").strip() or deterministic_color(a))
            if b:
                self.b_counts[b] += 1
                self.b_color.setdefault(b, row.get(b_clr, "").strip() or deterministic_color(b))
            if a and b:
                self.pairs[(a, b)] += 1
                self.n_cells += 1

    def agreement(self) -> float:
        if not self.n_cells:
            return 0.0
        same = sum(c for (a, b), c in self.pairs.items() if a == b)
        return 100 * same / self.n_cells


def abundance_chart(cmp: Comparison, label_a, label_b, title, out_path) -> None:
    labels = sorted(set(cmp.a_counts) | set(cmp.b_counts),
                    key=lambda l: -(cmp.a_counts.get(l, 0) + cmp.b_counts.get(l, 0)))
    a_prop = [100 * cmp.a_counts.get(l, 0) / cmp.n_cells for l in labels]
    b_prop = [100 * cmp.b_counts.get(l, 0) / cmp.n_cells for l in labels]
    y = range(len(labels)); h = 0.4
    fig, ax = plt.subplots(figsize=(9, max(5, 0.42 * len(labels))))
    ax.barh([i + h/2 for i in y], a_prop, height=h, color=COLUMN_A_SERIES, label=label_a)
    ax.barh([i - h/2 for i in y], b_prop, height=h, color=COLUMN_B_SERIES, label=label_b)
    ax.set_yticks(list(y)); ax.set_yticklabels(labels, fontsize=8); ax.invert_yaxis()
    ax.set_xlabel("% of cells"); ax.set_title(title)
    ax.legend(loc="lower right"); ax.grid(axis="x", alpha=0.3); fig.tight_layout()
    fig.savefig(out_path, dpi=150); plt.close(fig)


def sankey(cmp: Comparison, label_a, label_b, title, out_png, out_html,
           min_flow_fraction: float) -> None:
    pool = min_flow_fraction > 0
    thresh = min_flow_fraction * cmp.n_cells
    flows = defaultdict(int)
    for (a, b), c in cmp.pairs.items():
        flows[(a, b if (not pool or c >= thresh) else "(minor)")] += c
    left = sorted({a for (a, _) in flows}, key=lambda l: -sum(c for (x, _), c in flows.items() if x == l))
    right = sorted({b for (_, b) in flows}, key=lambda l: -sum(c for (_, x), c in flows.items() if x == l))
    left_nodes = [f"A: {l}" for l in left]
    right_nodes = [f"B: {r}" for r in right]
    idx = {n: i for i, n in enumerate(left_nodes + right_nodes)}

    def hexa(hh, alpha=0.4):
        hh = hh.lstrip("#")
        return f"rgba({int(hh[0:2],16)},{int(hh[2:4],16)},{int(hh[4:6],16)},{alpha})"

    node_colors = [cmp.a_color.get(l, "#888888") for l in left] + \
                  [("#bbbbbb" if r == "(minor)" else cmp.b_color.get(r, "#888888")) for r in right]
    src, tgt, val, lc = [], [], [], []
    for (a, b), c in sorted(flows.items(), key=lambda kv: -kv[1]):
        src.append(idx[f"A: {a}"]); tgt.append(idx[f"B: {b}"]); val.append(c)
        lc.append(hexa(cmp.a_color.get(a, "#888888")))
    note = "flows below the threshold pooled as (minor)" if pool else "all flows shown"
    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(label=left_nodes + right_nodes, color=node_colors,
                  pad=12, thickness=14, line=dict(width=0.4, color="#444")),
        link=dict(source=src, target=tgt, value=val, color=lc)))
    fig.update_layout(
        title_text=f"{title}  (A={label_a} → B={label_b}; n={cmp.n_cells:,}; {note})",
        font_size=11, width=1200, height=900 if pool else 1600)
    fig.write_image(out_png, scale=2)
    fig.write_html(out_html)


def main() -> None:
    p = argparse.ArgumentParser(description="Compare two cell-type columns per run.")
    p.add_argument("--bucket", required=True)
    p.add_argument("--prefix", required=True, help="S3 experiment prefix (contains run/slide dirs)")
    p.add_argument("--column-a", default="celltype_refit")
    p.add_argument("--column-b", default="celltype_norefit")
    p.add_argument("--label-a", default="Refit=TRUE")
    p.add_argument("--label-b", default="Refit=FALSE")
    p.add_argument("--min-flow-fraction", type=float, default=0.005,
                   help="Sankey flows below this fraction of cells are pooled into "
                        "(minor). Set 0 to draw every flow.")
    p.add_argument("--pool-all", action="store_true",
                   help="Pool all runs into a single comparison instead of per-run.")
    p.add_argument("--output-dir", default="compare_out")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    s3 = boto3.client("s3")
    groups = discover_slides(s3, args.bucket, args.prefix)
    if not groups:
        print(f"ERROR: no _metadata.csv found under s3://{args.bucket}/{args.prefix}",
              file=sys.stderr)
        sys.exit(1)

    if args.pool_all:
        merged = [(s, k) for slides in groups.values() for (s, k) in slides]
        groups = {"all": merged}

    for run, slides in groups.items():
        cmp = Comparison(args.column_a, args.column_b)
        for slide, key in sorted(slides):
            cmp.add_slide(s3, args.bucket, key)
        if cmp.n_cells == 0:
            print(f"[{run}] no comparable cells — skipped")
            continue
        name = short_run_name(run)
        print(f"[{name}] {len(slides)} slides, {cmp.n_cells:,} cells | "
              f"A types={len(cmp.a_counts)} B types={len(cmp.b_counts)} | "
              f"identical-label agreement={cmp.agreement():.1f}%")
        base = os.path.join(args.output_dir, name)
        abundance_chart(
            cmp, args.label_a, args.label_b,
            f"{name}: cell-type abundance — {args.label_a} vs {args.label_b}\n"
            f"({len(slides)} slides, n={cmp.n_cells:,} cells)",
            f"{base}_abundance.png")
        suffix = "" if args.min_flow_fraction > 0 else "_allflows"
        sankey(
            cmp, args.label_a, args.label_b, f"{name}",
            f"{base}_sankey{suffix}.png", f"{base}_sankey{suffix}.html",
            args.min_flow_fraction)
        print(f"  wrote {base}_abundance.png, {base}_sankey{suffix}.png/.html")


if __name__ == "__main__":
    main()
