#!/usr/bin/env python3
"""Compare two cell-type columns in per-slide metadata files.

Reads per-slide metadata CSVs from S3 and, for each group of slides, compares two
of their columns per cell:

  - a grouped abundance bar chart (proportion of cells per type, column A vs B)
  - a Sankey tracing each cell from its column-A type to its column-B type

Two metadata layouts are supported (pick with --metadata-suffix):

  - stitched `_metadata.csv` from generate-slide-metadata.py — friendly column
    names (celltype_refit / celltype_norefit) each paired with a <column>_color;
    laid out as `.../<run>/<slide>/_metadata.csv` (the default).
  - raw AtoMx flatFiles `<slide>_metadata_file.csv.gz` — the original AtoMx column
    names (e.g. RNA_RNA_Cell.Typing.InSituType.1_1_clusters) and no color sidecars
    (colors fall back to a stable per-label hash); gzipped and laid out as
    `.../flatFiles/<slide>/<slide>_metadata_file.csv.gz`.

Slides are grouped by their run folder (the directory just above each slide) by
default; --per-slide compares each slide on its own (the right granularity when the
columns are de-novo cluster labels, which InSituType assigns independently per
slide), and --pool-all pools everything into one comparison.

Usage:
    # Stitched metadata, per-run (default), Refit=TRUE vs Refit=FALSE:
    uv run --with boto3 --with matplotlib --with plotly --with kaleido \
        python scripts/compare-cell-typing.py \
        --bucket keene-cosmx-data \
        --prefix napari-stitched/CosMx-Victoria/<experiment> \
        --output-dir compare_out

    # Raw AtoMx flatFiles, per-slide, with the original InSituType column names:
    uv run --with boto3 --with matplotlib --with plotly --with kaleido \
        python scripts/compare-cell-typing.py \
        --bucket keene-cosmx-data \
        --prefix CosMx-Victoria/<experiment>/flatFiles/<inner>/flatFiles \
        --metadata-suffix _metadata_file.csv.gz --per-slide \
        --column-a RNA_RNA_Cell.Typing.InSituType.1_1_clusters \
        --column-b RNA_RNA_Cell.Typing.InSituType.No.Refit_1_clusters \
        --label-a "Refit=TRUE" --label-b "Refit=FALSE" \
        --output-dir compare_out
"""

import argparse
import csv
import gzip
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


def safe_name(name: str) -> str:
    """Filesystem-safe token for an output basename."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") or "group"


def discover_slides(s3, bucket: str, prefix: str, suffix: str) -> list:
    """Find every per-slide metadata file under the prefix.

    Matches keys ending in `suffix` (e.g. `_metadata.csv` for stitched output, or
    `_metadata_file.csv.gz` for raw AtoMx flatFiles). Returns a list of
    (run, slide, key): the slide is the directory holding the file and the run is
    the directory immediately above it; files sitting directly under the prefix are
    grouped under the prefix's last path part.
    """
    records = []
    paginator = s3.get_paginator("list_objects_v2")
    base = prefix.rstrip("/")
    for page in paginator.paginate(Bucket=bucket, Prefix=base + "/"):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if not key.endswith(suffix):
                continue
            parts = key.split("/")
            slide = parts[-2] if len(parts) >= 2 else base.rsplit("/", 1)[-1]
            run = parts[-3] if len(parts) >= 3 else base.rsplit("/", 1)[-1]
            records.append((run, slide, key))
    return records


def group_records(records: list, mode: str) -> dict:
    """Bucket (run, slide, key) records into named comparison groups.

    mode="run"   -> one group per run folder (readable run label)
    mode="slide" -> one group per slide (raw slide name; de-novo labels are
                    per-slide, so pooling slides would mix unrelated label spaces)
    mode="all"   -> a single pooled group
    """
    groups = defaultdict(list)
    for run, slide, key in records:
        if mode == "all":
            name = "all"
        elif mode == "slide":
            name = slide
        else:
            name = short_run_name(run)
        groups[name].append((slide, key))
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
        raw = s3.get_object(Bucket=bucket, Key=key)["Body"].read()
        if key.endswith(".gz"):
            raw = gzip.decompress(raw)
        reader = csv.DictReader(io.StringIO(raw.decode("utf-8")))
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
    p = argparse.ArgumentParser(description="Compare two cell-type columns per run/slide.")
    p.add_argument("--bucket", required=True)
    p.add_argument("--prefix", required=True, help="S3 experiment prefix (contains run/slide dirs)")
    p.add_argument("--metadata-suffix", default="_metadata.csv",
                   help="Match metadata objects ending in this suffix. Use "
                        "_metadata_file.csv.gz for raw AtoMx flatFiles. "
                        "Gzipped (.gz) files are decompressed automatically.")
    p.add_argument("--column-a", default="celltype_refit")
    p.add_argument("--column-b", default="celltype_norefit")
    p.add_argument("--label-a", default="Refit=TRUE")
    p.add_argument("--label-b", default="Refit=FALSE")
    p.add_argument("--min-flow-fraction", type=float, default=0.005,
                   help="Sankey flows below this fraction of cells are pooled into "
                        "(minor). Set 0 to draw every flow.")
    grouping = p.add_mutually_exclusive_group()
    grouping.add_argument("--per-slide", action="store_true",
                          help="Compare each slide on its own (right granularity when "
                               "the columns are de-novo cluster labels).")
    grouping.add_argument("--pool-all", action="store_true",
                          help="Pool all slides into a single comparison.")
    p.add_argument("--output-dir", default="compare_out")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    s3 = boto3.client("s3")
    records = discover_slides(s3, args.bucket, args.prefix, args.metadata_suffix)
    if not records:
        print(f"ERROR: no *{args.metadata_suffix} found under "
              f"s3://{args.bucket}/{args.prefix}", file=sys.stderr)
        sys.exit(1)

    mode = "slide" if args.per_slide else "all" if args.pool_all else "run"
    groups = group_records(records, mode)

    for group, slides in sorted(groups.items()):
        cmp = Comparison(args.column_a, args.column_b)
        for slide, key in sorted(slides):
            cmp.add_slide(s3, args.bucket, key)
        if cmp.n_cells == 0:
            print(f"[{group}] no comparable cells — skipped")
            continue
        print(f"[{group}] {len(slides)} slide(s), {cmp.n_cells:,} cells | "
              f"A types={len(cmp.a_counts)} B types={len(cmp.b_counts)} | "
              f"identical-label agreement={cmp.agreement():.1f}%")
        base = os.path.join(args.output_dir, safe_name(group))
        abundance_chart(
            cmp, args.label_a, args.label_b,
            f"{group}: cell-type abundance — {args.label_a} vs {args.label_b}\n"
            f"({len(slides)} slide(s), n={cmp.n_cells:,} cells)",
            f"{base}_abundance.png")
        suffix = "" if args.min_flow_fraction > 0 else "_allflows"
        sankey(
            cmp, args.label_a, args.label_b, f"{group}",
            f"{base}_sankey{suffix}.png", f"{base}_sankey{suffix}.html",
            args.min_flow_fraction)
        print(f"  wrote {base}_abundance.png, {base}_sankey{suffix}.png/.html")


if __name__ == "__main__":
    main()
