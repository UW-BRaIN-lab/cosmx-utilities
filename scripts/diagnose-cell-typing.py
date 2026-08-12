#!/usr/bin/env python3
"""Diagnose a CosMx InSituType run: are de-novo clusters real or low-signal?

Joins per-cell cluster labels (from a stitched _metadata.csv column) to the raw
exprMat flat files, then per run/study reports:

  * de-novo fraction and per-cluster stats (n, mean total counts, "flatness" =
    share of expression in the top-10 genes, and whether it looks like a
    low-signal catch-all vs a coherent, marker-defined population)
  * how many of the panel genes are ABOVE BACKGROUND in the data — the key
    number for whether gene pruning would help InSituType on a 6k panel
    (InSituType FAQ: on 6k panels, subset to genes above background / informative
    in the reference, ~3-5k genes)
  * a marker heatmap (top enriched genes x clusters)

Cluster labels come from the stitched _metadata.csv (napari-stitched), so it uses
the same run/slide names as the raw flat files. exprMat cells are matched to
metadata cells by (fov, cell index) parsed from the c_<slide>_<fov>_<cell> id.

Usage:
    uv run --with pandas --with numpy --with matplotlib python \
        scripts/diagnose-cell-typing.py \
        --bucket keene-cosmx-data \
        --cluster-prefix napari-stitched/CosMx-retina/CosMx-retina-brain-segmentation-test-4.1.26 \
        --expr-prefix   CosMx-retina/CosMx-retina-brain-segmentation-test-4.1.26 \
        --cluster-column celltype_norefit \
        --output-dir typing_diag
"""
import argparse
import io
import os
import re

import boto3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DENOVO = re.compile(r"^[a-z](_[0-9]+)?$")           # de-novo cluster id: a, b, d_1, d_2
CONTROL = re.compile(r"^(SystemControl|Negative|NegPrb|FalseCode)", re.I)


def s3_get(s3, bucket, key):
    return s3.get_object(Bucket=bucket, Key=key)["Body"].read()


def list_dirs(s3, bucket, prefix):
    resp = s3.list_objects_v2(Bucket=bucket, Prefix=prefix.rstrip("/") + "/", Delimiter="/")
    return [p["Prefix"].rstrip("/").rsplit("/", 1)[-1] for p in resp.get("CommonPrefixes", [])]


def load_clusters(s3, bucket, cluster_prefix, run, slide, column):
    """cell (fov, cell) -> cluster, from the stitched _metadata.csv."""
    key = f"{cluster_prefix.rstrip('/')}/{run}/{slide}/_metadata.csv"
    df = pd.read_csv(io.BytesIO(s3_get(s3, bucket, key)), usecols=["cell_ID", column])
    df = df.rename(columns={column: "cluster"}).dropna(subset=["cluster"])
    ids = df["cell_ID"].str.split("_")
    df["fov"] = ids.str[-2].astype(np.int32)
    df["cell_ID_int"] = ids.str[-1].astype(np.int32)
    return df[["fov", "cell_ID_int", "cluster"]].rename(columns={"cell_ID_int": "cell_ID"})


def accumulate(s3, bucket, expr_prefix, run, slide, clusters, sums, counts):
    key = f"{expr_prefix.rstrip('/')}/{run}/flatFiles/{slide}/{slide}_exprMat_file.csv.gz"
    reader = pd.read_csv(io.BytesIO(s3_get(s3, bucket, key)), compression="gzip", chunksize=20000)
    gene_cols = None
    for chunk in reader:
        if gene_cols is None:
            gene_cols = [c for c in chunk.columns if c not in ("fov", "cell_ID")]
        chunk["fov"] = chunk["fov"].astype(np.int32)
        chunk["cell_ID"] = chunk["cell_ID"].astype(np.int32)
        m = chunk.merge(clusters, on=["fov", "cell_ID"], how="inner")
        if m.empty:
            continue
        g = m.groupby("cluster")[gene_cols].sum()
        n = m.groupby("cluster").size()
        for cl in g.index:
            sums[cl] = sums.get(cl, 0) + g.loc[cl]
            counts[cl] = counts.get(cl, 0) + int(n.loc[cl])
    return gene_cols


def gene_informativeness(overall, means, counts, min_cluster_n=500):
    """Estimate the InSituType gene-pruning target (6k-panel FAQ).

    Two views:
      * above background: real genes whose pooled mean clears the negative-probe
        background (a weak floor — usually most genes pass).
      * informative: real genes that are a MARKER for at least one sizeable
        population (max cross-cluster enrichment >= 2x / 3x). This is the number
        that matters for typing — the rest add per-cell noise on a 6k panel.
    """
    idx = overall.index
    ctrl = [g for g in idx if CONTROL.match(g)]
    real = [g for g in idx if not CONTROL.match(g)]
    if not ctrl:
        return None
    bg = float(overall[ctrl].mean())
    big = [cl for cl in counts if counts[cl] >= min_cluster_n]
    # per-gene max enrichment across sizeable clusters
    mat = np.vstack([((means[cl] + 0.05) / (overall + 0.05)).values for cl in big])
    max_enr = pd.Series(mat.max(axis=0), index=idx)[real]
    real_mu = overall[real]
    return dict(background=bg, n_real=len(real), n_ctrl=len(ctrl), n_big=len(big),
                above_2x_bg=int((real_mu > 2 * bg).sum()),
                informative_2x=int((max_enr >= 2).sum()),
                informative_3x=int((max_enr >= 3).sum()))


def heatmap(means, counts, overall, run, out_path, top_clusters=12, per_cluster=3):
    order = sorted(counts, key=lambda c: -counts[c])[:top_clusters]
    genes = []
    for cl in order:
        enr = (means[cl] + 0.05) / (overall + 0.05)
        for g in enr[[x for x in enr.index if not CONTROL.match(x)]].sort_values(ascending=False).head(per_cluster).index:
            if g not in genes:
                genes.append(g)
    M = np.array([[means[cl][g] for cl in order] for g in genes])
    Z = (M - M.mean(axis=1, keepdims=True)) / (M.std(axis=1, keepdims=True) + 1e-9)
    fig, ax = plt.subplots(figsize=(max(8, 0.7 * len(order)), max(6, 0.28 * len(genes))))
    im = ax.imshow(Z, aspect="auto", cmap="magma")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{c}\nn={counts[c]:,}\n{means[c].sum():.0f}ct" for c in order],
                       fontsize=7, rotation=90)
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=6)
    ax.set_title(f"{run}: top marker expression (row z-score) per cluster\n"
                 f"(de-novo ids a/b/d_1/… ; 'ct' = mean total counts)")
    fig.colorbar(im, ax=ax, label="z-score", shrink=0.6)
    fig.tight_layout(); fig.savefig(out_path, dpi=150); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="Diagnose de-novo / low-signal clusters in a typing run.")
    p.add_argument("--bucket", required=True)
    p.add_argument("--cluster-prefix", required=True, help="stitched napari prefix (has run/slide/_metadata.csv)")
    p.add_argument("--expr-prefix", required=True, help="raw prefix (has run/flatFiles/slide/exprMat)")
    p.add_argument("--cluster-column", default="celltype_norefit")
    p.add_argument("--sample-slides", type=int, default=3)
    p.add_argument("--output-dir", default="typing_diag")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    s3 = boto3.client("s3")
    runs = list_dirs(s3, args.bucket, args.cluster_prefix)
    if not runs:
        raise SystemExit(f"No runs under s3://{args.bucket}/{args.cluster_prefix}")

    for run in runs:
        slides = sorted(list_dirs(s3, args.bucket, f"{args.cluster_prefix.rstrip('/')}/{run}"))
        sample = slides[:args.sample_slides]
        print(f"\n########## {run}  (sample: {', '.join(sample)}) ##########")
        sums, counts = {}, {}
        for slide in sample:
            cl = load_clusters(s3, args.bucket, args.cluster_prefix, run, slide, args.cluster_column)
            accumulate(s3, args.bucket, args.expr_prefix, run, slide, cl, sums, counts)
            print(f"  {slide}: {len(cl):,} cells")
        if not counts:
            print("  no cells matched — skipping"); continue

        total_n = sum(counts.values())
        overall = sum(sums.values()) / total_n
        means = {cl: sums[cl] / counts[cl] for cl in counts}
        clusters = sorted(counts, key=lambda c: -counts[c])
        denovo_n = sum(counts[c] for c in clusters if DENOVO.match(c))
        print(f"  de-novo fraction: {100*denovo_n/total_n:.1f}%  (n={total_n:,})")

        gb = gene_informativeness(overall, means, counts)
        if gb:
            print(f"  GENE PRUNING (bg={gb['background']:.3f} from {gb['n_ctrl']} controls, "
                  f"{gb['n_big']} sizeable clusters):")
            print(f"    of {gb['n_real']} real genes: {gb['above_2x_bg']} clear 2x background, "
                  f"but only {gb['informative_2x']} are a marker (>=2x in some population) "
                  f"and {gb['informative_3x']} at >=3x")
            print(f"    → InSituType pruning target ~= {gb['informative_2x']} informative genes")

        print(f"  {'cluster':<28}{'n':>10}{'counts':>9}{'top10%':>8}  top enriched genes")
        for cl in clusters[:14]:
            mu = means[cl]
            enr = (mu + 0.05) / (overall + 0.05)
            real = [g for g in enr.index if not CONTROL.match(g)]
            top = enr[real].sort_values(ascending=False).head(5)
            tot = float(mu.sum())
            top10 = float(mu.sort_values(ascending=False).head(10).sum()) / tot * 100 if tot else 0
            kind = "LOW-SIGNAL/flat" if (DENOVO.match(cl) and top.iloc[0] < 3 and tot < 1.2 * overall.sum()) else ""
            markers = ", ".join(f"{g}={enr[g]:.1f}x" for g in top.index)
            print(f"  {cl:<28}{counts[cl]:>10,}{tot:>9.0f}{top10:>7.0f}%  {markers} {kind}")

        heatmap(means, counts, overall, run,
                os.path.join(args.output_dir, f"{run[:40]}_marker_heatmap.png"))
        print(f"  wrote {run[:40]}_marker_heatmap.png")


if __name__ == "__main__":
    main()
