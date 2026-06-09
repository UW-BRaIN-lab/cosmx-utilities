#!/usr/bin/env Rscript
# Marker-gene heatmap renderer (the render half; compute half is
# pipeline/python/marker_pseudobulk.py). Adapted from the lab's InSituType
# marker_heatmap.R, keyed on Stage 3c Leiden clusters.
#
# Reads the small CSVs written by marker_pseudobulk.py and draws a ComplexHeatmap
# of z-scored pseudobulk expression: top markers (rows, split by the cluster they
# mark) x cluster|Region groups (columns, split by cluster, annotated by Region).
# Needs only ComplexHeatmap + circlize — no Seurat/AnnData — so it runs locally on
# a laptop that already has ComplexHeatmap, or in a small R container.
#
# Usage:
#   Rscript pipeline/R/marker_heatmap.R <input_dir> [output_dir]
# where <input_dir> holds marker_heatmap_zmatrix.csv + top_markers_per_cluster.csv.

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

TOP_TITLE <- "Top markers per Leiden cluster, pseudobulk by region (z-scored log-norm expression)"
# Tumor -> edge -> normal; only those present are used.
REGION_ORDER  <- c("Tumor bulk", "Infiltrating edge", "Contralateral uninvolved")
REGION_COLORS <- c("Tumor bulk" = "#E7298A",
                   "Infiltrating edge" = "#009E73",
                   "Contralateral uninvolved" = "#7570B3")

# Font sizes (pt) and per-row/col canvas size (inches). With many Leiden clusters
# the heatmap is large; bump these up if labels render too small (or render fewer
# markers with a smaller TOP_N in marker_pseudobulk.py).
ROW_FONTSIZE   <- 11   # marker gene names
COL_FONTSIZE   <- 12   # column (Region) labels
SPLIT_FONTSIZE <- 13   # cluster split titles
TITLE_FONTSIZE <- 16   # overall title
PER_ROW_IN     <- 0.28
PER_COL_IN     <- 0.55

args   <- commandArgs(trailingOnly = TRUE)
indir  <- if (length(args) >= 1) args[[1]] else "."
outdir <- if (length(args) >= 2) args[[2]] else indir

zmat_path <- file.path(indir, "marker_heatmap_zmatrix.csv")
mark_path <- file.path(indir, "top_markers_per_cluster.csv")
stopifnot("missing marker_heatmap_zmatrix.csv" = file.exists(zmat_path),
          "missing top_markers_per_cluster.csv" = file.exists(mark_path))

message("Reading ", zmat_path)
# check.names = FALSE keeps the "<cluster> | <Region>" column names intact.
pb_z <- as.matrix(read.csv(zmat_path, row.names = 1, check.names = FALSE))
markers <- read.csv(mark_path, stringsAsFactors = FALSE)
gene_to_cluster <- setNames(as.character(markers$cluster), markers$gene)

# --- derive column (cluster, Region) metadata from the column names -----------
col_cluster <- sub(" \\| .*", "", colnames(pb_z))
col_region  <- sub(".* \\| ", "", colnames(pb_z))

# Cluster ordering: numeric when leiden-like, else lexical.
uniq_cl <- unique(c(col_cluster, unname(gene_to_cluster)))
cl_levels <- tryCatch(as.character(sort(as.integer(uniq_cl))),
                      warning = function(w) sort(uniq_cl))
region_levels <- intersect(REGION_ORDER, unique(col_region))

col_cluster <- factor(col_cluster, levels = cl_levels)
col_region  <- factor(col_region, levels = region_levels)
row_cluster <- factor(gene_to_cluster[rownames(pb_z)], levels = cl_levels)

# --- palettes -----------------------------------------------------------------
cl_palette <- setNames(
  colorRampPalette(c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                     "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"))(length(cl_levels)),
  cl_levels)
region_colors <- REGION_COLORS[region_levels]
heatmap_col <- colorRamp2(c(-3, 0, 3), c("navy", "white", "firebrick"))

top_anno <- HeatmapAnnotation(
  Region = col_region,
  col = list(Region = region_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(Region = list(title = "Region", nrow = 1))
)
left_anno <- rowAnnotation(
  cluster = row_cluster,
  col = list(cluster = cl_palette),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  width = unit(3, "mm")
)

ht <- Heatmap(
  pb_z,
  name              = "z-score",
  col               = heatmap_col,
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  show_row_names    = TRUE,
  show_column_names = TRUE,
  column_labels     = as.character(col_region),
  column_names_rot  = 90,
  column_names_side = "top",
  column_names_gp   = gpar(fontsize = COL_FONTSIZE),
  column_split      = col_cluster,
  cluster_column_slices = FALSE,
  column_title_side = "bottom",
  column_title_rot  = 0,
  column_title_gp   = gpar(fontsize = SPLIT_FONTSIZE, fontface = "bold"),
  column_gap        = unit(1, "mm"),
  row_split         = row_cluster,
  cluster_row_slices = FALSE,
  row_title_side    = "left",
  row_title_rot     = 0,
  row_title_gp      = gpar(fontsize = SPLIT_FONTSIZE, fontface = "bold"),
  row_gap           = unit(0.7, "mm"),
  row_names_gp      = gpar(fontsize = ROW_FONTSIZE),
  row_names_side    = "right",
  top_annotation    = top_anno,
  left_annotation   = left_anno,
  border            = TRUE,
  heatmap_legend_param = list(title = "z-score", direction = "horizontal",
                              title_position = "topcenter")
)

n_genes <- nrow(pb_z)
height  <- max(10, n_genes * PER_ROW_IN + 5)
width   <- max(14, ncol(pb_z) * PER_COL_IN + 4)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(outdir, "marker_heatmap.pdf")
pdf(out_pdf, width = width, height = height)
draw(ht,
     heatmap_legend_side    = "bottom",
     annotation_legend_side = "bottom",
     merge_legend           = TRUE,
     column_title           = TOP_TITLE,
     column_title_gp        = gpar(fontsize = TITLE_FONTSIZE, fontface = "bold"))
invisible(dev.off())
message("Wrote: ", out_pdf)
