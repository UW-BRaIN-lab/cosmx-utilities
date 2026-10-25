# Stage 4 de-novo cluster annotations

InSituType's de-novo clusters (cells that don't match the GBmap reference) come out as
single letters `a`–`l`. These tables assign each letter a biological identity from that
run's de-novo marker heatmap (`marker_pseudobulk.py --clusters a,…,l`, rendered with
`R/marker_heatmap.R`). Apply them with `pipeline/python/annotate_denovo.py`, which
rewrites `obs.cell_type` so a letter becomes e.g. `a - MES/AC-like tumor` — **keeping the
letter prefix** so it stays clear the cluster was de-novo discovered, not reference-named.
Named GBmap reference types are left unchanged.

The letter labels are **arbitrary and run-specific** — `a` in one run is unrelated to `a`
in another. Each table stands alone.

Columns: `denovo_label`, `annotation` (what gets written to `cell_type`), `top_markers`,
`quality` (confidence note — clean / moderate / low-signal / low-confidence / mixed).

## The four Wenyu InSituType runs

All on the 2.33M-cell min50 cohort, GBmap reference, non-preemptible gpu-l40s slice.

| table | reference | profile update | notes |
|---|---|---|---|
| `stage4.csv` | Core, level 4 | rescale only | ~90% de-novo; baseline |
| `stage4_refit.csv` | Core, level 4 | rescale + refit | refit inflated rare immune types (Reg_T/DC3) |
| `stage4_ext_l3.csv` | Extended, level 3 | rescale + refit | refit inflated rare types (B_cell/RG) |
| `stage4_extl3_rescale.csv` | Extended, level 3 | rescale only | **keeper** — no fabricated labels; ~36% in a tumor-restricted stressed/low-signal cluster |

See memory `project_insitutype_letter_clusters` for the full run history and the
rescale-vs-refit / Core-vs-Extended findings.

## SORL1 3D cohort (Maddie), Brain_6k reference

`stage4_k30.csv`: 30 de novo clusters (`a`-`z`, `aa`-`ad`) from the K=30 run
(`update_reference=TRUE, rescale=TRUE, refit=FALSE`, 3,348-gene pruned panel, 793,380 cells), read from
marker lists, the Leiden cross-tab, per-cluster donor spread and carrier enrichment. Letters are specific to
that run. The 13 named Brain_6k types are unchanged. Confidence is in `quality`; several rows are
donor-specific stress states or unresolved low-information clusters and should not be read as cell types.

`stage4_k21.csv`: the 21 de novo clusters of the K=21 run (the recommended working resolution), same method.
Marker lists for every type in a fit come from `python/top_markers_from_profiles.py`.

`stage4_k21_seaad_markers.csv`: SEA-AD MTG marker evidence for the same K=21 run, one row per type (21 de novo +
13 named Brain_6k controls), beside the current annotation. Built by `python/marker_overlap_annotation_table.py`
from `marker_overlap_vs_reference.py`. Evidence to check in Napari, not a rename: it has no `keep`/`new_name`
columns. `support` = strong (concordance >= 0.95 and overlap p < 1e-3) / partial (p < 0.05) / none. This copy was
built from the top three matches per cluster printed by the marker-overlap run; regenerate from the full
`marker_overlap.csv` for exact `also_close` and `top_overlap_*`.
