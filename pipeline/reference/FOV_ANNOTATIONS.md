# Per-FOV tissue annotations (retina study)

`fov_annotations.csv` is the neuropathologist's FOV → tissue assignment for all 12 retina
slides, transcribed from `cosmx_annotations_8.1.25.pptx` (one table per donor). It is the
authoritative anatomical annotation; the `Region` column that arrives in the flat-file
metadata is an older, coarser AtoMx annotation and is **demonstrably stale on at least one
slide** — on `uwa7761eyes` it labels all 67 FOVs from 2–68 as `Retina`, where only 16–23,
34 and 35 are retina and 57 of them are ciliary body / cornea.

## Columns

| column | meaning |
|---|---|
| `slide_id` | pipeline slide id (matches `manifest_retina.csv`) |
| `fov` | per-slide integer FOV, as in `obs['fov']` |
| `region` | canonical tissue: `Retina`, `Optic nerve`, `Gray matter`, `White matter`, `Adjacent soft tissue`, or `EXCLUDE_not_in_atlas` |
| `mixed_adjacent` | 1 if the FOV straddles a boundary (the deck's parenthesised `(… + adjacent)` rows) |
| `exclude` | 1 if the SME ruled it out of analysis |

The parenthesised rows are recorded as a **flag, not a competing assignment** — they
overlap the primary rows by design, marking FOVs that contain two tissues.

`EXCLUDE_not_in_atlas` is 57 FOVs on `uwa7761eyes` (2–15, 24–33, 36–68): ciliary body and
cornea, which the SME excluded because no cell type in the HRCA + Monavarfeshani + Allen
reference covers them. Those 14,433 cells formed their own Leiden clusters (11 and 21)
rather than contaminating others, and are *deeper* than the real tissue (469 vs 302 median
counts) — epithelium is RNA-rich.

## Two vocabularies in the source deck

Ten donors are annotated only with the simple scheme above. **`UWA7753` (slide 26) and
`UWA7575` (slide 29) additionally carry a finer one** — `Retina + choroid`,
`Retina + choroid + sclera/soft tissue`, `Optic nerve + sheath`, `Soft tissue/sheath/sclera
only`, `Cortex`, `White matter + cortex`, `Cortex + leptomeninges`. Since it exists for only
2 of 12 donors it cannot be the cohort-wide `Region`, so this file uses the simple scheme
throughout. The detailed tables remain in the deck if a sub-analysis of those two wants them.

`uwa7761eyes` has two tables; this file uses the one marked **NEW** (slide 15), which covers
FOVs 1–215 and carries the exclusion. The older table (slide 16) covers only 40–215.

## ⚠ Open items for the SME

**28 FOVs are claimed by two primary rows.** Not silently resolved — the first assignment
wins in the CSV and every collision is listed below. Most are boundary cases of exactly the
kind the `(… + adjacent)` rows express, so the likely intent is "both".

| slide | FOVs | conflict |
|---|---|---|
| uwa7634eyes | 170–173, 178, 179 | White matter vs Gray matter |
| uwa7634eyes | 29, 36, 37, 38 | Optic nerve vs Adjacent soft tissue |
| uwa7634eyes | 135 | Retina vs Adjacent soft tissue |
| UWA7689eyes | 42, 43, 48, 49, 50, 51 | Optic nerve vs Adjacent soft tissue |
| UWA7697eyes | 33, 34 | Optic nerve vs Retina |
| UWA7697eyes | 80, 81 | Retina vs Adjacent soft tissue |
| uwa7509eyes | 94 | Optic nerve vs Adjacent soft tissue |
| uwa7509eyes | 104 | Optic nerve vs Retina |
| uwa7509eyes | 217 | White matter vs Gray matter |
| eyes7597 | 67 | Optic nerve vs Retina |

**5 unparseable FOV tokens** — transcription typos, left out of the CSV rather than guessed:

| slide | row | token | likely intent |
|---|---|---|---|
| UWA7740eyes | White matter | `179-174` | reversed range → `174-179` |
| uwa7753eyes | Optic nerve | `118-120131` | missing comma → `118-120, 131` |
| UWA7689eyes | Retina | `86-102226` | missing comma, but `226` exceeds the FOV count |
| eyes7517 | (Optic nerve + adjacent) | `819` | ambiguous — `8, 19`? `81, 9`? |
| eyes7597 | White matter | `188-197?` | the SME's own question mark |

## Regenerating

Transcribed by unpacking the deck and parsing the DrawingML tables. The CSV is committed
rather than re-derived at runtime: the source is a slide deck, the ranges needed manual
adjudication, and a 2,528-row table is reviewable in a diff.
