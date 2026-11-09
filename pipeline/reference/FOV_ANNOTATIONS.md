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

## SME rulings, 2026-11-10

Answers to the five asks sent after the 2026-11-09 meeting. The FOV questions (the 7
collisions and the 46 unruled FOVs) are **not** here: she takes those in Napari on Friday
2026-11-13 from about noon.

**Adjacent soft tissue on `uwa7761eyes` — EXCLUDED.** Her ruling was conditional: *"if the
soft tissue also has very low quality, it can be excluded."* Stage 21 (job 43071516) measured
it as the worst region of any donor in the cohort — median 101 panel counts against 191 for
other donors' adjacent soft tissue, 0.53x, rank 12 of 12, with 49.4% of its cells under 100
counts. The condition is met, so the 26 FOVs (1, 69-75, 79, 80, 86, 87, 93, 94, 100-103,
108-115) now carry `exclude=1`. **183 FOVs are excluded cohort-wide**: 57 ciliary/cornea not
in any atlas, 100 brain, 26 soft tissue.

That donor keeps 22 optic nerve and 10 retina FOVs, both of which measured above the cohort
average.

**Retina is almost all CENTRAL (near macula) — `uwa7761eyes` is the exception, peripheral.**
This kills the standing plan to collapse the fovea-enriched HRCA subtypes on the grounds that
the cohort is peripheral-only; the opposite is true. It also means the one donor whose retina
we kept after excluding its brain is anatomically different from the other eleven, which is a
confound for any retina comparison and is not a quality problem.

**Retinal granularity goes to Kathryn Pepple.** Our neuropathologist asked for her perspective
rather than giving one, so every retina merge is parked until then.

**Brain granularity — settled.** `Excitatory neuron` is sufficient, *without* the deep vs
superficial labelling. The brain blocks are cortex and subcortical white matter only, so
**hippocampus and the basal-ganglia spiny neurons are excluded**.

**Lymphocytes — out.** Not on our detection argument, which was the weaker one: she reviewed
the H&E on all of these cases and saw no concerning inflammation.

## ⚠ Open items for the SME

**28 FOVs were claimed by two primary rows — 21 RESOLVED by the SME 2026-11-06, 7 still
open.** The first assignment used to win by default; the resolved rows now carry her call,
with `mixed_adjacent=1` reserved for FOVs she described as genuinely holding both tissues.

| slide | FOVs | conflict | SME's ruling | applied as |
|---|---|---|---|---|
| uwa7634eyes | 29 | Optic nerve vs Adjacent soft tissue | mostly soft tissue, some optic nerve | **Adjacent soft tissue**, mixed |
| uwa7634eyes | 36, 37 | Optic nerve vs Adjacent soft tissue | soft tissue | **Adjacent soft tissue** |
| uwa7634eyes | 38 | Optic nerve vs Adjacent soft tissue | mostly soft tissue, some optic nerve | **Adjacent soft tissue**, mixed |
| uwa7634eyes | 170 | White matter vs Gray matter | pretty much all white matter | White matter, no longer mixed |
| uwa7634eyes | 171 | White matter vs Gray matter | white + gray | White matter, mixed (unchanged) |
| uwa7634eyes | 172, 173 | White matter vs Gray matter | white + gray | White matter, now mixed |
| UWA7689eyes | 42 | Optic nerve vs Adjacent soft tissue | optic nerve | Optic nerve (unchanged) |
| UWA7689eyes | 43, 51 | Optic nerve vs Adjacent soft tissue | optic nerve + soft tissue | Optic nerve, mixed (unchanged) |
| UWA7689eyes | 48 | Optic nerve vs Adjacent soft tissue | almost entirely optic nerve, 3–4 other cells | Optic nerve, no longer mixed |
| UWA7689eyes | 49, 50 | Optic nerve vs Adjacent soft tissue | optic nerve | Optic nerve, no longer mixed |
| UWA7697eyes | 33 | Optic nerve vs Retina | right at the optic nerve head + retina | Optic nerve, now mixed |
| UWA7697eyes | 34 | Optic nerve vs Retina | almost all optic nerve, a few retinal cells | Optic nerve (unchanged) |
| UWA7697eyes | 80, 81 | Retina vs Adjacent soft tissue | **choroid** | **Choroid** — a new region value |

`Choroid` is the sixth region in the vocabulary (with Retina, Optic nerve, White matter,
Gray matter, Adjacent soft tissue). Two FOVs is too few to support a cohort-wide statement,
but it does give the "RPE scattered into choroid" question a denominator on this slide, and it
is where the SME placed the unnamed deep population `a_deep_retinal`.

**Still open — 7 FOVs the SME did not reach.** Deferred to the next round of InSituTree
review; until then the first assignment still wins and `mixed_adjacent` is unchanged.

| slide | FOVs | conflict |
|---|---|---|
| uwa7634eyes | 178, 179 | White matter vs Gray matter — the rest of the 170–179 run |
| uwa7634eyes | 135 | Retina vs Adjacent soft tissue |
| uwa7509eyes | 94 | Optic nerve vs Adjacent soft tissue |
| uwa7509eyes | 104 | Optic nerve vs Retina |
| uwa7509eyes | 217 | White matter vs Gray matter |
| eyes7597 | 67 | Optic nerve vs Retina |

**5 unparseable FOV tokens — RESOLVED by the SME 2026-09-01.** All five are now in the CSV
(2,528 → 2,549 rows). Two of our guesses were wrong, which is why they were left out:

| slide | row | token | our guess | SME's answer | applied as |
|---|---|---|---|---|---|
| UWA7740eyes | White matter | `179-174` | `174-179` ✗ | **`169-174`** — typo from 179 | 6 new rows, White matter |
| uwa7753eyes | Optic nerve | `118-120131` | `118-120, 131` ✓ | `118-120, 131` | 4 new rows, Optic nerve |
| UWA7689eyes | Retina | `86-102226` | `226` exceeds FOV count ✗ | fuzzy labels, "`227` and either `226` or `228`, keep with retina" | 1 new row: **226**, Retina |
| eyes7517 | (Optic nerve + adjacent) | `819` | `8, 19`? `81, 9`? | **`8, 19`** | flag only — see below |
| eyes7597 | White matter | `188-197?` | her own question mark | **`188-197` confirmed** | 10 new rows, White matter |

Two of those needed a judgement beyond the token itself:

- **`UWA7689eyes` 226 vs 228.** The SME could not read which of the pair it was, and the note
  above that `226` exceeds the FOV count was **wrong**. The slide carries exactly 227 FOVs
  (verified from `cell_ID` in the Napari metadata: 226 present, 228 absent), so the pair is
  226 + 227. `227` was already Retina; only `226` needed adding.
- **`eyes7517` 8 and 19 got no new region.** That row is parenthesised — `(Optic nerve +
  adjacent)` — so by this file's own convention it is a `mixed_adjacent` flag, not a competing
  assignment. Both FOVs already had primary rows (8 = Optic nerve, 19 = Adjacent soft tissue);
  those are unchanged and `mixed_adjacent` went 0 → 1.

**46 FOVs are in the data but NOT in this sheet — 21,331 cells, found 2026-11-08.** Stage 21
reported 20,191 such cells on the QC cohort and named the worst slides; the FOV identities
below come from the stitched Napari files. These are not excluded and not mislabelled: a FOV
absent from the sheet keeps whatever `Region` the AtoMx export gave it, which is the stale
annotation this file exists to correct. They are simply unruled.

| slide | FOVs | n | cells |
|---|---|---|---|
| `UWA7689eyes` | 86–102 | 17 | 10,512 |
| `uwa7753eyes` | 110–113, 116–117, 193–194, 196–197 | 10 | 6,419 |
| `UWA7740eyes` | 36–40, 43–46 | 9 | 1,774 |
| `UWA7697eyes` | 45–46, 69, 134 | 4 | 905 |
| `eyes7597` | 57, 206 | 2 | 1,042 |
| `uwa7611eyes` | 82 | 1 | 285 |
| `uwa7509eyes` | 73 | 1 | 171 |
| `UWA7575eyes` | 49 | 1 | 168 |
| `eyes7517` | 81 | 1 | 55 |

The nearest labelled FOV on each side, which is what an SME needs to place them in Napari:

| slide | FOVs | below | above |
|---|---|---|---|
| `UWA7689eyes` | 86–102 | 84–85 Optic nerve | 103–107 Adjacent soft tissue |
| `uwa7753eyes` | 110–113 | 109 Adjacent soft tissue | 114 Retina |
| `uwa7753eyes` | 116–117 | 115 Retina | 118 Optic nerve |
| `uwa7753eyes` | 193–194 | 192 White matter | 195 White matter |
| `uwa7753eyes` | 196–197 | 195 White matter | 198 White matter |
| `UWA7740eyes` | 36–40 | 35 Adjacent soft tissue | 41 Retina |
| `UWA7740eyes` | 43–46 | 42 Retina | 47 Retina |
| `UWA7697eyes` | 45–46 | 44 Retina | 47 Optic nerve |
| `UWA7697eyes` | 69 | 68 Optic nerve | 70 Adjacent soft tissue |
| `UWA7697eyes` | 134 | 133 Adjacent soft tissue | 135 Retina |
| `eyes7597` | 57 | 56 Retina | 58 Adjacent soft tissue |
| `eyes7597` | 206 | 205 White matter | none — the sheet ends at 205 |
| `uwa7611eyes` | 82 | 81 Optic nerve | 83 Retina |
| `uwa7509eyes` | 73 | 72 Optic nerve | 74 Adjacent soft tissue |
| `UWA7575eyes` | 49 | 48 Optic nerve | 50 Adjacent soft tissue |
| `eyes7517` | 81 | 80 Retina | 82 Optic nerve |

⚠ **`UWA7689eyes` 86–102 looks like a dropped half of a resolved token, not a true gap.**
That slide's Retina row carried `86-102226`. The 2026-09-01 resolution read the tail as
`226`/`227` and added **only** 226 — the `86-102` half was never written. Retina on that
slide is 60–80 and 226–227, with 86–102 sitting in an unlabelled hole between `Optic nerve`
84–85 and `Adjacent soft tissue` 103–107. It is 10,512 cells, half of this entire list, and
it should be confirmed rather than assumed — the token says retina but the neighbourhood
does not corroborate it, so this is not a gap to fill by inference. `uwa7753eyes` 193–197 is the same shape: a hole
inside `White matter` 192, 195, 198–209.

Gap-fill sanity: every added FOV was **absent** from the CSV, never an overwrite, and each new
range abuts a compatible neighbour — `eyes7597` 188-197 runs straight into 198-199 which were
already White matter; `UWA7740eyes` 169-174 sits in a Gray matter surround, which is what makes
a distinct White matter band there plausible.

⚠ **Stage 1 must be re-run for these to reach the data.** The 21 new rows change region
assignments for cells that currently carry the stale AtoMx `Region`, and §05-style per-region
figures computed before this date do not include them.

## Regenerating

Transcribed by unpacking the deck and parsing the DrawingML tables. The CSV is committed
rather than re-derived at runtime: the source is a slide deck, the ranges needed manual
adjudication, and a 2,528-row table is reviewable in a diff.
