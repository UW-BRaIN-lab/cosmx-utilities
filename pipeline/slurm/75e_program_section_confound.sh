#!/bin/bash
# Is a de-novo letter's gene programme a cell state, or a slide-level handling artefact?
#
# 75d showed de-novo `t` is separated from every natively-called GBmap class by a heat-shock
# programme (HSPA1A/HSPA1B/HSPB1/DNAJB1/HSP90AA1/HSPH1, mean z +0.3..+1.6 in every t-> column
# against negative in every native). Heat-shock is exactly what warm ischaemia and slow fixation
# produce, so it has to be ruled out before `t` is named a tumour state.
#
# The unit is the TISSUE SECTION, not the slide. Two sections are mounted per CosMx slide and
# need not be from the same donor, so slide_id pools two pieces with different blocks and
# fixation histories — and if the letter's cells sit mostly in one piece, a "within-slide" gap is
# partly a between-piece comparison. On a fixture where one section per slide was hot and the
# letter was 80% concentrated there, slide-level grouping called it a CELL STATE (100% positive
# gaps, median +1.37) while section-level correctly called it handling. Sections are derived from
# the stage-1 cell ids ("<slide>_F<fov>_C<cell>") as contiguous FOV runs, the same notion
# tissue_section_gap.py uses. That check EARNED ITS KEEP on run 39972319: it found 56 of 57
# slides with a single run, because this cohort's FOV numbering is continuous 1-200 across both
# pieces — so the run heuristic collapses to the slide here and must not be used.
#
# The unit is therefore the FOV, which lies entirely within one piece by construction and so
# controls for the section and for position inside it. See python/diagnose_program_by_section.py.
#
# Submit (defaults to t + the heat-shock set):
#   sbatch pipeline/slurm/75e_program_slide_confound.sh
# Another letter or another programme:
#   LETTER=b GENES=CES3,A1BG,SNRPA1,GPX1 sbatch pipeline/slurm/75e_program_slide_confound.sh
# The amplicon question, against the malignant reference types rather than everyone:
#   LETTER=b GENES=CDK4,MDM2,NUP107,OS9,SLC35E3 COMPARE_TO=malignant NULL_SETS=200 \
#       TAG=amp12qnull sbatch pipeline/slurm/75e_program_section_confound.sh
#   LETTER=b GENES=EGFR,CDK4,MDM2,NUP107,OS9,SLC35E3 COMPARE_TO=malignant \
#       sbatch pipeline/slurm/75e_program_section_confound.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR  Kopah sub-dir with anchor/anchor_typing.h5 (default stage4_anchor_pruned).
#   INPUT_DIR   Kopah sub-dir with anchor/anchor_input.h5 (default stage4_anchor).
#   UNIT        fov (default) or fov-run. fov-run is kept for cohorts whose FOV numbering does
#               break between the two pieces; check its runs-per-slide line says 2 first.
#   MIN_CELLS   minimum cells of BOTH groups per unit (default 25).
#   LETTER      de-novo letter under test (default t).
#   GENES       comma-separated programme genes (default: the heat-shock set from 75d).
#   TAG         label folded into every output filename. Defaults to the Slurm job id, because
#               without it two runs that differ only in GENES write to the SAME Kopah key and the
#               second silently destroys the first — which is exactly what the control run did to
#               the amplicon run's CSVs. Pass something readable, e.g. TAG=amplicon.
#   NULL_SETS   draw this many MATCHED-NULL gene sets and print where the observed deficit falls
#               among them (default 0 = off). REQUIRED before any deficit goes in a figure or to
#               the PI: the score is library-size normalised, but normalisation cannot recover an
#               undetected transcript, so a shallower population shows a deficit on lightly
#               expressed genes of ANY biology. The de-novo letters run as low as 0.37x their
#               comparators' depth (75g), so this is the live confound. 200 sets resolves p to
#               about 0.005 and costs roughly 10-15 minutes.
#   NULL_SEED   seed for that draw (default 0), so a reported p is reproducible.
#               With NULL_SETS set, the run also writes <stem>_program_by_unit_null.csv: the mean
#               per-unit gap across the null sets, which is the properly matched per-unit CONTROL
#               for program_gap_by_group.py (donor-level excess). The old hand-picked 30-gene
#               control was matched to the COMBINED six-gene set and is the WRONG control once
#               the loci are split -- especially for the one-gene EGFR run.
#   NULL_NEIGHBOURS  candidate pool per programme gene (default 100).
#   SPANNING_SETS    run the BACKGROUND-DILUTION correction with this many spanning sets
#               (default 0 = off). REQUIRED FOR ANY AMPLICON CLAIM. 75q put CDK4/MDM2/OS9 in the
#               panel's bottom 0.25% for how their expression tracks the per-cell negprobe load,
#               so a background-heavier population shows a deficit in them with no copy number
#               involved. The matched null CANNOT control this (too few genes are that extreme to
#               match against); this lays null sets ACROSS the affinity range instead, fits
#               gap ~ affinity + abundance, and reports the residual. 300 is plenty. Check the
#               `bracketed` column is True -- False means the prediction is an extrapolation.
#               Costs roughly as much again as the same number of NULL_SETS, so for a long run
#               either set NULL_SETS=0 (the spanning fit controls abundance too) or move off ckpt
#               with:  --account=carrolllab --partition=cpu-g2-mem2x --qos=''
#               which is non-preemptible.
#   COMPARE_TO  restrict the comparison group to these cell types instead of every other cell.
#               'malignant' expands to the nine malignant Core-L4 columns. REQUIRED whenever the
#               reference itself stratifies on the programme: gbmap_level4_panel.csv puts all
#               nine malignant columns at ranks 1-9 of 54 on the amplicon genes (mean 1.93 vs
#               0.72), so an all-other-cells pool is mostly non-malignant and dilutes a
#               malignant-vs-malignant deficit away. Run 40050341/40050443 hit exactly that.

#SBATCH --job-name=cosmx-program-confound
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
# Only the programme genes' rows and per-cell column sums are materialised, so this never
# transposes or densifies the counts — far lighter than 75d.
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=pipeline/logs/program_confound_%j.out
#SBATCH --error=pipeline/logs/program_confound_%j.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/lmod.sh 2>/dev/null \
        || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module load apptainer
export PATH="${HOME}/bin:${PATH}"

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PIPELINE_DIR="${SLURM_SUBMIT_DIR}/pipeline"
else
    PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
fi

set -a
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/.env"
set +a

STAGE4="${STAGE4_DIR:-stage4_anchor_pruned}"
INPUT="${INPUT_DIR:-stage4_anchor}"
LETTER="${LETTER:-t}"
# Never let two gene sets collide on one key; the job id is unique even when TAG is forgotten.
TAG="${TAG:-${SLURM_JOB_ID:-local}}"
STEM="${LETTER}_${TAG}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_program_confound_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

# Echo the RESOLVED plan before any work: a stale checkout silently ignores new env knobs, and
# that has already cost a run that re-did the defaults while looking like it honoured the request.
echo "PLAN  letter=${LETTER}  tag=${TAG}  unit=${UNIT:-fov}  min_cells=${MIN_CELLS:-25}"
echo "      genes=${GENES:-<default heat-shock set>}"
echo "      compare_to=${COMPARE_TO:-<every other cell>}"
echo "      null_sets=${NULL_SETS:-0}  spanning_sets=${SPANNING_SETS:-0}"
if [[ "${SPANNING_SETS:-0}" -eq 0 ]]; then
    echo "      NOTE: no dilution correction. Any AMPLICON claim needs SPANNING_SETS>0."
fi

echo "Staging anchor counts + labels from Kopah (sections come from the cell ids)..."
s5cmd cp "${BASE}/${INPUT}/anchor/anchor_input.h5" "$WORK/anchor_input.h5"
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"

GENES_ARG=()
if [[ -n "${GENES:-}" ]]; then GENES_ARG=(--genes "$GENES"); fi
COMPARE_ARG=()
if [[ -n "${COMPARE_TO:-}" ]]; then COMPARE_ARG=(--compare-to "$COMPARE_TO"); fi
NULL_ARG=()
if [[ "${NULL_SETS:-0}" -gt 0 ]]; then
    NULL_ARG=(--null-sets "$NULL_SETS" --null-seed "${NULL_SEED:-0}"
              --null-neighbours "${NULL_NEIGHBOURS:-100}")
fi
SPAN_ARG=()
if [[ "${SPANNING_SETS:-0}" -gt 0 ]]; then
    SPAN_ARG=(--spanning-sets "$SPANNING_SETS")
fi

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python "${PIPELINE_DIR}/python/diagnose_program_by_section.py" \
        --counts-h5 "$WORK/anchor_input.h5" \
        --typing-h5 "$WORK/anchor_typing.h5" \
        --letter "$LETTER" \
        --unit "${UNIT:-fov}" \
        --min-section-cells "${MIN_CELLS:-25}" \
        ${GENES_ARG[@]+"${GENES_ARG[@]}"} \
        ${COMPARE_ARG[@]+"${COMPARE_ARG[@]}"} \
        ${NULL_ARG[@]+"${NULL_ARG[@]}"} \
        ${SPAN_ARG[@]+"${SPAN_ARG[@]}"} \
        --output-csv "$WORK/${STEM}_program_by_section.csv" \
        --per-unit-csv "$WORK/${STEM}_program_by_unit.csv" \
        --null-per-unit-csv "$WORK/${STEM}_program_by_unit_null.csv"

echo "Uploading the tables to Kopah (tag: ${TAG})..."
DEST="${BASE}/${STAGE4}/supervised_gbmap/program_confound"
s5cmd cp "$WORK/${STEM}_program_by_section.csv" "${DEST}/${STEM}_program_by_section.csv"

# The first two are written only when COMPARE_TO named more than one type; the null pair only
# when NULL_SETS was set.
for extra in "${STEM}_program_by_section_by_type.csv" "${STEM}_program_by_unit.csv" \
             "${STEM}_program_by_section_null.csv" \
             "${STEM}_program_by_section_null_draws.csv" \
             "${STEM}_program_by_unit_null.csv" \
             "${STEM}_program_by_section_dilution.csv"; do
    if [[ -f "$WORK/$extra" ]]; then
        s5cmd cp "$WORK/$extra" "${DEST}/${extra}"
    fi
done

echo "Done. The verdict is in the summary block printed above."
