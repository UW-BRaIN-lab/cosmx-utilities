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
#SBATCH --time=02:00:00
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

echo "Staging anchor counts + labels from Kopah (sections come from the cell ids)..."
s5cmd cp "${BASE}/${INPUT}/anchor/anchor_input.h5" "$WORK/anchor_input.h5"
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"

GENES_ARG=()
if [[ -n "${GENES:-}" ]]; then GENES_ARG=(--genes "$GENES"); fi
COMPARE_ARG=()
if [[ -n "${COMPARE_TO:-}" ]]; then COMPARE_ARG=(--compare-to "$COMPARE_TO"); fi

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
        --output-csv "$WORK/${STEM}_program_by_section.csv" \
        --per-unit-csv "$WORK/${STEM}_program_by_unit.csv"

echo "Uploading the tables to Kopah (tag: ${TAG})..."
DEST="${BASE}/${STAGE4}/supervised_gbmap/program_confound"
s5cmd cp "$WORK/${STEM}_program_by_section.csv" "${DEST}/${STEM}_program_by_section.csv"

# Both are written only when COMPARE_TO named more than one type.
for extra in "${STEM}_program_by_section_by_type.csv" "${STEM}_program_by_unit.csv"; do
    if [[ -f "$WORK/$extra" ]]; then
        s5cmd cp "$WORK/$extra" "${DEST}/${extra}"
    fi
done

echo "Done. The verdict is in the summary block printed above."
