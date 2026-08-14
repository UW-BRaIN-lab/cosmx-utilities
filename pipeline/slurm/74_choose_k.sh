#!/bin/bash
# De-novo cluster-number sweep for the anchor reference rebuild (choose_denovo_k.R).
# Re-runs InSituType's own AIC-based K selection across a WIDER range (default 15:35) and
# reports the whole curve, so we can confirm the optimum turns over instead of censoring at
# the range edge the way the pilot's 10:20 did (it picked K=20, the ceiling).
#
# Run it TWICE to answer "does pruning change the optimal K?":
#   full panel:   sbatch pipeline/slurm/74_choose_k.sh
#   pruned panel: KEEP_GENES=stage4_anchor/gene_selection/kept_genes.txt \
#                 K_SWEEP_OUT=k_sweep_pruned sbatch pipeline/slurm/74_choose_k.sh
#
# Required env (pipeline/.env): KOPAH_*, APPTAINER_INSITUTYPE.

#SBATCH --job-name=cosmx-choose-k
# CPU-only InSituType EM on subsets — heavier than gene selection, lighter than the anchor
# fit. ckpt is the default; if it gets preempted+requeued, override to a dedicated partition
# (e.g. --account=carrolllab --partition=cpu-g2-mem2x) since this does NOT checkpoint.
#SBATCH --account=glioblastoma
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=06:00:00
#SBATCH --output=pipeline/logs/choose_k_%j.out
#SBATCH --error=pipeline/logs/choose_k_%j.err

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

STAGE4="${STAGE4_DIR:-stage4_anchor}"
REFERENCE_BASENAME="${REFERENCE_BASENAME:-gbmap_level4_panel.csv}"
K_SWEEP_OUT="${K_SWEEP_OUT:-k_sweep}"
: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_choose_k_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging anchor input + reference from Kopah..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/anchor/anchor_input.h5" \
    "$WORK/anchor_input.h5"
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/reference/${REFERENCE_BASENAME}" \
    "$WORK/reference.csv"

# Optional pruned panel: KEEP_GENES is a Kopah key (under ${KOPAH_PREFIX}/) to kept_genes.txt.
KEEP_ARG=()
if [[ -n "${KEEP_GENES:-}" ]]; then
    s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${KEEP_GENES}" "$WORK/kept_genes.txt"
    KEEP_ARG=(--keep-genes "$WORK/kept_genes.txt")
    echo "Using pruned panel from ${KEEP_GENES}"
fi

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_INSITUTYPE" \
    Rscript "${PIPELINE_DIR}/R/choose_denovo_k.R" \
        --input "$WORK/anchor_input.h5" \
        --reference "$WORK/reference.csv" \
        --out-dir "$WORK/$K_SWEEP_OUT" \
        --n-clusts "${K_SWEEP_RANGE:-15:35}" \
        --subset-size "${K_SWEEP_SUBSET:-10000}" \
        --n-reps "${K_SWEEP_REPS:-2}" \
        --max-iters "${K_SWEEP_MAX_ITERS:-20}" \
        --max-cells "${K_SWEEP_MAX_CELLS:-500000}" \
        ${KEEP_ARG[@]+"${KEEP_ARG[@]}"}

echo "Uploading K-sweep outputs to Kopah (${STAGE4}/${K_SWEEP_OUT})..."
for f in k_sweep.csv k_sweep.png k_sweep_summary.txt; do
    s5cmd cp "$WORK/$K_SWEEP_OUT/$f" \
        "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE4}/${K_SWEEP_OUT}/$f"
done

echo "=== K-sweep summary ==="
cat "$WORK/$K_SWEEP_OUT/k_sweep_summary.txt"
echo "Done: de-novo K-sweep (${K_SWEEP_OUT})."
