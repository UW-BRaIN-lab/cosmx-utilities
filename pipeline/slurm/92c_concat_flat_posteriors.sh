#!/bin/bash
# Stage 5f (SHARDED): concatenate the per-slide flat-posterior CSVs (92b) into one cohort
# flat_posteriors.csv. All per-slide CSVs share identical columns (same flat_posteriors.R, same
# TOP_K, same profile set), so this just keeps one header and appends the rest. No container.
#
# Submit after the 92b array completes:
#   sbatch --dependency=afterok:<92b_jid> pipeline/slurm/92c_concat_flat_posteriors.sh
#
# Env knobs (KOPAH_* from pipeline/.env): STAGE5_DIR (default stage5_insitucnv).

#SBATCH --job-name=cosmx-flat-post-concat
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=pipeline/logs/flat_post_concat_%j.out
#SBATCH --error=pipeline/logs/flat_post_concat_%j.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/lmod.sh 2>/dev/null \
        || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
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

STAGE5="${STAGE5_DIR:-stage5_insitucnv}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_flat_post_concat_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/per_slide"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

echo "Staging per-slide flat posteriors from Kopah (${STAGE5}/posteriors/per_slide)..."
s5cmd cp "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/posteriors/per_slide/*" \
    "$WORK/per_slide/"

mapfile -t CSVS < <(find "$WORK/per_slide" -name '*.csv' | sort)
if (( ${#CSVS[@]} == 0 )); then
    echo "ERROR: no per-slide CSVs staged" >&2; exit 1
fi
echo "Concatenating ${#CSVS[@]} per-slide CSVs..."

OUT="$WORK/flat_posteriors.csv"
head -n 1 "${CSVS[0]}" > "$OUT"
for f in "${CSVS[@]}"; do
    tail -n +2 "$f" >> "$OUT"
done
rows=$(( $(wc -l < "$OUT") - 1 ))
echo "Cohort flat_posteriors.csv: ${rows} cells across ${#CSVS[@]} slides."

echo "Uploading cohort flat posteriors to Kopah (${STAGE5}/posteriors)..."
s5cmd cp "$OUT" \
    "s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}/${STAGE5}/posteriors/flat_posteriors.csv"

echo "Done: cohort flat posteriors assembled (${rows} cells)."
