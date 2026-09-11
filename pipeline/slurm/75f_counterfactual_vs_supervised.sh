#!/bin/bash
# Head-to-head: the de-novo-removed COUNTERFACTUAL against a FRESH SUPERVISED run.
#
# The PI asked what the 27 de-novo letters would be called if forced onto GBmap. There are two
# defensible readings and she wants both:
#   A  counterfactual   hold the k=27 fit fixed, drop its 27 de-novo columns, take the argmax of
#                       the 54 named ones (75c's forced_named_posteriors.csv). The profiles are
#                       whatever that semi-supervised fit converged on.
#   B  fresh supervised insitutype(n_clusts = 0) from the published GBmap reference, no de-novo
#                       clusters at any point, so the profile update and rescale adapt to ALL
#                       2.54M cells rather than the ~6.7% a semi-supervised fit leaves named.
#
# Neither is more correct. Where they agree, a letter's GBmap identity is robust to how the
# question is asked; where they diverge, it depends on whether the de-novo clusters were present
# while the profiles were fitted — worth knowing before those letters get named.
#
# PREREQUISITE — produce B first (~1-2h; supervised mode skips every clustering phase):
#   INPUT_DIR=stage4_anchor STAGE4_DIR=stage4_anchor_supervised N_CLUSTS=0 \
#     KEEP_GENES=stage4_anchor/gene_selection/kept_genes.txt \
#     sbatch pipeline/slurm/72_anchor_typing.sh
#
# Then this job:
#   sbatch pipeline/slurm/75f_counterfactual_vs_supervised.sh
#
# And, to get the letters-vs-B Sankey and cross-tabs in the same shape as the A ones:
#   FORCED_H5=stage4_anchor_supervised/anchor/anchor_typing.h5 \
#     sbatch pipeline/slurm/75b_denovo_vs_supervised.sh
#
# Env knobs (KOPAH_*, APPTAINER_RSC from pipeline/.env):
#   STAGE4_DIR       Kopah sub-dir with the k=27 fit (default stage4_anchor_pruned).
#   SUPERVISED_DIR   Kopah sub-dir with the N_CLUSTS=0 run (default stage4_anchor_supervised).

#SBATCH --job-name=cosmx-cf-vs-supervised
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=pipeline/logs/cf_vs_supervised_%j.out
#SBATCH --error=pipeline/logs/cf_vs_supervised_%j.err

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
SUPERVISED="${SUPERVISED_DIR:-stage4_anchor_supervised}"
ANNOTATIONS="${ANNOTATIONS:-reference/denovo_annotations/fullcohort_pruned_k27.csv}"
: "${APPTAINER_RSC:?must be set in pipeline/.env}"

WORK="${SLURM_TMPDIR:-/tmp}/cosmx_cf_vs_supervised_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/out"
trap 'rm -rf "$WORK"' EXIT

export AWS_ACCESS_KEY_ID="$KOPAH_ACCESS_KEY_ID"
export AWS_SECRET_ACCESS_KEY="$KOPAH_SECRET_ACCESS_KEY"
export S3_ENDPOINT_URL="$KOPAH_ENDPOINT_URL"

BASE="s3://${KOPAH_BUCKET}/${KOPAH_PREFIX}"

echo "Staging the k=27 fit, its counterfactual call, and the supervised run..."
s5cmd cp "${BASE}/${STAGE4}/anchor/anchor_typing.h5" "$WORK/anchor_typing.h5"
s5cmd cp "${BASE}/${STAGE4}/supervised_gbmap/forced_named_posteriors.csv" "$WORK/counterfactual.csv"
s5cmd cp "${BASE}/${SUPERVISED}/anchor/anchor_typing.h5" "$WORK/supervised.h5"

apptainer exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --bind "${WORK}:${WORK}" \
    "$APPTAINER_RSC" \
    python "${PIPELINE_DIR}/python/compare_forced_call_sources.py" \
        --typing-h5 "$WORK/anchor_typing.h5" \
        --source-a "$WORK/counterfactual.csv" \
        --source-b "$WORK/supervised.h5" \
        --annotations "${PIPELINE_DIR}/${ANNOTATIONS}" \
        --output-dir "$WORK/out"

echo "Uploading the comparison to Kopah..."
s5cmd cp "$WORK/out/*" "${BASE}/${STAGE4}/supervised_gbmap/counterfactual_vs_supervised/"

echo "Done. Agreement summary is printed above."
