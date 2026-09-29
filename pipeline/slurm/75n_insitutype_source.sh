#!/bin/bash
# What actually sets `clust`? Dump the INSTALLED InSituType source, not GitHub's.
#
# 75c measures that 1.13% of cells (28,672, every one of them named-labelled) carry a label that
# is NOT the argmax of the fit's own stored log-likelihoods. Per the published source that
# cannot happen: insitutypeML() does
#     clust <- colnames(logliks)[apply(logliks, 1, which.max)]
# the plain argmax, with the cohort-frequency adjustment already folded into the matrix it
# returns, and insitutype()'s only post-step pins cells passed as `anchors` — which this fit
# does not use (insitutype_typing.R passes neither `anchors` nor `cohort`).
#
# Two candidate explanations are already ruled out: anchors (not passed) and the frequency
# penalty on rare types (it predicts small types mismatch most; size and mismatch rate are
# uncorrelated at r = -0.16, and in a size-matched comparison of the 20 smallest named types the
# vascular ones run 92.6% against 12.2% for the rest). So this checks the remaining candidate:
# that the CONTAINER's version differs from the published source.
#
# What to look for in the output:
#   1. is `clust` derived from the SAME matrix the function returns?
#   2. is update_logliks_with_cohort_freqs applied BEFORE the argmax, and is its result the
#      copy that gets returned?
#   3. any post-hoc reassignment — a minimum-cluster-size rule, a refinement pass, a merge that
#      relabels — able to move a cell after it has been scored. This fit passes `refinement`.
#   4. the version, against the published one.
#
# Everything goes to stdout, so the answer is in the .out log — no Kopah round trip.
#
# MUST run under Slurm: on a Hyak login node the site apptainer config bind-mounts
# /var/run/slurm and /var/spool/slurmd, which exist only on compute nodes, and the container
# refuses to start.
#
# Submit:
#   sbatch pipeline/slurm/75n_insitutype_source.sh

#SBATCH --job-name=cosmx-insitutype-source
#SBATCH --account=glioblastoma-ckpt
#SBATCH --partition=ckpt
#SBATCH --qos=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=pipeline/logs/insitutype_source_%j.out
#SBATCH --error=pipeline/logs/insitutype_source_%j.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
    source /etc/profile.d/lmod.sh 2>/dev/null \
        || source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
fi
module load apptainer

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PIPELINE_DIR="${SLURM_SUBMIT_DIR}/pipeline"
else
    PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
fi

set -a
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/.env"
set +a

: "${APPTAINER_INSITUTYPE:?must be set in pipeline/.env}"

apptainer exec "$APPTAINER_INSITUTYPE" Rscript "${PIPELINE_DIR}/R/dump_insitutype_source.R"

echo "Done. The verdict is whether clust comes off the matrix that is returned."
