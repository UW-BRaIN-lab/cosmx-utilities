#!/bin/bash
# Load pipeline/.env WITHOUT overwriting anything the caller set on the command line.
#
# Sourced by every stage in place of a bare `set -a; source .env; set +a`. That bare form
# gives .env the last word, so `FOO=bar sbatch stage.sh` silently loses to a `FOO=` line in
# .env. It cost a 12-hour InSituTree fit on 2026-11-04: .env defines REFERENCE_BASENAME for
# stage 4b, stage 85 read the same name for a different object, and the fit ran against the
# atlas reference instead of the InSituTree profile matrix. Nothing errored at the time --
# it loaded a real file of the wrong shape and failed much later, pointing at the hierarchy.
#
# `FOO=bar sbatch stage.sh` now wins, and the NAMES of the overrides are echoed, so a log
# shows which values a run actually used rather than which ones were requested.
#
# TWO THINGS THIS DELIBERATELY DOES NOT DO, both learned on 2026-11-06:
#
# 1. It never prints a VALUE. sbatch inherits the submitting shell, so the caller set
#    includes credentials, and the first version wrote KOPAH_SECRET_ACCESS_KEY and
#    AWS_SOURCE_SECRET_ACCESS_KEY in clear text into every job log on a shared filesystem.
#
# 2. It never lets the AMBIENT ENVIRONMENT override the study's identity. `--export=ALL` is
#    sbatch's default, so a variable exported in the login shell arrives looking exactly
#    like a deliberate `FOO=bar sbatch`. A `kopah` helper in .bashrc exporting
#    KOPAH_PREFIX=cosmx-gbm therefore silently redirected a retina worktree's jobs at the
#    GBM study's data. Names matching PROTECTED_PATTERN are taken from .env and nowhere
#    else: a worktree IS a study, and which bucket/prefix/container/account it uses is a
#    property of the checkout, not of whoever's shell launched it. Per-run knobs
#    (INPUT_KEY, KEEP_GENES, N_CLUSTS, K_SWEEP_OUT, *_BASENAME, ...) override as before.
#
# Expects PIPELINE_DIR to be set. Sourced, not executed:
#     source "${PIPELINE_DIR}/slurm/_load_env.sh"

ENV_FILE="${PIPELINE_DIR}/.env"
if [[ ! -f "$ENV_FILE" ]]; then
    echo "ERROR: no such env file: $ENV_FILE" >&2
    exit 1
fi

# Only names .env actually assigns can be clobbered, so only those need protecting. Reading
# them from the file keeps this in step with .env automatically -- a per-stage list of knobs
# would go stale the first time someone added one.
_ENV_ASSIGNED_NAMES="$(sed -nE 's/^[[:space:]]*(export[[:space:]]+)?([A-Za-z_][A-Za-z0-9_]*)=.*/\2/p' \
    "$ENV_FILE" | sort -u)"

# Credentials, endpoints, containers and scheduler identity. These describe the environment
# a checkout belongs to, never a per-run choice, so .env is their only source. A pattern
# rather than a list, so adding a KOPAH_* or AWS_* name to .env protects it automatically.
PROTECTED_PATTERN='^(KOPAH_|AWS_|SOURCE_S3_|APPTAINER_|SLURM_)'

_CALLER_NAMES=()
_IGNORED_NAMES=()
for _name in $_ENV_ASSIGNED_NAMES; do
    # ${!_name+set} is true for "set but empty" too: passing FOO= is a deliberate override
    # to empty, and .env must not quietly refill it.
    [[ -n "${!_name+set}" ]] || continue
    if [[ "$_name" =~ $PROTECTED_PATTERN ]]; then
        _IGNORED_NAMES+=("$_name")
    else
        _CALLER_NAMES+=("$_name")
        declare "_CALLER_VALUE_${_name}=${!_name}"
    fi
done

set -a
# shellcheck disable=SC1090
source "$ENV_FILE"
set +a

if (( ${#_CALLER_NAMES[@]} )); then
    for _name in "${_CALLER_NAMES[@]}"; do
        _held="_CALLER_VALUE_${_name}"
        export "${_name}=${!_held}"
        unset "$_held"
    done
    # NAMES ONLY -- see the header. Never "${_name}=${!_name}".
    echo "Overrides kept over .env: ${_CALLER_NAMES[*]}"
fi

if (( ${#_IGNORED_NAMES[@]} )); then
    # Loud, because the inherited value may be pointing at another study's data and the job
    # is about to run anyway. Silence here is what made the GBM/retina mix-up invisible.
    echo "NOTE: ignoring inherited ${_IGNORED_NAMES[*]} -- .env wins for these."
fi
echo "Using KOPAH_PREFIX=${KOPAH_PREFIX:-<unset>} from ${ENV_FILE}"

unset _ENV_ASSIGNED_NAMES _CALLER_NAMES _IGNORED_NAMES _name _held PROTECTED_PATTERN
