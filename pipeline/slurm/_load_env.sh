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
# `FOO=bar sbatch stage.sh` now wins, and the override is ECHOED, so a log shows which
# values a run actually used rather than which ones were requested.
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

_CALLER_SET=()
for _name in $_ENV_ASSIGNED_NAMES; do
    # ${!_name+set} is true for "set but empty" too: passing FOO= is a deliberate override
    # to empty, and .env must not quietly refill it.
    if [[ -n "${!_name+set}" ]]; then
        _CALLER_SET+=("${_name}=${!_name}")
    fi
done

set -a
# shellcheck disable=SC1090
source "$ENV_FILE"
set +a

if (( ${#_CALLER_SET[@]} )); then
    for _override in "${_CALLER_SET[@]}"; do
        export "$_override"
    done
    echo "Command-line overrides kept over .env: ${_CALLER_SET[*]}"
fi

unset _ENV_ASSIGNED_NAMES _CALLER_SET _name _override
