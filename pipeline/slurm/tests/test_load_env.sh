#!/bin/bash
# Tests for slurm/_load_env.sh. Run directly:
#     bash pipeline/slurm/tests/test_load_env.sh
#
# The two behaviours under test are the two ways the loader has gone wrong in production,
# so each has a named case: a protected name must not be overridable from the ambient
# shell (the retina worktree ran against the GBM prefix), and no VALUE may ever be printed
# (job logs on a shared filesystem carried live secret keys).

set -uo pipefail

LOADER="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../_load_env.sh"
FAILURES=0

# Each case gets its own PIPELINE_DIR so one case's .env cannot leak into another.
make_pipeline_dir() {
    local dir
    dir="$(mktemp -d)"
    mkdir -p "${dir}/slurm"
    cp "$LOADER" "${dir}/slurm/_load_env.sh"
    cat > "${dir}/.env" <<'ENV'
KOPAH_PREFIX=cosmx-retina
KOPAH_SECRET_ACCESS_KEY=env-file-secret
KEEP_GENES=3390
REFERENCE_BASENAME=from_env.csv
ENV
    echo "$dir"
}

# Runs the loader with the given `NAME=value` assignments exported first, as sbatch's
# --export=ALL does, and echoes the resulting values of the names asked for.
run_loader() {
    local pipeline_dir="$1"; shift
    local assignments="$1"; shift
    env -i PATH="$PATH" HOME="$HOME" bash -c "
        ${assignments}
        PIPELINE_DIR='${pipeline_dir}'
        source '${pipeline_dir}/slurm/_load_env.sh'
        for name in $*; do printf '%s=%s\n' \"\$name\" \"\${!name-<unset>}\"; done
    " 2>&1
}

check() {
    local name="$1" expected="$2" actual="$3"
    if [[ "$actual" == *"$expected"* ]]; then
        echo "  PASS  ${name}"
    else
        echo "  FAIL  ${name}"
        echo "        expected to find: ${expected}"
        echo "        in output:"
        sed 's/^/          /' <<< "$actual"
        FAILURES=$((FAILURES + 1))
    fi
}

refute() {
    local name="$1" forbidden="$2" actual="$3"
    if [[ "$actual" != *"$forbidden"* ]]; then
        echo "  PASS  ${name}"
    else
        echo "  FAIL  ${name}"
        echo "        output must NOT contain: ${forbidden}"
        echo "        in output:"
        sed 's/^/          /' <<< "$actual"
        FAILURES=$((FAILURES + 1))
    fi
}

echo "An inherited KOPAH_PREFIX must lose to .env"
DIR="$(make_pipeline_dir)"
OUT="$(run_loader "$DIR" "export KOPAH_PREFIX=cosmx-gbm" KOPAH_PREFIX)"
check "the study's own prefix wins" "KOPAH_PREFIX=cosmx-retina" "$OUT"
check "and the override is announced" "ignoring inherited KOPAH_PREFIX" "$OUT"
check "every run records the prefix it used" "Using KOPAH_PREFIX=cosmx-retina" "$OUT"
rm -rf "$DIR"

echo "A per-run knob must still override .env"
DIR="$(make_pipeline_dir)"
OUT="$(run_loader "$DIR" "export KEEP_GENES=1200; export REFERENCE_BASENAME=tree.csv" \
    KEEP_GENES REFERENCE_BASENAME)"
check "KEEP_GENES overrides" "KEEP_GENES=1200" "$OUT"
check "REFERENCE_BASENAME overrides" "REFERENCE_BASENAME=tree.csv" "$OUT"
check "the names are listed" "Overrides kept over .env:" "$OUT"
rm -rf "$DIR"

echo "An override to EMPTY is deliberate and must not be refilled"
DIR="$(make_pipeline_dir)"
OUT="$(run_loader "$DIR" "export KEEP_GENES=" KEEP_GENES)"
check "KEEP_GENES stays empty" "KEEP_GENES=" "$OUT"
refute "it is not refilled from .env" "KEEP_GENES=3390" "$OUT"
rm -rf "$DIR"

echo "No VALUE may ever be printed"
DIR="$(make_pipeline_dir)"
OUT="$(run_loader "$DIR" \
    "export KOPAH_SECRET_ACCESS_KEY=live-secret-from-the-login-shell; export KEEP_GENES=1200" \
    KEEP_GENES)"
refute "not the inherited secret" "live-secret-from-the-login-shell" "$OUT"
refute "not the .env secret" "env-file-secret" "$OUT"
check "only the name" "Overrides kept over .env: KEEP_GENES" "$OUT"
# A value cannot leak without an '=' on the loader's own line, whatever the name is.
ANNOUNCED="$(grep '^Overrides kept over .env:' <<< "$OUT")"
refute "the announcement carries no values at all" "=" "${ANNOUNCED#Overrides kept over .env:}"
rm -rf "$DIR"

echo "A missing .env is fatal, not silent"
DIR="$(mktemp -d)"
mkdir -p "${DIR}/slurm"
cp "$LOADER" "${DIR}/slurm/_load_env.sh"
OUT="$(run_loader "$DIR" "true" KOPAH_PREFIX)"
STATUS=$?
check "it says which file is missing" "no such env file" "$OUT"
if (( STATUS != 0 )); then
    echo "  PASS  and exits non-zero"
else
    echo "  FAIL  and exits non-zero (got ${STATUS})"
    FAILURES=$((FAILURES + 1))
fi
rm -rf "$DIR"

echo
if (( FAILURES )); then
    echo "${FAILURES} failure(s)"
    exit 1
fi
echo "all cases pass"
