#!/bin/bash
# Tests for slurm/_load_env.sh: the caller's command line must beat pipeline/.env.
#
# Run:  bash pipeline/slurm/tests/test_load_env.sh
#
# The property under test is the one that cost a 12-hour fit: an override that is silently
# ignored. Every case here asserts on the VALUE a stage would actually use, not on whether
# the load succeeded -- the failure mode is a successful load of the wrong thing.

set -euo pipefail

HELPER="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/_load_env.sh"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

mkdir -p "$WORK/pipeline/slurm"
cp "$HELPER" "$WORK/pipeline/slurm/_load_env.sh"

cat > "$WORK/pipeline/.env" <<'ENVFILE'
# comment line, and a blank one follow

REFERENCE_BASENAME=from_env.csv
KOPAH_BUCKET=brainlabkg
export EXPORTED_NAME=from_env
EMPTY_IN_ENV=
  INDENTED_NAME=from_env
ENVFILE

cat > "$WORK/stage.sh" <<'STAGE'
set -euo pipefail
PIPELINE_DIR="$1"
# shellcheck disable=SC1091
source "${PIPELINE_DIR}/slurm/_load_env.sh"
echo "REFERENCE_BASENAME=${REFERENCE_BASENAME}"
echo "KOPAH_BUCKET=${KOPAH_BUCKET}"
echo "EXPORTED_NAME=${EXPORTED_NAME}"
echo "INDENTED_NAME=${INDENTED_NAME}"
echo "EMPTY_IN_ENV=[${EMPTY_IN_ENV}]"
STAGE

fail() { echo "FAIL: $1" >&2; exit 1; }
assert_contains() {
    grep -qF -- "$2" <<< "$1" || fail "$3"$'\n'"--- got ---"$'\n'"$1"
}

echo "test_load_env"

out="$(bash "$WORK/stage.sh" "$WORK/pipeline")"
assert_contains "$out" "REFERENCE_BASENAME=from_env.csv" ".env should supply unset names"
assert_contains "$out" "INDENTED_NAME=from_env" "indented assignments must be seen"
[[ "$out" != *"overrides kept"* ]] || fail "nothing was overridden; should not announce one"
echo "  .env supplies what the caller did not set"

out="$(REFERENCE_BASENAME=from_cli.csv EXPORTED_NAME=from_cli \
       bash "$WORK/stage.sh" "$WORK/pipeline")"
assert_contains "$out" "REFERENCE_BASENAME=from_cli.csv" "caller must beat .env"
assert_contains "$out" "EXPORTED_NAME=from_cli" "caller must beat an 'export NAME=' line"
assert_contains "$out" "KOPAH_BUCKET=brainlabkg" ".env must still fill the rest"
assert_contains "$out" "overrides kept over .env" "an override must be announced in the log"
echo "  the command line wins, and says so"

# Passing FOO= is a deliberate override to empty; .env must not quietly refill it.
out="$(INDENTED_NAME= bash "$WORK/stage.sh" "$WORK/pipeline")"
assert_contains "$out" "INDENTED_NAME=" "an explicitly empty override must survive"
[[ "$out" != *"INDENTED_NAME=from_env"* ]] || fail ".env refilled an explicitly emptied name"
echo "  an explicitly empty override is honoured"

if (cd "$WORK" && PIPELINE_DIR="$WORK/absent" bash -c \
        'set -euo pipefail; source "'"$WORK"'/pipeline/slurm/_load_env.sh"') 2>/dev/null; then
    fail "a missing .env must be fatal"
fi
echo "  a missing .env is fatal"

echo "All tests passed."
