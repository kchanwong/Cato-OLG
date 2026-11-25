#!/usr/bin/env bash
# Runs jobs inside a built image the way AWS Batch does, to check the image
# itself works before it is published: packages load without recompiling, a
# baseline and a full-size reform job finish inside Batch's limits, and a bad
# job fails.
#
#   docker build -t olg:test . && test/smoke_image.sh olg:test
#
# Limits mirror the Batch job definition in aws-platform
# (stacks/backburner/stack.py): 8 GB of memory and 900 seconds per attempt.
set -euo pipefail

image=${1:?usage: test/smoke_image.sh IMAGE}
out=$(mktemp -d)
trap 'rm -rf "$out"' EXIT

run_job() {  # run_job NAME JSON: runs one job, results in $out/NAME
  mkdir -p "$out/$1"
  echo "::group::$1"
  timeout 900 docker run --rm --memory 8g -v "$out/$1:/data" "$image" "$2" 2>&1 | tee "$out/$1.log"
  echo "::endgroup::"
}

fail() { echo "FAIL: $*" >&2; exit 1; }

run_job baseline '{"projection":{"n_years":10}}'
[[ -s $out/baseline/score.json && -s $out/baseline/score.xlsx ]] || fail "baseline wrote no results"
[[ $(jq -c 'keys' "$out/baseline/score.json") == '["Baseline"]' ]] || fail "baseline results have the wrong sections"
# Every package should already be compiled into the image. Recompiling here
# would repeat at the start of every Batch job and eat into its time limit.
if grep -q "Precompiling" "$out/baseline.log"; then fail "packages were recompiled when the job started"; fi

# The same kind of job the dashboard submits, at full model size.
run_job reform '{"label":"FRA 70","reform_year":2034,"params":{"J_retire":70}}'
[[ $(jq -c 'keys' "$out/reform/score.json") == '["Baseline","Dynamic","Static"]' ]] \
  || fail "reform results are missing a score"

# A bad job must exit non-zero so Batch marks it FAILED.
if docker run --rm "$image" '{"benefit_fn":"nope"}' >/dev/null 2>&1; then
  fail "a job with an unknown benefit rule succeeded"
fi

echo "image smoke test passed"
