#!/usr/bin/env bash
# Run the validation a release is held to, on the machine it is run on.
#
# CI runs the network-marked corpus nightly on a hosted runner with no GPU.
# This runs the same steps, and the live RCSB fetch check and a CUDA smoke
# campaign besides, on a workstation such as aailab01, one stage after
# another so that only one CUDA process runs at a time. Each stage writes
# its own log, JUnit XML where it is pytest, and a row in registry.csv; a
# stage that fails is recorded and the next one runs.
#
# Usage, from a checkout, with the environment activated or named:
#
#   scripts/run_validation_campaign.sh --output ~/validation-2026-10-15/run-001 \
#       [--env ~/.conda/envs/fastmdx-validation-v260] [--stages fetch,corpus]
#
# The output must not exist: no stage reuses or overwrites one.
set -uo pipefail

OUTPUT=""
ENV_PREFIX=""
STAGES="info,fetch-unit,fetch-live,openmm,network,corpus,smoke,suite"
REPEATS=20

while [[ $# -gt 0 ]]; do
  case "$1" in
    --output) OUTPUT="$2"; shift 2 ;;
    --env) ENV_PREFIX="$2"; shift 2 ;;
    --stages) STAGES="$2"; shift 2 ;;
    --repeats) REPEATS="$2"; shift 2 ;;
    -h|--help) sed -n '2,17p' "$0"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

[[ -n "$OUTPUT" ]] || { echo "--output is required" >&2; exit 2; }
if [[ -e "$OUTPUT" ]]; then
  echo "$OUTPUT exists; a campaign is never written over another." >&2
  exit 2
fi
if [[ -n "$ENV_PREFIX" ]]; then
  export PATH="$ENV_PREFIX/bin:$PATH"
fi
cd "$(dirname "$0")/.." || exit 2
mkdir -p "$OUTPUT"
OUTPUT="$(cd "$OUTPUT" && pwd)"
REGISTRY="$OUTPUT/registry.csv"
echo "stage,status,exit_code,started_at,ended_at,seconds,log,junit" > "$REGISTRY"

{
  echo "commit: $(git rev-parse HEAD)"
  echo "describe: $(git describe --tags --always --dirty)"
  echo "host: $(hostname)"
  echo "python: $(command -v python)"
  nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader 2>/dev/null \
    | sed 's/^/gpu: /'
} > "$OUTPUT/identity.txt"

gpu_busy() {
  nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null | grep -q .
}

stage() {
  local name="$1"; shift
  [[ ",$STAGES," == *",$name,"* ]] || return 0
  local log="$OUTPUT/$name.log" junit="" started ended code status t0
  if [[ "$*" == *pytest* ]]; then
    junit="$OUTPUT/$name.junit.xml"
    set -- "$@" --junitxml="$junit"
  fi
  # Another process on the card would make every timing here, and the
  # CUDA stages' memory, someone else's.
  while gpu_busy; do
    echo "$(date -u +%FT%TZ) $name waits for the GPU to be free" | tee -a "$OUTPUT/waits.log"
    sleep 300
  done
  started="$(date -u +%FT%TZ)"; t0=$SECONDS
  echo "== $started $name: $*" | tee -a "$OUTPUT/campaign.log"
  "$@" > "$log" 2>&1
  code=$?
  ended="$(date -u +%FT%TZ)"
  status=$([[ $code -eq 0 ]] && echo PASS || echo FAIL)
  echo "== $ended $name: $status (exit $code, $((SECONDS - t0)) s)" | tee -a "$OUTPUT/campaign.log"
  echo "$name,$status,$code,$started,$ended,$((SECONDS - t0)),$log,$junit" >> "$REGISTRY"
}

stage info        fastmdx info --json
stage fetch-unit  python -m pytest tests/test_a_fetch_is_tried_again.py -p no:cacheprovider
stage fetch-live  python scripts/check_rcsb_fetch.py --output "$OUTPUT/fetch-live" --repeats "$REPEATS"
stage openmm      python -m pytest tests/test_seeding_end_to_end.py -p no:cacheprovider \
                    --timeout-method=signal --timeout=600
stage network     python -m pytest -m network tests --ignore=tests/validation -p no:cacheprovider \
                    --timeout-method=signal --timeout=600
stage corpus      python -m pytest -m network tests/validation -p no:cacheprovider \
                    --timeout-method=signal --timeout=2400
stage smoke       python scripts/run_pdb_smoke_campaign.py --input-list examples/pdb_list.txt \
                    --output-root "$OUTPUT/smoke" --platform CUDA --max-setup-atoms 250000 \
                    --continue-on-error
stage suite       python -m pytest -p no:cacheprovider

echo "Registry: $REGISTRY"
column -s, -t < "$REGISTRY" | cut -c1-120
! grep -q ',FAIL,' "$REGISTRY"
