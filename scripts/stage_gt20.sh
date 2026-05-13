#!/usr/bin/env bash
set -euo pipefail

FIXTURE="hela_iccoff_gt20peps"
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LOG_DIR="${REPO_ROOT}/bench_out/logs"
LOG_FILE="${LOG_DIR}/stage_${FIXTURE}.log"

mkdir -p "${LOG_DIR}"

cd "${REPO_ROOT}"
echo "Staging ${FIXTURE} -> bench_out/staged"
echo "Log: ${LOG_FILE}"

uv run --group bench python -m bench.stage_fixture "${FIXTURE}" \
  > "${LOG_FILE}" 2>&1

echo "Done. Staged inputs under bench_out/staged/${FIXTURE}/"
echo "Run with:"
echo "  uv run --group bench python -m bench.wandb_bench \\"
echo "    --fixtures-dir bench_out/staged ${FIXTURE}"
