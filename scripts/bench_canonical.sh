#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${REPO_ROOT}"

# Run every fixture tagged "canonical" against wandb. Pass extra args through
# (e.g. --notes, --fixtures-dir bench_out/staged, --dry-run).
exec uv run --group bench python -m bench.wandb_bench --tag canonical "$@"
