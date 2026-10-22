#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
rvb_env="${RVB_ENV_DIR:-/workspace/peps-environment}"
mkdir -p "$rvb_env"
export UV_CACHE_DIR="$rvb_env/uv-cache"
if [[ ! -x "$rvb_env/.venv/bin/python" ]]; then
  uv venv "$rvb_env/.venv" --python 3.12
fi
uv pip install --python "$rvb_env/.venv/bin/python" --require-hashes -r requirements.lock
uv pip install --python "$rvb_env/.venv/bin/python" --require-hashes -r requirements-paper.lock
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
export MPLCONFIGDIR="$rvb_env/matplotlib" XDG_CACHE_HOME="$rvb_env/cache"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
"$rvb_env/.venv/bin/python" -m pytest -q
printf '%s\n' "Core environment verified. Run: bash run.sh --profile reference --output results/reference"
