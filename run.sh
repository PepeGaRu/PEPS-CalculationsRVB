#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export OPENBLAS_NUM_THREADS="${RVB_NUM_THREADS:-2}"
export OMP_NUM_THREADS="${RVB_NUM_THREADS:-2}"
export NUMBA_NUM_THREADS="${RVB_NUM_THREADS:-2}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/workspace/peps-environment/matplotlib}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/workspace/peps-environment/cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
rvb_python="${RVB_PYTHON:-/workspace/peps-environment/.venv/bin/python}"
if [[ ! -x "$rvb_python" ]]; then
  printf '%s\n' 'Python environment missing. Follow the installation commands in README.md.' >&2
  exit 1
fi
exec "$rvb_python" -m rvb "$@"
