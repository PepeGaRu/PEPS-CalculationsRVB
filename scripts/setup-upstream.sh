#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
rvb_env="${RVB_ENV_DIR:-/workspace/peps-environment}"
export UV_CACHE_DIR="$rvb_env/uv-cache"
if [[ ! -x "$rvb_env/.venv/bin/python" ]]; then
  printf '%s\n' 'Run bash scripts/setup.sh first.' >&2
  exit 1
fi
uv pip install --python "$rvb_env/.venv/bin/python" --require-hashes \
  --index-url https://download.pytorch.org/whl/cpu -r requirements-upstream.lock
rvb_upstream="$rvb_env/upstream/peps-torch"
rvb_commit=239a2e717a567df0ece80a5bce5443a6f78d048b
if [[ ! -e "$rvb_upstream" ]]; then
  mkdir -p "$(dirname "$rvb_upstream")"
  git clone --filter=blob:none --no-checkout https://github.com/jurajHasik/peps-torch.git "$rvb_upstream"
  git -C "$rvb_upstream" checkout --detach "$rvb_commit"
fi
if [[ "$(git -C "$rvb_upstream" rev-parse HEAD)" != "$rvb_commit" ]]; then
  printf '%s\n' 'Existing upstream source has a different commit; preserve it and choose a fresh RVB_ENV_DIR.' >&2
  exit 1
fi
git -C "$rvb_upstream" diff --quiet HEAD --
"$rvb_env/.venv/bin/python" -c 'import torch; assert torch.__version__ == "2.6.0+cpu"; assert torch.version.cuda is None; print("CPU PyTorch:", torch.__version__)'
printf '%s\n' "Pinned peps-torch source ready at $rvb_upstream"
