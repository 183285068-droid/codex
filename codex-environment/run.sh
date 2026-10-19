#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"
if [[ ! -x .venv/bin/python ]]; then
  echo 'Please run bash setup.sh first.' >&2
  exit 1
fi
export PYTHONUTF8=1
export MPLBACKEND=Agg
export MPLCONFIGDIR="$project_dir/.cache/matplotlib"
export PYTHONPYCACHEPREFIX="$project_dir/.cache/python"
exec "$project_dir/.venv/bin/python" "$@"
