#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"
python_command="${CODEX_PRIMARY_RUNTIME_PYTHON:-python3}"
"$python_command" -c 'import sys; assert sys.version_info >= (3, 12), "Python 3.12 or newer is required"'
if [[ ! -x .venv/bin/python ]]; then
  "$python_command" -m venv --system-site-packages .venv
fi
if ! .venv/bin/python scripts/check_env.py --versions-only; then
  .venv/bin/python -m pip install --disable-pip-version-check --timeout 15 --retries 0 -r requirements.txt
fi
./run.sh scripts/check_env.py
