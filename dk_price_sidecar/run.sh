#!/usr/bin/env bash
# Launch the DraftKings price sidecar (issue #102).
#
#   dk_price_sidecar/run.sh                       # http://127.0.0.1:8095
#   DK_SIDECAR_PORT=9095 dk_price_sidecar/run.sh
#
# Needs a Python with Playwright and a real Google Chrome install. The repo's
# default python3 has neither; DK_SIDECAR_PYTHON overrides the interpreter.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PY="${DK_SIDECAR_PYTHON:-/Library/Frameworks/Python.framework/Versions/3.12/bin/python3}"
"$PY" -c "import playwright" 2>/dev/null || {
  echo "error: $PY has no playwright — pip install playwright, or set DK_SIDECAR_PYTHON" >&2
  exit 1
}
cd "$REPO_ROOT"
exec "$PY" -m dk_price_sidecar.server
