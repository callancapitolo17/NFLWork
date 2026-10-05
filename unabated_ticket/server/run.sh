#!/usr/bin/env bash
# Launch the Unabated Ticket server runner (the headless Edges scan, and the
# closing fairs it saves for the Bet Tracker's CLV) from the repo root.
#
#   ./unabated_ticket/server/run.sh        # http://127.0.0.1:8095/edges.json
#
# Needs Node 18+ on PATH (Homebrew's /opt/homebrew/bin on the Mac). No npm install.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"
exec node unabated_ticket/server/runner.js "$@"
