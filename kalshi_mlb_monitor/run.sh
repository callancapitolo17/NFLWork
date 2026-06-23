#!/usr/bin/env bash
# Launch the Kalshi MLB Bots Monitor (read-only dashboard over both bots).
#
#   ./run.sh                  # serves on http://127.0.0.1:8092
#   KALSHI_MLB_MONITOR_PORT=9000 ./run.sh
#
# Reads the LIVE bot DBs under ~/NFLWork by default (override with
# KALSHI_MLB_MONITOR_ROOT). It never writes to them.
set -euo pipefail

# Repo root = two levels up from this script (kalshi_mlb_monitor/ -> repo).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Default the data root to the canonical checkout so we always read live data,
# even when this code is launched from a git worktree.
export KALSHI_MLB_MONITOR_ROOT="${KALSHI_MLB_MONITOR_ROOT:-$HOME/NFLWork}"

cd "$REPO_ROOT"
echo "Kalshi MLB Bots Monitor — data root: $KALSHI_MLB_MONITOR_ROOT"
echo "Serving on http://127.0.0.1:${KALSHI_MLB_MONITOR_PORT:-8092}"
exec python3 -m kalshi_mlb_monitor.app
