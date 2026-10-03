#!/usr/bin/env bash
# Launch the NFL fecta pricer from the repo root.
#
#   nfl_specials/run.sh          # http://127.0.0.1:8096
#
# Start dk_price_sidecar/run.sh first (in its own terminal) for DraftKings
# prices; without it the board still prices trifectas at FanDuel + BetMGM.
# Credentials: WAGERZON[_X]_USERNAME / _PASSWORD in bet_logger/.env.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."
exec python3 -m nfl_specials.app "$@"
