#!/usr/bin/env bash
# Launch the SGP compare page from the repo root.
#
#   sgp_compare/run.sh          # http://127.0.0.1:8097
#
# Needs: curl_cffi, duckdb, scipy. No credentials (ProphetX and BetMGM are read anonymously).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."
exec python3 -m sgp_compare.app "$@"
