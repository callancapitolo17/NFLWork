#!/usr/bin/env bash
# Launch the SGP compare page from the repo root.
#
#   sgp_compare/run.sh          # http://127.0.0.1:8097
#
# Needs: curl_cffi, duckdb, scipy, cryptography. ProphetX needs no credentials; Kalshi
# reads KALSHI_API_KEY_ID / KALSHI_PRIVATE_KEY_PATH / KALSHI_USER_ID (env or kalshi_mlb_rfq/.env).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."
exec python3 -m sgp_compare.app "$@"
