#!/usr/bin/env bash
# Launch the Unabated Ticket bets service from the repo root.
#
#   ./unabated_ticket/bets_service/run.sh        # http://127.0.0.1:8094
#   BETS_SERVICE_PORT=9000 ./unabated_ticket/bets_service/run.sh
#
# Credentials: unabated_ticket/bets_service/.env, else the environment, else
# kalshi_draft/.env in the main checkout (see .env.example). Uses the
# kalshi_draft venv when present (it has duckdb + cryptography), else python3.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
MAIN_ROOT="${REPO_ROOT%%/.worktrees/*}"
PYTHON="$MAIN_ROOT/kalshi_draft/venv/bin/python3"
if [ ! -x "$PYTHON" ]; then PYTHON=python3; fi
cd "$REPO_ROOT"
exec "$PYTHON" -m unabated_ticket.bets_service.service "$@"
