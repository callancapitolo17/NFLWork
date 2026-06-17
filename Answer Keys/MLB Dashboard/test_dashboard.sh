#!/bin/bash
# Run the MLB dashboard from the current worktree against seeded, isolated data.
# Usage:
#   ./test_dashboard.sh                 # seed if needed, render-only, serve :8093
#   ./test_dashboard.sh --mlb           # re-run MLB.R against seeded DBs, then serve
#   ./test_dashboard.sh --pipeline      # full run.py mlb (live re-scrape), then serve
#   ./test_dashboard.sh --refresh ...   # re-clone seed data first
set -euo pipefail

DASH_DIR="$(cd "$(dirname "$0")" && pwd)"
WT_ROOT="$(git -C "$DASH_DIR" rev-parse --show-toplevel)"
PORT="${MLB_DASHBOARD_PORT:-8093}"
export MLB_DASHBOARD_PORT="$PORT"

MODE="render"
REFRESH=""
for arg in "$@"; do
  case "$arg" in
    --refresh)  REFRESH="--refresh" ;;
    --mlb)      MODE="mlb" ;;
    --pipeline) MODE="pipeline" ;;
    --render)   MODE="render" ;;
    *) echo "Unknown arg: $arg" >&2; exit 1 ;;
  esac
done

echo "== Seeding worktree data =="
bash "$DASH_DIR/seed_test_data.sh" $REFRESH

cd "$WT_ROOT"
case "$MODE" in
  pipeline)
    echo "== Full pipeline (run.py mlb) =="
    python3 "Answer Keys/run.py" mlb || echo "pipeline had errors, continuing"
    ;;
  mlb)
    echo "== MLB.R only =="
    Rscript "Answer Keys/MLB Answer Key/MLB.R" || echo "MLB.R had errors, continuing"
    ;;
esac

echo "== Rendering dashboard =="
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R" || echo "render had errors, continuing"

echo "== Starting test server on :$PORT =="
(sleep 2 && open "http://localhost:$PORT") &
python3 "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
