#!/bin/bash
# Seed the current worktree with CoW clones of main's live DuckDB files so the
# worktree can run the MLB pipeline + dashboard in full isolation. Safe to run
# while the live dashboard / RFQ bot are running: cp is a byte-level read of
# main's files and never blocks or mutates them.
set -euo pipefail

MAIN_ROOT="$HOME/NFLWork"
WT_ROOT="$(git rev-parse --show-toplevel)"

if [ "$WT_ROOT" = "$MAIN_ROOT" ]; then
  echo "Refusing to seed: you are in the main repo, not a worktree." >&2
  exit 1
fi

# DBs to clone: <relative path under repo root>
DBS=(
  "Answer Keys/pbp.duckdb"
  "Answer Keys/mlb.duckdb"
  "Answer Keys/mlb_mm.duckdb"
  "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
  "dk_odds/dk.duckdb"
  "fd_odds/fd.duckdb"
  "wagerzon_odds/wagerzon.duckdb"
  "hoop88_odds/hoop88.duckdb"
  "bfa_odds/bfa.duckdb"
  "bookmaker_odds/bookmaker.duckdb"
  "bet105_odds/bet105.duckdb"
)

REFRESH=0
[ "${1:-}" = "--refresh" ] && REFRESH=1

for rel in "${DBS[@]}"; do
  src="$MAIN_ROOT/$rel"
  dst="$WT_ROOT/$rel"
  if [ ! -f "$src" ]; then
    echo "  skip (no source): $rel"
    continue
  fi
  if [ -f "$dst" ] && [ "$REFRESH" -eq 0 ]; then
    echo "  keep (exists):   $rel"
    continue
  fi
  mkdir -p "$(dirname "$dst")"
  # CoW clone on APFS; falls back to a normal copy on non-APFS volumes.
  cp -c "$src" "$dst" 2>/dev/null || cp "$src" "$dst"
  echo "  seeded:          $rel"
done

echo "Seed complete in: $WT_ROOT"
