#!/usr/bin/env python3
"""BetMGM MLB SGP Scraper — thin shim invoking the betmgm library.

Reads target_lines from MLB_DB (default: mlb_mm.duckdb; bot overrides via
MLB_SGP_DB_PATH env var). Calls mlb_sgp.betmgm.price_sgps() with the periods
configured via MLB_SGP_PERIODS env var (default: FG,F5). Writes PricedRow
results back to MLB_DB via mlb_sgp.db.upsert_priced_rows.

Same shim contract as scraper_novig_sgp.py: load_target_lines → price_sgps →
clear_source + upsert_priced_rows. The orchestrator's optional
``parallelism``/``fetchers`` knobs (added by the SGP speedup refactor for the
other books) are forward-compatible here — this shim passes neither, so it
keeps the current serial-per-cycle behavior the bot's cadence loop expects.
"""
import os
import sys
from pathlib import Path

# Resolve repo root dynamically (works from worktrees too)
_THIS_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _THIS_DIR.parent
if ".worktrees" in str(_REPO_ROOT):
    _REPO_ROOT = Path(str(_REPO_ROOT).split(".worktrees")[0].rstrip("/"))
_ANSWER_KEYS = _REPO_ROOT / "Answer Keys"

# canonical_match (used by mlb_sgp.betmgm._match_events) lives in Answer Keys.
if str(_ANSWER_KEYS) not in sys.path:
    sys.path.insert(0, str(_ANSWER_KEYS))


def main():
    if str(_REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(_REPO_ROOT))

    from mlb_sgp import db
    from mlb_sgp._shared import load_target_lines

    db_path = str(db.MLB_DB)
    db.ensure_table(db_path)

    targets = load_target_lines(db_path)
    if not targets:
        print(f"  No target lines in {db_path} — nothing to scrape.")
        return 0

    periods_raw = os.environ.get("MLB_SGP_PERIODS", "FG,F5").split(",")
    periods = tuple(p.strip() for p in periods_raw if p.strip())

    from mlb_sgp import betmgm
    print(f"  MGM shim: {len(targets)} target lines, periods={periods}")
    try:
        rows = betmgm.price_sgps(targets, periods=periods, verbose=False)
    except Exception as e:
        # Hard failure (e.g. fixtures/markets fetch blew up): leave the previous
        # cycle's rows in place rather than wiping the source. The downstream
        # fetch_time freshness gate filters anything stale.
        print(f"  MGM shim: price_sgps failed ({e}) — preserving last cycle's rows")
        return 1
    print(f"  MGM shim: priced {len(rows)} rows")

    db.clear_source("betmgm_direct", db_path=db_path)
    db.upsert_priced_rows(rows, db_path=db_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
