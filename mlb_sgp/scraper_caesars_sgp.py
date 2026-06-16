#!/usr/bin/env python3
"""Caesars MLB SGP Scraper — thin shim invoking the caesars library.

Reads target_lines from MLB_DB (bot overrides via MLB_SGP_DB_PATH). Calls
mlb_sgp.caesars.price_sgps() for the periods in MLB_SGP_PERIODS (default
FG,F5). Writes PricedRow results back via mlb_sgp.db.upsert_priced_rows.

Same shim contract as the other books: load_target_lines → price_sgps →
clear_source + upsert_priced_rows. Caesars needs a one-time AWS-WAF token mint
(headless Playwright, cached ~5 min by CaesarsClient); cycles where the token
can't be validated produce no rows rather than bad data.
"""
import os
import sys
from pathlib import Path

_THIS_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _THIS_DIR.parent
if ".worktrees" in str(_REPO_ROOT):
    _REPO_ROOT = Path(str(_REPO_ROOT).split(".worktrees")[0].rstrip("/"))
_ANSWER_KEYS = _REPO_ROOT / "Answer Keys"
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

    from mlb_sgp import caesars
    print(f"  CZR shim: {len(targets)} target lines, periods={periods}")
    rows = caesars.price_sgps(targets, periods=periods, verbose=False)
    print(f"  CZR shim: priced {len(rows)} rows")

    db.clear_source("caesars_direct", db_path=db_path)
    db.upsert_priced_rows(rows, db_path=db_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
