"""SGP scrape orchestration for the Kalshi MLB RFQ bot.

Owns the cycle that:
  1. Enumerates Kalshi MVE markets per open MLB game
  2. Writes target lines (one row per game x (spread, total)) to bot DB
  3. Spawns the 4 scraper subprocesses with MLB_SGP_DB_PATH redirect
  4. Reads back priced SGP odds into the bot's _SGP_ODDS_CACHE

This module is invoked from main_loop on the SGP cadence tick.
"""
from __future__ import annotations
from datetime import datetime, timezone
from pathlib import Path

import duckdb


def should_scrape(last_fetch_time: datetime | None,
                   now: datetime,
                   min_interval_sec: int) -> bool:
    """True if we should scrape this tick. Guards against tight cadences
    after crash-recovery restarts that hit an already-fresh DB.

    Both `last_fetch_time` and `now` are normalized to UTC when naive."""
    if last_fetch_time is None:
        return True
    if last_fetch_time.tzinfo is None:
        last_fetch_time = last_fetch_time.replace(tzinfo=timezone.utc)
    if now.tzinfo is None:
        now = now.replace(tzinfo=timezone.utc)
    age = (now - last_fetch_time).total_seconds()
    return age > min_interval_sec


def latest_sgp_fetch_time(bot_market_db: str) -> datetime | None:
    """Read MAX(fetch_time) from mlb_sgp_odds in bot_market_db.
    Returns None for missing DB / missing table / empty table."""
    if not Path(bot_market_db).exists():
        return None
    con = duckdb.connect(bot_market_db, read_only=True)
    try:
        tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
        if "mlb_sgp_odds" not in tables:
            return None
        row = con.execute("SELECT MAX(fetch_time) FROM mlb_sgp_odds").fetchone()
        return row[0] if row and row[0] is not None else None
    finally:
        con.close()
