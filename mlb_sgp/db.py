#!/usr/bin/env python3
"""
DuckDB helpers for mlb_sgp_odds table.

Stores SGP (Same Game Parlay) odds from multiple sources (Pikkit, DraftKings direct)
in the shared MLB database. Downstream, mlb_correlated_parlay.R can join these
against sample-based fair odds for cross-validation.
"""

import duckdb
import random
import time
from pathlib import Path
from datetime import datetime

# Resolve repo root dynamically — works from main repo or worktrees.
_THIS_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _THIS_DIR.parent
if ".worktrees" in str(_REPO_ROOT):
    _REPO_ROOT = Path(str(_REPO_ROOT).split(".worktrees")[0].rstrip("/"))

MLB_DB = _REPO_ROOT / "Answer Keys" / "mlb.duckdb"

# DuckDB allows only one writer per file. When the four SGP scrapers run in
# parallel they will occasionally collide at connect() — this helper retries
# with exponential backoff + jitter so the contention is invisible to callers.
def _connect_with_retry(db_path, *, read_only=False,
                        max_attempts=10, base_delay=0.1, max_delay=1.5):
    last_err = None
    for attempt in range(max_attempts):
        try:
            return duckdb.connect(db_path, read_only=read_only)
        except duckdb.IOException as e:
            msg = str(e).lower()
            if "lock" not in msg and "in use" not in msg:
                raise
            last_err = e
            if attempt == max_attempts - 1:
                break
            delay = min(base_delay * (2 ** attempt), max_delay)
            time.sleep(delay + random.uniform(0, 0.05))
    raise last_err


CREATE_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS mlb_sgp_odds (
    game_id       VARCHAR,
    combo         VARCHAR,
    period        VARCHAR,
    bookmaker     VARCHAR,
    sgp_decimal   DOUBLE,
    sgp_american  INTEGER,
    fetch_time    TIMESTAMP,
    source        VARCHAR
);
"""


def ensure_table(db_path: str = None):
    """Create the mlb_sgp_odds table if it doesn't exist."""
    db_path = db_path or str(MLB_DB)
    con = _connect_with_retry(db_path)
    try:
        con.execute(CREATE_TABLE_SQL)
    finally:
        con.close()


def clear_source(source: str, db_path: str = None):
    """Delete all rows for a given source (e.g. 'draftkings_direct').

    Called at the start of each scraper run so only fresh prices exist.
    If the scraper fails to price a combo, there's simply no row for it
    rather than a stale price from a previous run.
    """
    db_path = db_path or str(MLB_DB)
    con = _connect_with_retry(db_path)
    try:
        ensure_table(db_path)
        con.execute("DELETE FROM mlb_sgp_odds WHERE source = ?", [source])
    finally:
        con.close()


def upsert_sgp_odds(rows: list[dict], db_path: str = None):
    """
    Insert SGP odds rows, replacing any existing rows for the same
    (game_id, combo, bookmaker, source) combination.

    Each row dict should have keys:
        game_id, combo, period, bookmaker, sgp_decimal, sgp_american, source

    fetch_time is set automatically to now().
    """
    if not rows:
        return

    db_path = db_path or str(MLB_DB)
    con = _connect_with_retry(db_path)
    try:
        ensure_table(db_path)

        # Batch delete stale rows
        keys = [(r["game_id"], r["combo"], r["bookmaker"], r["source"]) for r in rows]
        con.execute("""
            DELETE FROM mlb_sgp_odds
            WHERE (game_id, combo, bookmaker, source) IN (
                SELECT * FROM (VALUES """ +
            ",".join(["(?, ?, ?, ?)"] * len(keys)) +
            """))""",
            [v for k in keys for v in k],
        )

        # Batch insert fresh rows
        now = datetime.now()
        values = []
        for row in rows:
            values.extend([
                row["game_id"], row["combo"], row["period"], row["bookmaker"],
                row["sgp_decimal"], row["sgp_american"], now, row["source"],
            ])

        placeholders = ",".join(["(?, ?, ?, ?, ?, ?, ?, ?)"] * len(rows))
        con.execute(f"""
            INSERT INTO mlb_sgp_odds
                (game_id, combo, period, bookmaker, sgp_decimal, sgp_american, fetch_time, source)
            VALUES {placeholders}
        """, values)
    finally:
        con.close()


def get_sgp_odds(game_id: str = None, max_age_minutes: int = 15, db_path: str = None):
    """
    Read SGP odds from the table, optionally filtered by game_id.

    Returns list of dicts with all columns.
    """
    db_path = db_path or str(MLB_DB)
    con = duckdb.connect(db_path, read_only=True)
    try:
        tables = con.execute("SHOW TABLES").fetchall()
        if not any(t[0] == "mlb_sgp_odds" for t in tables):
            return []

        query = """
            SELECT game_id, combo, period, bookmaker,
                   sgp_decimal, sgp_american, fetch_time, source
            FROM mlb_sgp_odds
            WHERE fetch_time > now() - INTERVAL ? MINUTE
        """
        params = [max_age_minutes]

        if game_id:
            query += " AND game_id = ?"
            params.append(game_id)

        query += " ORDER BY game_id, combo, bookmaker"
        rows = con.execute(query, params).fetchall()
        cols = ["game_id", "combo", "period", "bookmaker",
                "sgp_decimal", "sgp_american", "fetch_time", "source"]
        return [dict(zip(cols, row)) for row in rows]
    finally:
        con.close()
