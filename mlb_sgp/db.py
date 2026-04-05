#!/usr/bin/env python3
"""
DuckDB helpers for mlb_sgp_odds table.

Stores SGP (Same Game Parlay) odds from multiple sources (Pikkit, DraftKings direct)
in the shared MLB database. Downstream, mlb_correlated_parlay.R can join these
against sample-based fair odds for cross-validation.
"""

import duckdb
from pathlib import Path
from datetime import datetime

# Default path to the shared MLB database.
# Use absolute path so this works from worktrees too.
_REPO_ROOT = Path("/Users/callancapitolo/NFLWork")
MLB_DB = _REPO_ROOT / "Answer Keys" / "mlb.duckdb"

CREATE_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS mlb_sgp_odds (
    game_id       VARCHAR,    -- Odds API event ID (matches mlb_parlay_opportunities)
    combo         VARCHAR,    -- e.g., "Home Spread + Over" (matches mlb_parlay_opportunities)
    period        VARCHAR,    -- "FG" or "F5"
    bookmaker     VARCHAR,    -- "draftkings", "fanduel", "prophetx", "novig"
    sgp_decimal   DOUBLE,
    sgp_american  INTEGER,
    fetch_time    TIMESTAMP,
    source        VARCHAR     -- "pikkit" or "draftkings_direct"
);
"""


def ensure_table(db_path: str = None):
    """Create the mlb_sgp_odds table if it doesn't exist."""
    db_path = db_path or str(MLB_DB)
    con = duckdb.connect(db_path)
    try:
        con.execute(CREATE_TABLE_SQL)
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
    con = duckdb.connect(db_path)
    try:
        ensure_table(db_path)

        # Delete stale rows for these combos so we always have fresh prices
        for row in rows:
            con.execute("""
                DELETE FROM mlb_sgp_odds
                WHERE game_id = ? AND combo = ? AND bookmaker = ? AND source = ?
            """, [row["game_id"], row["combo"], row["bookmaker"], row["source"]])

        # Insert fresh rows
        now = datetime.now()
        for row in rows:
            con.execute("""
                INSERT INTO mlb_sgp_odds
                    (game_id, combo, period, bookmaker, sgp_decimal, sgp_american, fetch_time, source)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, [
                row["game_id"],
                row["combo"],
                row["period"],
                row["bookmaker"],
                row["sgp_decimal"],
                row["sgp_american"],
                now,
                row["source"],
            ])
    finally:
        con.close()


def get_sgp_odds(game_id: str = None, max_age_minutes: int = 15, db_path: str = None):
    """
    Read SGP odds from the table, optionally filtered by game_id.

    Args:
        game_id: Filter to a specific game (None = all games)
        max_age_minutes: Only return rows fetched within this window
        db_path: Override database path

    Returns:
        List of dicts with all columns
    """
    db_path = db_path or str(MLB_DB)
    con = duckdb.connect(db_path, read_only=True)
    try:
        # Table might not exist yet on first run
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
        result = con.execute(query, params).fetchdf()
        return result.to_dict("records")
    finally:
        con.close()
