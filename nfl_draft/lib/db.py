"""DuckDB connection management, schema, and TZ boundary helpers."""

import duckdb
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

DB_PATH = Path(__file__).resolve().parent.parent / "nfl_draft.duckdb"


@contextmanager
def write_connection() -> Iterator[duckdb.DuckDBPyConnection]:
    """Short-lived write connection. Lock held for milliseconds."""
    con = duckdb.connect(str(DB_PATH))
    try:
        yield con
    finally:
        con.close()


@contextmanager
def read_connection() -> Iterator[duckdb.DuckDBPyConnection]:
    """Short-lived read connection per Dash callback.

    Opened as read/write (no read_only=True) so it can coexist with other
    read/write connections in the same process. DuckDB disallows mixing
    read-only and read/write handles on the same file, and Dash callbacks
    are strictly SELECT-only, so a R/W handle is safe here. Cross-process
    contention with cron writers is bounded by the short cron window.
    """
    con = duckdb.connect(str(DB_PATH))
    try:
        yield con
    finally:
        con.close()


def local_to_utc_iso(local_ts: datetime) -> str:
    """Convert naive local datetime to UTC ISO-8601 string (for Kalshi API)."""
    aware_local = local_ts.astimezone()  # interprets naive as system local
    return aware_local.astimezone(timezone.utc).isoformat()


def utc_iso_to_local(iso_str: str) -> datetime:
    """Convert Kalshi ISO-8601 (with Z or offset) to naive local datetime."""
    # Handle Z suffix (Kalshi style) by replacing with +00:00
    s = iso_str.replace("Z", "+00:00")
    aware_utc = datetime.fromisoformat(s)
    return aware_utc.astimezone().replace(tzinfo=None)


SCHEMA = [
    """
    CREATE TABLE IF NOT EXISTS draft_markets (
        market_id TEXT PRIMARY KEY,
        market_type TEXT NOT NULL,
        subject_player TEXT,
        subject_team TEXT,
        position TEXT,
        pick_number INT,
        range_low INT,
        range_high INT,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS draft_odds (
        market_id TEXT,
        book TEXT,
        american_odds INT,
        implied_prob DOUBLE,
        devig_prob DOUBLE,
        fetched_at TIMESTAMP
    )
    """,
    "CREATE INDEX IF NOT EXISTS idx_draft_odds_market_book_time ON draft_odds (market_id, book, fetched_at DESC)",
    """
    CREATE TABLE IF NOT EXISTS kalshi_trades (
        trade_id TEXT PRIMARY KEY,
        ticker TEXT,
        side TEXT,
        price_cents INT,
        count INT,
        notional_usd DOUBLE,
        traded_at TIMESTAMP,
        fetched_at TIMESTAMP
    )
    """,
    "CREATE INDEX IF NOT EXISTS idx_kalshi_trades_ticker_time ON kalshi_trades (ticker, traded_at DESC)",
    """
    CREATE TABLE IF NOT EXISTS kalshi_poll_state (
        series_ticker TEXT PRIMARY KEY,
        last_traded_at TIMESTAMP,
        last_polled_at TIMESTAMP
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS draft_bets (
        bet_id TEXT PRIMARY KEY,
        market_id TEXT,
        book TEXT,
        side TEXT,
        american_odds INT,
        stake_usd DOUBLE,
        taken_at TIMESTAMP,
        note TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS players (
        canonical_name TEXT PRIMARY KEY,
        position TEXT,
        college TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS player_aliases (
        alias TEXT PRIMARY KEY,
        canonical_name TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS teams (
        canonical_abbr TEXT PRIMARY KEY,
        full_name TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS team_aliases (
        alias TEXT PRIMARY KEY,
        canonical_abbr TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS market_map (
        book TEXT,
        book_label TEXT,
        book_subject TEXT,
        market_id TEXT,
        PRIMARY KEY (book, book_label, book_subject)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS draft_odds_unmapped (
        book TEXT,
        book_label TEXT,
        book_subject TEXT,
        american_odds INT,
        fetched_at TIMESTAMP,
        reason TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS draft_odds_unmapped_players (
        book TEXT,
        book_player_name TEXT,
        american_odds INT,
        fetched_at TIMESTAMP
    )
    """,
]


def init_schema() -> None:
    """Create all tables (idempotent)."""
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    with write_connection() as con:
        for stmt in SCHEMA:
            con.execute(stmt)
