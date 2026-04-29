"""DuckDB schema + helpers for the Kalshi MLB RFQ bot."""

import time
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from random import random

import duckdb

from kalshi_mlb_rfq.config import DB_PATH

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS combo_cache (
    leg_set_hash        VARCHAR PRIMARY KEY,
    collection_ticker   VARCHAR NOT NULL,
    combo_market_ticker VARCHAR NOT NULL,
    combo_event_ticker  VARCHAR NOT NULL,
    legs_json           VARCHAR NOT NULL,
    game_id             VARCHAR NOT NULL,
    created_at          TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX IF NOT EXISTS idx_combo_cache_game ON combo_cache(game_id);

CREATE TABLE IF NOT EXISTS live_rfqs (
    rfq_id                  VARCHAR PRIMARY KEY,
    combo_market_ticker     VARCHAR NOT NULL,
    leg_set_hash            VARCHAR NOT NULL,
    game_id                 VARCHAR NOT NULL,
    blended_fair_at_submit  DOUBLE,
    kalshi_ref_at_submit    DOUBLE,
    edge_at_submit          DOUBLE,
    status                  VARCHAR NOT NULL,
    submitted_at            TIMESTAMP NOT NULL,
    closed_at               TIMESTAMP,
    cancellation_reason     VARCHAR
);

CREATE TABLE IF NOT EXISTS quote_log (
    quote_id              VARCHAR PRIMARY KEY,
    rfq_id                VARCHAR NOT NULL,
    combo_market_ticker   VARCHAR,
    creator_id            VARCHAR,
    yes_bid_dollars       DOUBLE,
    no_bid_dollars        DOUBLE,
    blended_fair_at_eval  DOUBLE,
    post_fee_ev_pct       DOUBLE,
    decision              VARCHAR NOT NULL,
    reason_detail         VARCHAR,
    observed_at           TIMESTAMP NOT NULL
);

CREATE TABLE IF NOT EXISTS fills (
    fill_id              VARCHAR PRIMARY KEY,
    quote_id             VARCHAR NOT NULL,
    rfq_id               VARCHAR NOT NULL,
    combo_market_ticker  VARCHAR NOT NULL,
    game_id              VARCHAR NOT NULL,
    side                 VARCHAR NOT NULL,
    contracts            DOUBLE NOT NULL,
    price_dollars        DOUBLE NOT NULL,
    fee_dollars          DOUBLE NOT NULL,
    blended_fair_at_fill DOUBLE NOT NULL,
    expected_ev_dollars  DOUBLE NOT NULL,
    filled_at            TIMESTAMP NOT NULL,
    raw_response         VARCHAR
);
CREATE INDEX IF NOT EXISTS idx_fills_game ON fills(game_id);
CREATE INDEX IF NOT EXISTS idx_fills_filled_at ON fills(filled_at);

CREATE TABLE IF NOT EXISTS positions (
    combo_market_ticker  VARCHAR PRIMARY KEY,
    game_id              VARCHAR NOT NULL,
    net_contracts        DOUBLE NOT NULL,
    weighted_price       DOUBLE NOT NULL,
    legs_json            VARCHAR NOT NULL,
    updated_at           TIMESTAMP NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_positions_game ON positions(game_id);

CREATE TABLE IF NOT EXISTS sessions (
    session_id   VARCHAR PRIMARY KEY,
    started_at   TIMESTAMP NOT NULL,
    ended_at     TIMESTAMP,
    pid          INTEGER,
    bot_version  VARCHAR,
    dry_run      BOOLEAN NOT NULL,
    notes        VARCHAR
);

CREATE TABLE IF NOT EXISTS combo_cooldown (
    leg_set_hash  VARCHAR PRIMARY KEY,
    game_id       VARCHAR NOT NULL,
    cooled_until  TIMESTAMP NOT NULL,
    reason        VARCHAR
);

CREATE TABLE IF NOT EXISTS reference_lines (
    rfq_id     VARCHAR PRIMARY KEY,
    lines_json VARCHAR NOT NULL,
    snapped_at TIMESTAMP NOT NULL
);
"""


@contextmanager
def connect(read_only: bool = False, retries: int = 10):
    """DuckDB connection with retry-with-backoff for write contention."""
    last_err = None
    for attempt in range(retries):
        try:
            con = duckdb.connect(str(DB_PATH), read_only=read_only)
            try:
                yield con
            finally:
                con.close()
            return
        except duckdb.IOException as e:
            last_err = e
            time.sleep(0.05 * (2 ** attempt) + random() * 0.05)
    raise last_err


def init_database():
    """Apply schema migrations idempotently."""
    with connect() as con:
        con.execute(SCHEMA_SQL)


def start_session(pid: int, dry_run: bool, version: str) -> str:
    sid = str(uuid.uuid4())
    with connect() as con:
        con.execute(
            "INSERT INTO sessions (session_id, started_at, pid, bot_version, dry_run) "
            "VALUES (?, ?, ?, ?, ?)",
            [sid, datetime.now(timezone.utc), pid, version, dry_run],
        )
    return sid


def end_session(session_id: str, notes: str | None = None):
    with connect() as con:
        con.execute(
            "UPDATE sessions SET ended_at=?, notes=COALESCE(?, notes) "
            "WHERE session_id=?",
            [datetime.now(timezone.utc), notes, session_id],
        )
