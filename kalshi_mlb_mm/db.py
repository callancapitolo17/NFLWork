"""DuckDB schema + helpers for the maker bot."""
import time
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from random import random

import duckdb

from kalshi_mlb_mm.config import DB_PATH

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS seen_rfqs (
    rfq_id          VARCHAR PRIMARY KEY,
    market_ticker   VARCHAR,
    in_scope        BOOLEAN,
    game_id         VARCHAR,
    legs_json       VARCHAR,
    first_seen_at   TIMESTAMPTZ NOT NULL,
    last_decision   VARCHAR,
    creator_id      VARCHAR
);
CREATE TABLE IF NOT EXISTS live_quotes (
    quote_id            VARCHAR PRIMARY KEY,
    rfq_id              VARCHAR NOT NULL,
    combo_market_ticker VARCHAR NOT NULL,
    game_id             VARCHAR NOT NULL,
    yes_bid             DOUBLE,
    no_bid              DOUBLE,
    model_fair          DOUBLE,
    book_fair           DOUBLE,
    blended_fair        DOUBLE,
    status              VARCHAR NOT NULL,
    submitted_at        TIMESTAMPTZ NOT NULL,
    closed_at           TIMESTAMPTZ
);
CREATE TABLE IF NOT EXISTS quote_decisions (
    decision_id   VARCHAR PRIMARY KEY,
    rfq_id        VARCHAR,
    quote_id      VARCHAR,
    combo_market_ticker VARCHAR,
    game_id       VARCHAR,
    decision      VARCHAR NOT NULL,
    reason        VARCHAR,
    model_fair    DOUBLE,
    book_fair     DOUBLE,
    blended_fair  DOUBLE,
    yes_bid       DOUBLE,
    no_bid        DOUBLE,
    observed_at   TIMESTAMPTZ NOT NULL
);
CREATE TABLE IF NOT EXISTS fills (
    fill_id              VARCHAR PRIMARY KEY,
    quote_id             VARCHAR NOT NULL,
    rfq_id               VARCHAR NOT NULL,
    combo_market_ticker  VARCHAR NOT NULL,
    game_id              VARCHAR NOT NULL,
    side_held            VARCHAR NOT NULL,
    contracts            DOUBLE NOT NULL,
    price                DOUBLE NOT NULL,
    fee                  DOUBLE NOT NULL,
    model_fair_at_quote  DOUBLE,
    book_fair_at_quote   DOUBLE,
    blended_fair_at_quote DOUBLE,
    fair_at_confirm      DOUBLE,
    realized_pnl         DOUBLE,
    filled_at            TIMESTAMPTZ NOT NULL,
    reconciled           BOOLEAN NOT NULL DEFAULT FALSE
);
CREATE TABLE IF NOT EXISTS positions (
    combo_market_ticker VARCHAR NOT NULL,
    side                VARCHAR NOT NULL,
    game_id             VARCHAR NOT NULL,
    net_contracts       DOUBLE NOT NULL,
    weighted_price      DOUBLE NOT NULL,
    updated_at          TIMESTAMPTZ NOT NULL,
    PRIMARY KEY (combo_market_ticker, side)
);
CREATE TABLE IF NOT EXISTS sessions (
    session_id VARCHAR PRIMARY KEY,
    started_at TIMESTAMPTZ NOT NULL,
    ended_at   TIMESTAMPTZ,
    pid        INTEGER,
    dry_run    BOOLEAN NOT NULL
);
CREATE TABLE IF NOT EXISTS combo_cooldown (
    combo_market_ticker VARCHAR PRIMARY KEY,
    cooled_until        TIMESTAMPTZ NOT NULL
);
-- Per-game exposure ledger: maps each fill to EVERY game its combo touches
-- (one row per game). A cross-game combo lands one row per game so the
-- per-game exposure cap counts its full stake against each game (correlated
-- risk). `fills` stays one-row-per-combo, so the daily cap and P&L never
-- double-count. See main._today_fills_by_game / _fill_game_ids.
CREATE TABLE IF NOT EXISTS fill_games (
    fill_id  VARCHAR NOT NULL,
    game_id  VARCHAR NOT NULL,
    PRIMARY KEY (fill_id, game_id)
);
CREATE INDEX IF NOT EXISTS idx_quote_decisions_observed_at
    ON quote_decisions(observed_at);
CREATE INDEX IF NOT EXISTS idx_fills_reconciled
    ON fills(reconciled);
"""

# Idempotent column-adds for DBs created before the H4/N5 hardening pass.
# DuckDB accepts "ADD COLUMN IF NOT EXISTS" — safe to run on every startup.
#
# O-2 (issue #18): TIMESTAMPTZ migration. Pre-migration DBs stored naive-LOCAL
# timestamps (DuckDB converts an aware bind into the session TimeZone and
# drops the offset on a naive TIMESTAMP column). ALTER ... SET DATA TYPE
# TIMESTAMPTZ reinterprets the naive value in the session TimeZone — i.e.
# local — so the true instant is preserved (verified empirically). Re-running
# on an already-TIMESTAMPTZ column is a no-op cast, so this stays idempotent.
MIGRATE_SQL = """
ALTER TABLE fills ADD COLUMN IF NOT EXISTS reconciled BOOLEAN DEFAULT FALSE;
ALTER TABLE seen_rfqs ADD COLUMN IF NOT EXISTS creator_id VARCHAR;
DROP INDEX IF EXISTS idx_quote_decisions_observed_at;
DROP INDEX IF EXISTS idx_fills_reconciled;
ALTER TABLE seen_rfqs ALTER COLUMN first_seen_at SET DATA TYPE TIMESTAMPTZ;
ALTER TABLE live_quotes ALTER COLUMN submitted_at SET DATA TYPE TIMESTAMPTZ;
ALTER TABLE live_quotes ALTER COLUMN closed_at SET DATA TYPE TIMESTAMPTZ;
ALTER TABLE quote_decisions ALTER COLUMN observed_at SET DATA TYPE TIMESTAMPTZ;
ALTER TABLE fills ALTER COLUMN filled_at SET DATA TYPE TIMESTAMPTZ;
ALTER TABLE positions ALTER COLUMN updated_at SET DATA TYPE TIMESTAMPTZ;
ALTER TABLE sessions ALTER COLUMN started_at SET DATA TYPE TIMESTAMPTZ;
ALTER TABLE sessions ALTER COLUMN ended_at SET DATA TYPE TIMESTAMPTZ;
ALTER TABLE combo_cooldown ALTER COLUMN cooled_until SET DATA TYPE TIMESTAMPTZ;
CREATE INDEX IF NOT EXISTS idx_quote_decisions_observed_at
    ON quote_decisions(observed_at);
CREATE INDEX IF NOT EXISTS idx_fills_reconciled
    ON fills(reconciled);
"""


@contextmanager
def connect(read_only: bool = False, retries: int = 10):
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
    with connect() as con:
        con.execute(SCHEMA_SQL)
        # Run idempotent migrations for older DBs.
        for stmt in MIGRATE_SQL.strip().split(";"):
            s = stmt.strip()
            if s:
                con.execute(s)


def start_session(pid: int, dry_run: bool) -> str:
    sid = str(uuid.uuid4())
    with connect() as con:
        con.execute(
            "INSERT INTO sessions (session_id, started_at, pid, dry_run) VALUES (?,?,?,?)",
            [sid, datetime.now(timezone.utc), pid, dry_run],
        )
    return sid


def end_session(session_id: str):
    with connect() as con:
        con.execute("UPDATE sessions SET ended_at=? WHERE session_id=?",
                    [datetime.now(timezone.utc), session_id])
