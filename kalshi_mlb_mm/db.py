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
    closed_at           TIMESTAMPTZ,
    -- R-1 (issue #22): worst-case dollars at risk if this quote fills, frozen
    -- at quote time (main._worst_fill_exposure_usd). Summed over open quotes
    -- by the per-game and daily cap gates. NULL (pre-migration rows) is
    -- counted at the per-fill cap by readers — never treated as zero.
    worst_exposure_usd  DOUBLE,
    -- #17 singles-veto baseline: raw Kalshi odds of every leg at quote time,
    -- JSON {leg market_ticker: {yes_bid, yes_ask}} in dollars. The confirm
    -- tick voids the accept if any leg moved vs this snapshot. NULL
    -- (pre-migration rows) fails CLOSED — the accept voids.
    leg_prices_json     VARCHAR
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
-- Open-quote exposure ledger (R-1, issue #22): maps each live quote to EVERY
-- game its combo touches, written in the same transaction as the live_quotes
-- insert. The per-game cap gate sums OPEN quotes' worst_exposure_usd through
-- this map so a cross-game quote counts against each game it touches — same
-- attribution rule as fill_games. Rows for closed quotes are inert (readers
-- join on live_quotes.status='open').
CREATE TABLE IF NOT EXISTS quote_games (
    quote_id VARCHAR NOT NULL,
    game_id  VARCHAR NOT NULL,
    PRIMARY KEY (quote_id, game_id)
);
-- Settlement audit (issue #12 / audit M-1): one row per settled combo market,
-- written by settlement.settlement_sweep_tick in the same transaction as the
-- fills.realized_pnl updates. raw_payload preserves Kalshi's /markets/{ticker}
-- response verbatim — the verification record for the settlement response
-- shape (status/result field vocabulary is unverified until real data lands).
CREATE TABLE IF NOT EXISTS settlements (
    combo_market_ticker VARCHAR PRIMARY KEY,
    result              VARCHAR NOT NULL,
    settled_at          TIMESTAMPTZ,
    raw_payload         VARCHAR,
    recorded_at         TIMESTAMPTZ NOT NULL
);
-- Expired-quote outcome labels (expiry_outcome.expiry_outcome_sweep_tick):
-- one row per terminal unfilled quote, written once after the grace window
-- elapses. Doubles as the dedupe marker (PK quote_id) and the queryable
-- label store for the margin-tuning ratio: competitor_traded_in_window
-- (a competing maker won the RFQ) vs no_trade (the creator never executed).
-- window_start/window_end are the quote's UTC resting bounds.
CREATE TABLE IF NOT EXISTS quote_expiry_outcomes (
    quote_id            VARCHAR PRIMARY KEY,
    combo_market_ticker VARCHAR NOT NULL,
    quote_status        VARCHAR NOT NULL,
    label               VARCHAR NOT NULL,
    n_trades_in_window  INTEGER NOT NULL,
    first_trade_after_close_sec DOUBLE,
    window_start        TIMESTAMPTZ NOT NULL,
    window_end          TIMESTAMPTZ NOT NULL,
    recorded_at         TIMESTAMPTZ NOT NULL
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
ALTER TABLE live_quotes ADD COLUMN IF NOT EXISTS worst_exposure_usd DOUBLE;
ALTER TABLE live_quotes ADD COLUMN IF NOT EXISTS leg_prices_json VARCHAR;
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
