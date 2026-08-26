"""DuckDB persistence for the YRFI maker.

Owns kalshi_rfi/kalshi_rfi.duckdb (this bot is its only writer). Tables:
  snapshots   — one row per (cycle x game): fairs, Kalshi touch, decision.
                This IS the research dataset; shadow mode fills it for free.
  orders      — every place/cancel/replace attempt (shadow included).
  fills       — our executed trades (trade_id-deduped upsert).
  settlements — realized P&L per settled ticker.
All writes append except fills (INSERT OR IGNORE on trade_id).
"""
import json
import logging
from datetime import datetime, timezone

import duckdb

from kalshi_rfi import config

log = logging.getLogger(__name__)

_SCHEMA = """
CREATE TABLE IF NOT EXISTS snapshots (
    ts               TIMESTAMPTZ,
    ticker           VARCHAR,
    home_team        VARCHAR,
    away_team        VARCHAR,
    commence_utc     TIMESTAMP,
    sec_to_start     DOUBLE,
    n_books          INTEGER,
    book_fairs       VARCHAR,      -- JSON {book: devigged P(YRFI)}
    sigma_z          DOUBLE,
    consensus_yes    DOUBLE,
    kalshi_yes_bid   INTEGER,
    kalshi_yes_ask   INTEGER,
    action           VARCHAR,      -- quote | no_quote
    reason           VARCHAR,
    our_price_cents  INTEGER,
    our_count        INTEGER
);
CREATE TABLE IF NOT EXISTS orders (
    ts               TIMESTAMPTZ,
    ticker           VARCHAR,
    order_id         VARCHAR,
    action           VARCHAR,      -- place | cancel
    price_cents      INTEGER,
    count            INTEGER,
    mode             VARCHAR,      -- shadow | live
    reason           VARCHAR
);
CREATE TABLE IF NOT EXISTS fills (
    ts               TIMESTAMPTZ,
    trade_id         VARCHAR PRIMARY KEY,
    order_id         VARCHAR,
    ticker           VARCHAR,
    price_cents      INTEGER,
    count            DOUBLE,
    fee_usd          DOUBLE
);
CREATE TABLE IF NOT EXISTS settlements (
    ts               TIMESTAMPTZ,
    ticker           VARCHAR,
    pnl_usd          DOUBLE
);
"""


def connect(db_path: str | None = None) -> duckdb.DuckDBPyConnection:
    con = duckdb.connect(str(db_path or config.DB_PATH))
    con.execute(_SCHEMA)
    return con


def _now():
    return datetime.now(timezone.utc)


def log_snapshot(con, game, fair_result, decision, sec_to_start: float):
    """fair_result may be None (no consensus this cycle)."""
    try:
        con.execute(
            "INSERT INTO snapshots (ts, ticker, home_team, away_team, "
            "commence_utc, sec_to_start, n_books, book_fairs, sigma_z, "
            "consensus_yes, kalshi_yes_bid, kalshi_yes_ask, action, reason, "
            "our_price_cents, our_count) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [_now(), game.ticker, game.home_team, game.away_team,
             game.commence_utc, sec_to_start,
             len(fair_result.book_fairs) if fair_result else 0,
             json.dumps({b: round(f, 4)
                         for b, f in fair_result.book_fairs.items()})
             if fair_result else None,
             fair_result.sigma_z if fair_result else None,
             fair_result.consensus_yes if fair_result else None,
             game.yes_bid_cents, game.yes_ask_cents,
             decision.action, decision.reason,
             decision.price_cents, decision.count])
    except duckdb.Error as e:
        log.error("storage: snapshot insert failed for %s: %s", game.ticker, e)


def log_order(con, ticker: str, order_id: str | None, action: str,
              price_cents: int | None, count: int | None, mode: str,
              reason: str):
    try:
        con.execute(
            "INSERT INTO orders (ts, ticker, order_id, action, price_cents, "
            "count, mode, reason) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            [_now(), ticker, order_id, action, price_cents, count, mode,
             reason])
    except duckdb.Error as e:
        log.error("storage: order insert failed for %s: %s", ticker, e)


def log_fill(con, trade_id: str, order_id: str, ticker: str,
             price_cents: int, count: float, fee_usd: float):
    try:
        con.execute(
            "INSERT OR IGNORE INTO fills (ts, trade_id, order_id, ticker, "
            "price_cents, count, fee_usd) VALUES (?, ?, ?, ?, ?, ?, ?)",
            [_now(), trade_id, order_id, ticker, price_cents, count, fee_usd])
    except duckdb.Error as e:
        log.error("storage: fill insert failed for %s: %s", ticker, e)


def load_recent_fills(con, hours: float = 48.0) -> list[tuple]:
    """(trade_id, ticker, price_cents, count, ts) for recent fills — startup
    state hydration, so a restart doesn't forget filled exposure (caps) or
    orphan unsettled fills from settlement matching. 48h comfortably covers
    an ET trading day plus settlement lag."""
    try:
        return con.execute(
            "SELECT trade_id, ticker, price_cents, count, ts FROM fills "
            "WHERE ts >= now() - (? * INTERVAL 1 HOUR)", [hours]).fetchall()
    except duckdb.Error as e:
        log.error("storage: recent-fills load failed: %s", e)
        return []


def load_recent_placed_orders(con, hours: float = 48.0) -> list[tuple]:
    """(order_id, ticker) for recent live 'place' rows — startup rebuild of
    the order→ticker map, so fills on a PREVIOUS run's orders still
    attribute (poll_fills skips any order_id it doesn't know)."""
    try:
        return con.execute(
            "SELECT DISTINCT order_id, ticker FROM orders "
            "WHERE action = 'place' AND mode = 'live' "
            "AND order_id IS NOT NULL "
            "AND ts >= now() - (? * INTERVAL 1 HOUR)", [hours]).fetchall()
    except duckdb.Error as e:
        log.error("storage: recent-orders load failed: %s", e)
        return []


def load_settled_tickers(con) -> set:
    try:
        return {r[0] for r in con.execute(
            "SELECT DISTINCT ticker FROM settlements").fetchall()}
    except duckdb.Error as e:
        log.error("storage: settled-tickers load failed: %s", e)
        return set()


def log_settlement(con, ticker: str, pnl_usd: float):
    try:
        con.execute(
            "INSERT INTO settlements (ts, ticker, pnl_usd) VALUES (?, ?, ?)",
            [_now(), ticker, pnl_usd])
    except duckdb.Error as e:
        log.error("storage: settlement insert failed for %s: %s", ticker, e)
