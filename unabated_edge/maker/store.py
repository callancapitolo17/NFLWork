"""Maker sibling DB (unabated_edge_maker.duckdb): quote decisions, fills,
ledger snapshots, settlements. Separate file so maker state can be reset
without touching the capture history (three-DB discipline)."""
import json
import logging

from unabated_edge import config
from unabated_edge.storage import connect

log = logging.getLogger("unabated_edge")


def init():
    with connect(config.MAKER_DB_PATH) as c:
        c.execute("""CREATE TABLE IF NOT EXISTS maker_quotes(
            ts TIMESTAMPTZ, sport VARCHAR, event_id BIGINT, market_ticker VARCHAR,
            side VARCHAR, action VARCHAR, price DOUBLE, size DOUBLE,
            fair DOUBLE, margin DOUBLE, alt BOOLEAN, reason VARCHAR, order_id VARCHAR)""")
        c.execute("""CREATE TABLE IF NOT EXISTS maker_fills(
            trade_id VARCHAR PRIMARY KEY, ts TIMESTAMPTZ, sport VARCHAR, event_id BIGINT,
            order_id VARCHAR, market_ticker VARCHAR, side VARCHAR, price DOUBLE,
            contracts DOUBLE, fee DOUBLE, ledger_worst_after DOUBLE)""")
        c.execute("""CREATE TABLE IF NOT EXISTS ledger_snapshots(
            ts TIMESTAMPTZ, sport VARCHAR, event_id BIGINT, worst_case DOUBLE,
            pnl_grid JSON, quotes_live INTEGER)""")
        c.execute("""CREATE TABLE IF NOT EXISTS maker_pnl(
            market_ticker VARCHAR PRIMARY KEY, ts TIMESTAMPTZ, sport VARCHAR,
            settled_pnl DOUBLE)""")
        # Single-row latch for the daily-loss halt (finding #75/F6): persists
        # across a restart so a losing run can't silently un-halt just
        # because the process bounced. id is always 1 -- one latch, upserted.
        c.execute("""CREATE TABLE IF NOT EXISTS maker_halt_state(
            id INTEGER PRIMARY KEY, trading_day DATE, halted BOOLEAN)""")


def log_quote(ts, sport, event_id, ticker, side, action, price, size, fair, margin, alt, reason, order_id):
    with connect(config.MAKER_DB_PATH) as c:
        c.execute("INSERT INTO maker_quotes VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                  [ts, sport, event_id, ticker, side, action, price, size, fair, margin, alt, reason, order_id])


def log_fill(ts, sport, event_id, order_id, ticker, side, price, contracts, fee, worst_after, trade_id):
    with connect(config.MAKER_DB_PATH) as c:
        c.execute("INSERT OR IGNORE INTO maker_fills VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                  [trade_id, ts, sport, event_id, order_id, ticker, side, price, contracts, fee, worst_after])


def log_ledger(ts, sport, event_id, worst_case, grid, quotes_live):
    with connect(config.MAKER_DB_PATH) as c:
        c.execute("INSERT INTO ledger_snapshots VALUES (?,?,?,?,?,?)",
                  [ts, sport, event_id, worst_case, json.dumps(grid), quotes_live])


def log_settlement(ts, sport, ticker, pnl):
    with connect(config.MAKER_DB_PATH) as c:
        c.execute("INSERT OR IGNORE INTO maker_pnl VALUES (?,?,?,?)", [ticker, ts, sport, pnl])


def set_halt_latch(trading_day, halted):
    with connect(config.MAKER_DB_PATH) as c:
        c.execute("""INSERT INTO maker_halt_state VALUES (1, ?, ?)
                     ON CONFLICT (id) DO UPDATE SET
                     trading_day = excluded.trading_day, halted = excluded.halted""",
                  [trading_day, halted])


def get_halt_latch():
    """(trading_day, halted) for the persisted latch, or (None, False) if
    never written (fresh DB / first run)."""
    with connect(config.MAKER_DB_PATH) as c:
        row = c.execute("SELECT trading_day, halted FROM maker_halt_state WHERE id = 1").fetchone()
    return (row[0], row[1]) if row else (None, False)
