"""
Legacy Kalshi dashboard query helpers + legacy writer shim.

Historically this module opened raw ``duckdb.connect()`` handles. As of the
NFL Draft Portal work, all DuckDB access goes through the short-lived
context-managed helpers in ``nfl_draft.lib.db`` so cross-process lock
contention is bounded.

- The dashboard read helpers below now use ``read_connection()`` under the hood
  so every callback opens and closes its own short-lived handle. This is the
  path that matters for lock-contention safety: the Dash process reads while
  the cron writer holds the lock.
- ``get_connection()`` remains as a thin wrapper because the legacy writer
  pipeline (``fetcher.py`` / ``consensus.py`` / ``edge_detector.py``) still
  uses the older ``con = get_connection(); ... ; con.close()`` pattern. Those
  files live in the lint allowlist so the new-code discipline still holds.
"""

import duckdb

from nfl_draft.lib import db as _nfl_db
from nfl_draft.lib.db import read_connection


# The canonical DB location lives in nfl_draft.lib.db. We expose it here as a
# convenience for legacy callers (fetcher.py / consensus.py / edge_detector.py)
# that used to import DB_PATH from this module. Tests that want to redirect
# the database location should monkeypatch ``nfl_draft.lib.db.DB_PATH``, not
# ``kalshi_draft.db.DB_PATH`` — this constant is only captured at import time.
DB_PATH = _nfl_db.DB_PATH


def get_connection(read_only: bool = False):
    """Thin shim over ``duckdb.connect`` for legacy writer callers.

    Looks up ``nfl_draft.lib.db.DB_PATH`` at call time so tests that
    monkeypatch that module's ``DB_PATH`` take effect here automatically.

    New code should use ``nfl_draft.lib.db.write_connection`` /
    ``read_connection`` context managers instead. This function exists only
    to keep the legacy ``kalshi_draft/fetcher.py``, ``consensus.py``, and
    ``edge_detector.py`` callsites working without touching those files in
    the NFL Draft Portal work.
    """
    return duckdb.connect(str(_nfl_db.DB_PATH), read_only=read_only)


def init_schema():
    """DEPRECATED: schema now managed by nfl_draft/lib/db.py.init_schema()."""
    raise RuntimeError("Use nfl_draft.lib.db.init_schema() instead")


def get_latest_odds():
    """Get the most recent odds snapshot."""
    try:
        with read_connection() as con:
            return con.execute("""
                SELECT * FROM kalshi_odds
                WHERE fetch_time = (SELECT MAX(fetch_time) FROM kalshi_odds)
                ORDER BY series_ticker, last_price DESC
            """).fetchdf()
    except Exception:
        return None


def get_price_history(tickers=None, days=30):
    """Get price history for time-series charts."""
    try:
        with read_connection() as con:
            if tickers:
                placeholders = ",".join(["?" for _ in tickers])
                return con.execute(f"""
                    SELECT fetch_time, ticker, candidate, series_ticker,
                           last_price, yes_bid, yes_ask, volume
                    FROM kalshi_odds
                    WHERE ticker IN ({placeholders})
                      AND fetch_time >= CURRENT_TIMESTAMP - INTERVAL '{days} days'
                    ORDER BY fetch_time
                """, tickers).fetchdf()
            return con.execute(f"""
                SELECT fetch_time, ticker, candidate, series_ticker,
                       last_price, yes_bid, yes_ask, volume
                FROM kalshi_odds
                WHERE fetch_time >= CURRENT_TIMESTAMP - INTERVAL '{days} days'
                ORDER BY fetch_time
            """).fetchdf()
    except Exception:
        return None


def get_all_series():
    """List all discovered draft series."""
    try:
        with read_connection() as con:
            return con.execute(
                "SELECT * FROM draft_series ORDER BY series_ticker"
            ).fetchdf()
    except Exception:
        return None


def get_portfolio():
    """Get latest positions and resting orders."""
    try:
        with read_connection() as con:
            positions = con.execute("""
                SELECT p.*, m.title as market_title, m.subtitle as market_name
                FROM positions p
                LEFT JOIN market_info m ON p.ticker = m.ticker
                WHERE p.fetch_time = (SELECT MAX(fetch_time) FROM positions)
                ORDER BY ABS(p.market_exposure) DESC
            """).fetchdf()

            orders = con.execute("""
                SELECT r.*, m.title as market_title, m.subtitle as market_name
                FROM resting_orders r
                LEFT JOIN market_info m ON r.ticker = m.ticker
                WHERE r.fetch_time = (SELECT MAX(fetch_time) FROM resting_orders)
                ORDER BY r.ticker
            """).fetchdf()

            return positions, orders
    except Exception:
        import pandas as pd
        return pd.DataFrame(), pd.DataFrame()


def get_position_changes():
    """Compare current vs previous position snapshots."""
    try:
        with read_connection() as con:
            fetch_times = con.execute("""
                SELECT DISTINCT fetch_time FROM positions
                ORDER BY fetch_time DESC LIMIT 2
            """).fetchall()

            if len(fetch_times) < 2:
                return None

            current_time = fetch_times[0][0]
            prev_time = fetch_times[1][0]

            return con.execute("""
                WITH current AS (
                    SELECT ticker, position, market_exposure
                    FROM positions WHERE fetch_time = ?
                ),
                previous AS (
                    SELECT ticker, position, market_exposure
                    FROM positions WHERE fetch_time = ?
                )
                SELECT
                    COALESCE(c.ticker, p.ticker) as ticker,
                    COALESCE(m.subtitle, COALESCE(c.ticker, p.ticker)) as market_name,
                    COALESCE(c.position, 0) as current_pos,
                    COALESCE(p.position, 0) as prev_pos,
                    COALESCE(c.position, 0) - COALESCE(p.position, 0) as pos_change,
                    COALESCE(c.market_exposure, 0) - COALESCE(p.market_exposure, 0) as exp_change,
                    CASE
                        WHEN c.ticker IS NULL THEN 'CLOSED'
                        WHEN p.ticker IS NULL THEN 'NEW'
                        ELSE 'CHANGED'
                    END as change_type
                FROM current c
                FULL OUTER JOIN previous p ON c.ticker = p.ticker
                LEFT JOIN market_info m ON COALESCE(c.ticker, p.ticker) = m.ticker
                WHERE c.position != p.position
                   OR c.market_exposure != p.market_exposure
                   OR c.ticker IS NULL
                   OR p.ticker IS NULL
                ORDER BY ABS(COALESCE(c.position, 0) - COALESCE(p.position, 0)) DESC
            """, [current_time, prev_time]).fetchdf()
    except Exception:
        return None


def get_latest_edges():
    """Get the most recently computed edges."""
    try:
        with read_connection() as con:
            return con.execute("""
                SELECT * FROM detected_edges
                WHERE fetch_time = (SELECT MAX(fetch_time) FROM detected_edges)
                ORDER BY ABS(implied_edge) DESC
            """).fetchdf()
    except Exception:
        return None


def get_latest_consensus():
    """Get the most recent consensus board."""
    try:
        with read_connection() as con:
            return con.execute("""
                SELECT * FROM consensus_board
                WHERE fetch_time = (SELECT MAX(fetch_time) FROM consensus_board)
                ORDER BY rank
            """).fetchdf()
    except Exception:
        return None


def get_snapshot_count():
    """Get number of historical snapshots."""
    try:
        with read_connection() as con:
            return con.execute("""
                SELECT COUNT(DISTINCT fetch_time) as snapshots,
                       MIN(fetch_time) as first_fetch,
                       MAX(fetch_time) as last_fetch
                FROM kalshi_odds
            """).fetchone()
    except Exception:
        return (0, None, None)
