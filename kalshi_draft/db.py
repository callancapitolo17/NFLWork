"""
DuckDB schema and query helpers for NFL Draft dashboard.
"""

import duckdb
from pathlib import Path
from datetime import datetime, timezone

DB_PATH = Path(__file__).resolve().parent.parent / "nfl_draft" / "nfl_draft.duckdb"


def get_connection(read_only=False):
    """Get a DuckDB connection."""
    return duckdb.connect(str(DB_PATH), read_only=read_only)


def init_schema():
    """DEPRECATED: schema now managed by nfl_draft/lib/db.py.init_schema()."""
    raise RuntimeError("Use nfl_draft.lib.db.init_schema() instead")


def get_latest_odds():
    """Get the most recent odds snapshot."""
    con = get_connection(read_only=True)
    try:
        df = con.execute("""
            SELECT * FROM kalshi_odds
            WHERE fetch_time = (SELECT MAX(fetch_time) FROM kalshi_odds)
            ORDER BY series_ticker, last_price DESC
        """).fetchdf()
        return df
    except Exception:
        return None
    finally:
        con.close()


def get_price_history(tickers=None, days=30):
    """Get price history for time-series charts."""
    con = get_connection(read_only=True)
    try:
        if tickers:
            placeholders = ",".join(["?" for _ in tickers])
            df = con.execute(f"""
                SELECT fetch_time, ticker, candidate, series_ticker,
                       last_price, yes_bid, yes_ask, volume
                FROM kalshi_odds
                WHERE ticker IN ({placeholders})
                  AND fetch_time >= CURRENT_TIMESTAMP - INTERVAL '{days} days'
                ORDER BY fetch_time
            """, tickers).fetchdf()
        else:
            df = con.execute(f"""
                SELECT fetch_time, ticker, candidate, series_ticker,
                       last_price, yes_bid, yes_ask, volume
                FROM kalshi_odds
                WHERE fetch_time >= CURRENT_TIMESTAMP - INTERVAL '{days} days'
                ORDER BY fetch_time
            """).fetchdf()
        return df
    except Exception:
        return None
    finally:
        con.close()


def get_all_series():
    """List all discovered draft series."""
    con = get_connection(read_only=True)
    try:
        return con.execute("SELECT * FROM draft_series ORDER BY series_ticker").fetchdf()
    except Exception:
        return None
    finally:
        con.close()


def get_portfolio():
    """Get latest positions and resting orders."""
    con = get_connection(read_only=True)
    try:
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
    finally:
        con.close()


def get_position_changes():
    """Compare current vs previous position snapshots."""
    con = get_connection(read_only=True)
    try:
        fetch_times = con.execute("""
            SELECT DISTINCT fetch_time FROM positions
            ORDER BY fetch_time DESC LIMIT 2
        """).fetchall()

        if len(fetch_times) < 2:
            return None

        current_time = fetch_times[0][0]
        prev_time = fetch_times[1][0]

        changes = con.execute("""
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

        return changes
    except Exception:
        return None
    finally:
        con.close()


def get_latest_edges():
    """Get the most recently computed edges."""
    con = get_connection(read_only=True)
    try:
        return con.execute("""
            SELECT * FROM detected_edges
            WHERE fetch_time = (SELECT MAX(fetch_time) FROM detected_edges)
            ORDER BY ABS(implied_edge) DESC
        """).fetchdf()
    except Exception:
        return None
    finally:
        con.close()


def get_latest_consensus():
    """Get the most recent consensus board."""
    con = get_connection(read_only=True)
    try:
        return con.execute("""
            SELECT * FROM consensus_board
            WHERE fetch_time = (SELECT MAX(fetch_time) FROM consensus_board)
            ORDER BY rank
        """).fetchdf()
    except Exception:
        return None
    finally:
        con.close()


def get_snapshot_count():
    """Get number of historical snapshots."""
    con = get_connection(read_only=True)
    try:
        result = con.execute("""
            SELECT COUNT(DISTINCT fetch_time) as snapshots,
                   MIN(fetch_time) as first_fetch,
                   MAX(fetch_time) as last_fetch
            FROM kalshi_odds
        """).fetchone()
        return result
    except Exception:
        return (0, None, None)
    finally:
        con.close()
