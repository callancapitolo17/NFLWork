"""
DuckDB state management for the market maker and taker.
Tracks positions, fills, resting orders, quote history, and takes.

Uses a shared write connection (opened lazily, closed on shutdown) to avoid
the overhead of opening/closing 1000+ connections per quote cycle. The bot
is single-threaded, so no lock contention on the shared connection.
"""

import math
import time as _time
import duckdb
from datetime import datetime, timezone
from pathlib import Path

from config import MM_DB_PATH, CBB_DB_PATH


# ---------------------------------------------------------------------------
# Shared write connection
# ---------------------------------------------------------------------------

_write_conn = None


def get_write_conn():
    """Get or create the shared write connection to kalshi_mm.duckdb."""
    global _write_conn
    if _write_conn is None:
        _write_conn = duckdb.connect(str(MM_DB_PATH))
    return _write_conn


def _reset_write_conn():
    """Close bad connection so next get_write_conn() reconnects."""
    global _write_conn
    try:
        if _write_conn is not None:
            _write_conn.close()
    except Exception:
        pass
    _write_conn = None


def close_write_conn():
    """Close shared connection (call on shutdown to flush WAL)."""
    global _write_conn
    if _write_conn is not None:
        try:
            _write_conn.close()
        except Exception:
            pass
        _write_conn = None


# ---------------------------------------------------------------------------
# Retry helper (startup only — shared conn eliminates need in hot path)
# ---------------------------------------------------------------------------

def _with_retry(fn, max_retries=3, base_delay=0.1):
    """Execute a DB write function with retry on DuckDB lock conflict.

    Only used during startup (init_database, init_taker_tables) where a
    dying previous process may still hold the lock.
    """
    for attempt in range(max_retries):
        try:
            return fn()
        except duckdb.IOException as e:
            if attempt < max_retries - 1 and "lock" in str(e).lower():
                _time.sleep(base_delay * (2 ** attempt))
            else:
                raise


# ---------------------------------------------------------------------------
# Schema / init (use their own connections — run before shared conn)
# ---------------------------------------------------------------------------

def _tables_exist(required):
    """Check if all required tables exist using a read-only connection."""
    if not Path(str(MM_DB_PATH)).exists():
        return False
    conn = duckdb.connect(str(MM_DB_PATH), read_only=True)
    try:
        existing = {r[0] for r in conn.execute(
            "SELECT table_name FROM information_schema.tables"
        ).fetchall()}
        return all(t in existing for t in required)
    finally:
        conn.close()


def _migrate_positions():
    """Add missing columns to positions table (Kelly sizing migrations)."""
    if not Path(str(MM_DB_PATH)).exists():
        return
    conn = duckdb.connect(str(MM_DB_PATH), read_only=True)
    try:
        cols = {r[0] for r in conn.execute(
            "SELECT column_name FROM information_schema.columns "
            "WHERE table_name = 'positions'"
        ).fetchall()}
    finally:
        conn.close()

    migrations = []
    if "fair_prob" not in cols:
        migrations.append("ALTER TABLE positions ADD COLUMN fair_prob FLOAT DEFAULT 0.5")
    if "contract_team" not in cols:
        migrations.append("ALTER TABLE positions ADD COLUMN contract_team VARCHAR")

    if migrations:
        def _add_cols():
            c = duckdb.connect(str(MM_DB_PATH))
            try:
                for sql in migrations:
                    c.execute(sql)
            finally:
                c.close()
        _with_retry(_add_cols)


def init_database():
    """Create all market maker tables if they don't exist."""
    required = ["positions", "fills", "resting_orders", "quote_log", "sessions", "reference_lines"]
    if _tables_exist(required):
        _migrate_positions()
        return
    def _init():
        conn = duckdb.connect(str(MM_DB_PATH))
        try:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS positions (
                    ticker VARCHAR PRIMARY KEY,
                    event_ticker VARCHAR,
                    home_team VARCHAR,
                    away_team VARCHAR,
                    market_type VARCHAR,
                    line_value FLOAT,
                    net_yes INTEGER DEFAULT 0,
                    avg_entry_price FLOAT DEFAULT 0,
                    fair_prob FLOAT DEFAULT 0.5,
                    contract_team VARCHAR,
                    updated_at TIMESTAMP
                )
            """)
            conn.execute("""
                CREATE TABLE IF NOT EXISTS fills (
                    fill_id VARCHAR PRIMARY KEY,
                    ticker VARCHAR,
                    side VARCHAR,
                    action VARCHAR,
                    price INTEGER,
                    count INTEGER,
                    fee_cents FLOAT,
                    filled_at TIMESTAMP,
                    order_id VARCHAR
                )
            """)
            conn.execute("""
                CREATE TABLE IF NOT EXISTS resting_orders (
                    order_id VARCHAR PRIMARY KEY,
                    ticker VARCHAR,
                    side VARCHAR,
                    action VARCHAR,
                    price INTEGER,
                    remaining_count INTEGER,
                    created_at TIMESTAMP,
                    last_amended_at TIMESTAMP
                )
            """)
            conn.execute("""
                CREATE TABLE IF NOT EXISTS quote_log (
                    logged_at TIMESTAMP,
                    ticker VARCHAR,
                    fair_prob FLOAT,
                    bid_yes INTEGER,
                    ask_yes INTEGER,
                    spread_cents INTEGER,
                    net_position INTEGER,
                    skew_applied INTEGER,
                    prediction_age_sec FLOAT,
                    book_bid INTEGER,
                    book_ask INTEGER
                )
            """)
            conn.execute("""
                CREATE TABLE IF NOT EXISTS sessions (
                    session_id VARCHAR PRIMARY KEY,
                    started_at TIMESTAMP,
                    ended_at TIMESTAMP,
                    total_fills INTEGER DEFAULT 0,
                    realized_pnl FLOAT DEFAULT 0,
                    fees_paid FLOAT DEFAULT 0,
                    markets_quoted INTEGER DEFAULT 0
                )
            """)
            conn.execute("""
                CREATE TABLE IF NOT EXISTS reference_lines (
                    home_team VARCHAR,
                    away_team VARCHAR,
                    market VARCHAR,
                    line_value FLOAT,
                    snapshot_at TIMESTAMP
                )
            """)
            # Indexes for pruning performance
            conn.execute("CREATE INDEX IF NOT EXISTS idx_quote_log_time ON quote_log(logged_at)")
            # Migrate: add book columns if missing (existing DBs won't have them)
            try:
                conn.execute("ALTER TABLE quote_log ADD COLUMN book_bid INTEGER")
            except Exception:
                pass
            try:
                conn.execute("ALTER TABLE quote_log ADD COLUMN book_ask INTEGER")
            except Exception:
                pass
        finally:
            conn.close()
    _with_retry(_init, max_retries=10, base_delay=0.5)


# ---------------------------------------------------------------------------
# Read helpers (use read-only connections — no contention)
# ---------------------------------------------------------------------------

def _safe_float(val):
    """Return val if it's a valid float, else None. Guards against NaN from pandas."""
    if val is None:
        return None
    try:
        if math.isnan(val):
            return None
    except (TypeError, ValueError):
        return None
    return val


def get_prediction_timestamp():
    """Check when predictions were last updated (lightweight — 1 row read)."""
    db_path = str(CBB_DB_PATH)
    if not Path(db_path).exists():
        return None
    try:
        conn = duckdb.connect(db_path, read_only=True)
        try:
            meta = conn.execute(
                "SELECT updated_at FROM cbb_prediction_meta LIMIT 1"
            ).fetchone()
            return meta[0] if meta else None
        finally:
            conn.close()
    except Exception:
        return None


def load_predictions():
    """Load raw predictions from the answer key's DuckDB."""
    db_path = str(CBB_DB_PATH)
    if not Path(db_path).exists():
        print(f"Warning: CBB database not found at {db_path}")
        return [], None

    conn = duckdb.connect(db_path, read_only=True)
    try:
        # Check if table exists
        tables = [r[0] for r in conn.execute(
            "SELECT table_name FROM information_schema.tables WHERE table_name = 'cbb_raw_predictions'"
        ).fetchall()]
        if not tables:
            print("Warning: cbb_raw_predictions table not found. Run the pipeline first.")
            return [], None

        predictions = conn.execute(
            "SELECT * FROM cbb_raw_predictions"
        ).fetchdf()

        # Get prediction timestamp
        meta = conn.execute(
            "SELECT updated_at FROM cbb_prediction_meta LIMIT 1"
        ).fetchone()
        updated_at = meta[0] if meta else None

        # Sanitize NaN values from pandas — DuckDB NULL → pandas NaN
        records = predictions.to_dict("records")
        for rec in records:
            for key in ("prob_side1", "prob_side2", "prob_tie", "line_value"):
                if key in rec:
                    rec[key] = _safe_float(rec[key])

        return records, updated_at
    finally:
        conn.close()


def get_position(ticker):
    """Get current position for a ticker."""
    conn = get_write_conn()
    row = conn.execute(
        "SELECT net_yes, avg_entry_price FROM positions WHERE ticker = ?",
        [ticker]
    ).fetchone()
    return {"net_yes": row[0], "avg_entry_price": row[1]} if row else {"net_yes": 0, "avg_entry_price": 0}


def get_all_positions():
    """Get all open positions."""
    conn = get_write_conn()
    return conn.execute(
        "SELECT * FROM positions WHERE net_yes != 0"
    ).fetchdf().to_dict("records")


# ---------------------------------------------------------------------------
# Write functions (use shared connection)
# ---------------------------------------------------------------------------

def update_position(ticker, side, price, count, event_ticker=None,
                    home_team=None, away_team=None, market_type=None,
                    line_value=None, fair_prob=None, contract_team=None):
    """Update position after a fill."""
    conn = get_write_conn()
    now = datetime.now(timezone.utc)
    current = conn.execute(
        "SELECT net_yes, avg_entry_price FROM positions WHERE ticker = ?",
        [ticker]
    ).fetchone()

    if current:
        net_yes, avg_price = current
    else:
        net_yes, avg_price = 0, 0.0

    if side == "yes":
        delta = count
    else:
        delta = -count

    new_net = net_yes + delta

    if delta > 0 and net_yes >= 0:
        total_cost = avg_price * net_yes + price * delta
        new_avg = total_cost / new_net if new_net > 0 else 0
    elif delta < 0 and net_yes <= 0:
        total_cost = avg_price * abs(net_yes) + price * abs(delta)
        new_avg = total_cost / abs(new_net) if new_net != 0 else 0
    else:
        if new_net == 0:
            new_avg = 0
        elif (net_yes > 0 and new_net < 0) or (net_yes < 0 and new_net > 0):
            new_avg = price
        else:
            new_avg = avg_price

    conn.execute("""
        INSERT INTO positions (ticker, event_ticker, home_team, away_team,
                               market_type, line_value, net_yes, avg_entry_price,
                               fair_prob, contract_team, updated_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT (ticker) DO UPDATE SET
            net_yes = ?, avg_entry_price = ?,
            fair_prob = COALESCE(?, fair_prob, 0.5),
            line_value = COALESCE(?, line_value),
            contract_team = COALESCE(?, contract_team),
            updated_at = ?
    """, [ticker, event_ticker, home_team, away_team, market_type, line_value,
          new_net, new_avg, fair_prob or 0.5, contract_team, now,
          new_net, new_avg, fair_prob, line_value, contract_team, now])


def record_fill(fill_id, ticker, side, action, price, count, fee_cents, order_id):
    """Record a fill in the audit log."""
    conn = get_write_conn()
    conn.execute("""
        INSERT INTO fills (fill_id, ticker, side, action, price, count, fee_cents, filled_at, order_id)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT (fill_id) DO NOTHING
    """, [fill_id, ticker, side, action, price, count, fee_cents,
          datetime.now(timezone.utc), order_id])


# ---------------------------------------------------------------------------
# Quote logging (buffered + batched)
# ---------------------------------------------------------------------------

_quote_log_buffer = []
_quote_log_flush_count = 0


def log_quote(ticker, fair_prob, bid_yes, ask_yes, net_position, skew, prediction_age,
              book_bid=0, book_ask=0):
    """Append a quote to the batch buffer. Call flush_quote_log() to write."""
    spread = (ask_yes - bid_yes) if bid_yes and ask_yes else 0
    _quote_log_buffer.append((
        datetime.now(timezone.utc), ticker, fair_prob, bid_yes, ask_yes,
        spread, net_position, skew, prediction_age, book_bid, book_ask
    ))


def flush_quote_log():
    """Write buffered quote logs in one batch, prune old entries periodically."""
    global _quote_log_flush_count
    if not _quote_log_buffer:
        return
    try:
        conn = get_write_conn()
        conn.executemany("""
            INSERT INTO quote_log (logged_at, ticker, fair_prob, bid_yes, ask_yes,
                                   spread_cents, net_position, skew_applied, prediction_age_sec,
                                   book_bid, book_ask)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, _quote_log_buffer)
        _quote_log_buffer.clear()

        # Prune old entries every ~100 flushes (~16 min) instead of every cycle
        _quote_log_flush_count += 1
        if _quote_log_flush_count >= 100:
            conn.execute("""
                DELETE FROM quote_log
                WHERE logged_at < current_timestamp - INTERVAL '24 hours'
            """)
            _quote_log_flush_count = 0
    except duckdb.Error as e:
        print(f"  Warning: flush_quote_log DB error — resetting connection: {e}")
        _reset_write_conn()
        # Don't clear buffer — retry next cycle
    except Exception as e:
        print(f"  Warning: flush_quote_log failed: {e}")


# ---------------------------------------------------------------------------
# Resting order tracking (individual + batch)
# ---------------------------------------------------------------------------

def save_resting_order(order_id, ticker, side, action, price, count):
    """Track a resting order we placed."""
    conn = get_write_conn()
    now = datetime.now(timezone.utc)
    conn.execute("""
        INSERT INTO resting_orders (order_id, ticker, side, action, price,
                                    remaining_count, created_at, last_amended_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT (order_id) DO UPDATE SET
            price = ?, remaining_count = ?, last_amended_at = ?
    """, [order_id, ticker, side, action, price, count, now, now,
          price, count, now])


def batch_save_resting_orders(order_specs):
    """Save multiple resting orders in one batch write.

    Args:
        order_specs: list of dicts with keys: order_id, ticker, side, action, price, count
    """
    if not order_specs:
        return
    conn = get_write_conn()
    now = datetime.now(timezone.utc)
    conn.executemany("""
        INSERT INTO resting_orders (order_id, ticker, side, action, price,
                                    remaining_count, created_at, last_amended_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT (order_id) DO UPDATE SET
            price = excluded.price, remaining_count = excluded.remaining_count,
            last_amended_at = excluded.last_amended_at
    """, [(s["order_id"], s["ticker"], s["side"], s["action"], s["price"],
           s["count"], now, now) for s in order_specs])


def remove_resting_order(order_id):
    """Remove a resting order (cancelled or fully filled)."""
    conn = get_write_conn()
    conn.execute("DELETE FROM resting_orders WHERE order_id = ?", [order_id])


def batch_remove_resting_orders(order_ids):
    """Remove multiple resting orders in one query."""
    if not order_ids:
        return
    conn = get_write_conn()
    placeholders = ",".join(["?"] * len(order_ids))
    conn.execute(
        f"DELETE FROM resting_orders WHERE order_id IN ({placeholders})",
        list(order_ids)
    )


def clear_all_resting_orders():
    """Flush all resting order records from DB (startup cleanup)."""
    conn = get_write_conn()
    conn.execute("DELETE FROM resting_orders")


def get_resting_orders():
    """Get all our resting orders."""
    conn = get_write_conn()
    return conn.execute("SELECT * FROM resting_orders").fetchdf().to_dict("records")


# ---------------------------------------------------------------------------
# Reference lines
# ---------------------------------------------------------------------------

def save_reference_lines(lines):
    """Save reference lines from Bookmaker/Bet105 for line-move detection."""
    conn = get_write_conn()
    conn.execute("BEGIN TRANSACTION")
    try:
        conn.execute("DELETE FROM reference_lines")
        now = datetime.now(timezone.utc)
        for line in lines:
            conn.execute("""
                INSERT INTO reference_lines (home_team, away_team, market, line_value, snapshot_at)
                VALUES (?, ?, ?, ?, ?)
            """, [line["home_team"], line["away_team"], line["market"],
                  line["line_value"], now])
        conn.execute("COMMIT")
    except Exception:
        try:
            conn.execute("ROLLBACK")
        except Exception:
            _reset_write_conn()
        raise


def get_reference_lines():
    """Get stored reference lines."""
    conn = get_write_conn()
    return conn.execute("SELECT * FROM reference_lines").fetchdf().to_dict("records")


# ---------------------------------------------------------------------------
# Exposure
# ---------------------------------------------------------------------------

def compute_total_exposure():
    """Compute total dollars at risk across all positions AND resting orders."""
    conn = get_write_conn()
    filled = conn.execute("""
        SELECT COALESCE(SUM(
            avg_entry_price * ABS(net_yes) / 100.0
        ), 0) FROM positions WHERE net_yes != 0
    """).fetchone()[0]

    resting = conn.execute("""
        SELECT COALESCE(SUM(
            price * remaining_count / 100.0
        ), 0) FROM resting_orders
    """).fetchone()[0]

    return filled + resting


# ---------------------------------------------------------------------------
# Sessions
# ---------------------------------------------------------------------------

def start_session(session_id):
    """Record session start."""
    conn = get_write_conn()
    conn.execute("""
        INSERT INTO sessions (session_id, started_at)
        VALUES (?, ?)
    """, [session_id, datetime.now(timezone.utc)])


def end_session(session_id, total_fills, realized_pnl, fees_paid, markets_quoted):
    """Record session end."""
    conn = get_write_conn()
    conn.execute("""
        UPDATE sessions SET ended_at = ?, total_fills = ?,
            realized_pnl = ?, fees_paid = ?, markets_quoted = ?
        WHERE session_id = ?
    """, [datetime.now(timezone.utc), total_fills, realized_pnl,
          fees_paid, markets_quoted, session_id])


# ---------------------------------------------------------------------------
# Taker tables
# ---------------------------------------------------------------------------

def init_taker_tables():
    """Create taker-specific tables if they don't exist."""
    if _tables_exist(["take_log"]):
        return
    def _init():
        conn = duckdb.connect(str(MM_DB_PATH))
        try:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS take_log (
                    take_id VARCHAR PRIMARY KEY,
                    ticker VARCHAR,
                    event_ticker VARCHAR,
                    side VARCHAR,
                    price INTEGER,
                    fair_cents FLOAT,
                    ev_pct FLOAT,
                    fee_cents FLOAT,
                    count_requested INTEGER,
                    count_filled INTEGER,
                    order_id VARCHAR,
                    taken_at TIMESTAMP
                )
            """)
        finally:
            conn.close()
    _with_retry(_init, max_retries=10, base_delay=0.5)


def log_take(take_id, ticker, event_ticker, side, price, fair_cents, ev_pct,
             fee_cents, count_requested, count_filled, order_id):
    """Record a take attempt for analysis."""
    conn = get_write_conn()
    conn.execute("""
        INSERT INTO take_log (take_id, ticker, event_ticker, side, price,
            fair_cents, ev_pct, fee_cents, count_requested, count_filled,
            order_id, taken_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT (take_id) DO NOTHING
    """, [take_id, ticker, event_ticker, side, price, fair_cents, ev_pct,
          fee_cents, count_requested, count_filled, order_id,
          datetime.now(timezone.utc)])
