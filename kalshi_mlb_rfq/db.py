"""DuckDB schema + helpers for the Kalshi MLB RFQ bot."""

import shutil
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
    quote_id                         VARCHAR PRIMARY KEY,
    rfq_id                           VARCHAR NOT NULL,
    combo_market_ticker              VARCHAR,
    creator_id                       VARCHAR,
    yes_bid_dollars                  DOUBLE,
    no_bid_dollars                   DOUBLE,
    blended_fair_at_eval             DOUBLE,
    post_fee_ev_pct                  DOUBLE,
    decision                         VARCHAR NOT NULL,
    reason_detail                    VARCHAR,
    observed_at                      TIMESTAMP NOT NULL,
    competitor_count                 INTEGER,
    best_competitor_no_bid_dollars   DOUBLE,
    accept_response_body             VARCHAR,
    rfq_terminal_status              VARCHAR,
    quote_first_seen_at              TIMESTAMP,
    accept_attempted_at              TIMESTAMP,
    accept_response_at               TIMESTAMP
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
    combo_market_ticker  VARCHAR NOT NULL,
    side                 VARCHAR NOT NULL,
    game_id              VARCHAR NOT NULL,
    net_contracts        DOUBLE NOT NULL,
    weighted_price       DOUBLE NOT NULL,
    legs_json            VARCHAR NOT NULL,
    updated_at           TIMESTAMP NOT NULL,
    PRIMARY KEY (combo_market_ticker, side)
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
    leg_set_hash  VARCHAR NOT NULL,
    side          VARCHAR NOT NULL,
    game_id       VARCHAR NOT NULL,
    cooled_until  TIMESTAMP NOT NULL,
    reason        VARCHAR,
    PRIMARY KEY (leg_set_hash, side)
);
"""

MIGRATE_SQL = """
-- v1 quote_log diag columns (walk diagnostics)
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS competitor_count                 INTEGER;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS best_competitor_no_bid_dollars   DOUBLE;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS accept_response_body             VARCHAR;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS rfq_terminal_status              VARCHAR;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS quote_first_seen_at              TIMESTAMP;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS accept_attempted_at              TIMESTAMP;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS accept_response_at               TIMESTAMP;

-- v2 symmetric-side columns (NO-side accepts: 2026-05-22)
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS chosen_side          VARCHAR;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS ev_yes_pct           DOUBLE;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS ev_no_pct            DOUBLE;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS hedge_added          BOOLEAN;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS hedge_original_side  VARCHAR;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS hedge_original_price DOUBLE;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS hedge_new_price      DOUBLE;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS hedge_current_fair   DOUBLE;
ALTER TABLE quote_log ADD COLUMN IF NOT EXISTS hedge_projected_net  DOUBLE;
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


def _needs_v2_side_migration() -> bool:
    """True iff `positions` exists without a `side` column (= pre-v2 DB)."""
    if not DB_PATH.exists():
        return False
    with connect(read_only=True) as con:
        try:
            cols = {row[0] for row in con.execute(
                "SELECT column_name FROM information_schema.columns "
                "WHERE table_name='positions'"
            ).fetchall()}
        except duckdb.Error:
            return False
    if not cols:
        return False  # No positions table yet (fresh schema path)
    return "side" not in cols


def _backup_db_once(suffix: str) -> None:
    """Copy the .duckdb file to a timestamped .bak before a destructive migration.
    Bot must be stopped — DuckDB WAL files belong to the live process. We rely on
    the operational convention that init_database() only runs at startup."""
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    backup_path = DB_PATH.parent / f"{DB_PATH.name}.bak.{suffix}.{ts}"
    shutil.copy2(DB_PATH, backup_path)
    print(f"[db] Backed up {DB_PATH.name} to {backup_path.name} before {suffix}",
          flush=True)


def _migrate_v2_side_columns(con) -> None:
    """Add `side` to positions/combo_cooldown PKs via temp-table swap.

    Old PK (combo_market_ticker) → new PK (combo_market_ticker, side).
    Old PK (leg_set_hash)        → new PK (leg_set_hash, side).
    Existing rows backfill to side='yes' (correct historical interpretation —
    bot has only ever taken YES). Wrapped in a transaction so a mid-migration
    failure rolls back instead of leaving us with half-swapped tables."""
    con.execute("BEGIN TRANSACTION")
    try:
        # Defensive cleanup of leftover *_new tables from a prior failed run.
        con.execute("DROP TABLE IF EXISTS positions_new")
        con.execute("DROP TABLE IF EXISTS combo_cooldown_new")

        # --- positions ---
        con.execute("""
            CREATE TABLE positions_new (
                combo_market_ticker  VARCHAR NOT NULL,
                side                 VARCHAR NOT NULL,
                game_id              VARCHAR NOT NULL,
                net_contracts        DOUBLE NOT NULL,
                weighted_price       DOUBLE NOT NULL,
                legs_json            VARCHAR NOT NULL,
                updated_at           TIMESTAMP NOT NULL,
                PRIMARY KEY (combo_market_ticker, side)
            )
        """)
        con.execute("""
            INSERT INTO positions_new
                (combo_market_ticker, side, game_id, net_contracts,
                 weighted_price, legs_json, updated_at)
            SELECT combo_market_ticker, 'yes', game_id, net_contracts,
                   weighted_price, legs_json, updated_at
            FROM positions
        """)
        con.execute("DROP TABLE positions")
        con.execute("ALTER TABLE positions_new RENAME TO positions")
        con.execute("CREATE INDEX IF NOT EXISTS idx_positions_game "
                    "ON positions(game_id)")

        # --- combo_cooldown ---
        con.execute("""
            CREATE TABLE combo_cooldown_new (
                leg_set_hash  VARCHAR NOT NULL,
                side          VARCHAR NOT NULL,
                game_id       VARCHAR NOT NULL,
                cooled_until  TIMESTAMP NOT NULL,
                reason        VARCHAR,
                PRIMARY KEY (leg_set_hash, side)
            )
        """)
        con.execute("""
            INSERT INTO combo_cooldown_new
                (leg_set_hash, side, game_id, cooled_until, reason)
            SELECT leg_set_hash, 'yes', game_id, cooled_until, reason
            FROM combo_cooldown
        """)
        con.execute("DROP TABLE combo_cooldown")
        con.execute("ALTER TABLE combo_cooldown_new RENAME TO combo_cooldown")

        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise


def init_database():
    """Apply schema migrations idempotently.

    Order:
      1. If a pre-v2 DB exists (`positions` without `side`), back up first.
      2. Run SCHEMA_SQL — no-op for existing tables (CREATE IF NOT EXISTS).
      3. Run MIGRATE_SQL — additive ALTERs, idempotent via IF NOT EXISTS.
      4. Run v2 side-column migration if needed (destructive PK swap)."""
    needs_v2 = _needs_v2_side_migration()
    if needs_v2:
        _backup_db_once("v2_side_columns")

    with connect() as con:
        con.execute(SCHEMA_SQL)
        con.execute(MIGRATE_SQL)
        if needs_v2:
            _migrate_v2_side_columns(con)


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
