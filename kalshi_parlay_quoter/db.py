"""Decision log (DuckDB). One writer: the main thread flushes buffered rows
in one transaction per FLUSH_SEC, so worker threads never touch the
connection. Insert-only (INSERT OR IGNORE: the first row for a key wins).

Tables
  rfq_decisions  one row per STORED RFQ: every confirmed parlay, plus a
                 sample of the rest (main.py)
  rfq_legs       one row per leg of each stored RFQ that was priced
  minute_counts  (minute, key, n): every frame and every decision counted,
                 stored or not, so the totals never depend on the sample
  door_counts    written before 2026-11-13 only; kept in old files, unused
"""
import json
import threading

import duckdb

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS rfq_decisions (
    rfq_id                  VARCHAR PRIMARY KEY,
    seen_at                 TIMESTAMPTZ NOT NULL,
    created_ts              TIMESTAMPTZ,
    combo_ticker            VARCHAR,
    creator_id              VARCHAR,
    contracts               DOUBLE,
    target_cost_usd         DOUBLE,
    n_legs                  INTEGER,
    leagues                 VARCHAR,
    n_games                 INTEGER,
    same_game               BOOLEAN,
    decision                VARCHAR NOT NULL,
    reason                  VARCHAR,
    reason_detail           VARCHAR,
    parlay_fair             DOUBLE,
    parlay_fair_lo          DOUBLE,
    parlay_fair_hi          DOUBLE,
    max_leg_sigma           DOUBLE,
    legs_kalshi_ask_product DOUBLE,
    earliest_kickoff_ms     BIGINT,
    hours_to_kickoff        DOUBLE,
    pricer_ms               DOUBLE,
    priced_at_ms            BIGINT,
    queue_ms                DOUBLE
);
CREATE INDEX IF NOT EXISTS idx_decisions_seen ON rfq_decisions(seen_at);
CREATE TABLE IF NOT EXISTS rfq_legs (
    rfq_id          VARCHAR NOT NULL,
    leg_idx         INTEGER NOT NULL,
    market_ticker   VARCHAR,
    side            VARCHAR,
    series          VARCHAR,
    league          VARCHAR,
    market          VARCHAR,
    game_key        VARCHAR,
    ok              BOOLEAN,
    reason          VARCHAR,
    fair            DOUBLE,
    fair_lo         DOUBLE,
    fair_hi         DOUBLE,
    bracketed       BOOLEAN,
    sigma           DOUBLE,
    c               DOUBLE,
    tau             DOUBLE,
    line            JSON,
    kickoff_ms      BIGINT,
    fair_age_ms     BIGINT,
    kalshi_status   VARCHAR,
    kalshi_yes_bid  DOUBLE,
    kalshi_yes_ask  DOUBLE,
    kalshi_no_ask   DOUBLE,
    PRIMARY KEY (rfq_id, leg_idx)
);
-- Added after the first VM deploy (2026-11-07): the runner's per-leg detail.
ALTER TABLE rfq_legs ADD COLUMN IF NOT EXISTS detail VARCHAR;
-- 2026-11-12: when the runner priced the leg and how long the quoter reused it.
ALTER TABLE rfq_legs ADD COLUMN IF NOT EXISTS priced_at_ms BIGINT;
ALTER TABLE rfq_legs ADD COLUMN IF NOT EXISTS cache_age_ms DOUBLE;
-- 2026-11-13: priced through the Kalshi-checked /rfq_price.json (else off
-- the leg table alone, kalshi_* columns null).
ALTER TABLE rfq_decisions ADD COLUMN IF NOT EXISTS confirmed BOOLEAN;
CREATE TABLE IF NOT EXISTS minute_counts (
    minute  TIMESTAMPTZ NOT NULL,
    key     VARCHAR NOT NULL,
    n       BIGINT NOT NULL,
    PRIMARY KEY (minute, key)
);
"""

DECISION_COLS = [
    "rfq_id", "seen_at", "created_ts", "combo_ticker", "creator_id",
    "contracts", "target_cost_usd", "n_legs", "leagues", "n_games",
    "same_game", "decision", "reason", "reason_detail", "parlay_fair",
    "parlay_fair_lo", "parlay_fair_hi", "max_leg_sigma",
    "legs_kalshi_ask_product", "earliest_kickoff_ms", "hours_to_kickoff",
    "pricer_ms", "priced_at_ms", "queue_ms", "confirmed"]
LEG_COLS = [
    "rfq_id", "leg_idx", "market_ticker", "side", "series", "league", "market",
    "game_key", "ok", "reason", "fair", "fair_lo", "fair_hi", "bracketed",
    "sigma", "c", "tau", "line", "kickoff_ms", "fair_age_ms", "kalshi_status",
    "kalshi_yes_bid", "kalshi_yes_ask", "kalshi_no_ask", "detail",
    "priced_at_ms", "cache_age_ms"]
COUNT_COLS = ["minute", "key", "n"]


CHUNK_ROWS = 500


def _insert_many(con, table, cols, rows):
    """One multi-row INSERT per CHUNK_ROWS rows, not executemany, which
    DuckDB runs row by row (~4-12 ms a row with the PK check; ~12x slower)."""
    one = "(" + ",".join("?" * len(cols)) + ")"
    for i in range(0, len(rows), CHUNK_ROWS):
        chunk = rows[i:i + CHUNK_ROWS]
        con.execute(
            f"INSERT OR IGNORE INTO {table} ({','.join(cols)}) VALUES "
            + ",".join([one] * len(chunk)),
            [v for r in chunk for v in r])


def _int(v):
    try:
        return None if v is None else int(v)
    except (TypeError, ValueError):
        return None


class DecisionLog:
    def __init__(self, path):
        self._con = duckdb.connect(str(path))
        self._con.execute(SCHEMA_SQL)
        self._lock = threading.Lock()
        self._decisions: list[tuple] = []
        self._legs: list[tuple] = []
        self._counts: list[tuple] = []

    def add(self, row: dict, leg_rows: list[dict]) -> None:
        d = tuple(row.get(c) for c in DECISION_COLS)
        legs = []
        for lr in leg_rows:
            lr = dict(lr)
            lr["line"] = None if lr.get("line") is None else json.dumps(lr["line"])
            lr["kickoff_ms"] = _int(lr.get("kickoff_ms"))
            lr["fair_age_ms"] = _int(lr.get("fair_age_ms"))
            lr["priced_at_ms"] = _int(lr.get("priced_at_ms"))
            b = lr.get("bracketed")
            lr["bracketed"] = None if b is None else bool(b)
            legs.append(tuple(lr.get(c) for c in LEG_COLS))
        with self._lock:
            self._decisions.append(d)
            self._legs.extend(legs)

    def add_counts(self, minute, counts: dict[str, int]) -> None:
        """One minute's counts; zero counts are not stored."""
        with self._lock:
            self._counts.extend((minute, key, n) for key, n in counts.items() if n)

    def flush(self) -> int:
        with self._lock:
            d, l, o = self._decisions, self._legs, self._counts
            self._decisions, self._legs, self._counts = [], [], []
        if not (d or l or o):
            return 0
        self._con.execute("BEGIN")
        try:
            if d:
                _insert_many(self._con, "rfq_decisions", DECISION_COLS, d)
            if l:
                _insert_many(self._con, "rfq_legs", LEG_COLS, l)
            if o:
                _insert_many(self._con, "minute_counts", COUNT_COLS, o)
            self._con.execute("COMMIT")
        except Exception:
            self._con.execute("ROLLBACK")
            raise
        return len(d)

    def checkpoint(self) -> None:
        """Fold the write-ahead log into the file (once a minute, so the WAL
        on the VM's data disk stays small)."""
        self._con.execute("CHECKPOINT")

    def query(self, sql, params=None):
        return self._con.execute(sql, params or []).fetchall()

    def close(self) -> None:
        self.flush()
        self._con.close()
