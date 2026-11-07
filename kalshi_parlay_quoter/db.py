"""Decision log (DuckDB). One writer: the main thread flushes buffered rows
in one transaction per FLUSH_SEC, so pricer worker threads never touch the
connection.

Tables
  rfq_decisions  one row per in-league RFQ id (first sighting wins)
  rfq_legs       one row per leg of each logged RFQ
  door_counts    per-minute counts of RFQs that got no row (other leagues,
                 below the size floor, malformed) plus feed mode, so every
                 RFQ seen is accounted for
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
    priced_at_ms            BIGINT
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
CREATE TABLE IF NOT EXISTS door_counts (
    minute          TIMESTAMPTZ PRIMARY KEY,
    logged          BIGINT,
    other_league    BIGINT,
    below_floor     BIGINT,
    too_many_legs   BIGINT,
    malformed       BIGINT,
    feed_mode       VARCHAR
);
"""

DECISION_COLS = [
    "rfq_id", "seen_at", "created_ts", "combo_ticker", "creator_id",
    "contracts", "target_cost_usd", "n_legs", "leagues", "n_games",
    "same_game", "decision", "reason", "reason_detail", "parlay_fair",
    "parlay_fair_lo", "parlay_fair_hi", "max_leg_sigma",
    "legs_kalshi_ask_product", "earliest_kickoff_ms", "hours_to_kickoff",
    "pricer_ms", "priced_at_ms"]
LEG_COLS = [
    "rfq_id", "leg_idx", "market_ticker", "side", "series", "league", "market",
    "game_key", "ok", "reason", "fair", "fair_lo", "fair_hi", "bracketed",
    "sigma", "c", "tau", "line", "kickoff_ms", "fair_age_ms", "kalshi_status",
    "kalshi_yes_bid", "kalshi_yes_ask", "kalshi_no_ask"]
DOOR_COLS = ["minute", "logged", "other_league", "below_floor",
             "too_many_legs", "malformed", "feed_mode"]


def _ins(table, cols):
    return (f"INSERT OR IGNORE INTO {table} ({','.join(cols)}) "
            f"VALUES ({','.join('?' * len(cols))})")


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
        self._doors: list[tuple] = []

    def add(self, row: dict, leg_rows: list[dict]) -> None:
        d = tuple(row.get(c) for c in DECISION_COLS)
        legs = []
        for lr in leg_rows:
            lr = dict(lr)
            lr["line"] = None if lr.get("line") is None else json.dumps(lr["line"])
            lr["kickoff_ms"] = _int(lr.get("kickoff_ms"))
            lr["fair_age_ms"] = _int(lr.get("fair_age_ms"))
            b = lr.get("bracketed")
            lr["bracketed"] = None if b is None else bool(b)
            legs.append(tuple(lr.get(c) for c in LEG_COLS))
        with self._lock:
            self._decisions.append(d)
            self._legs.extend(legs)

    def add_door(self, counts: dict) -> None:
        with self._lock:
            self._doors.append(tuple(counts.get(c) for c in DOOR_COLS))

    def flush(self) -> int:
        with self._lock:
            d, l, o = self._decisions, self._legs, self._doors
            self._decisions, self._legs, self._doors = [], [], []
        if not (d or l or o):
            return 0
        self._con.execute("BEGIN")
        try:
            if d:
                self._con.executemany(_ins("rfq_decisions", DECISION_COLS), d)
            if l:
                self._con.executemany(_ins("rfq_legs", LEG_COLS), l)
            if o:
                self._con.executemany(_ins("door_counts", DOOR_COLS), o)
            self._con.execute("COMMIT")
        except Exception:
            self._con.execute("ROLLBACK")
            raise
        return len(d)

    def query(self, sql, params=None):
        return self._con.execute(sql, params or []).fetchall()

    def close(self) -> None:
        self.flush()
        self._con.close()
