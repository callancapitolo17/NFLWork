"""DuckDB log of the automated taker's decisions (automated_taker_plan.md, T1).

Owns unabated_ticket/bets_service/taker.duckdb, apart from bets.duckdb so the
taker's log never touches Cal's bet history (this service is its only writer).
Table:
  taker_decisions  APPEND, one row per decision the server runner
               (server/runner.js + server/taker.js) POSTs to
               /taker_decisions.json. The runner sends a row only when that
               line's decision changed (reason, ticker, live ask, Unabated
               price or contracts), so the table is the history of what the
               bot would have done and why, not one row per 30 s pass. T1
               places nothing: reason "bet" is a would-bet. Never pruned
               (a few thousand rows a day at most); `decision` holds the
               whole JSON the runner sent, the other columns are for queries.

The DuckDB connection is shared between the HTTP handler threads behind one
lock (DuckDB connections are not thread-safe).
"""
import json
import threading
from datetime import datetime
from pathlib import Path

import duckdb

# taker.js REASONS: a reason outside it is a runner the service does not know.
DECISION_REASONS = ("stale_feed", "no_ticker", "no_book", "price_moved", "edge_short", "not_best",
                    "cap", "size_declined", "size_zero", "bet")

_SCHEMA = """
CREATE TABLE IF NOT EXISTS taker_decisions (
    decided_at      TIMESTAMPTZ NOT NULL,
    received_at     TIMESTAMPTZ NOT NULL,
    row_key         VARCHAR NOT NULL,     -- the Edges line key (feed.js lineKeyOf / altLineKeyOf)
    league          VARCHAR,
    event_id        BIGINT,
    event_start     TIMESTAMPTZ,
    reason          VARCHAR NOT NULL,     -- DECISION_REASONS
    ticker          VARCHAR,              -- the Kalshi market, NULL on no_ticker / stale_feed
    contract_side   VARCHAR,              -- yes | no
    unabated_price  INTEGER,              -- Unabated's shown American price
    live_all_in     DOUBLE,               -- Kalshi ask + taker fee, as a probability
    live_edge       DOUBLE,               -- fair / live_all_in - 1
    contracts       INTEGER,              -- what the bot would buy (reason bet / size_zero)
    decision        VARCHAR NOT NULL      -- JSON, the whole decision as the runner sent it
);
"""

_INSERT = """
INSERT INTO taker_decisions (decided_at, received_at, row_key, league, event_id, event_start, reason, ticker,
                             contract_side, unabated_price, live_all_in, live_edge, contracts, decision)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""


class TakerStore:
    def __init__(self, db_path: str | Path):
        self._lock = threading.Lock()
        self._con = duckdb.connect(str(db_path))
        self._con.execute(_SCHEMA)

    def close(self) -> None:
        with self._lock:
            self._con.close()

    def append_decisions(self, rows: list[dict], received_at: datetime) -> int:
        """APPEND the validated rows ({decision, decided_at, event_start}, from
        service.validate_taker_decisions); returns how many."""
        values = []
        for row in rows:
            decision = row["decision"]
            values.append((
                row["decided_at"], received_at, decision["rowKey"], decision.get("league"), decision.get("eventId"),
                row["event_start"], decision["reason"], decision.get("ticker"), decision.get("side"),
                decision.get("unabatedPrice"), decision.get("liveAllIn"), decision.get("liveEdge"),
                decision.get("contracts"), json.dumps(decision, separators=(",", ":")),
            ))
        with self._lock:
            if values:
                self._con.executemany(_INSERT, values)
        return len(values)

    def count_decisions(self) -> int:
        with self._lock:
            return self._con.execute("SELECT count(*) FROM taker_decisions").fetchone()[0]
