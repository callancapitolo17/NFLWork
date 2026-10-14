"""DuckDB log of the automated taker's decisions (automated_taker_plan.md, T1).

Owns unabated_ticket/bets_service/taker.duckdb, apart from bets.duckdb so the
taker's log never touches Cal's bet history (this service is its only writer).
Tables:
  taker_decisions  APPEND, one row per decision the server runner
               (server/runner.js + server/taker.js) POSTs to
               /taker_decisions.json. The runner sends a row only when that
               line's decision changed (reason, ticker, live ask, Unabated
               price or contracts), so the table is the history of what the
               bot would have done and why, not one row per 30 s pass. T1
               places nothing: reason "bet" is a would-bet. Never pruned
               (a few thousand rows a day at most); `decision` holds the
               whole JSON the runner sent, the other columns are for queries.
               "bet" means the bot sent (or tried to send) an order;
               `decision.order` carries what came back.
  taker_orders  one row per order the bot sent to Kalshi, keyed by its
               `client_order_id` ("ubot-…", the tag that marks a bot bet).
               INSERTED before the order is sent (so a crash mid-send still
               counts it against the ceiling), then UPDATED as Kalshi
               answers, its fills are read, and its market settles. The $500
               ceiling is the sum of `open_cost` over rows not yet settled
               (open_exposure). Never pruned (a handful a day).

The DuckDB connection is shared between the HTTP handler threads behind one
lock (DuckDB connections are not thread-safe).
"""
import json
import threading
from datetime import datetime, timezone
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
CREATE TABLE IF NOT EXISTS taker_orders (
    client_order_id  VARCHAR PRIMARY KEY,  -- "ubot-<hash>", deterministic per (ticker, side, pass)
    placed_at        TIMESTAMPTZ NOT NULL,
    row_key          VARCHAR NOT NULL,
    ticker           VARCHAR NOT NULL,
    contract_side    VARCHAR NOT NULL,     -- yes | no
    ask_price        DOUBLE NOT NULL,      -- the limit: the ask the decision was checked at
    all_in_price     DOUBLE NOT NULL,      -- ask + taker fee per contract
    count_requested  INTEGER NOT NULL,
    event_start      TIMESTAMPTZ,
    status           VARCHAR NOT NULL,     -- sending | sent | filled | partial | unfilled | rejected | unknown
    kalshi_order_id  VARCHAR,
    count_filled     DOUBLE,               -- NULL until the fills are read
    cost_filled      DOUBLE,               -- dollars paid incl. fees, NULL until read
    fills_final      BOOLEAN NOT NULL DEFAULT FALSE,  -- a fills read at least FILLS_SETTLE_SEC after placing
    open_cost        DOUBLE NOT NULL,      -- what counts against the ceiling while the market is open
    error            VARCHAR,
    settled_at       TIMESTAMPTZ           -- when its market was seen settled; frees open_cost
);
"""

_INSERT = """
INSERT INTO taker_decisions (decided_at, received_at, row_key, league, event_id, event_start, reason, ticker,
                             contract_side, unabated_price, live_all_in, live_edge, contracts, decision)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""


def _utc_from_ms(epoch_ms: int | None) -> datetime | None:
    return None if epoch_ms is None else datetime.fromtimestamp(epoch_ms / 1000, tz=timezone.utc)


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

    def insert_order(self, order: dict) -> None:
        """INSERT a new order row (status "sending", open_cost the full requested cost)."""
        with self._lock:
            self._con.execute(
                "INSERT INTO taker_orders (client_order_id, placed_at, row_key, ticker, contract_side, ask_price, "
                "all_in_price, count_requested, event_start, status, open_cost) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [order["client_order_id"], order["placed_at"], order["row_key"], order["ticker"], order["contract_side"],
                 order["ask_price"], order["all_in_price"], order["count_requested"], order["event_start"],
                 "sending", order["open_cost"]])

    def has_order(self, client_order_id: str) -> bool:
        with self._lock:
            return self._con.execute("SELECT count(*) FROM taker_orders WHERE client_order_id = ?",
                                     [client_order_id]).fetchone()[0] > 0

    def update_order(self, client_order_id: str, fields: dict) -> None:
        """UPDATE the named columns of one order (column names are this module's own, never caller text)."""
        allowed = {"status", "kalshi_order_id", "count_filled", "cost_filled", "fills_final", "open_cost", "error",
                   "settled_at"}
        unknown = set(fields) - allowed
        if unknown:
            raise ValueError(f"update_order: unknown columns {sorted(unknown)}")
        assignments = ", ".join(f"{column} = ?" for column in fields)
        with self._lock:
            self._con.execute(f"UPDATE taker_orders SET {assignments} WHERE client_order_id = ?",
                              [*fields.values(), client_order_id])

    def open_orders(self) -> list[dict]:
        """Every order whose market has not been seen settled, oldest first."""
        # Times leave DuckDB as epoch ms: a TIMESTAMPTZ fetched as a datetime needs pytz, which the service lacks.
        with self._lock:
            rows = self._con.execute(
                "SELECT client_order_id, epoch_ms(placed_at), ticker, contract_side, status, kalshi_order_id, "
                "fills_final, open_cost, epoch_ms(event_start), all_in_price, count_requested "
                "FROM taker_orders WHERE settled_at IS NULL ORDER BY placed_at").fetchall()
        return [{
            "client_order_id": row[0], "placed_at": _utc_from_ms(row[1]), "ticker": row[2], "contract_side": row[3],
            "status": row[4], "kalshi_order_id": row[5], "fills_final": row[6], "open_cost": row[7],
            "event_start": _utc_from_ms(row[8]), "all_in_price": row[9], "count_requested": row[10],
        } for row in rows]

    def open_exposure(self) -> float:
        """Dollars the bot has in markets not yet settled: the ceiling's measure."""
        with self._lock:
            return float(self._con.execute(
                "SELECT coalesce(sum(open_cost), 0) FROM taker_orders WHERE settled_at IS NULL").fetchone()[0])

    def count_decisions(self) -> int:
        with self._lock:
            return self._con.execute("SELECT count(*) FROM taker_decisions").fetchone()[0]
