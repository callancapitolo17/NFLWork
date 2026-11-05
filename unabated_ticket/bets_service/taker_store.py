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
               answers, the rest is cancelled, its fills are read, and its
               market settles (taker_orders.py has the status vocabulary).
               The $500 ceiling is the sum of `open_cost` over rows not yet
               settled (open_exposure). Never pruned (a handful a day).
               `subaccount` (added 2026-10-16 by ALTER on open) is the Kalshi
               subaccount the bot's key is locked to; NULL on the orders the
               bot placed from the primary account before that.
  novig_taker_orders  one row per order the bot sent to Novig (added
               2026-10-16, novig_taker_desk.py), keyed by its `client_id` (a
               UUID; Novig does not dedupe it, so the desk never re-sends).
               INSERTED before the send, then UPDATED as Novig answers, its
               fills are read, and its market settles. Its own $850 ceiling is
               the sum of `open_cost` over rows not yet settled. The row also
               keeps the board line the bot priced (league, teams, market,
               side, number), which the bot's subaccount read
               (sources/novig_subaccount.py) gives its bets: Novig's own
               catalog names neither the board's teams nor its periods.

The DuckDB connection is shared between the HTTP handler threads behind one
lock (DuckDB connections are not thread-safe).
"""
import json
import threading
from datetime import datetime, timezone
from pathlib import Path

import duckdb

# taker.js REASONS: a reason outside it is a runner the service does not know.
DECISION_REASONS = ("stale_feed", "no_ticker", "no_outcome", "no_book", "price_moved", "edge_short", "not_best",
                    "cap", "size_declined", "size_zero", "market_mismatch", "bet")

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
    event_id         BIGINT NOT NULL,      -- Unabated's event: the game an add must see earlier fills on
    ticker           VARCHAR NOT NULL,
    contract_side    VARCHAR NOT NULL,     -- yes | no
    ask_price        DOUBLE NOT NULL,      -- the limit: the ask the decision was checked at
    all_in_price     DOUBLE NOT NULL,      -- ask + taker fee per contract
    count_requested  INTEGER NOT NULL,
    event_start      TIMESTAMPTZ,
    status           VARCHAR NOT NULL,     -- sending | unknown | sent | resting | filled | partial | unfilled | rejected
    kalshi_order_id  VARCHAR,
    exchange_index   INTEGER,              -- Kalshi's shard for the order, when its reply named one; routes the cancel
    cancelled_at     TIMESTAMPTZ,          -- when the cancel of the rest was confirmed; fills are final 60 s on
    count_reported   DOUBLE,               -- contracts Kalshi's create reply said filled; the fills are final only once they show them all
    count_filled     DOUBLE,               -- NULL until the fills are read
    cost_filled      DOUBLE,               -- dollars paid incl. fees, NULL until read
    fills_final      BOOLEAN NOT NULL DEFAULT FALSE,  -- a fills read at least FILLS_SETTLE_SEC after the cancel
    open_cost        DOUBLE NOT NULL,      -- what counts against the ceiling while the market is open
    error            VARCHAR,
    settled_at       TIMESTAMPTZ,          -- when its market was seen settled; frees open_cost
    subaccount       INTEGER               -- the Kalshi subaccount it traded in; NULL = primary (before 2026-10-16)
);
ALTER TABLE taker_orders ADD COLUMN IF NOT EXISTS subaccount INTEGER;
CREATE TABLE IF NOT EXISTS novig_taker_orders (
    client_id        VARCHAR PRIMARY KEY,  -- the UUID the order carried (clientId)
    placed_at        TIMESTAMPTZ NOT NULL,
    row_key          VARCHAR NOT NULL,
    event_id         BIGINT NOT NULL,      -- Unabated's event: the game an add must see earlier fills on
    event_start      TIMESTAMPTZ NOT NULL,
    outcome_id       VARCHAR NOT NULL,     -- the Novig outcome bought
    market_id        VARCHAR NOT NULL,
    price            DOUBLE NOT NULL,      -- the limit: the ask the decision was checked at (no fee pregame)
    qty_requested    INTEGER NOT NULL,     -- 1c contracts (100 = $1 payout)
    status           VARCHAR NOT NULL,     -- sending | unknown | sent | resting | filled | partial | unfilled | rejected
    novig_order_id   VARCHAR,
    final_at         TIMESTAMPTZ,          -- when Novig showed the order final (FILLED / CANCELED / REJECTED)
    qty_filled       INTEGER,              -- NULL until the fills are read
    cost_filled      DOUBLE,               -- dollars paid incl. fees, NULL until read
    fills_final      BOOLEAN NOT NULL DEFAULT FALSE,  -- a fills read at least FILLS_SETTLE_SEC after final_at
    open_cost        DOUBLE NOT NULL,      -- what counts against the ceiling while the market is open
    error            VARCHAR,
    settled_at       TIMESTAMPTZ,          -- when its market was seen settled; frees open_cost
    league           VARCHAR NOT NULL,     -- the board line, as a bet record names it (feed.js league path)
    bet_type         VARCHAR NOT NULL,     -- spread | total | moneyline
    period           VARCHAR NOT NULL,     -- FG | 1H | ...
    side             VARCHAR NOT NULL,     -- away | home | over | under
    points           DOUBLE,
    away_team        VARCHAR,
    home_team        VARCHAR,
    away_team_id     VARCHAR,              -- Unabated's team ids
    home_team_id     VARCHAR
);
"""

NOVIG_ORDER_COLUMNS = ("client_id", "placed_at", "row_key", "event_id", "event_start", "outcome_id", "market_id", "price",
                       "qty_requested", "status", "open_cost", "league", "bet_type", "period", "side", "points",
                       "away_team", "home_team", "away_team_id", "home_team_id")
NOVIG_UPDATABLE_COLUMNS = ("status", "novig_order_id", "final_at", "qty_filled", "cost_filled", "fills_final", "open_cost",
                           "error", "settled_at")
# What open_novig_orders reads back; times as epoch ms (see open_orders).
_NOVIG_OPEN_COLUMNS = ("client_id", "placed_at", "event_id", "event_start", "outcome_id", "market_id", "price",
                       "qty_requested", "status", "novig_order_id", "final_at", "qty_filled", "fills_final", "open_cost")
_NOVIG_TIME_COLUMNS = ("placed_at", "event_start", "final_at")
NOVIG_CONTEXT_COLUMNS = ("league", "bet_type", "period", "side", "points", "away_team", "home_team", "away_team_id",
                         "home_team_id", "event_start", "event_id", "market_id", "outcome_id", "price", "settled_at")
_NOVIG_CONTEXT_TIME_COLUMNS = ("event_start", "settled_at")

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
        """INSERT a new order row (status "sending", open_cost the full requested cost, the
        subaccount it trades in)."""
        with self._lock:
            self._con.execute(
                "INSERT INTO taker_orders (client_order_id, placed_at, row_key, event_id, ticker, contract_side, "
                "ask_price, all_in_price, count_requested, event_start, status, open_cost, subaccount) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [order["client_order_id"], order["placed_at"], order["row_key"], order["event_id"], order["ticker"],
                 order["contract_side"], order["ask_price"], order["all_in_price"], order["count_requested"],
                 order["event_start"], "sending", order["open_cost"], order["subaccount"]])

    def has_order(self, client_order_id: str) -> bool:
        with self._lock:
            return self._con.execute("SELECT count(*) FROM taker_orders WHERE client_order_id = ?",
                                     [client_order_id]).fetchone()[0] > 0

    def update_order(self, client_order_id: str, fields: dict) -> None:
        """UPDATE the named columns of one order (column names are this module's own, never caller text)."""
        allowed = {"status", "kalshi_order_id", "exchange_index", "cancelled_at", "count_reported", "count_filled",
                   "cost_filled", "fills_final", "open_cost", "error", "settled_at"}
        unknown = set(fields) - allowed
        if unknown:
            raise ValueError(f"update_order: unknown columns {sorted(unknown)}")
        assignments = ", ".join(f"{column} = ?" for column in fields)
        with self._lock:
            self._con.execute(f"UPDATE taker_orders SET {assignments} WHERE client_order_id = ?",
                              [*fields.values(), client_order_id])

    def open_orders(self) -> list[dict]:
        """Every order whose market has not been seen settled, oldest first
        (`subaccount` 0 for the orders placed from the primary account)."""
        # Times leave DuckDB as epoch ms: a TIMESTAMPTZ fetched as a datetime needs pytz, which the service lacks.
        with self._lock:
            rows = self._con.execute(
                "SELECT client_order_id, epoch_ms(placed_at), ticker, contract_side, status, kalshi_order_id, "
                "fills_final, open_cost, epoch_ms(event_start), all_in_price, count_requested, event_id, "
                "exchange_index, epoch_ms(cancelled_at), count_filled, count_reported, coalesce(subaccount, 0) "
                "FROM taker_orders WHERE settled_at IS NULL ORDER BY placed_at").fetchall()
        return [{
            "client_order_id": row[0], "placed_at": _utc_from_ms(row[1]), "ticker": row[2], "contract_side": row[3],
            "status": row[4], "kalshi_order_id": row[5], "fills_final": row[6], "open_cost": row[7],
            "event_start": _utc_from_ms(row[8]), "all_in_price": row[9], "count_requested": row[10],
            "event_id": row[11], "exchange_index": row[12], "cancelled_at": _utc_from_ms(row[13]), "count_filled": row[14],
            "count_reported": row[15], "subaccount": row[16],
        } for row in rows]

    def open_exposure(self) -> float:
        """Dollars the bot has in markets not yet settled: the ceiling's measure."""
        with self._lock:
            return float(self._con.execute(
                "SELECT coalesce(sum(open_cost), 0) FROM taker_orders WHERE settled_at IS NULL").fetchone()[0])

    def count_decisions(self) -> int:
        with self._lock:
            return self._con.execute("SELECT count(*) FROM taker_decisions").fetchone()[0]

    # ---- Novig (novig_taker_desk.py) ----

    def insert_novig_order(self, order: dict) -> None:
        """INSERT a new Novig order row (status "sending", open_cost the full requested cost)."""
        values = [order[column] for column in NOVIG_ORDER_COLUMNS]
        placeholders = ", ".join("?" for _ in NOVIG_ORDER_COLUMNS)
        with self._lock:
            self._con.execute(f"INSERT INTO novig_taker_orders ({', '.join(NOVIG_ORDER_COLUMNS)}) VALUES ({placeholders})",
                              values)

    def has_novig_order(self, client_id: str) -> bool:
        with self._lock:
            return self._con.execute("SELECT count(*) FROM novig_taker_orders WHERE client_id = ?",
                                     [client_id]).fetchone()[0] > 0

    def update_novig_order(self, client_id: str, fields: dict) -> None:
        """UPDATE the named columns of one Novig order (this module's column names only)."""
        unknown = set(fields) - set(NOVIG_UPDATABLE_COLUMNS)
        if unknown:
            raise ValueError(f"update_novig_order: unknown columns {sorted(unknown)}")
        assignments = ", ".join(f"{column} = ?" for column in fields)
        with self._lock:
            self._con.execute(f"UPDATE novig_taker_orders SET {assignments} WHERE client_id = ?",
                              [*fields.values(), client_id])

    def open_novig_orders(self) -> list[dict]:
        """Every Novig order whose market has not been seen settled, oldest first."""
        selected = ", ".join(f"epoch_ms({column})" if column in _NOVIG_TIME_COLUMNS else column
                             for column in _NOVIG_OPEN_COLUMNS)
        with self._lock:
            rows = self._con.execute(f"SELECT {selected} FROM novig_taker_orders WHERE settled_at IS NULL "
                                     f"ORDER BY placed_at").fetchall()
        return [{column: (_utc_from_ms(value) if column in _NOVIG_TIME_COLUMNS else value)
                 for column, value in zip(_NOVIG_OPEN_COLUMNS, row)} for row in rows]

    def novig_open_exposure(self) -> float:
        """Dollars the bot has on Novig in markets not yet settled: the Novig ceiling's measure."""
        with self._lock:
            return float(self._con.execute(
                "SELECT coalesce(sum(open_cost), 0) FROM novig_taker_orders WHERE settled_at IS NULL").fetchone()[0])

    def novig_order_contexts(self) -> dict[str, dict]:
        """Novig order id -> the board line the bot priced it on (NOVIG_CONTEXT_COLUMNS), for
        every order Novig gave an id. The bot subaccount's read names its bets with these;
        `settled_at` (when the desk first saw its market SETTLED) is the bet's fixed closedAt."""
        selected = ", ".join(f"epoch_ms({column})" if column in _NOVIG_CONTEXT_TIME_COLUMNS else column
                             for column in NOVIG_CONTEXT_COLUMNS)
        with self._lock:
            rows = self._con.execute(f"SELECT novig_order_id, {selected} FROM novig_taker_orders "
                                     f"WHERE novig_order_id IS NOT NULL").fetchall()
        contexts = {}
        for row in rows:
            context = dict(zip(NOVIG_CONTEXT_COLUMNS, row[1:]))
            for column in _NOVIG_CONTEXT_TIME_COLUMNS:
                context[column] = _utc_from_ms(context[column])
            contexts[row[0]] = context
        return contexts

    # ---- the tracker's Bot tab (bot_report.py); read only ----

    def bot_activity(self, window_start: datetime | None, window_end: datetime | None) -> dict[str, dict]:
        """Per venue ("kalshi", "novig"): the bot's orders placed in [window_start,
        window_end) (either bound None = open) and, from the decision log, how
        many distinct lines took each decision then. A decision is logged only
        when it changes (taker.js decisionSignature), so `reasons` counts lines
        that reached that decision in the window, not how many passes saw it.
        A bet decision carries the desk's answer (`order.status`): one the
        desk refused (the ceiling, the stop file, ...) is counted under it."""
        # Only bet decisions carry an order, so the joins read the JSON of those rows alone: the
        # lock is the order desks' too, and a long read would hold up an order.
        with self._lock:
            return {
                "kalshi": {**self._order_totals("taker_orders", "client_order_id", "$.order.clientOrderId",
                                                "count_requested", "count_filled", window_start, window_end),
                           "reasons": self._decision_reasons("kalshi", window_start, window_end)},
                "novig": {**self._order_totals("novig_taker_orders", "client_id", "$.order.clientId",
                                               "qty_requested", "qty_filled", window_start, window_end),
                          "reasons": self._decision_reasons("novig", window_start, window_end)},
            }

    def _order_totals(self, table: str, id_column: str, decision_id_path: str, requested_column: str,
                      filled_column: str, window_start: datetime | None, window_end: datetime | None) -> dict:
        # Table and column names are this module's own literals, never caller text.
        where, params = _window_clause("o.placed_at", window_start, window_end)
        row = self._con.execute(
            f"SELECT count(*), count(*) FILTER (WHERE o.{filled_column} > 0), "
            f"coalesce(sum(o.{requested_column}) FILTER (WHERE o.{filled_column} IS NOT NULL), 0), "
            f"coalesce(sum(o.{filled_column}), 0), coalesce(sum(o.cost_filled), 0), "
            f"sum(o.cost_filled * d.live_edge) FILTER (WHERE o.cost_filled > 0 AND d.live_edge IS NOT NULL), "
            f"sum(o.cost_filled) FILTER (WHERE o.cost_filled > 0 AND d.live_edge IS NOT NULL) "
            f"FROM {table} o LEFT JOIN ("
            f"  SELECT json_extract_string(decision, '{decision_id_path}') AS order_id, any_value(live_edge) AS live_edge "
            f"  FROM taker_decisions WHERE reason = 'bet' AND json_extract_string(decision, '{decision_id_path}') IS NOT NULL GROUP BY 1"
            f") d ON d.order_id = o.{id_column} WHERE {where}", params).fetchone()
        orders, filled_orders, requested_read, filled, cost, edge_dollars, edge_cost = row
        return {
            "orders": orders, "filledOrders": filled_orders,
            # Fill rate over the orders whose fills have been read (a NULL fill is not yet known, never 0).
            "fillRate": float(filled) / float(requested_read) if requested_read else None,
            "staked": float(cost),
            "avgEdgeAtFill": float(edge_dollars) / float(edge_cost) if edge_cost else None,
        }

    def _decision_reasons(self, venue: str, window_start: datetime | None, window_end: datetime | None) -> list[dict]:
        where, params = _window_clause("decided_at", window_start, window_end)
        # Decisions logged before the Novig taker (2026-10-16) carry no venue: they are all Kalshi's.
        rows = self._con.execute(
            "SELECT reason, json_extract_string(decision, '$.order.status') AS order_status, count(DISTINCT row_key) "
            "FROM taker_decisions WHERE coalesce(json_extract_string(decision, '$.venue'), 'kalshi') = ? "
            f"AND {where} GROUP BY 1, 2 ORDER BY 3 DESC", [venue, *params]).fetchall()
        return [{"reason": reason, "orderStatus": order_status, "lines": lines} for reason, order_status, lines in rows]

    def bot_open_orders(self) -> dict[str, dict]:
        """Per venue: orders whose market is not yet seen settled and still hold
        cost against the ceiling (open_cost > 0), and when the bot last placed one."""
        with self._lock:
            out = {}
            for venue, table in (("kalshi", "taker_orders"), ("novig", "novig_taker_orders")):
                count, last_ms = self._con.execute(
                    f"SELECT count(*) FILTER (WHERE settled_at IS NULL AND open_cost > 0), epoch_ms(max(placed_at)) "
                    f"FROM {table}").fetchone()
                out[venue] = {"openOrders": count, "lastOrderAt": _utc_from_ms(last_ms)}
            return out

    def recent_bot_orders(self, limit: int) -> list[dict]:
        """The newest `limit` bot orders on either venue, with the line the
        decision that sent it named (eventName, period, betType, sideLabel,
        points; null when that decision is not in the log)."""
        with self._lock:
            rows = self._con.execute(
                "WITH lines AS ("
                "  SELECT coalesce(json_extract_string(decision, '$.order.clientOrderId'), "
                "                  json_extract_string(decision, '$.order.clientId')) AS order_id, "
                "         any_value(decision) AS decision, any_value(live_edge) AS live_edge "
                "  FROM taker_decisions WHERE reason = 'bet' AND json_extract_string(decision, '$.order') IS NOT NULL GROUP BY 1"
                "), orders AS ("
                "  SELECT 'kalshi' AS venue, client_order_id AS order_id, placed_at, ticker AS contract, "
                "         all_in_price AS price, count_requested AS requested, count_filled AS filled, cost_filled, status, error "
                "  FROM taker_orders "
                "  UNION ALL "
                "  SELECT 'novig', client_id, placed_at, outcome_id, price, qty_requested, qty_filled, cost_filled, status, error "
                "  FROM novig_taker_orders"
                ") "
                "SELECT o.venue, o.order_id, epoch_ms(o.placed_at), o.contract, o.price, o.requested, o.filled, o.cost_filled, "
                "       o.status, o.error, l.live_edge, json_extract_string(l.decision, '$.eventName'), "
                "       json_extract_string(l.decision, '$.league'), json_extract_string(l.decision, '$.period'), "
                "       json_extract_string(l.decision, '$.betType'), json_extract_string(l.decision, '$.sideLabel'), "
                "       TRY_CAST(json_extract(l.decision, '$.points') AS DOUBLE) "
                "FROM orders o LEFT JOIN lines l ON l.order_id = o.order_id ORDER BY o.placed_at DESC LIMIT ?",
                [limit]).fetchall()
        keys = ("venue", "orderId", "placedAt", "contract", "price", "requested", "filled", "cost", "status", "error",
                "liveEdge", "eventName", "league", "period", "betType", "sideLabel", "points")
        orders = []
        for row in rows:
            order = dict(zip(keys, row))
            order["placedAt"] = _utc_from_ms(order["placedAt"])
            orders.append(order)
        return orders


def _window_clause(column: str, window_start: datetime | None, window_end: datetime | None) -> tuple[str, list]:
    """A SQL condition on `column` (this module's own name) for [start, end), and its parameters."""
    parts, params = ["TRUE"], []
    if window_start is not None:
        parts.append(f"{column} >= ?")
        params.append(window_start)
    if window_end is not None:
        parts.append(f"{column} < ?")
        params.append(window_end)
    return " AND ".join(parts), params
