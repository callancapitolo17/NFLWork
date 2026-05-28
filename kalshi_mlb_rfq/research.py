"""Research firehose for the Kalshi MLB RFQ bot.

Captures the full RFQ lifecycle (candidates considered, per-book fairs, gate
internals, Kelly math, position changes) that the trading path computes and
discards. Writes to a SEPARATE DuckDB file (its own write lock — zero
contention with the trading DB). Events are buffered in-process and flushed
in one batched transaction per main-loop tick.

SAFETY INVARIANT: nothing here may ever raise into the trading loop. emit()
only appends to a list; flush() swallows every exception.
"""

import json
import logging
import time
import uuid
from datetime import datetime, timezone
from random import random

import duckdb

from kalshi_mlb_rfq import config

log = logging.getLogger("kalshi_mlb_rfq.research")

_BUFFER: list[dict] = []
_SESSION_ID: str | None = None
_LAST_FLUSH_WARN = 0.0
FLUSH_WARN_INTERVAL_SEC = 60.0   # rate-limit failure warnings (flush runs ~2x/sec)

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS events (
    event_id     VARCHAR PRIMARY KEY,
    session_id   VARCHAR,
    event_type   VARCHAR NOT NULL,
    ts           TIMESTAMPTZ NOT NULL,
    game_id      VARCHAR,
    combo_ticker VARCHAR,
    rfq_id       VARCHAR,
    quote_id     VARCHAR,
    payload      JSON
);
CREATE INDEX IF NOT EXISTS idx_events_type_ts ON events(event_type, ts);
CREATE INDEX IF NOT EXISTS idx_events_game    ON events(game_id);
CREATE INDEX IF NOT EXISTS idx_events_combo   ON events(combo_ticker);
"""


def _connect(read_only: bool = False, retries: int = 5):
    """Open the research DB with retry-backoff (mirrors db.connect)."""
    last_err = None
    for attempt in range(retries):
        try:
            return duckdb.connect(str(config.RESEARCH_DB_PATH),
                                  read_only=read_only)
        except duckdb.IOException as e:
            last_err = e
            time.sleep(0.05 * (2 ** attempt) + random() * 0.05)
    raise last_err


def init_research_db() -> None:
    """Idempotently create the research DB + events table."""
    con = _connect()
    try:
        con.execute(SCHEMA_SQL)
    finally:
        con.close()


def set_session(session_id: str) -> None:
    global _SESSION_ID
    _SESSION_ID = session_id


def emit(event_type: str, *, game_id: str | None = None,
         combo_ticker: str | None = None, rfq_id: str | None = None,
         quote_id: str | None = None, **payload) -> None:
    """Append one research event to the buffer. O(1). NEVER raises.

    `payload` holds event-specific fields; it is JSON-serialized at flush
    time. A non-serializable value is coerced to its repr so emit cannot
    throw into the trading loop.
    """
    try:
        _BUFFER.append({
            "event_id": str(uuid.uuid4()),
            "session_id": _SESSION_ID,
            "event_type": event_type,
            "ts": datetime.now(timezone.utc),
            "game_id": game_id,
            "combo_ticker": combo_ticker,
            "rfq_id": rfq_id,
            "quote_id": quote_id,
            "payload": payload,
        })
        cap = config.RESEARCH_BUFFER_MAX
        if len(_BUFFER) > cap:
            del _BUFFER[:len(_BUFFER) - cap]   # drop oldest
    except Exception:  # pragma: no cover — emit must never break trading
        pass


def _json_default(o):
    return repr(o)


def flush() -> None:
    """Write the buffered events to the research DB in one batched
    transaction, then clear the buffer. Swallows all errors — research
    logging must never break trading. On write failure the buffer is
    retained (capped by emit) so the next flush can retry. The failure
    warning is rate-limited (flush runs ~2x/sec).
    """
    if not _BUFFER:
        return
    batch = list(_BUFFER)
    rows = []
    for e in batch:
        rows.append([
            e["event_id"], e["session_id"], e["event_type"], e["ts"],
            e["game_id"], e["combo_ticker"], e["rfq_id"], e["quote_id"],
            json.dumps(e["payload"], default=_json_default),
        ])
    try:
        con = _connect()
        try:
            con.executemany(
                "INSERT INTO events (event_id, session_id, event_type, ts, "
                "game_id, combo_ticker, rfq_id, quote_id, payload) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?) ON CONFLICT DO NOTHING",
                rows,
            )
        finally:
            con.close()
        del _BUFFER[:len(batch)]   # only clear what we successfully wrote
    except Exception as exc:
        global _LAST_FLUSH_WARN
        now = time.time()
        if now - _LAST_FLUSH_WARN >= FLUSH_WARN_INTERVAL_SEC:
            log.warning("research flush failed (%d events retained): %s",
                        len(batch), exc)
            _LAST_FLUSH_WARN = now


def prune_research(days: int = 90) -> int:
    """Delete events older than `days`. Returns rows deleted.

    NOT scheduled anywhere in v1 — call it from a cron line when the DB
    gets chunky. Retention is otherwise unbounded by design.
    """
    con = _connect()
    try:
        before = con.execute("SELECT count(*) FROM events").fetchone()[0]
        con.execute(
            "DELETE FROM events WHERE ts < now() - "
            "(INTERVAL 1 DAY * ?)", [days])
        after = con.execute("SELECT count(*) FROM events").fetchone()[0]
        return before - after
    finally:
        con.close()
