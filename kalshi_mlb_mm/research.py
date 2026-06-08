"""Research firehose for the Kalshi MLB MM (maker) bot.

Captures the full maker lifecycle (RFQs considered, quotes priced, decisions,
fills, reconcile outcomes, halts, circuit breakers) that the trading path
computes and discards. Writes to a SEPARATE DuckDB file (its own write lock —
zero contention with the trading DB). Events are buffered in-process and
flushed in one batched transaction per main-loop tick.

SAFETY INVARIANT: nothing here may ever raise into the trading loop. emit()
only appends to a list; flush() swallows every exception.

Mirrors kalshi_mlb_rfq/research.py exactly in structure and contract.
"""

import json
import logging
import threading
import time
import uuid
from datetime import datetime, timezone
from random import random

import duckdb

from kalshi_mlb_mm import config

log = logging.getLogger("kalshi_mlb_mm.research")

_BUFFER: list[dict] = []
_BUFFER_LOCK = threading.Lock()
_SESSION_ID: str | None = None
_LAST_FLUSH_WARN = 0.0
FLUSH_WARN_INTERVAL_SEC = float(config.RESEARCH_FLUSH_WARN_RATE_LIMIT_SEC)

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS events (
    event_id     VARCHAR PRIMARY KEY,
    session_id   VARCHAR,
    event_type   VARCHAR NOT NULL,
    ts           TIMESTAMPTZ NOT NULL,
    ticker       VARCHAR,
    rfq_id       VARCHAR,
    quote_id     VARCHAR,
    payload      JSON
);
CREATE INDEX IF NOT EXISTS idx_events_type ON events(event_type);
CREATE INDEX IF NOT EXISTS idx_events_ts   ON events(ts);
"""


def _connect(read_only: bool = False, retries: int = 5):
    """Open the research DB with retry-backoff (mirrors taker research._connect)."""
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
    """Idempotently create the research DB + events table. Fail-fast on bad path."""
    con = _connect()
    try:
        con.execute(SCHEMA_SQL)
    finally:
        con.close()


def set_session_id(sid: str) -> None:
    """Store the session ID so all subsequent emit() calls stamp it."""
    global _SESSION_ID
    _SESSION_ID = sid


def emit(event_type: str, *, ticker: str | None = None,
         rfq_id: str | None = None, quote_id: str | None = None,
         payload: dict | None = None, **kwargs) -> None:
    """Append one research event to the buffer. O(1). NEVER raises.

    `payload` holds event-specific fields; extra kwargs are merged into it.
    JSON-serialization happens at flush time. A non-serializable value is
    coerced to its repr() so emit cannot throw into the trading loop.
    """
    try:
        merged = {}
        if payload:
            merged.update(payload)
        if kwargs:
            merged.update(kwargs)
        _BUFFER.append({
            "event_id": str(uuid.uuid4()),
            "session_id": _SESSION_ID,
            "event_type": event_type,
            "ts": datetime.now(timezone.utc),
            "ticker": ticker,
            "rfq_id": rfq_id,
            "quote_id": quote_id,
            "payload": merged,
        })
        cap = config.RESEARCH_BUFFER_MAX
        if len(_BUFFER) > cap:
            del _BUFFER[:len(_BUFFER) - cap]   # drop oldest
    except Exception:  # pragma: no cover — emit must never break trading
        pass


def _json_default(o):
    return repr(o)


def flush() -> int:
    """Write buffered events to the research DB in one batched transaction.
    Returns count flushed. Swallows all errors — research logging must never
    break trading. On write failure the buffer is retained (capped by emit)
    so the next flush can retry. The failure warning is rate-limited (flush
    runs every tick). Uses retries=1 (fail-fast: a locked research DB can't
    stall the trading tick).
    """
    if not _BUFFER:
        return 0
    try:
        batch = list(_BUFFER)
        rows = [
            [e["event_id"], e["session_id"], e["event_type"], e["ts"],
             e["ticker"], e["rfq_id"], e["quote_id"],
             json.dumps(e["payload"], default=_json_default)]
            for e in batch
        ]
        con = _connect(retries=1)
        try:
            con.executemany(
                "INSERT INTO events (event_id, session_id, event_type, ts, "
                "ticker, rfq_id, quote_id, payload) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?) ON CONFLICT DO NOTHING",
                rows,
            )
        finally:
            con.close()
        del _BUFFER[:len(batch)]   # only clear what we successfully wrote
        return len(batch)
    except Exception as exc:
        global _LAST_FLUSH_WARN
        now = time.time()
        if now - _LAST_FLUSH_WARN >= FLUSH_WARN_INTERVAL_SEC:
            log.warning("research flush failed (%d events retained): %s",
                        len(_BUFFER), exc)
            _LAST_FLUSH_WARN = now
        return 0


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
