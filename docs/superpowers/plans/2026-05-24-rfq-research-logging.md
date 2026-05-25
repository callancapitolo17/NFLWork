# RFQ Research Logging Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the Kalshi MLB RFQ taker bot end-to-end observability of the RFQ lifecycle — capture every candidate's pricing, gate, and sizing internals that are currently discarded — plus structured rotating logs, without adding write-lock pressure to the trading DB.

**Architecture:** A new sibling `kalshi_mlb_rfq_research.duckdb` (Approach 2b), written only by the bot, batched once per main-loop tick via a `research.py` event buffer (`emit`/`flush`). One wide `events` table with a JSON `payload` column. Operational `print()` is replaced with Python `logging` + a size-capped `RotatingFileHandler`. Research logging failures can never propagate into the trading loop.

**Tech Stack:** Python 3.11, DuckDB, pytest. Package `kalshi_mlb_rfq` (imported as `kalshi_mlb_rfq.<module>`); tests in `kalshi_mlb_rfq/tests/`, run with `pytest` from the worktree root.

**Spec:** `docs/superpowers/specs/2026-05-23-rfq-research-logging-design.md`

---

## ⚠️ Churn warning (read before starting)

`kalshi_mlb_rfq/main.py` is under heavy active development — `main` advanced ~20 commits (NO-side accepts, two-RFQ per-side Kelly) during the spec session alone. Phases 0, 1, and 4 touch **stable, self-contained** code and should be done and merged first. Phases 2–3 add `emit()` calls into the churning hot path (`_enumerate_and_score_all_games`, `_evaluate_quote`, gate/Kelly helpers); do them **last and quickly**, and before applying each, re-read the target function in case its lines shifted. The line numbers below are accurate as of `main` @ `5bd5e6c`.

---

## File Structure

**New files:**
- `kalshi_mlb_rfq/log_setup.py` — one function, `setup_logging()`, configures the root logger (stdout + rotating file). One responsibility: logging configuration.
- `kalshi_mlb_rfq/research.py` — the research event store: `init_research_db()`, `set_session()`, `emit()`, `flush()`, `prune_research()`, and a module-private buffer + connection helper. One responsibility: capture/persist research events off the trading lock.
- `kalshi_mlb_rfq/tests/test_research.py` — unit tests for `research.py`.
- `kalshi_mlb_rfq/tests/test_log_setup.py` — unit tests for `log_setup.py`.

**Modified files:**
- `kalshi_mlb_rfq/config.py` — new knobs (logging + research).
- `kalshi_mlb_rfq/main.py` — wire `setup_logging()`, `init_research_db()`, `set_session()`, `flush()`; replace `print()`; add `emit()` hooks (Phases 2–3).
- `kalshi_mlb_rfq/README.md`, repo `CLAUDE.md`, `.gitignore` — docs + ignore rules (Phase 4).

---

# PHASE 0 — Operational logging (gap C)

Self-contained. No dependency on the research store. Ship/merge first.

### Task 1: Logging config knobs

**Files:**
- Modify: `kalshi_mlb_rfq/config.py` (append near the end, after line ~104)
- Test: `kalshi_mlb_rfq/tests/test_config.py`

- [ ] **Step 1: Write the failing test** — append to `test_config.py`:

```python
def test_logging_knobs_have_defaults():
    assert config_mod.LOG_MAX_BYTES == 50 * 1024 * 1024
    assert config_mod.LOG_BACKUP_COUNT == 5
    assert config_mod.LOG_LEVEL == "INFO"
```

- [ ] **Step 2: Run test, verify it fails**

Run: `pytest kalshi_mlb_rfq/tests/test_config.py::test_logging_knobs_have_defaults -v`
Expected: FAIL — `AttributeError: module ... has no attribute 'LOG_MAX_BYTES'`

- [ ] **Step 3: Add the knobs** — append to `config.py` (after the `MIN_BOOK_COUNT_FOR_BLEND` line):

```python
# Logging (operational log rotation)
LOG_LEVEL = _get("LOG_LEVEL", "INFO")
LOG_MAX_BYTES = int(_get("LOG_MAX_BYTES", str(50 * 1024 * 1024)))  # 50 MB
LOG_BACKUP_COUNT = int(_get("LOG_BACKUP_COUNT", "5"))
```

- [ ] **Step 4: Run test, verify it passes**

Run: `pytest kalshi_mlb_rfq/tests/test_config.py -v`
Expected: PASS (all config tests)

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/config.py kalshi_mlb_rfq/tests/test_config.py
git commit -m "feat(rfq-logging): add log-rotation config knobs"
```

### Task 2: `log_setup.py`

**Files:**
- Create: `kalshi_mlb_rfq/log_setup.py`
- Test: `kalshi_mlb_rfq/tests/test_log_setup.py`

- [ ] **Step 1: Write the failing test** — create `kalshi_mlb_rfq/tests/test_log_setup.py`:

```python
"""Tests for kalshi_mlb_rfq/log_setup.py."""
import logging
from logging.handlers import RotatingFileHandler

from kalshi_mlb_rfq.log_setup import setup_logging


def test_setup_logging_installs_rotating_and_stream_handlers(tmp_path):
    log_file = tmp_path / "bot.log"
    logger = setup_logging(log_path=log_file, max_bytes=1024,
                           backup_count=3, level="DEBUG")
    handler_types = {type(h) for h in logger.handlers}
    assert RotatingFileHandler in handler_types
    assert logging.StreamHandler in handler_types
    rfh = next(h for h in logger.handlers
               if isinstance(h, RotatingFileHandler))
    assert rfh.maxBytes == 1024
    assert rfh.backupCount == 3
    assert logger.level == logging.DEBUG


def test_setup_logging_is_idempotent(tmp_path):
    """Calling twice must not double the handlers (no duplicate log lines)."""
    log_file = tmp_path / "bot.log"
    setup_logging(log_path=log_file)
    logger = setup_logging(log_path=log_file)
    rfh_count = sum(1 for h in logger.handlers
                    if isinstance(h, RotatingFileHandler))
    assert rfh_count == 1
```

- [ ] **Step 2: Run test, verify it fails**

Run: `pytest kalshi_mlb_rfq/tests/test_log_setup.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'kalshi_mlb_rfq.log_setup'`

- [ ] **Step 3: Create `kalshi_mlb_rfq/log_setup.py`:**

```python
"""Central logging configuration for the Kalshi MLB RFQ bot.

Replaces ad-hoc print() with the stdlib logging module: levels (filterable
without code changes) + a size-capped RotatingFileHandler (so bot.log can
never grow unbounded). Call setup_logging() once at process start.
"""

import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path

from kalshi_mlb_rfq import config

_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"


def setup_logging(log_path: Path | None = None,
                  max_bytes: int | None = None,
                  backup_count: int | None = None,
                  level: str | None = None) -> logging.Logger:
    """Configure the root logger with a stdout handler + a rotating file
    handler. Idempotent: repeated calls do not stack handlers.

    Defaults come from config (LOG_PATH, LOG_MAX_BYTES, LOG_BACKUP_COUNT,
    LOG_LEVEL); args override for tests.
    """
    log_path = Path(log_path) if log_path else config.LOG_PATH
    max_bytes = max_bytes if max_bytes is not None else config.LOG_MAX_BYTES
    backup_count = (backup_count if backup_count is not None
                    else config.LOG_BACKUP_COUNT)
    level = level or config.LOG_LEVEL

    root = logging.getLogger()
    root.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Idempotency: drop our previously-installed handlers before re-adding.
    for h in list(root.handlers):
        if getattr(h, "_rfq_managed", False):
            root.removeHandler(h)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter(_FORMAT)

    file_handler = RotatingFileHandler(
        log_path, maxBytes=max_bytes, backupCount=backup_count)
    file_handler.setFormatter(fmt)
    file_handler._rfq_managed = True
    root.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(fmt)
    stream_handler._rfq_managed = True
    root.addHandler(stream_handler)

    return root
```

- [ ] **Step 4: Run test, verify it passes**

Run: `pytest kalshi_mlb_rfq/tests/test_log_setup.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/log_setup.py kalshi_mlb_rfq/tests/test_log_setup.py
git commit -m "feat(rfq-logging): add setup_logging with rotating file handler"
```

### Task 3: Replace `print()` with logging in `main.py`

This is a mechanical refactor (no behavior change). There are ~40 `print(..., flush=True)` calls plus `notify.py`'s file writes.

**Files:**
- Modify: `kalshi_mlb_rfq/main.py`
- Modify: `kalshi_mlb_rfq/notify.py`

- [ ] **Step 1: Add a module logger + setup call to `main.py`.** Near the top of `main.py` (after the existing imports), add:

```python
import logging
from kalshi_mlb_rfq.log_setup import setup_logging

log = logging.getLogger("kalshi_mlb_rfq")
```

In `main_loop()` (line ~1492), make `setup_logging()` the **first** statement, before `db.init_database()`:

```python
def main_loop(dry_run: bool):
    setup_logging()
    db.init_database()
    ...
```

- [ ] **Step 2: Replace each `print(...)` with the matching log level.** Apply this mapping (drop the `, flush=True`):
  - Normal beats → `log.info(...)`: session banner, `cache_refresh:` counts, `sgp_cache_refresh:`, `sgp_cycle:`, `rfq_refresh:` counts, `[HB] ... alive`, `startup:` messages, `shutdown:` messages, `=== shutdown complete ===`.
  - Recoverable errors → `log.warning(...)`: `rfq_refresh error`, `quote_poll error`, `risk_sweep error`, `sgp_cycle error`, `cache_refresh: gave up ...`, `add ... failed`, `drop ... failed`, `shutdown: cancel ... failed`.

  Example transform (line ~1541):

```python
# before:
print(f"  rfq_refresh: {len(candidates)} candidates ({time.time()-t_ref:.1f}s)", flush=True)
# after:
log.info("rfq_refresh: %d candidates (%.1fs)", len(candidates), time.time()-t_ref)
```

```python
# before:
print(f"  rfq_refresh error: {e}", flush=True)
# after:
log.warning("rfq_refresh error: %s", e)
```

- [ ] **Step 3: Route `notify.py` file writes through logging.** In `kalshi_mlb_rfq/notify.py`, replace `_append_log` with a logger call (keep the webhook unchanged):

```python
import logging
log = logging.getLogger("kalshi_mlb_rfq.notify")

def _append_log(line: str):
    log.info(line)
```

(Remove the now-unused `LOG_PATH` import and the file-open code.)

- [ ] **Step 4: Verify no behavior change.** Confirm no stray prints remain and the module imports:

Run: `grep -n "print(" kalshi_mlb_rfq/main.py kalshi_mlb_rfq/notify.py`
Expected: no matches (or only inside comments/docstrings).
Run: `python -c "import kalshi_mlb_rfq.main"`
Expected: no error.
Run: `pytest kalshi_mlb_rfq/tests/ -v`
Expected: PASS (existing tests unaffected).

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/main.py kalshi_mlb_rfq/notify.py
git commit -m "refactor(rfq-logging): replace print() with structured logging"
```

---

# PHASE 1 — Research store plumbing (Section 1 of spec)

Creates the firehose infrastructure and proves it writes safely without emitting real events yet.

### Task 4: Research config knobs

**Files:**
- Modify: `kalshi_mlb_rfq/config.py`
- Test: `kalshi_mlb_rfq/tests/test_config.py`

- [ ] **Step 1: Write the failing test** — append to `test_config.py`:

```python
def test_research_knobs_have_defaults():
    assert str(config_mod.RESEARCH_DB_PATH).endswith(
        "kalshi_mlb_rfq_research.duckdb")
    assert config_mod.RESEARCH_CANDIDATE_SAMPLING == 1.0
    assert config_mod.RESEARCH_BUFFER_MAX == 50000
```

- [ ] **Step 2: Run test, verify it fails**

Run: `pytest kalshi_mlb_rfq/tests/test_config.py::test_research_knobs_have_defaults -v`
Expected: FAIL — `AttributeError: ... 'RESEARCH_DB_PATH'`

- [ ] **Step 3: Add the knobs** — append to `config.py`:

```python
# Research firehose (separate sibling DB, off the trading write lock)
RESEARCH_DB_PATH = Path(_get("RESEARCH_DB_PATH",
                             str(PKG_DIR / "kalshi_mlb_rfq_research.duckdb")))
# Fraction of candidate_evaluated events to keep per tick (1.0 = full movie).
RESEARCH_CANDIDATE_SAMPLING = float(_get("RESEARCH_CANDIDATE_SAMPLING", "1.0"))
# Buffer cap: if flush keeps failing, drop oldest beyond this to bound memory.
RESEARCH_BUFFER_MAX = int(_get("RESEARCH_BUFFER_MAX", "50000"))
```

- [ ] **Step 4: Run test, verify it passes**

Run: `pytest kalshi_mlb_rfq/tests/test_config.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/config.py kalshi_mlb_rfq/tests/test_config.py
git commit -m "feat(rfq-logging): add research firehose config knobs"
```

### Task 5: `research.py` — schema + `init_research_db()`

**Files:**
- Create: `kalshi_mlb_rfq/research.py`
- Test: `kalshi_mlb_rfq/tests/test_research.py`

- [ ] **Step 1: Write the failing test** — create `kalshi_mlb_rfq/tests/test_research.py`:

```python
"""Tests for kalshi_mlb_rfq/research.py — the research event firehose."""
import duckdb
import pytest

import kalshi_mlb_rfq.research as research


@pytest.fixture(autouse=True)
def isolated_research_db(tmp_path, monkeypatch):
    """Point research at a temp DB and reset module state per test."""
    db = tmp_path / "research.duckdb"
    monkeypatch.setattr(research.config, "RESEARCH_DB_PATH", db)
    research._BUFFER.clear()
    research._SESSION_ID = None
    yield db


def test_init_creates_events_table(isolated_research_db):
    research.init_research_db()
    con = duckdb.connect(str(isolated_research_db), read_only=True)
    cols = {r[0] for r in con.execute(
        "SELECT column_name FROM information_schema.columns "
        "WHERE table_name='events'").fetchall()}
    con.close()
    assert {"event_id", "session_id", "event_type", "ts",
            "game_id", "combo_ticker", "rfq_id", "quote_id",
            "payload"}.issubset(cols)


def test_init_is_idempotent(isolated_research_db):
    research.init_research_db()
    research.init_research_db()  # must not raise
```

- [ ] **Step 2: Run test, verify it fails**

Run: `pytest kalshi_mlb_rfq/tests/test_research.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'kalshi_mlb_rfq.research'`

- [ ] **Step 3: Create `kalshi_mlb_rfq/research.py`** (this task only adds the schema + init + connect helper; emit/flush come next):

```python
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
```

- [ ] **Step 4: Run test, verify it passes**

Run: `pytest kalshi_mlb_rfq/tests/test_research.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/research.py kalshi_mlb_rfq/tests/test_research.py
git commit -m "feat(rfq-logging): research DB schema + init_research_db"
```

### Task 6: `set_session()` + `emit()`

**Files:**
- Modify: `kalshi_mlb_rfq/research.py`
- Test: `kalshi_mlb_rfq/tests/test_research.py`

- [ ] **Step 1: Write the failing test** — append to `test_research.py`:

```python
def test_emit_appends_to_buffer(isolated_research_db):
    research.set_session("sess-1")
    research.emit("candidate_evaluated", game_id="g1",
                  combo_ticker="KX-X", model_fair=0.4, edge=0.03)
    assert len(research._BUFFER) == 1
    ev = research._BUFFER[0]
    assert ev["event_type"] == "candidate_evaluated"
    assert ev["session_id"] == "sess-1"
    assert ev["game_id"] == "g1"
    payload = ev["payload"]  # dict at buffer stage
    assert payload["model_fair"] == 0.4 and payload["edge"] == 0.03


def test_emit_never_raises_on_bad_payload(isolated_research_db):
    # Non-JSON-serializable payload must not raise out of emit.
    class Weird:
        pass
    research.emit("candidate_evaluated", obj=Weird())  # must not raise
    assert len(research._BUFFER) == 1


def test_emit_respects_buffer_cap(isolated_research_db, monkeypatch):
    monkeypatch.setattr(research.config, "RESEARCH_BUFFER_MAX", 3)
    for i in range(5):
        research.emit("candidate_evaluated", i=i)
    assert len(research._BUFFER) == 3
    # Oldest dropped: the surviving events are the last three (i=2,3,4).
    assert [e["payload"]["i"] for e in research._BUFFER] == [2, 3, 4]
```

- [ ] **Step 2: Run test, verify it fails**

Run: `pytest kalshi_mlb_rfq/tests/test_research.py -k "emit or session" -v`
Expected: FAIL — `AttributeError: module ... has no attribute 'set_session'`

- [ ] **Step 3: Implement `set_session` + `emit`** — append to `research.py`:

```python
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
```

Note: `emit` stores `payload` as a dict in the buffer; `flush` (next task) handles JSON serialization, with a fallback for non-serializable values.

- [ ] **Step 4: Run test, verify it passes**

Run: `pytest kalshi_mlb_rfq/tests/test_research.py -k "emit or session" -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/research.py kalshi_mlb_rfq/tests/test_research.py
git commit -m "feat(rfq-logging): research emit() + session stamping + buffer cap"
```

### Task 7: `flush()` — batched write, error-swallowing

**Files:**
- Modify: `kalshi_mlb_rfq/research.py`
- Test: `kalshi_mlb_rfq/tests/test_research.py`

- [ ] **Step 1: Write the failing test** — append to `test_research.py`:

```python
def test_flush_writes_batch_and_clears(isolated_research_db):
    research.init_research_db()
    research.set_session("sess-1")
    research.emit("candidate_evaluated", game_id="g1", model_fair=0.4)
    research.emit("gate_evaluated", game_id="g1", gate="tipoff")
    research.flush()
    assert research._BUFFER == []          # buffer cleared
    con = duckdb.connect(str(isolated_research_db), read_only=True)
    n = con.execute("SELECT count(*) FROM events").fetchone()[0]
    fair = con.execute(
        "SELECT payload->>'model_fair' FROM events "
        "WHERE event_type='candidate_evaluated'").fetchone()[0]
    con.close()
    assert n == 2
    assert float(fair) == 0.4


def test_flush_swallows_db_error_and_keeps_running(isolated_research_db,
                                                   monkeypatch, caplog):
    research.set_session("sess-1")
    research.emit("candidate_evaluated", game_id="g1")

    def boom(*a, **k):
        raise duckdb.IOException("locked")
    monkeypatch.setattr(research, "_connect", boom)

    research.flush()   # MUST NOT raise
    assert "research flush failed" in caplog.text.lower()


def test_flush_on_empty_buffer_is_noop(isolated_research_db):
    research.flush()   # must not raise, must not create rows
```

- [ ] **Step 2: Run test, verify it fails**

Run: `pytest kalshi_mlb_rfq/tests/test_research.py -k flush -v`
Expected: FAIL — `AttributeError: ... 'flush'`

- [ ] **Step 3: Implement `flush`** — append to `research.py`:

```python
def _json_default(o):
    return repr(o)


def flush() -> None:
    """Write the buffered events to the research DB in one batched
    transaction, then clear the buffer. Swallows all errors — research
    logging must never break trading. On write failure the buffer is
    retained (capped by emit) so the next flush can retry.
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
        log.warning("research flush failed (%d events retained): %s",
                    len(batch), exc)
```

- [ ] **Step 4: Run test, verify it passes**

Run: `pytest kalshi_mlb_rfq/tests/test_research.py -k flush -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/research.py kalshi_mlb_rfq/tests/test_research.py
git commit -m "feat(rfq-logging): research flush() — batched, error-swallowing"
```

### Task 8: `prune_research()` (written, unscheduled)

**Files:**
- Modify: `kalshi_mlb_rfq/research.py`
- Test: `kalshi_mlb_rfq/tests/test_research.py`

- [ ] **Step 1: Write the failing test** — append to `test_research.py`:

```python
def test_prune_deletes_only_old_rows(isolated_research_db):
    research.init_research_db()
    con = duckdb.connect(str(isolated_research_db))
    con.execute(
        "INSERT INTO events (event_id, event_type, ts, payload) VALUES "
        "('old', 't', now() - INTERVAL 100 DAY, '{}'), "
        "('new', 't', now() - INTERVAL 1 DAY, '{}')")
    con.close()
    deleted = research.prune_research(days=90)
    assert deleted == 1
    con = duckdb.connect(str(isolated_research_db), read_only=True)
    ids = {r[0] for r in con.execute("SELECT event_id FROM events").fetchall()}
    con.close()
    assert ids == {"new"}
```

- [ ] **Step 2: Run test, verify it fails**

Run: `pytest kalshi_mlb_rfq/tests/test_research.py -k prune -v`
Expected: FAIL — `AttributeError: ... 'prune_research'`

- [ ] **Step 3: Implement `prune_research`** — append to `research.py`:

```python
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
```

- [ ] **Step 4: Run test, verify it passes**

Run: `pytest kalshi_mlb_rfq/tests/test_research.py -k prune -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/research.py kalshi_mlb_rfq/tests/test_research.py
git commit -m "feat(rfq-logging): prune_research helper (unscheduled)"
```

### Task 9: Wire research store into `main_loop`

**Files:**
- Modify: `kalshi_mlb_rfq/main.py` (`main_loop`, lines ~1491–1613)

- [ ] **Step 1: Add the import.** Near the other `kalshi_mlb_rfq` imports in `main.py`:

```python
from kalshi_mlb_rfq import research
```

- [ ] **Step 2: Init + session at startup.** In `main_loop`, right after `sid = db.start_session(...)` (line ~1493):

```python
    research.init_research_db()
    research.set_session(sid)
```

- [ ] **Step 3: Flush once per tick.** In the `while _running.is_set():` loop, immediately before `time.sleep(0.5)` (line ~1588):

```python
            research.flush()
            time.sleep(0.5)
```

- [ ] **Step 4: Final flush on shutdown.** In the `finally:` block (after `db.end_session(sid)`, line ~1612):

```python
        research.flush()   # persist the last tick's buffered events
        db.end_session(sid)
        log.info("=== shutdown complete ===")
```

- [ ] **Step 5: Verify import + tests still green.**

Run: `python -c "import kalshi_mlb_rfq.main"`
Expected: no error.
Run: `pytest kalshi_mlb_rfq/tests/ -v`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add kalshi_mlb_rfq/main.py
git commit -m "feat(rfq-logging): wire research init/flush into main loop"
```

**→ Merge checkpoint A:** Phases 0–1 are a coherent, mergeable unit (real logging + a working firehose that emits nothing yet). Consider an executive-engineer review + merge here before the churn-sensitive Phases 2–3.

---

# PHASE 2 — `candidate_evaluated` (gaps A + B)

⚠️ Churn-sensitive. Re-read `_enumerate_and_score_all_games` (currently lines ~1419–1484) before editing.

### Task 10: Emit `candidate_evaluated` at every enumeration outcome

**Files:**
- Modify: `kalshi_mlb_rfq/main.py` (`_enumerate_and_score_all_games`)

The per-book `books` dict, `model`, `blended`, `kalshi_ref`, and the per-side Kelly tuple are all already in scope in the loop — no refactor of `_load_book_fairs` is needed. Add an `emit` at each `continue` (rejection) and at the accept point. Sampling guard wraps the whole emit.

- [ ] **Step 1: Add a sampling helper** near the top of `main.py` (after the `research` import):

```python
import random as _random

def _research_sample() -> bool:
    """True if this candidate event should be logged, per the sampling knob."""
    s = config.RESEARCH_CANDIDATE_SAMPLING
    return s >= 1.0 or _random.random() < s
```

- [ ] **Step 2: Instrument the rejection + accept points.** Rewrite the body of the `for cand in combo_enumerator.enumerate_2leg(...)` loop so each exit emits first. Replace lines ~1465–1482 with:

```python
            typed = [_leg_dict_to_typed(dict(l), game_id) for l in cand.legs]
            spread_line = _spread_line_from_legs([dict(l) for l in cand.legs])
            total_line = _total_line_from_legs([dict(l) for l in cand.legs])

            def _emit_cand(outcome, *, model=None, books=None,
                           blended=None, kalshi_ref=None, kelly=None):
                if not _research_sample():
                    return
                research.emit(
                    "candidate_evaluated",
                    game_id=game_id, combo_ticker=cand.legs[0]["market_ticker"],
                    leg_set_hash=cand.leg_set_hash,
                    spread_line=spread_line, total_line=total_line,
                    outcome=outcome, model_fair=model, book_fairs=books,
                    n_books=(len(books) if books else 0),
                    blended_fair=blended, kalshi_ref=kalshi_ref,
                    kelly_yes_n=(kelly[0] if kelly else None),
                    kelly_no_n=(kelly[1] if kelly else None),
                    worst_yes_ask=(kelly[2] if kelly else None),
                    worst_no_ask=(kelly[3] if kelly else None),
                )

            if any(l is None for l in typed):
                _emit_cand("rejected_no_mapping")
                continue
            model = fair_value.model_fair(samples, typed)
            books = _load_book_fairs(game_id, spread_line, total_line)
            blended = fair_value.blend(model, books)
            if blended is None:
                _emit_cand("rejected_no_book_data", model=model, books=books)
                continue
            if not (config.MIN_FAIR_PROB <= blended <= config.MAX_FAIR_PROB):
                _emit_cand("rejected_fair_oob", model=model, books=books,
                           blended=blended)
                continue
            kalshi_ref = _kalshi_last_price(cand.legs[0]["market_ticker"])
            yes_n, no_n, yes_ask, no_ask = _kelly_size_for_candidate(
                game_id, typed, samples, blended)
            _emit_cand("submitted", model=model, books=books, blended=blended,
                       kalshi_ref=kalshi_ref,
                       kelly=(yes_n, no_n, yes_ask, no_ask))
            candidates_all.append(cand)
            fair_scores[cand.leg_set_hash] = (blended, kalshi_ref)
            kelly_sizes[cand.leg_set_hash] = (yes_n, no_n, yes_ask, no_ask)
```

(Note: `outcome="submitted"` here means "passed enumeration and was ranked." The further `rejected_below_topN` distinction happens in `_refresh_rfqs`; capturing it is deferred — `live_rfqs` already records which candidates became RFQs, so the join `candidate_evaluated.outcome='submitted'` minus `live_rfqs` gives the below-topN set.)

- [ ] **Step 3: Verify import + tests.**

Run: `python -c "import kalshi_mlb_rfq.main"`
Expected: no error.
Run: `pytest kalshi_mlb_rfq/tests/ -v`
Expected: PASS.

- [ ] **Step 4: Manual smoke (optional, dry-run, NOT in worktree).** After merge to `main`, on the live host run one dry cycle and confirm `candidate_evaluated` rows appear:

```bash
python -c "import duckdb; print(duckdb.connect('kalshi_mlb_rfq/kalshi_mlb_rfq_research.duckdb', read_only=True).execute(\"SELECT outcome, count(*) FROM events WHERE event_type='candidate_evaluated' GROUP BY 1\").fetchall())"
```

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/main.py
git commit -m "feat(rfq-logging): emit candidate_evaluated (gaps A+B)"
```

---

# PHASE 3 — Remaining lifecycle events

⚠️ Churn-sensitive. These hook into `_all_per_accept_gates_pass`, `_kelly_size_for_candidate`, and `_evaluate_quote`. Re-read each before editing.

### Task 11: `gate_evaluated` — headroom on each gate

**Files:**
- Modify: `kalshi_mlb_rfq/main.py` (`_all_per_accept_gates_pass` ~529–604, `_evaluate_quote` call site ~725)

`_all_per_accept_gates_pass` returns `(bool, label)`. The simplest non-invasive capture: emit a `gate_evaluated` event at the point it returns, with the decision label and the few headroom numbers cheaply available (it already computes `today_fills`, the fill-ratio window, the cooldown map).

- [ ] **Step 1:** At the end of `_all_per_accept_gates_pass`, before each `return`, the function knows the failing gate via the label. Add a single emit just before the final `return True, "passed"` and convert early `return False, X` into a local helper. Replace the function's `return False, "<label>"` statements by assigning `decision` and `goto`-style fallthrough is not Pythonic — instead, wrap the body so every exit routes through one emit. Minimal approach: capture `combo_meta` + label at the *call site* in `_evaluate_quote` instead, where `passed, decision` are already returned (line ~725):

```python
        passed, decision = _all_per_accept_gates_pass(
            quote, fair, {"leg_set_hash": leg_set_hash, "game_id": game_id,
                          "legs_json": legs_json})
        research.emit("gate_evaluated", game_id=game_id,
                      combo_ticker=combo_market_ticker, rfq_id=rfq_id,
                      quote_id=quote.get("id"),
                      decision=decision, passed=passed,
                      blended_fair=fair)
        if not passed:
            _log_quote_decision(quote, fair, decision)
            return
```

(Per-gate *headroom* numbers — e.g. spend vs cap — are a deeper refactor of `_all_per_accept_gates_pass` to return a dict; deferred to a follow-up. The decision label already tells you which gate fired, which covers the primary "why" question.)

- [ ] **Step 2: Verify + commit**

Run: `python -c "import kalshi_mlb_rfq.main"` → no error; `pytest kalshi_mlb_rfq/tests/ -v` → PASS

```bash
git add kalshi_mlb_rfq/main.py
git commit -m "feat(rfq-logging): emit gate_evaluated at quote eval"
```

### Task 12: `kelly_sized`

**Files:**
- Modify: `kalshi_mlb_rfq/main.py` (`_kelly_size_for_candidate` ~744–791)

- [ ] **Step 1:** After computing `yes_n, no_n, yes_ask, no_ask` (line ~790, before `return`), add:

```python
    research.emit("kelly_sized", game_id=game_id,
                  blended_fair=blended_fair,
                  kelly_yes_n=yes_n, kelly_no_n=no_n,
                  worst_yes_ask=yes_ask, worst_no_ask=no_ask,
                  n_existing_positions=len(existing),
                  bankroll=config.BANKROLL,
                  kelly_fraction=config.KELLY_FRACTION)
    return yes_n, no_n, yes_ask, no_ask
```

- [ ] **Step 2: Verify + commit**

Run: `pytest kalshi_mlb_rfq/tests/ -v` → PASS

```bash
git add kalshi_mlb_rfq/main.py
git commit -m "feat(rfq-logging): emit kelly_sized at candidate sizing"
```

### Task 13: `position_snapshot` on fill

**Files:**
- Modify: `kalshi_mlb_rfq/main.py` (`_evaluate_quote`, the fill/positions write block)

- [ ] **Step 1:** In `_evaluate_quote`, immediately after the `INSERT INTO positions ... ON CONFLICT ...` upsert and before/after the `_log_quote_decision(..., "accepted", ...)` call, read back the post-fill position and emit:

```python
        with db.connect(read_only=True) as con:
            after = con.execute(
                "SELECT net_contracts, weighted_price FROM positions "
                "WHERE combo_market_ticker=? AND side=?",
                [combo_market_ticker, chosen_side]).fetchone()
        research.emit("position_snapshot", game_id=game_id,
                      combo_ticker=combo_market_ticker, rfq_id=rfq_id,
                      quote_id=quote.get("id"), side=chosen_side,
                      contracts_added=actual,
                      net_contracts_after=(after[0] if after else None),
                      weighted_price_after=(after[1] if after else None))
```

(`chosen_side` is the NO-side-aware variable already in scope at the accept block; if the local is named differently in current code, use that name.)

- [ ] **Step 2: Verify + commit**

Run: `pytest kalshi_mlb_rfq/tests/ -v` → PASS

```bash
git add kalshi_mlb_rfq/main.py
git commit -m "feat(rfq-logging): emit position_snapshot on fill"
```

### Task 14: `quote_priced`, `walk_diagnosed`, `rfq_submit_failed`

**Files:**
- Modify: `kalshi_mlb_rfq/main.py` (`_evaluate_quote` accept-fail path; `_refresh_rfqs` add-failure ~line where `add ... failed` is logged)

- [ ] **Step 1: `rfq_submit_failed`** — at the `add ... failed` warning in `_refresh_rfqs`:

```python
            except Exception as e:
                log.warning("add %s failed: %s", c.leg_set_hash[:8], e)
                research.emit("rfq_submit_failed",
                              game_id=getattr(c, "game_id", None),
                              combo_ticker=c.legs[0]["market_ticker"],
                              leg_set_hash=c.leg_set_hash, error=str(e))
```

- [ ] **Step 2: `walk_diagnosed`** — in `_evaluate_quote` where the accept response is `None`/non-2xx and `diag` already holds `accept_response_body` + `rfq_terminal_status`, add:

```python
            research.emit("walk_diagnosed", game_id=game_id,
                          combo_ticker=combo_market_ticker, rfq_id=rfq_id,
                          quote_id=quote.get("id"),
                          accept_response_body=diag.get("accept_response_body"),
                          rfq_terminal_status=diag.get("rfq_terminal_status"))
```

- [ ] **Step 3: `quote_priced`** — at the top of the `with ACCEPT_LOCK:` block right after `fair = _fresh_blended_fair(...)` succeeds, emit the components. NOTE: `_fresh_blended_fair` currently discards the per-book dict; to capture per-book at quote time, change it to also return `books` (small refactor) OR re-call `_load_book_fairs` here. Minimal version (re-call):

```python
        # quote_priced: re-derive components for the research log (cheap; cache hit)
        try:
            _qp_legs = json.loads(legs_json)
            _qp_sl = _spread_line_from_legs(_qp_legs)
            _qp_tl = _total_line_from_legs(_qp_legs)
            research.emit("quote_priced", game_id=game_id,
                          combo_ticker=combo_market_ticker, rfq_id=rfq_id,
                          quote_id=quote.get("id"), blended_fair=fair,
                          book_fairs=_load_book_fairs(game_id, _qp_sl, _qp_tl))
        except Exception:
            pass
```

- [ ] **Step 4: Verify + commit**

Run: `pytest kalshi_mlb_rfq/tests/ -v` → PASS

```bash
git add kalshi_mlb_rfq/main.py
git commit -m "feat(rfq-logging): emit quote_priced, walk_diagnosed, rfq_submit_failed"
```

---

# PHASE 4 — Analysis surface + docs

### Task 15: Canned analysis queries

**Files:**
- Create: `kalshi_mlb_rfq/research_queries.sql`

- [ ] **Step 1: Create `kalshi_mlb_rfq/research_queries.sql`** with the three example queries:

```sql
-- Research analysis queries. Run from the kalshi_mlb_rfq/ dir:
--   duckdb -init research_queries.sql
ATTACH 'kalshi_mlb_rfq.duckdb'          AS state    (READ_ONLY);
ATTACH 'kalshi_mlb_rfq_research.duckdb' AS research (READ_ONLY);

-- 1) RFQ JOURNEY: every event for one combo, in time order.
--    Replace :combo with a combo_ticker.
-- SELECT ts, event_type, payload FROM research.events
-- WHERE combo_ticker = :combo ORDER BY ts;

-- 2) MISSED EDGES: candidates we priced but did not submit, with edge.
SELECT ts, game_id, combo_ticker,
       payload->>'blended_fair'  AS blended_fair,
       payload->>'kalshi_ref'    AS kalshi_ref,
       payload->>'outcome'       AS outcome
FROM research.events
WHERE event_type = 'candidate_evaluated'
  AND payload->>'outcome' <> 'submitted'
ORDER BY ts DESC
LIMIT 200;

-- 3) GATE BREAKDOWN: which gate rejected accepts, last 24h.
SELECT payload->>'decision' AS decision, count(*) AS n
FROM research.events
WHERE event_type = 'gate_evaluated'
  AND ts > now() - INTERVAL 1 DAY
GROUP BY 1 ORDER BY n DESC;
```

- [ ] **Step 2: Verify queries parse** (after some events exist):

Run: `cd kalshi_mlb_rfq && duckdb -init research_queries.sql -c "SELECT 1"`
Expected: no SQL error (queries are valid; result sets may be empty pre-data).

- [ ] **Step 3: Commit**

```bash
git add kalshi_mlb_rfq/research_queries.sql
git commit -m "feat(rfq-logging): canned research analysis queries"
```

### Task 16: Docs + gitignore

**Files:**
- Modify: `kalshi_mlb_rfq/README.md`, repo `CLAUDE.md`, `.gitignore`

- [ ] **Step 1: gitignore.** Confirm `*.duckdb` is ignored (it covers the new research DB). Add rotated-log coverage:

Run: `grep -nE "bot\.log|\*\.duckdb" .gitignore`
If `bot.log*` is not covered, add to `.gitignore`:

```
kalshi_mlb_rfq/bot.log
kalshi_mlb_rfq/bot.log.*
```

- [ ] **Step 2: README.** Add an "Observability / research logging" section to `kalshi_mlb_rfq/README.md` documenting: the three sibling DBs; the `events` table + event catalog; the `ATTACH` query pattern (point to `research_queries.sql`); config knobs (`RESEARCH_DB_PATH`, `RESEARCH_CANDIDATE_SAMPLING`, `RESEARCH_BUFFER_MAX`, `LOG_*`); the unscheduled `prune_research()`; and the no-auto-prune retention note.

- [ ] **Step 3: CLAUDE.md.** Extend the `kalshi_mlb_rfq` bullet in the repo `CLAUDE.md` to mention the third sibling research DB and structured rotating logging.

- [ ] **Step 4: Commit**

```bash
git add kalshi_mlb_rfq/README.md CLAUDE.md .gitignore
git commit -m "docs(rfq-logging): observability section + gitignore log rotation"
```

---

## Pre-merge (per repo rules)

- Executive-engineer review of `git diff main..HEAD` (data integrity, resource safety, dead code, log hygiene, no secrets).
- Run full bot test suite: `pytest kalshi_mlb_rfq/tests/ -v`.
- The live bot is restarted **from the main repo cwd** after merge (never run from the worktree — see restart-gotchas memory).
- Explicit user approval before merging to `main`.
- After merge: add a memory note on the research-logging architecture.

---

## Self-Review

**Spec coverage:**
- Section 1 (research store) → Tasks 4–9. ✅
- Section 2 (event catalog) → Tasks 10–14. ✅ (`rejected_below_topN` partially deferred — covered via join, noted in Task 10; per-gate headroom deferred — noted in Task 11.)
- Section 3 (operational logging) → Tasks 1–3. ✅
- Section 4 (analysis) → Task 15; retention `prune_research` unscheduled → Task 8. ✅
- Section 5 (rollout/testing/VC/docs) → phase structure, merge checkpoints, Task 16. ✅
- Safety invariant (flush never breaks trading) → Task 7 test `test_flush_swallows_db_error_and_keeps_running`. ✅

**Placeholder scan:** No TBD/TODO in executable steps. Two explicit scope deferrals (below-topN distinction, per-gate headroom) are documented with the reason + the workaround, not left vague.

**Type consistency:** `emit(event_type, *, game_id, combo_ticker, rfq_id, quote_id, **payload)` signature used identically in Tasks 6, 10–14. `_BUFFER`/`_SESSION_ID` names match across tests and impl. `flush()`/`init_research_db()`/`prune_research()` names consistent. `config.RESEARCH_*` names match Task 4 definitions.

**Known fragility:** Phases 2–3 line numbers are vs `main` @ `5bd5e6c`; the churn warning instructs re-reading each target function before editing. Variable names at hot-path call sites (`chosen_side`, `diag`) must be confirmed against current code at edit time.
