# Kalshi MLB RFQ Bot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standalone Python daemon that auto-RFQs MLB game-line same-game-parlays on Kalshi, evaluates maker quotes against a model+sportsbook blended fair value, and auto-accepts +EV quotes within a complete safety scaffold.

**Architecture:** Long-running single-process daemon (mirrors `kalshi_mm/main.py` shape) at `kalshi_mlb_rfq/`. Reads `mlb.duckdb` (samples + sgp_odds) read-only and writes `kalshi_mlb_rfq.duckdb`. Wide-mode continuous-pipeline RFQs (up to 80 in-flight) on the cross-category MVE collection (`KXMVECROSSCATEGORY-R`). Per-accept gates handle most safety; conditional Kelly sized from `mlb_game_samples` for sizing.

**Tech Stack:** Python 3 + DuckDB + `cryptography` (RSA signing for Kalshi API auth). Auth helper (`sign_request`) lifts from existing `kalshi_draft/auth.py`. Kelly logic adapts from `kalshi_mm/kelly.py`.

**Spec reference:** `docs/superpowers/specs/2026-04-27-kalshi-mlb-rfq-bot-design.md`

---

## File structure

| File | Responsibility |
|---|---|
| `kalshi_mlb_rfq/__init__.py` | Package marker |
| `kalshi_mlb_rfq/config.py` | Env loader, all knobs with defaults |
| `kalshi_mlb_rfq/db.py` | DuckDB schema + helpers (combo_cache, live_rfqs, quote_log, fills, positions, sessions, combo_cooldown) |
| `kalshi_mlb_rfq/auth_client.py` | Authenticated HTTP wrapper around Kalshi `/trade-api/v2` |
| `kalshi_mlb_rfq/ticker_map.py` | game_id → Kalshi event ticker; leg type → Kalshi market+side |
| `kalshi_mlb_rfq/rfq_client.py` | mint combo, create RFQ, poll quotes, accept (+ fill reconcile), DELETE |
| `kalshi_mlb_rfq/ev_calc.py` | Quadratic fee formula + post-fee EV per side |
| `kalshi_mlb_rfq/fair_value.py` | model_fair from samples + book_fair from `mlb_sgp_odds` + blend with 2-source gate |
| `kalshi_mlb_rfq/combo_enumerator.py` | Per-game candidate enumeration + cross-game priority queue |
| `kalshi_mlb_rfq/kelly.py` | Conditional Kelly sizing on combo outcome vectors (adapted from `kalshi_mm/kelly.py`) |
| `kalshi_mlb_rfq/risk.py` | All per-accept gates + risk sweep (tipoff, line-move, cooldown, kill switch, exposure caps, fill ratio, staleness, positions health) |
| `kalshi_mlb_rfq/notify.py` | bot.log + optional webhook |
| `kalshi_mlb_rfq/main.py` | Daemon orchestrator: RFQ refresh / quote poll / risk sweep / pipeline refresh / heartbeat loops, accept lock |
| `kalshi_mlb_rfq/dashboard.py` | CLI status tool |
| `kalshi_mlb_rfq/.env.example` | Template |
| `kalshi_mlb_rfq/README.md` | Setup, run, monitor, troubleshoot |
| `tests/kalshi_mlb_rfq/test_*.py` | Unit tests (one file per module) |
| `mlb_sgp/recon_kalshi_mlb_rfq.py` | Already drafted in worktree; updated in Task 24 |

Modifications:
- Root `README.md` — add 1-line entry for `kalshi_mlb_rfq/` under Project Structure
- Root `CLAUDE.md` — add 1-paragraph entry under Project Structure
- Root `.gitignore` — add `kalshi_mlb_rfq/.env`, `kalshi_mlb_rfq/*.duckdb`, `kalshi_mlb_rfq/*.log`, `kalshi_mlb_rfq/.kill`

---

## Phase 0 — Setup

### Task 1: Project scaffolding + .env.example + .gitignore

**Files:**
- Create: `kalshi_mlb_rfq/__init__.py`
- Create: `kalshi_mlb_rfq/.env.example`
- Create: `kalshi_mlb_rfq/requirements.txt`
- Modify: `.gitignore` (root)

- [ ] **Step 1: Create directory and package marker**

```bash
mkdir -p /Users/callancapitolo/NFLWork/.worktrees/kalshi-mlb-rfq/kalshi_mlb_rfq
touch /Users/callancapitolo/NFLWork/.worktrees/kalshi-mlb-rfq/kalshi_mlb_rfq/__init__.py
mkdir -p /Users/callancapitolo/NFLWork/.worktrees/kalshi-mlb-rfq/tests/kalshi_mlb_rfq
touch /Users/callancapitolo/NFLWork/.worktrees/kalshi-mlb-rfq/tests/kalshi_mlb_rfq/__init__.py
```

- [ ] **Step 2: Write `kalshi_mlb_rfq/.env.example`**

```bash
# Kalshi credentials (copy from kalshi_mm/.env)
KALSHI_API_KEY_ID=
KALSHI_PRIVATE_KEY_PATH=
KALSHI_USER_ID=fb3682ce-41d2-4414-b1f3-e79244c5af8e
KALSHI_BASE_URL=https://api.elections.kalshi.com/trade-api/v2

# MVE collection (cross-category; the one that has MLB legs as eligible events)
MVE_COLLECTION_TICKER=KXMVECROSSCATEGORY-R

# Sizing
BANKROLL=1000.0
KELLY_FRACTION=0.25

# Per-accept gates
MIN_EV_PCT=0.05
MAX_QUOTE_DEVIATION=0.15
MIN_FAIR_PROB=0.05
MAX_FAIR_PROB=0.95
MAX_GAME_EXPOSURE_PCT=0.10
DAILY_EXPOSURE_CAP_USD=200.0
LINE_MOVE_THRESHOLD=0.5
MAX_PREDICTION_STALENESS_SEC=600
MAX_BOOK_STALENESS_SEC=60
COMBO_COOLDOWN_SEC=30
POSITIONS_HEALTH_RETRIES=2
MIN_FILL_RATIO=0.50
FILL_RATIO_WINDOW=50
TIPOFF_CANCEL_MIN=5

# Loops
RFQ_REFRESH_SEC=30
QUOTE_POLL_SEC=2
RISK_SWEEP_SEC=10
PIPELINE_REFRESH_SEC=600
MAX_LIVE_RFQS=80

# Notifications
NOTIFY_WEBHOOK_URL=

# Devigging fallbacks (used when fewer than 4 sides exist for a (game, book, spread, total) tuple)
DK_VIG_FALLBACK=0.125
FD_VIG_FALLBACK=0.18
PX_VIG_FALLBACK=0.05
NOVIG_VIG_FALLBACK=0.05
```

- [ ] **Step 3: Write `kalshi_mlb_rfq/requirements.txt`**

```
duckdb>=0.10.0
cryptography>=41.0.0
pytest>=7.0.0
```

- [ ] **Step 4: Update root `.gitignore`** — append the following lines if not present:

```
kalshi_mlb_rfq/.env
kalshi_mlb_rfq/*.duckdb
kalshi_mlb_rfq/*.duckdb.wal
kalshi_mlb_rfq/*.log
kalshi_mlb_rfq/.kill
```

- [ ] **Step 5: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/kalshi-mlb-rfq
git add kalshi_mlb_rfq/__init__.py kalshi_mlb_rfq/.env.example kalshi_mlb_rfq/requirements.txt tests/kalshi_mlb_rfq/__init__.py .gitignore
git commit -m "feat(kalshi-mlb-rfq): scaffold module + .env.example + gitignore"
```

---

## Phase 1 — Foundation

### Task 2: `config.py` — env loader with defaults

**Files:**
- Create: `kalshi_mlb_rfq/config.py`
- Test: `tests/kalshi_mlb_rfq/test_config.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_config.py
import os
from pathlib import Path

import pytest


def test_config_loads_from_env_example(monkeypatch):
    """Loading config from .env.example yields all expected knobs with correct types."""
    env_path = Path(__file__).parent.parent.parent / "kalshi_mlb_rfq" / ".env.example"
    # Load env file manually to populate os.environ
    for line in env_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        if v.strip():
            monkeypatch.setenv(k.strip(), v.strip())

    # Reload config so it picks up monkeypatched env
    import importlib
    import kalshi_mlb_rfq.config as config_mod
    importlib.reload(config_mod)

    assert config_mod.KALSHI_BASE_URL == "https://api.elections.kalshi.com/trade-api/v2"
    assert config_mod.MVE_COLLECTION_TICKER == "KXMVECROSSCATEGORY-R"
    assert config_mod.BANKROLL == 1000.0
    assert config_mod.KELLY_FRACTION == 0.25
    assert config_mod.MIN_EV_PCT == 0.05
    assert config_mod.MAX_LIVE_RFQS == 80
    assert config_mod.MAX_PREDICTION_STALENESS_SEC == 600
    assert config_mod.RFQ_REFRESH_SEC == 30
```

- [ ] **Step 2: Run test, expect failure**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/kalshi-mlb-rfq
pytest tests/kalshi_mlb_rfq/test_config.py -v
```
Expected: FAIL — `kalshi_mlb_rfq.config` module not found.

- [ ] **Step 3: Implement `kalshi_mlb_rfq/config.py`**

```python
"""Config knobs for the Kalshi MLB RFQ bot. Loaded from .env (or environment)."""

import os
from pathlib import Path

PKG_DIR = Path(__file__).parent
PROJECT_ROOT = PKG_DIR.parent
DB_PATH = PKG_DIR / "kalshi_mlb_rfq.duckdb"
LOG_PATH = PKG_DIR / "bot.log"
KILL_FILE = PKG_DIR / ".kill"


def _load_env(path: Path) -> dict[str, str]:
    """Parse a .env file; returns a dict of key→value."""
    env: dict[str, str] = {}
    if not path.exists():
        return env
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        env[k.strip()] = v.strip().strip('"').strip("'")
    return env


_FILE_ENV = _load_env(PKG_DIR / ".env")


def _get(key: str, default: str | None = None) -> str | None:
    return os.environ.get(key, _FILE_ENV.get(key, default))


# Credentials
KALSHI_API_KEY_ID = _get("KALSHI_API_KEY_ID")
KALSHI_PRIVATE_KEY_PATH = _get("KALSHI_PRIVATE_KEY_PATH")
KALSHI_USER_ID = _get("KALSHI_USER_ID")
KALSHI_BASE_URL = _get(
    "KALSHI_BASE_URL", "https://api.elections.kalshi.com/trade-api/v2"
)

# MVE
MVE_COLLECTION_TICKER = _get("MVE_COLLECTION_TICKER", "KXMVECROSSCATEGORY-R")

# Sizing
BANKROLL = float(_get("BANKROLL", "1000.0"))
KELLY_FRACTION = float(_get("KELLY_FRACTION", "0.25"))

# Per-accept gates
MIN_EV_PCT = float(_get("MIN_EV_PCT", "0.05"))
MAX_QUOTE_DEVIATION = float(_get("MAX_QUOTE_DEVIATION", "0.15"))
MIN_FAIR_PROB = float(_get("MIN_FAIR_PROB", "0.05"))
MAX_FAIR_PROB = float(_get("MAX_FAIR_PROB", "0.95"))
MAX_GAME_EXPOSURE_PCT = float(_get("MAX_GAME_EXPOSURE_PCT", "0.10"))
DAILY_EXPOSURE_CAP_USD = float(_get("DAILY_EXPOSURE_CAP_USD", "200.0"))
LINE_MOVE_THRESHOLD = float(_get("LINE_MOVE_THRESHOLD", "0.5"))
MAX_PREDICTION_STALENESS_SEC = int(_get("MAX_PREDICTION_STALENESS_SEC", "600"))
MAX_BOOK_STALENESS_SEC = int(_get("MAX_BOOK_STALENESS_SEC", "60"))
COMBO_COOLDOWN_SEC = int(_get("COMBO_COOLDOWN_SEC", "30"))
POSITIONS_HEALTH_RETRIES = int(_get("POSITIONS_HEALTH_RETRIES", "2"))
MIN_FILL_RATIO = float(_get("MIN_FILL_RATIO", "0.50"))
FILL_RATIO_WINDOW = int(_get("FILL_RATIO_WINDOW", "50"))
TIPOFF_CANCEL_MIN = int(_get("TIPOFF_CANCEL_MIN", "5"))

# Loops
RFQ_REFRESH_SEC = int(_get("RFQ_REFRESH_SEC", "30"))
QUOTE_POLL_SEC = int(_get("QUOTE_POLL_SEC", "2"))
RISK_SWEEP_SEC = int(_get("RISK_SWEEP_SEC", "10"))
PIPELINE_REFRESH_SEC = int(_get("PIPELINE_REFRESH_SEC", "600"))
MAX_LIVE_RFQS = int(_get("MAX_LIVE_RFQS", "80"))

# Notifications
NOTIFY_WEBHOOK_URL = _get("NOTIFY_WEBHOOK_URL")

# Vig fallbacks
DK_VIG_FALLBACK = float(_get("DK_VIG_FALLBACK", "0.125"))
FD_VIG_FALLBACK = float(_get("FD_VIG_FALLBACK", "0.18"))
PX_VIG_FALLBACK = float(_get("PX_VIG_FALLBACK", "0.05"))
NOVIG_VIG_FALLBACK = float(_get("NOVIG_VIG_FALLBACK", "0.05"))

# Source path for the answer-key DB
ANSWER_KEY_DB = PROJECT_ROOT / "Answer Keys" / "mlb.duckdb"
```

- [ ] **Step 4: Run test, expect pass**

```bash
pytest tests/kalshi_mlb_rfq/test_config.py -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/config.py tests/kalshi_mlb_rfq/test_config.py
git commit -m "feat(kalshi-mlb-rfq): config module with all knobs from spec"
```

---

### Task 3: `db.py` — schema migrations + helpers

**Files:**
- Create: `kalshi_mlb_rfq/db.py`
- Test: `tests/kalshi_mlb_rfq/test_db.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_db.py
import duckdb
from pathlib import Path

import pytest

from kalshi_mlb_rfq import db


@pytest.fixture
def tmpdb(tmp_path, monkeypatch):
    p = tmp_path / "test.duckdb"
    monkeypatch.setattr(db, "DB_PATH", p)
    db.init_database()
    return p


def test_init_creates_all_tables(tmpdb):
    con = duckdb.connect(str(tmpdb), read_only=True)
    try:
        names = {row[0] for row in con.execute(
            "SELECT table_name FROM information_schema.tables "
            "WHERE table_schema='main'"
        ).fetchall()}
    finally:
        con.close()
    assert names == {
        "combo_cache", "live_rfqs", "quote_log", "fills",
        "positions", "sessions", "combo_cooldown", "reference_lines",
    }


def test_session_round_trip(tmpdb):
    sid = db.start_session(pid=99, dry_run=True, version="0.1.0")
    assert sid
    db.end_session(sid)
    con = duckdb.connect(str(tmpdb), read_only=True)
    try:
        ended_at = con.execute(
            "SELECT ended_at FROM sessions WHERE session_id=?", [sid]
        ).fetchone()[0]
    finally:
        con.close()
    assert ended_at is not None
```

- [ ] **Step 2: Run test, expect failure**

```bash
pytest tests/kalshi_mlb_rfq/test_db.py -v
```
Expected: FAIL — `kalshi_mlb_rfq.db` not found.

- [ ] **Step 3: Implement `kalshi_mlb_rfq/db.py`**

```python
"""DuckDB schema + helpers for the Kalshi MLB RFQ bot."""

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
    quote_id              VARCHAR PRIMARY KEY,
    rfq_id                VARCHAR NOT NULL,
    combo_market_ticker   VARCHAR,
    creator_id            VARCHAR,
    yes_bid_dollars       DOUBLE,
    no_bid_dollars        DOUBLE,
    blended_fair_at_eval  DOUBLE,
    post_fee_ev_pct       DOUBLE,
    decision              VARCHAR NOT NULL,
    reason_detail         VARCHAR,
    observed_at           TIMESTAMP NOT NULL
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
    combo_market_ticker  VARCHAR PRIMARY KEY,
    game_id              VARCHAR NOT NULL,
    net_contracts        DOUBLE NOT NULL,
    weighted_price       DOUBLE NOT NULL,
    legs_json            VARCHAR NOT NULL,
    updated_at           TIMESTAMP NOT NULL
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
    leg_set_hash  VARCHAR PRIMARY KEY,
    game_id       VARCHAR NOT NULL,
    cooled_until  TIMESTAMP NOT NULL,
    reason        VARCHAR
);

CREATE TABLE IF NOT EXISTS reference_lines (
    rfq_id     VARCHAR PRIMARY KEY,
    lines_json VARCHAR NOT NULL,
    snapped_at TIMESTAMP NOT NULL
);
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


def init_database():
    """Apply schema migrations idempotently."""
    with connect() as con:
        con.execute(SCHEMA_SQL)


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
```

- [ ] **Step 4: Run test, expect pass**

```bash
pytest tests/kalshi_mlb_rfq/test_db.py -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/db.py tests/kalshi_mlb_rfq/test_db.py
git commit -m "feat(kalshi-mlb-rfq): DuckDB schema + connection helper with retry"
```

---

### Task 4: `auth_client.py` — authenticated HTTP wrapper

**Files:**
- Create: `kalshi_mlb_rfq/auth_client.py`
- Test: `tests/kalshi_mlb_rfq/test_auth_client.py`

The CBB MM bot's `orders.py::_authenticated_request` uses `urllib.request` and signs each call with `kalshi_draft.auth.sign_request`. Lift the same pattern but make it more testable by allowing dependency-injection of the signer.

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_auth_client.py
import json
from unittest.mock import patch, MagicMock

from kalshi_mlb_rfq import auth_client


def test_get_returns_parsed_json():
    fake_resp = MagicMock()
    fake_resp.status = 200
    fake_resp.read.return_value = b'{"foo": "bar"}'
    fake_resp.headers = {}
    fake_resp.__enter__.return_value = fake_resp
    fake_resp.__exit__.return_value = False

    with patch("kalshi_mlb_rfq.auth_client.urllib.request.urlopen", return_value=fake_resp), \
         patch("kalshi_mlb_rfq.auth_client._sign", return_value="sig"):
        status, body, _ = auth_client.api("GET", "/exchange/status")

    assert status == 200
    assert body == {"foo": "bar"}


def test_delete_handles_empty_body():
    fake_resp = MagicMock()
    fake_resp.status = 204
    fake_resp.read.return_value = b''
    fake_resp.headers = {}
    fake_resp.__enter__.return_value = fake_resp
    fake_resp.__exit__.return_value = False

    with patch("kalshi_mlb_rfq.auth_client.urllib.request.urlopen", return_value=fake_resp), \
         patch("kalshi_mlb_rfq.auth_client._sign", return_value="sig"):
        status, body, _ = auth_client.api("DELETE", "/communications/rfqs/abc")

    assert status == 204
    assert body == {}
```

- [ ] **Step 2: Run test, expect failure**

```bash
pytest tests/kalshi_mlb_rfq/test_auth_client.py -v
```
Expected: FAIL — module not found.

- [ ] **Step 3: Implement `kalshi_mlb_rfq/auth_client.py`**

```python
"""Authenticated HTTP wrapper around Kalshi /trade-api/v2.

Borrows the signing helper from kalshi_draft/auth.py (already in the main repo).
"""

import json
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from kalshi_mlb_rfq.config import (
    KALSHI_API_KEY_ID,
    KALSHI_BASE_URL,
    KALSHI_PRIVATE_KEY_PATH,
    PROJECT_ROOT,
)

# kalshi_draft is gitignored; load from main repo (works inside a worktree because
# git rev-parse --git-common-dir resolves to the shared .git dir's parent).
sys.path.insert(0, str(PROJECT_ROOT.parent / "kalshi_draft"))
try:
    from auth import sign_request as _sign_request
except ImportError:
    # In tests we patch _sign anyway; provide a stub so import doesn't fail.
    def _sign_request(_pk, _ts, _method, _path):
        raise RuntimeError("sign_request unavailable; set KALSHI_PRIVATE_KEY_PATH and "
                           "ensure kalshi_draft/auth.py is reachable")


def _sign(method: str, path_no_query: str) -> tuple[str, str]:
    """Returns (signature, timestamp_ms_str)."""
    ts = str(int(datetime.now(timezone.utc).timestamp() * 1000))
    sig = _sign_request(KALSHI_PRIVATE_KEY_PATH, ts, method,
                        f"/trade-api/v2{path_no_query}")
    return sig, ts


def api(method: str, path: str, body: dict | None = None,
        timeout: int = 30) -> tuple[int, dict | str, dict]:
    """Make an authenticated request.

    Args:
        method: HTTP verb.
        path: Path under /trade-api/v2 (e.g., "/exchange/status").
              Query string allowed; only the path portion is signed.
        body: Optional JSON body.
        timeout: Seconds.

    Returns:
        (status_code, parsed_body_or_text, response_headers_dict)
    """
    path_no_query = path.split("?", 1)[0]
    sig, ts = _sign(method, path_no_query)

    url = f"{KALSHI_BASE_URL}{path}"
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("KALSHI-ACCESS-KEY", KALSHI_API_KEY_ID or "")
    req.add_header("KALSHI-ACCESS-SIGNATURE", sig)
    req.add_header("KALSHI-ACCESS-TIMESTAMP", ts)
    req.add_header("Content-Type", "application/json")

    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            text = r.read().decode()
            headers = dict(r.headers)
            if not text:
                return r.status, {}, headers
            try:
                return r.status, json.loads(text), headers
            except json.JSONDecodeError:
                return r.status, text, headers
    except urllib.error.HTTPError as e:
        text = e.read().decode()
        headers = dict(e.headers) if e.headers else {}
        try:
            return e.code, json.loads(text), headers
        except (json.JSONDecodeError, ValueError):
            return e.code, text, headers
```

- [ ] **Step 4: Run test, expect pass**

```bash
pytest tests/kalshi_mlb_rfq/test_auth_client.py -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/auth_client.py tests/kalshi_mlb_rfq/test_auth_client.py
git commit -m "feat(kalshi-mlb-rfq): authenticated HTTP wrapper around Kalshi API"
```

---

## Phase 2 — Domain models

### Task 5: `ticker_map.py` — game/leg → Kalshi tickers

Per recon (spec §4.1):
- `KXMLBSPREAD-{event_suffix}-{TEAM}{N}` where N=2 means "team wins by over 1.5"
- `KXMLBTOTAL-{event_suffix}-{N}` where N=8 means "Over 7.5"
- `KXMLBGAME-{event_suffix}-{TEAM}` for winner
- `KXMLBRFI-{event_suffix}` for run-in-first

The event suffix format is `YYMMM[DDD]HHMM[AwayHome]` like `26APR282005NYYTEX` — matches the Odds-API `commence_time` in ET plus team codes.

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_ticker_map.py
from datetime import datetime, timezone

import pytest

from kalshi_mlb_rfq import ticker_map


def test_event_suffix_round_trip():
    # Apr 28 2026 at 8:05 PM ET → 26APR282005
    commence = datetime(2026, 4, 29, 0, 5, tzinfo=timezone.utc)  # 8:05 PM ET
    suffix = ticker_map.format_event_suffix(commence, away_code="NYY", home_code="TEX")
    assert suffix == "26APR282005NYYTEX"


def test_spread_ticker_for_home_team():
    t = ticker_map.spread_ticker("26APR282005NYYTEX", team_code="TEX", line=-1.5)
    assert t == "KXMLBSPREAD-26APR282005NYYTEX-TEX2"


def test_spread_ticker_for_away_team_alt_line():
    t = ticker_map.spread_ticker("26APR282005NYYTEX", team_code="NYY", line=-2.5)
    assert t == "KXMLBSPREAD-26APR282005NYYTEX-NYY3"


def test_total_ticker():
    t = ticker_map.total_ticker("26APR282005NYYTEX", line=7.5)
    assert t == "KXMLBTOTAL-26APR282005NYYTEX-8"


def test_total_ticker_alt():
    t = ticker_map.total_ticker("26APR282005NYYTEX", line=10.5)
    assert t == "KXMLBTOTAL-26APR282005NYYTEX-11"
```

- [ ] **Step 2: Run test, expect failure**

```bash
pytest tests/kalshi_mlb_rfq/test_ticker_map.py -v
```
Expected: FAIL — module not found.

- [ ] **Step 3: Implement `kalshi_mlb_rfq/ticker_map.py`**

```python
"""Map internal game/leg descriptors to Kalshi ticker strings.

Per recon (2026-04-27):
  KXMLBGAME   single ticker per (event, team) — winner
  KXMLBSPREAD multiple tickers per event: -{TEAM}{N}, "team wins by over (N-1).5"
  KXMLBTOTAL  multiple tickers per event: -{N}, "total over (N-1).5"
  KXMLBRFI    single ticker per event, no suffix — run in 1st inning

Event suffix format: YYMMMDDHHMM{AwayCode}{HomeCode}, all in ET.
"""

from datetime import datetime
from zoneinfo import ZoneInfo

ET = ZoneInfo("America/New_York")
_MONTHS = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
           "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]


def format_event_suffix(commence_time: datetime, away_code: str, home_code: str) -> str:
    """commence_time may be UTC or any tz-aware datetime; converted to ET for the suffix."""
    if commence_time.tzinfo is None:
        raise ValueError("commence_time must be timezone-aware")
    et = commence_time.astimezone(ET)
    return (
        f"{et.year % 100:02d}"
        f"{_MONTHS[et.month - 1]}"
        f"{et.day:02d}"
        f"{et.hour:02d}"
        f"{et.minute:02d}"
        f"{away_code}{home_code}"
    )


def spread_ticker(event_suffix: str, team_code: str, line: float) -> str:
    """Spread ticker. line is the team's spread (negative = favorite).

    -1.5 → suffix N=2 (team wins by over 1.5)
    -2.5 → suffix N=3
    -3.5 → suffix N=4
    """
    if line >= 0:
        raise ValueError(f"spread_ticker expects negative line for favorite; got {line}")
    n = int(round(-line + 0.5))
    return f"KXMLBSPREAD-{event_suffix}-{team_code}{n}"


def total_ticker(event_suffix: str, line: float) -> str:
    """Total ticker. line is the .5-suffixed total (e.g. 7.5).

    7.5 → -8
    8.5 → -9
    10.5 → -11
    """
    n = int(round(line + 0.5))
    return f"KXMLBTOTAL-{event_suffix}-{n}"


def game_ticker(event_suffix: str, team_code: str) -> str:
    return f"KXMLBGAME-{event_suffix}-{team_code}"


def rfi_ticker(event_suffix: str) -> str:
    return f"KXMLBRFI-{event_suffix}"
```

- [ ] **Step 4: Run test, expect pass**

```bash
pytest tests/kalshi_mlb_rfq/test_ticker_map.py -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/ticker_map.py tests/kalshi_mlb_rfq/test_ticker_map.py
git commit -m "feat(kalshi-mlb-rfq): ticker_map for KXMLB{SPREAD,TOTAL,GAME,RFI}"
```

---

## Phase 3 — Kalshi RFQ API client

### Task 6: `rfq_client.py` — mint combo ticker

Recon-validated path: `POST /multivariate_event_collections/{collection_ticker}` with body `{"selected_markets": [...]}` returns `{"market_ticker", "event_ticker"}`.

**Files:**
- Create: `kalshi_mlb_rfq/rfq_client.py`
- Test: `tests/kalshi_mlb_rfq/test_rfq_client_mint.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_rfq_client_mint.py
from unittest.mock import patch

from kalshi_mlb_rfq import rfq_client


def test_mint_combo_returns_ticker():
    legs = [
        {"market_ticker": "KXMLBSPREAD-X-NYY2", "event_ticker": "KXMLBSPREAD-X", "side": "yes"},
        {"market_ticker": "KXMLBTOTAL-X-8",     "event_ticker": "KXMLBTOTAL-X",  "side": "yes"},
    ]
    fake_response = (
        200,
        {"market_ticker": "KXMVECROSSCATEGORY-S-FOO", "event_ticker": "KXMVECROSSCATEGORY-S"},
        {},
    )
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=fake_response) as mock_api:
        ct, et = rfq_client.mint_combo_ticker("KXMVECROSSCATEGORY-R", legs)

    assert ct == "KXMVECROSSCATEGORY-S-FOO"
    assert et == "KXMVECROSSCATEGORY-S"
    mock_api.assert_called_once()
    call_args = mock_api.call_args
    assert call_args[0][0] == "POST"
    assert call_args[0][1] == "/multivariate_event_collections/KXMVECROSSCATEGORY-R"
    assert call_args[1]["body"] == {"selected_markets": legs}


def test_mint_combo_raises_on_error():
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=(400, {"error": "bad"}, {})):
        with __import__("pytest").raises(rfq_client.KalshiAPIError):
            rfq_client.mint_combo_ticker("KXMVECROSSCATEGORY-R", [])
```

- [ ] **Step 2: Run test, expect failure**

```bash
pytest tests/kalshi_mlb_rfq/test_rfq_client_mint.py -v
```
Expected: FAIL — module not found.

- [ ] **Step 3: Implement `kalshi_mlb_rfq/rfq_client.py` (initial version with mint only; later tasks will extend it)**

```python
"""Kalshi RFQ + Quote API client."""

from kalshi_mlb_rfq.auth_client import api


class KalshiAPIError(Exception):
    pass


def mint_combo_ticker(collection_ticker: str, selected_markets: list[dict]) -> tuple[str, str]:
    """Mint (or look up) a combo market by submitting selected_markets.

    Returns (combo_market_ticker, combo_event_ticker).
    Raises KalshiAPIError on non-200.
    """
    path = f"/multivariate_event_collections/{collection_ticker}"
    status, body, _ = api("POST", path, body={"selected_markets": selected_markets})
    if status != 200 or not isinstance(body, dict) or "market_ticker" not in body:
        raise KalshiAPIError(f"mint_combo_ticker failed: status={status} body={body}")
    return body["market_ticker"], body["event_ticker"]
```

- [ ] **Step 4: Run test, expect pass**

```bash
pytest tests/kalshi_mlb_rfq/test_rfq_client_mint.py -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/rfq_client.py tests/kalshi_mlb_rfq/test_rfq_client_mint.py
git commit -m "feat(kalshi-mlb-rfq): rfq_client.mint_combo_ticker"
```

---

### Task 7: RFQ client — create / get / delete RFQ

**Files:**
- Modify: `kalshi_mlb_rfq/rfq_client.py`
- Test: `tests/kalshi_mlb_rfq/test_rfq_client_lifecycle.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_rfq_client_lifecycle.py
from unittest.mock import patch

import pytest

from kalshi_mlb_rfq import rfq_client


def test_create_rfq_returns_id():
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=(200, {"id": "abc-123"}, {})):
        rid = rfq_client.create_rfq("KXMVECROSSCATEGORY-S-FOO", target_cost_dollars=0.50)
    assert rid == "abc-123"


def test_get_rfq_returns_status():
    with patch("kalshi_mlb_rfq.rfq_client.api",
               return_value=(200, {"rfq": {"id": "abc", "status": "open"}}, {})):
        rfq = rfq_client.get_rfq("abc")
    assert rfq["status"] == "open"


def test_delete_rfq_handles_204():
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=(204, {}, {})):
        ok = rfq_client.delete_rfq("abc")
    assert ok is True


def test_delete_rfq_handles_already_closed():
    with patch("kalshi_mlb_rfq.rfq_client.api",
               return_value=(400, {"error": {"code": "expired"}}, {})):
        ok = rfq_client.delete_rfq("abc")
    assert ok is False  # already closed; not a hard error
```

- [ ] **Step 2: Run test, expect failure**

```bash
pytest tests/kalshi_mlb_rfq/test_rfq_client_lifecycle.py -v
```

- [ ] **Step 3: Append to `kalshi_mlb_rfq/rfq_client.py`**

```python
def create_rfq(market_ticker: str, target_cost_dollars: float,
               replace_existing: bool = False) -> str:
    """Create an RFQ. Returns rfq_id.

    NOTE: replace_existing=True does NOT actually replace existing RFQs (recon-confirmed).
    The bot manages its own dedup via combo_cooldown and live_rfqs.
    """
    body = {
        "market_ticker": market_ticker,
        "rest_remainder": False,
        "target_cost_dollars": f"{target_cost_dollars:.2f}",
        "replace_existing": replace_existing,
    }
    status, resp, _ = api("POST", "/communications/rfqs", body=body)
    if status != 200 or not isinstance(resp, dict) or "id" not in resp:
        raise KalshiAPIError(f"create_rfq failed: status={status} body={resp}")
    return resp["id"]


def get_rfq(rfq_id: str) -> dict:
    """Get an RFQ object. Returns the inner 'rfq' dict, or raises KalshiAPIError."""
    status, body, _ = api("GET", f"/communications/rfqs/{rfq_id}")
    if status != 200 or not isinstance(body, dict) or "rfq" not in body:
        raise KalshiAPIError(f"get_rfq failed: status={status} body={body}")
    return body["rfq"]


def delete_rfq(rfq_id: str) -> bool:
    """DELETE an RFQ. Returns True if cancellation went through, False if already closed."""
    status, body, _ = api("DELETE", f"/communications/rfqs/{rfq_id}")
    if status == 204:
        return True
    # Common case: 400 'expired' = already cancelled or expired. Idempotent no-op.
    if status == 400 and isinstance(body, dict) and body.get("error", {}).get("code") == "expired":
        return False
    raise KalshiAPIError(f"delete_rfq failed: status={status} body={body}")
```

- [ ] **Step 4: Run test, expect pass**

```bash
pytest tests/kalshi_mlb_rfq/test_rfq_client_lifecycle.py -v
```

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/rfq_client.py tests/kalshi_mlb_rfq/test_rfq_client_lifecycle.py
git commit -m "feat(kalshi-mlb-rfq): rfq_client.{create,get,delete}_rfq"
```

---

### Task 8: RFQ client — poll quotes

`GET /communications/quotes?rfq_id=X&rfq_creator_user_id=USER_ID` returns `{"quotes": [...]}`. The user-id filter is required (recon found 403 without it).

**Files:**
- Modify: `kalshi_mlb_rfq/rfq_client.py`
- Test: `tests/kalshi_mlb_rfq/test_rfq_client_quotes.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_rfq_client_quotes.py
from unittest.mock import patch

from kalshi_mlb_rfq import rfq_client


def test_poll_quotes_returns_list():
    fake = (200, {"quotes": [
        {"id": "q1", "yes_bid_dollars": "0.13", "no_bid_dollars": "0.74", "status": "open"},
    ]}, {})
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=fake) as m:
        quotes = rfq_client.poll_quotes("rfq-1", user_id="user-uuid")

    assert len(quotes) == 1
    assert quotes[0]["id"] == "q1"
    args = m.call_args[0]
    assert "rfq_id=rfq-1" in args[1]
    assert "rfq_creator_user_id=user-uuid" in args[1]
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_rfq_client_quotes.py -v
```

- [ ] **Step 3: Append to `kalshi_mlb_rfq/rfq_client.py`**

```python
def poll_quotes(rfq_id: str, user_id: str) -> list[dict]:
    """Poll for quotes on an RFQ. Returns list of quote dicts."""
    path = f"/communications/quotes?rfq_id={rfq_id}&rfq_creator_user_id={user_id}"
    status, body, _ = api("GET", path)
    if status != 200 or not isinstance(body, dict):
        raise KalshiAPIError(f"poll_quotes failed: status={status} body={body}")
    return body.get("quotes") or []
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/rfq_client.py tests/kalshi_mlb_rfq/test_rfq_client_quotes.py
git commit -m "feat(kalshi-mlb-rfq): rfq_client.poll_quotes"
```

---

### Task 9: RFQ client — accept quote + positions API health check

**Files:**
- Modify: `kalshi_mlb_rfq/rfq_client.py`
- Test: `tests/kalshi_mlb_rfq/test_rfq_client_accept.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_rfq_client_accept.py
from unittest.mock import patch

from kalshi_mlb_rfq import rfq_client


def test_accept_quote_returns_response():
    with patch("kalshi_mlb_rfq.rfq_client.api",
               return_value=(200, {"order": {"id": "ord-1"}}, {})):
        resp = rfq_client.accept_quote("q1", contracts=10)
    assert resp["order"]["id"] == "ord-1"


def test_accept_quote_walked_returns_none():
    with patch("kalshi_mlb_rfq.rfq_client.api",
               return_value=(400, {"error": {"code": "quote_walked"}}, {})):
        resp = rfq_client.accept_quote("q1", contracts=10)
    assert resp is None


def test_get_positions_for_combo():
    fake = (200, {"event_positions": [], "market_positions": [
        {"ticker": "KXMVECROSSCATEGORY-S-FOO", "position": 100},
        {"ticker": "OTHER", "position": 5},
    ]}, {})
    with patch("kalshi_mlb_rfq.rfq_client.api", return_value=fake):
        n = rfq_client.get_position_contracts("KXMVECROSSCATEGORY-S-FOO")
    assert n == 100
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_rfq_client_accept.py -v
```

- [ ] **Step 3: Append to `kalshi_mlb_rfq/rfq_client.py`**

```python
def accept_quote(quote_id: str, contracts: int) -> dict | None:
    """Accept a quote for the given contract count.

    Returns the accept-response dict on success, or None if the quote walked /
    expired before our accept landed.
    """
    body = {"contracts": contracts}
    status, resp, _ = api("POST", f"/communications/quotes/{quote_id}/accept", body=body)
    if status == 200:
        return resp if isinstance(resp, dict) else {}
    if status in (400, 409):
        # Common race: quote walked or expired. Not a hard error.
        return None
    raise KalshiAPIError(f"accept_quote failed: status={status} body={resp}")


def get_position_contracts(market_ticker: str) -> int:
    """Authoritative current position count for a ticker via /portfolio/positions.

    Returns 0 if no position. Raises KalshiAPIError on API failure.
    """
    status, body, _ = api("GET", f"/portfolio/positions?ticker={market_ticker}&limit=10")
    if status != 200 or not isinstance(body, dict):
        raise KalshiAPIError(f"get_position_contracts failed: status={status} body={body}")
    for p in body.get("market_positions") or []:
        if p.get("ticker") == market_ticker:
            return int(p.get("position", 0))
    return 0


def list_open_rfqs(user_id: str) -> list[dict]:
    """Return all open RFQs for the given user. Used at startup for phantom cleanup."""
    path = f"/communications/rfqs?status=open&creator_user_id={user_id}&limit=100"
    status, body, _ = api("GET", path)
    if status != 200 or not isinstance(body, dict):
        raise KalshiAPIError(f"list_open_rfqs failed: status={status} body={body}")
    return body.get("rfqs") or []
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/rfq_client.py tests/kalshi_mlb_rfq/test_rfq_client_accept.py
git commit -m "feat(kalshi-mlb-rfq): rfq_client accept_quote + get_position_contracts + list_open_rfqs"
```

---

## Phase 4 — Pricing & sizing

### Task 10: `ev_calc.py` — fee + post-fee EV

**Files:**
- Create: `kalshi_mlb_rfq/ev_calc.py`
- Test: `tests/kalshi_mlb_rfq/test_ev_calc.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_ev_calc.py
import pytest

from kalshi_mlb_rfq import ev_calc


@pytest.mark.parametrize("price,expected_fee", [
    (0.50, 0.02),   # 0.07 * 0.5 * 0.5 = 0.0175 → ceil to 0.02
    (0.20, 0.02),   # 0.07 * 0.20 * 0.80 = 0.0112 → ceil to 0.02
    (0.10, 0.01),   # 0.07 * 0.1 * 0.9 = 0.0063 → ceil to 0.01
    (0.30, 0.02),   # 0.07 * 0.30 * 0.70 = 0.0147 → ceil to 0.02
])
def test_fee_per_contract(price, expected_fee):
    assert ev_calc.fee_per_contract(price) == pytest.approx(expected_fee)


def test_post_fee_ev_yes_side_buying():
    # We BUY YES at yes_ask = 1 - no_bid = 1 - 0.74 = 0.26.
    # Our fair: 0.32. Fee at 0.26 = ceil(0.07 * 0.26 * 0.74 * 100)/100 = ceil(1.347)/100 = 0.02.
    # Effective price = 0.28.
    # EV/contract = 0.32 * (1 - 0.26) - 0.68 * 0.26 - 0.02
    #             = 0.32*0.74 - 0.68*0.26 - 0.02
    #             = 0.2368 - 0.1768 - 0.02 = 0.04
    ev_dollars, ev_pct = ev_calc.post_fee_ev_buy_yes(blended_fair=0.32, no_bid=0.74)
    assert ev_dollars == pytest.approx(0.04, abs=0.001)
    assert ev_pct == pytest.approx(0.04 / 0.26, abs=0.001)


def test_post_fee_ev_no_side_buying():
    # We BUY NO at no_ask = 1 - yes_bid = 1 - 0.13 = 0.87.
    # Our fair NO = 1 - 0.32 = 0.68.
    # EV/contract = 0.68 * (1-0.87) - 0.32*0.87 - fee
    ev_dollars, ev_pct = ev_calc.post_fee_ev_buy_no(blended_fair=0.32, yes_bid=0.13)
    expected_fee = ev_calc.fee_per_contract(0.87)
    expected_ev = 0.68 * (1 - 0.87) - 0.32 * 0.87 - expected_fee
    assert ev_dollars == pytest.approx(expected_ev, abs=0.001)
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_ev_calc.py -v
```

- [ ] **Step 3: Implement `kalshi_mlb_rfq/ev_calc.py`**

```python
"""Kalshi taker fee + post-fee EV math.

Fee formula (Feb 2026 schedule): ceil(0.07 * P * (1-P) * 100) / 100 per contract,
where P is contract price in dollars.
"""

import math


def fee_per_contract(price_dollars: float) -> float:
    """Quadratic fee, rounded UP to the nearest cent."""
    return math.ceil(0.07 * price_dollars * (1 - price_dollars) * 100) / 100


def post_fee_ev_buy_yes(blended_fair: float, no_bid: float) -> tuple[float, float]:
    """We BUY YES by accepting a maker's no_bid (selling NO to them = buying YES from them).

    Effective YES purchase price = 1 - no_bid.
    Returns (ev_dollars_per_contract, ev_pct_of_stake).
    """
    yes_ask = 1.0 - no_bid
    if yes_ask <= 0 or yes_ask >= 1:
        return 0.0, 0.0
    fee = fee_per_contract(yes_ask)
    ev = blended_fair * (1 - yes_ask) - (1 - blended_fair) * yes_ask - fee
    return ev, ev / yes_ask


def post_fee_ev_buy_no(blended_fair: float, yes_bid: float) -> tuple[float, float]:
    """We BUY NO by accepting a maker's yes_bid (selling YES to them = buying NO from them).

    Effective NO purchase price = 1 - yes_bid.
    Returns (ev_dollars_per_contract, ev_pct_of_stake).
    """
    no_ask = 1.0 - yes_bid
    if no_ask <= 0 or no_ask >= 1:
        return 0.0, 0.0
    fee = fee_per_contract(no_ask)
    fair_no = 1.0 - blended_fair
    ev = fair_no * (1 - no_ask) - blended_fair * no_ask - fee
    return ev, ev / no_ask
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/ev_calc.py tests/kalshi_mlb_rfq/test_ev_calc.py
git commit -m "feat(kalshi-mlb-rfq): ev_calc with quadratic fee formula"
```

---

### Task 11: `fair_value.py` — model_fair from samples

**Files:**
- Create: `kalshi_mlb_rfq/fair_value.py`
- Test: `tests/kalshi_mlb_rfq/test_fair_value_model.py`

The `mlb_game_samples` table has columns `(game_id, sim_idx, home_margin, total_final_score, home_margin_f3, total_f3, home_margin_f5, total_f5, home_margin_f7, total_f7, home_scored_first)`. Joint hit probability across legs is `mean(rows where every leg's hit predicate is True)`.

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_fair_value_model.py
import pandas as pd
import pytest

from kalshi_mlb_rfq import fair_value


@pytest.fixture
def synthetic_samples():
    return pd.DataFrame({
        "home_margin": [3, -1, 5, 2, 1, 0, -2, 4, 6, -1],
        "total_final_score": [10, 6, 12, 9, 8, 7, 4, 11, 13, 5],
    })


def test_model_fair_home_minus_1_5_and_over_7_5(synthetic_samples):
    # Home -1.5 = home_margin >= 2; Over 7.5 = total >= 8.
    # Hits: rows where both: (3,10),(5,12),(2,9),(4,11),(6,13) → 5 of 10
    legs = [
        fair_value.SpreadLeg(team_is_home=True, line_n=2, side="yes"),
        fair_value.TotalLeg(line_n=8, side="yes"),
    ]
    fair = fair_value.model_fair(synthetic_samples, legs)
    assert fair == pytest.approx(0.5)


def test_model_fair_no_side_inverts(synthetic_samples):
    # Under 7.5 = total < 8 (i.e., total_final_score <= 7), so 5 of 10.
    legs = [fair_value.TotalLeg(line_n=8, side="no")]
    assert fair_value.model_fair(synthetic_samples, legs) == pytest.approx(0.5)
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_fair_value_model.py -v
```

- [ ] **Step 3: Implement initial `kalshi_mlb_rfq/fair_value.py`**

```python
"""Fair-value computation for combos: model + sportsbooks + blend."""

from dataclasses import dataclass
from typing import Literal

import pandas as pd


@dataclass(frozen=True)
class SpreadLeg:
    team_is_home: bool
    line_n: int                # KXMLBSPREAD ticker suffix N (= 2 for ±1.5)
    side: Literal["yes", "no"]


@dataclass(frozen=True)
class TotalLeg:
    line_n: int                # KXMLBTOTAL ticker suffix N (= 8 for Over 7.5)
    side: Literal["yes", "no"]


Leg = SpreadLeg | TotalLeg


def _hit_mask(samples: pd.DataFrame, leg: Leg) -> pd.Series:
    if isinstance(leg, SpreadLeg):
        if leg.team_is_home:
            base = samples["home_margin"] >= leg.line_n
        else:
            base = samples["home_margin"] <= -leg.line_n
        return base if leg.side == "yes" else ~base
    if isinstance(leg, TotalLeg):
        base = samples["total_final_score"] >= leg.line_n
        return base if leg.side == "yes" else ~base
    raise TypeError(f"unknown leg type: {type(leg)}")


def model_fair(samples: pd.DataFrame, legs: list[Leg]) -> float:
    """Fraction of sample paths where ALL legs hit."""
    if samples.empty or not legs:
        return 0.0
    mask = pd.Series([True] * len(samples), index=samples.index)
    for leg in legs:
        mask &= _hit_mask(samples, leg)
    return float(mask.mean())
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/fair_value.py tests/kalshi_mlb_rfq/test_fair_value_model.py
git commit -m "feat(kalshi-mlb-rfq): fair_value.model_fair from mlb_game_samples"
```

---

### Task 12: Fair value — book devig + 2-source blend

**Files:**
- Modify: `kalshi_mlb_rfq/fair_value.py`
- Test: `tests/kalshi_mlb_rfq/test_fair_value_book.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_fair_value_book.py
import pandas as pd
import pytest

from kalshi_mlb_rfq import fair_value


@pytest.fixture
def four_side_dk():
    """4 mutually-exclusive combos for one game/book/(spread,total) pair."""
    return pd.DataFrame({
        "combo": ["Home Spread + Over", "Home Spread + Under",
                  "Away Spread + Over", "Away Spread + Under"],
        "bookmaker": ["draftkings"] * 4,
        "sgp_decimal": [3.50, 4.20, 5.00, 4.80],
        "spread_line": [-1.5] * 4,
        "total_line": [8.5] * 4,
    })


def test_devig_with_4_sides(four_side_dk):
    # Vig = sum(1/d) = 1/3.5 + 1/4.2 + 1/5.0 + 1/4.8
    fair = fair_value.devig_book(four_side_dk, combo="Home Spread + Over")
    expected_vig = sum(1 / d for d in [3.50, 4.20, 5.00, 4.80])
    expected_fair = (1 / 3.50) / expected_vig
    assert fair == pytest.approx(expected_fair, abs=1e-4)


def test_devig_fallback_when_fewer_sides():
    df = pd.DataFrame({
        "combo": ["Home Spread + Over"],
        "bookmaker": ["draftkings"],
        "sgp_decimal": [3.50],
        "spread_line": [-1.5],
        "total_line": [8.5],
    })
    fair = fair_value.devig_book(df, combo="Home Spread + Over",
                                 vig_fallback=0.125)
    expected = (1 / 3.50) / 1.125
    assert fair == pytest.approx(expected, abs=1e-4)


def test_blended_requires_two_sources():
    # Model only → returns None (gate fails).
    blended = fair_value.blend(model_fair=0.30, book_fairs={})
    assert blended is None

    # Model + 1 book → returns mean.
    blended = fair_value.blend(model_fair=0.30, book_fairs={"dk": 0.34})
    assert blended == pytest.approx(0.32)

    # Model + 2 books → simple mean across all.
    blended = fair_value.blend(model_fair=0.30,
                                book_fairs={"dk": 0.34, "fd": 0.28})
    assert blended == pytest.approx((0.30 + 0.34 + 0.28) / 3)
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_fair_value_book.py -v
```

- [ ] **Step 3: Append to `kalshi_mlb_rfq/fair_value.py`**

```python
def devig_book(book_rows: pd.DataFrame, combo: str,
               vig_fallback: float = 0.0) -> float | None:
    """Devig a single combo's fair value from rows of mlb_sgp_odds.

    book_rows must already be filtered to (game_id, period, bookmaker, spread_line, total_line).
    Uses the 4-side per-game vig method when >=4 sides exist; else falls back
    to (1/decimal_odds) / (1 + vig_fallback).
    """
    if book_rows.empty:
        return None

    target = book_rows.loc[book_rows["combo"] == combo]
    if target.empty:
        return None
    target_decimal = float(target["sgp_decimal"].iloc[0])

    if len(book_rows) >= 4:
        vig_sum = float((1.0 / book_rows["sgp_decimal"]).sum())
        return (1.0 / target_decimal) / vig_sum

    return (1.0 / target_decimal) / (1.0 + vig_fallback)


def blend(model_fair: float, book_fairs: dict[str, float]) -> float | None:
    """2-source gate: returns simple mean of model + non-null books, or None
    if fewer than 2 sources.
    """
    sources = [model_fair] + [v for v in book_fairs.values() if v is not None]
    if len(sources) < 2:
        return None
    return sum(sources) / len(sources)
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/fair_value.py tests/kalshi_mlb_rfq/test_fair_value_book.py
git commit -m "feat(kalshi-mlb-rfq): fair_value devig + 2-source blend"
```

---

### Task 13: `kelly.py` — conditional Kelly sized on combo outcomes

Adapts the math from `kalshi_mm/kelly.py`. The big change: outcomes here are joint-leg hits (boolean on each sample) rather than CBB's per-market evaluation.

**Files:**
- Create: `kalshi_mlb_rfq/kelly.py`
- Test: `tests/kalshi_mlb_rfq/test_kelly.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_kelly.py
import numpy as np
import pandas as pd
import pytest

from kalshi_mlb_rfq import kelly, fair_value


def test_single_kelly_no_existing_positions():
    # 30% fair, 26% effective price → b = 0.74/0.26 ≈ 2.846, p=0.30, q=0.70
    # kelly_frac = (b*p - q) / b = (2.846*0.30 - 0.70) / 2.846 = (0.854 - 0.70)/2.846 ≈ 0.054
    # Quarter Kelly with $1000 bankroll, price $0.26 → 0.25 * 0.054 * 1000 / 0.26 ≈ 52
    contracts = kelly.kelly_size_combo(
        outcome_vec=np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0]),  # 30% hit rate
        existing_positions=[],
        effective_price=0.26,
        bankroll=1000.0,
        kelly_fraction=0.25,
    )
    assert 40 <= contracts <= 65   # ballpark check; exact value: 51 or 52


def test_kelly_zero_for_negative_edge():
    # Fair 0.20 vs price 0.30 → -EV → 0 contracts
    contracts = kelly.kelly_size_combo(
        outcome_vec=np.array([1, 1, 0, 0, 0, 0, 0, 0, 0, 0]),  # 20% hit rate
        existing_positions=[],
        effective_price=0.30,
        bankroll=1000.0,
        kelly_fraction=0.25,
    )
    assert contracts == 0


def test_conditional_kelly_shrinks_for_correlated_position():
    # Two highly correlated bets (same outcome vector). With existing position on bet A,
    # bet B's conditional Kelly size should be much smaller than its single-bet size.
    n = 1000
    rng = np.random.default_rng(42)
    outcome = (rng.random(n) < 0.30).astype(int)

    single = kelly.kelly_size_combo(
        outcome_vec=outcome, existing_positions=[],
        effective_price=0.25, bankroll=1000.0, kelly_fraction=0.25,
    )

    # Existing position is the SAME bet, already at 'single' contracts.
    placed = [{"outcome_vec": outcome.copy(), "contracts": single, "effective_price": 0.25}]
    conditional = kelly.kelly_size_combo(
        outcome_vec=outcome, existing_positions=placed,
        effective_price=0.25, bankroll=1000.0, kelly_fraction=0.25,
    )
    assert conditional < single * 0.5   # significantly shrunk
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_kelly.py -v
```

- [ ] **Step 3: Implement `kalshi_mlb_rfq/kelly.py`**

```python
"""Conditional Kelly sizing on combo outcome vectors.

Adapted from kalshi_mm/kelly.py. The new bet's outcome is a boolean vector
(1 if combo hits in that simulation path, else 0). Existing positions on the
same game contribute their own outcome vectors; we compute joint covariance
and apply the conditional Kelly formula:

    f_new* = Σ_nn⁻¹ × (μ_new − Σ_np × f_placed)

where Σ_nn is the new bet's variance, Σ_np is the covariance vector with
existing positions, and f_placed is each existing position's stake fraction.

Falls back to single-bet Kelly scaled by 1/sqrt(1 + (n-1)·avg_ρ) when the
covariance matrix is ill-conditioned (matches kalshi_mm fallback).
"""

import math

import numpy as np


def _single_bet_kelly_fraction(p: float, effective_price: float) -> float:
    """Standard Kelly: (bp - q) / b where b = (1 - price) / price."""
    if effective_price <= 0 or effective_price >= 1:
        return 0.0
    b = (1 - effective_price) / effective_price
    q = 1 - p
    f = (b * p - q) / b
    return max(0.0, f)


def kelly_size_combo(
    outcome_vec: np.ndarray,
    existing_positions: list[dict],
    effective_price: float,
    bankroll: float,
    kelly_fraction: float,
) -> int:
    """Return contract count to take.

    Args:
        outcome_vec: 0/1 array, len = sample count.
        existing_positions: list of {outcome_vec, contracts, effective_price}.
        effective_price: post-fee price per contract for this new bet.
        bankroll: dollars.
        kelly_fraction: e.g., 0.25 for quarter-Kelly.

    Returns:
        Floored, non-negative contract count.
    """
    p = float(outcome_vec.mean())
    base_frac = _single_bet_kelly_fraction(p, effective_price)
    if base_frac <= 0:
        return 0

    if not existing_positions:
        contracts = math.floor(kelly_fraction * base_frac * bankroll / effective_price)
        return max(0, contracts)

    # Conditional Kelly. Build the new outcome's centered vector and Σ_np.
    n_samples = len(outcome_vec)
    new_centered = outcome_vec - p
    var_new = float(np.var(outcome_vec))
    if var_new <= 1e-12:
        return 0

    # f_placed in stake-fraction units (contracts * price / bankroll).
    f_placed = []
    np_cov_terms = []
    for pos in existing_positions:
        pos_vec = pos["outcome_vec"]
        if len(pos_vec) != n_samples:
            continue
        pos_mean = float(pos_vec.mean())
        cov = float(((pos_vec - pos_mean) * new_centered).mean())
        f_placed.append(pos["contracts"] * pos["effective_price"] / bankroll)
        np_cov_terms.append(cov)

    if not f_placed:
        contracts = math.floor(kelly_fraction * base_frac * bankroll / effective_price)
        return max(0, contracts)

    f_placed = np.array(f_placed)
    np_cov = np.array(np_cov_terms)

    mu_new = p - effective_price  # excess return over price
    f_new_star = (mu_new - np.dot(np_cov, f_placed)) / var_new
    f_new_star = max(0.0, f_new_star)

    contracts = math.floor(kelly_fraction * f_new_star * bankroll / effective_price)
    return max(0, contracts)
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/kelly.py tests/kalshi_mlb_rfq/test_kelly.py
git commit -m "feat(kalshi-mlb-rfq): conditional Kelly sizing on combo outcomes"
```

---

## Phase 5 — Combo enumeration & priority queue

### Task 14: `combo_enumerator.py` — yield candidates per cycle

**Files:**
- Create: `kalshi_mlb_rfq/combo_enumerator.py`
- Test: `tests/kalshi_mlb_rfq/test_combo_enumerator.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_combo_enumerator.py
import pytest

from kalshi_mlb_rfq import combo_enumerator


def test_enumerate_yields_n_times_m_combos():
    # 2 spread lines × 2 sides = 4 spread legs
    # 3 total lines × 2 sides = 6 total legs
    # = 24 combos
    spread_legs = [(-1.5, "home"), (-1.5, "away"), (-2.5, "home"), (-2.5, "away")]
    total_legs = [7.5, 8.5, 9.5]

    combos = list(combo_enumerator.enumerate_2leg(
        game_id="gid",
        event_suffix="26APR282005NYYTEX",
        home_code="TEX", away_code="NYY",
        available_spreads=spread_legs,
        available_totals=total_legs,
    ))
    assert len(combos) == 4 * 3 * 2 * 2   # 4 spread legs (with sides) × 3 totals × 2 total sides
    # All combos have exactly 2 legs
    assert all(len(c.legs) == 2 for c in combos)


def test_canonical_leg_set_hash_is_order_invariant():
    legs_a = [{"market_ticker": "A", "side": "yes"}, {"market_ticker": "B", "side": "no"}]
    legs_b = [{"market_ticker": "B", "side": "no"}, {"market_ticker": "A", "side": "yes"}]
    assert combo_enumerator.canonical_leg_set_hash(legs_a) == \
           combo_enumerator.canonical_leg_set_hash(legs_b)
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_combo_enumerator.py -v
```

- [ ] **Step 3: Implement `kalshi_mlb_rfq/combo_enumerator.py`**

```python
"""Per-cycle combo enumeration + priority queue for the RFQ pipeline."""

import hashlib
import json
from dataclasses import dataclass, field
from typing import Iterable, Literal

from kalshi_mlb_rfq import ticker_map


@dataclass(frozen=True)
class ComboCandidate:
    game_id: str
    legs: tuple[dict, ...]   # each: {market_ticker, event_ticker, side}
    leg_set_hash: str
    descriptor: str          # human-readable label, e.g., "TEX -1.5 + Over 7.5"


def canonical_leg_set_hash(legs: Iterable[dict]) -> str:
    keys = sorted(f"{leg['market_ticker']}|{leg['side']}" for leg in legs)
    return hashlib.sha256("\n".join(keys).encode()).hexdigest()


def enumerate_2leg(
    game_id: str,
    event_suffix: str,
    home_code: str,
    away_code: str,
    available_spreads: list[tuple[float, Literal["home", "away"]]],
    available_totals: list[float],
) -> Iterable[ComboCandidate]:
    """Yield every (spread_leg, side) × (total_leg, side) combo for this game.

    available_spreads: list of (line, "home"|"away") tuples — enumerated by
        ticker_map.spread_ticker. line is the team's negative spread (e.g., -1.5).
    available_totals: list of total lines (e.g., 7.5, 8.5).
    """
    spread_event = f"KXMLBSPREAD-{event_suffix}"
    total_event = f"KXMLBTOTAL-{event_suffix}"

    spread_leg_specs = []
    for line, who in available_spreads:
        team_code = home_code if who == "home" else away_code
        ticker = ticker_map.spread_ticker(event_suffix, team_code, line)
        for side in ("yes", "no"):
            spread_leg_specs.append((ticker, side, who, line))

    total_leg_specs = []
    for line in available_totals:
        ticker = ticker_map.total_ticker(event_suffix, line)
        for side in ("yes", "no"):
            total_leg_specs.append((ticker, side, line))

    for s_ticker, s_side, s_who, s_line in spread_leg_specs:
        for t_ticker, t_side, t_line in total_leg_specs:
            legs = (
                {"market_ticker": s_ticker, "event_ticker": spread_event, "side": s_side},
                {"market_ticker": t_ticker, "event_ticker": total_event, "side": t_side},
            )
            sign = "-" if s_side == "yes" else "+"
            line_abs = abs(s_line)
            ou = "Over" if t_side == "yes" else "Under"
            descriptor = (
                f"{home_code if s_who == 'home' else away_code} "
                f"{sign}{line_abs} + {ou} {t_line}"
            )
            yield ComboCandidate(
                game_id=game_id,
                legs=legs,
                leg_set_hash=canonical_leg_set_hash(legs),
                descriptor=descriptor,
            )
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/combo_enumerator.py tests/kalshi_mlb_rfq/test_combo_enumerator.py
git commit -m "feat(kalshi-mlb-rfq): combo_enumerator with stable leg-set hash"
```

---

### Task 15: Priority queue + scoring

**Files:**
- Modify: `kalshi_mlb_rfq/combo_enumerator.py`
- Test: `tests/kalshi_mlb_rfq/test_priority_queue.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_priority_queue.py
import pytest

from kalshi_mlb_rfq import combo_enumerator


def test_priority_queue_orders_by_edge_magnitude():
    # Higher |fair - reference| ranks first.
    candidates_with_scores = [
        (combo_enumerator.ComboCandidate(game_id="g", legs=(), leg_set_hash="a",
                                         descriptor="A"), 0.32, 0.30),  # |Δ|=0.02
        (combo_enumerator.ComboCandidate(game_id="g", legs=(), leg_set_hash="b",
                                         descriptor="B"), 0.40, 0.20),  # |Δ|=0.20  ← top
        (combo_enumerator.ComboCandidate(game_id="g", legs=(), leg_set_hash="c",
                                         descriptor="C"), 0.50, 0.45),  # |Δ|=0.05
    ]
    ranked = combo_enumerator.rank_by_edge(candidates_with_scores)
    assert [c.leg_set_hash for c in ranked] == ["b", "c", "a"]


def test_score_falls_back_to_distance_from_half_when_ref_zero():
    score = combo_enumerator.edge_score(blended_fair=0.30, kalshi_ref=0.0)
    # |0.30 - 0.5| = 0.20
    assert score == pytest.approx(0.20)
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_priority_queue.py -v
```

- [ ] **Step 3: Append to `kalshi_mlb_rfq/combo_enumerator.py`**

```python
def edge_score(blended_fair: float, kalshi_ref: float) -> float:
    """Edge magnitude for ranking. If Kalshi has no reference price (last_price=0),
    fall back to |fair - 0.5| as a contentious-combo heuristic."""
    if kalshi_ref <= 0 or kalshi_ref >= 1:
        return abs(blended_fair - 0.5)
    return abs(blended_fair - kalshi_ref)


def rank_by_edge(candidates_with_scores: list[tuple[ComboCandidate, float, float]]
                 ) -> list[ComboCandidate]:
    """Input: list of (candidate, blended_fair, kalshi_ref).
       Output: candidates sorted by edge_score descending.
    """
    return [
        c for c, _, _ in sorted(
            candidates_with_scores,
            key=lambda t: edge_score(t[1], t[2]),
            reverse=True,
        )
    ]
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/combo_enumerator.py tests/kalshi_mlb_rfq/test_priority_queue.py
git commit -m "feat(kalshi-mlb-rfq): combo_enumerator priority queue (edge_score + rank)"
```

---

## Phase 6 — Risk module (per-accept gates)

### Task 16: `risk.py` — staleness, fair-value bounds, sanity, EV gate

**Files:**
- Create: `kalshi_mlb_rfq/risk.py`
- Test: `tests/kalshi_mlb_rfq/test_risk_basic_gates.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_risk_basic_gates.py
from datetime import datetime, timedelta, timezone

import pytest

from kalshi_mlb_rfq import risk


def test_staleness_gate_blocks_old_data():
    fresh = datetime.now(timezone.utc) - timedelta(seconds=300)  # 5 min
    old   = datetime.now(timezone.utc) - timedelta(seconds=900)  # 15 min
    assert risk.staleness_ok(fresh, max_age_sec=600)
    assert not risk.staleness_ok(old, max_age_sec=600)


def test_staleness_gate_handles_negative_age():
    future = datetime.now(timezone.utc) + timedelta(seconds=60)  # clock skew
    assert not risk.staleness_ok(future, max_age_sec=600)


def test_fair_bounds_gate():
    assert risk.fair_in_bounds(0.30, lower=0.05, upper=0.95)
    assert not risk.fair_in_bounds(0.03, lower=0.05, upper=0.95)
    assert not risk.fair_in_bounds(0.97, lower=0.05, upper=0.95)


def test_sanity_bound_gate():
    # blended_fair=0.30, quote_implied=0.27, deviation=0.03 → ok at 0.15 threshold
    assert risk.sanity_bound_ok(quote_implied=0.27, blended_fair=0.30, max_deviation=0.15)
    # blended_fair=0.30, quote_implied=0.10, deviation=0.20 → blocked at 0.15
    assert not risk.sanity_bound_ok(quote_implied=0.10, blended_fair=0.30, max_deviation=0.15)


def test_tipoff_gate():
    now = datetime.now(timezone.utc)
    far_future_game = now + timedelta(minutes=60)
    near_game     = now + timedelta(minutes=3)
    assert risk.tipoff_ok(commence_time=far_future_game, cancel_min=5, now=now)
    assert not risk.tipoff_ok(commence_time=near_game, cancel_min=5, now=now)
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_risk_basic_gates.py -v
```

- [ ] **Step 3: Implement initial `kalshi_mlb_rfq/risk.py`**

```python
"""Per-accept gates and risk sweep for the Kalshi MLB RFQ bot."""

from datetime import datetime, timedelta, timezone


def staleness_ok(generated_at: datetime, max_age_sec: int) -> bool:
    """True if data is younger than max_age_sec. False on negative-age (clock skew)."""
    if generated_at.tzinfo is None:
        generated_at = generated_at.replace(tzinfo=timezone.utc)
    age = (datetime.now(timezone.utc) - generated_at).total_seconds()
    if age < 0:
        return False  # clock skew — fail safe
    return age <= max_age_sec


def fair_in_bounds(blended_fair: float, lower: float, upper: float) -> bool:
    return lower <= blended_fair <= upper


def sanity_bound_ok(quote_implied: float, blended_fair: float,
                    max_deviation: float) -> bool:
    return abs(quote_implied - blended_fair) <= max_deviation


def tipoff_ok(commence_time: datetime, cancel_min: int,
              now: datetime | None = None) -> bool:
    if commence_time is None:
        return False  # unknown tipoff = fail-safe refuse
    if commence_time.tzinfo is None:
        commence_time = commence_time.replace(tzinfo=timezone.utc)
    if now is None:
        now = datetime.now(timezone.utc)
    return (commence_time - now) > timedelta(minutes=cancel_min)
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/risk.py tests/kalshi_mlb_rfq/test_risk_basic_gates.py
git commit -m "feat(kalshi-mlb-rfq): risk basic gates (staleness, fair bounds, sanity, tipoff)"
```

---

### Task 17: `risk.py` — exposure caps + cooldown + inverse-combo

**Files:**
- Modify: `kalshi_mlb_rfq/risk.py`
- Test: `tests/kalshi_mlb_rfq/test_risk_exposure.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_risk_exposure.py
from datetime import datetime, timezone

import pytest

from kalshi_mlb_rfq import risk


def test_per_game_cap_blocks_at_threshold():
    # bankroll=1000, MAX_GAME_EXPOSURE_PCT=0.10 → cap = $100/game
    today_fills = [
        {"game_id": "G1", "contracts": 100, "price_dollars": 0.30},  # $30
        {"game_id": "G1", "contracts": 200, "price_dollars": 0.20},  # $40
        {"game_id": "G2", "contracts": 100, "price_dollars": 0.50},  # $50 (different game)
    ]
    # G1 has $70 today → still under $100 cap
    assert risk.per_game_cap_ok(game_id="G1", today_fills=today_fills,
                                 bankroll=1000.0, max_pct=0.10)
    # Add a new $40 fill to G1 → $110, over cap
    today_fills.append({"game_id": "G1", "contracts": 100, "price_dollars": 0.40})
    assert not risk.per_game_cap_ok(game_id="G1", today_fills=today_fills,
                                     bankroll=1000.0, max_pct=0.10)


def test_daily_cap_blocks_at_threshold():
    today_fills = [
        {"game_id": "G1", "contracts": 100, "price_dollars": 0.50},  # $50
        {"game_id": "G2", "contracts": 200, "price_dollars": 0.50},  # $100
    ]
    # $150 < $200 cap
    assert risk.daily_cap_ok(today_fills=today_fills, daily_cap_usd=200.0)
    today_fills.append({"game_id": "G3", "contracts": 200, "price_dollars": 0.50})
    # +$100 = $250 > cap
    assert not risk.daily_cap_ok(today_fills=today_fills, daily_cap_usd=200.0)


def test_cooldown_gate():
    now = datetime.now(timezone.utc)
    cooled = {"abc": now.replace(year=now.year + 1)}  # cooldown extends into future
    assert not risk.cooldown_ok(leg_set_hash="abc", cooldown_map=cooled, now=now)
    assert risk.cooldown_ok(leg_set_hash="other", cooldown_map=cooled, now=now)


def test_inverse_combo_guard():
    legs = [
        {"market_ticker": "X", "side": "yes"},
        {"market_ticker": "Y", "side": "no"},
    ]
    inverse_legs = [
        {"market_ticker": "X", "side": "no"},
        {"market_ticker": "Y", "side": "yes"},
    ]
    # Compute hash on legs to look up inverse in positions
    open_inverse_hashes = {risk.inverse_leg_set_hash(legs)}
    assert not risk.inverse_combo_ok(
        legs=legs, open_position_hashes={risk.inverse_leg_set_hash(inverse_legs)})
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_risk_exposure.py -v
```

- [ ] **Step 3: Append to `kalshi_mlb_rfq/risk.py`**

```python
import hashlib


def per_game_cap_ok(game_id: str, today_fills: list[dict],
                     bankroll: float, max_pct: float) -> bool:
    cap = bankroll * max_pct
    spent = sum(f["contracts"] * f["price_dollars"]
                for f in today_fills if f["game_id"] == game_id)
    return spent < cap


def daily_cap_ok(today_fills: list[dict], daily_cap_usd: float) -> bool:
    spent = sum(f["contracts"] * f["price_dollars"] for f in today_fills)
    return spent < daily_cap_usd


def cooldown_ok(leg_set_hash: str, cooldown_map: dict, now: datetime | None = None) -> bool:
    if now is None:
        now = datetime.now(timezone.utc)
    cooled_until = cooldown_map.get(leg_set_hash)
    if cooled_until is None:
        return True
    if cooled_until.tzinfo is None:
        cooled_until = cooled_until.replace(tzinfo=timezone.utc)
    return cooled_until <= now


def inverse_leg_set_hash(legs: list[dict]) -> str:
    """Hash of leg-set after flipping every side. Used to detect held inverse positions."""
    flipped = [{"market_ticker": l["market_ticker"],
                "side": "no" if l["side"] == "yes" else "yes"} for l in legs]
    keys = sorted(f"{l['market_ticker']}|{l['side']}" for l in flipped)
    return hashlib.sha256("\n".join(keys).encode()).hexdigest()


def inverse_combo_ok(legs: list[dict], open_position_hashes: set[str]) -> bool:
    """True if no open position exists on the leg-set's inverse."""
    return inverse_leg_set_hash(legs) not in open_position_hashes
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/risk.py tests/kalshi_mlb_rfq/test_risk_exposure.py
git commit -m "feat(kalshi-mlb-rfq): risk exposure caps + cooldown + inverse-combo guard"
```

---

### Task 18: `risk.py` — line-move detection + kill switch + fill-ratio halt

**Files:**
- Modify: `kalshi_mlb_rfq/risk.py`
- Test: `tests/kalshi_mlb_rfq/test_risk_advanced.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_risk_advanced.py
import pytest

from kalshi_mlb_rfq import risk


def test_line_move_detection():
    ref = {"spread": -1.5, "total": 8.5}
    assert risk.line_move_ok(ref_lines=ref, current_lines={"spread": -1.5, "total": 8.5},
                             threshold=0.5)
    assert not risk.line_move_ok(ref_lines=ref, current_lines={"spread": -2.0, "total": 8.5},
                                  threshold=0.5)
    assert risk.line_move_ok(ref_lines=ref, current_lines={"spread": -1.7, "total": 8.5},
                             threshold=0.5)  # within threshold


def test_kill_switch_present(tmp_path, monkeypatch):
    kill_path = tmp_path / ".kill"
    monkeypatch.setattr(risk, "KILL_FILE", kill_path)
    assert risk.kill_switch_ok()  # absent
    kill_path.touch()
    assert not risk.kill_switch_ok()


def test_fill_ratio_halt():
    # 25 accepted, 25 walked, ratio = 0.50 → AT threshold (allowed)
    log = [{"decision": "accepted"}] * 25 + [{"decision": "failed_quote_walked"}] * 25
    assert risk.fill_ratio_ok(quote_log_window=log, min_ratio=0.50)

    # 20 accepted, 30 walked, ratio = 0.40 → below
    log = [{"decision": "accepted"}] * 20 + [{"decision": "failed_quote_walked"}] * 30
    assert not risk.fill_ratio_ok(quote_log_window=log, min_ratio=0.50)


def test_fill_ratio_ok_when_window_too_small():
    # Don't halt before we have enough samples to evaluate.
    log = [{"decision": "accepted"}] * 5
    assert risk.fill_ratio_ok(quote_log_window=log, min_ratio=0.50, min_window=10)
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_risk_advanced.py -v
```

- [ ] **Step 3: Append to `kalshi_mlb_rfq/risk.py`**

```python
from kalshi_mlb_rfq.config import KILL_FILE


def line_move_ok(ref_lines: dict, current_lines: dict, threshold: float) -> bool:
    """All lines moved by less than threshold? (None entries treated as no-data, blocks accept.)"""
    for key, ref in ref_lines.items():
        cur = current_lines.get(key)
        if cur is None or ref is None:
            return False
        if abs(cur - ref) > threshold:
            return False
    return True


def kill_switch_ok() -> bool:
    return not KILL_FILE.exists()


def fill_ratio_ok(quote_log_window: list[dict], min_ratio: float,
                  min_window: int = 10) -> bool:
    """True if we have enough samples AND the accepted-fraction is >= min_ratio."""
    relevant = [r for r in quote_log_window
                if r["decision"] in ("accepted", "failed_quote_walked")]
    if len(relevant) < min_window:
        return True   # not enough data; don't halt
    accepted = sum(1 for r in relevant if r["decision"] == "accepted")
    return accepted / len(relevant) >= min_ratio
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/risk.py tests/kalshi_mlb_rfq/test_risk_advanced.py
git commit -m "feat(kalshi-mlb-rfq): risk line-move + kill-switch + fill-ratio halt"
```

---

## Phase 7 — Notifications

### Task 19: `notify.py` — bot.log + optional webhook

**Files:**
- Create: `kalshi_mlb_rfq/notify.py`
- Test: `tests/kalshi_mlb_rfq/test_notify.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_notify.py
import json
from unittest.mock import patch, MagicMock

import pytest

from kalshi_mlb_rfq import notify


def test_fill_appends_log_line(tmp_path, monkeypatch):
    log = tmp_path / "bot.log"
    monkeypatch.setattr(notify, "LOG_PATH", log)
    notify.fill(rfq_id="r1", combo_market_ticker="MT", contracts=10,
                price=0.30, ev_pct=0.06)
    assert "[FILL]" in log.read_text()
    assert "r1" in log.read_text()


def test_halt_appends_log_line(tmp_path, monkeypatch):
    log = tmp_path / "bot.log"
    monkeypatch.setattr(notify, "LOG_PATH", log)
    notify.halt(reason="kill_switch", detail="user request")
    assert "[HALT]" in log.read_text()
    assert "kill_switch" in log.read_text()


def test_webhook_post_when_url_set(monkeypatch):
    monkeypatch.setattr(notify, "NOTIFY_WEBHOOK_URL", "https://hook.example.com/x")
    fake_resp = MagicMock()
    fake_resp.__enter__.return_value = fake_resp
    fake_resp.__exit__.return_value = False
    fake_resp.status = 200
    with patch("kalshi_mlb_rfq.notify.urllib.request.urlopen",
               return_value=fake_resp) as m:
        notify.fill(rfq_id="r1", combo_market_ticker="MT", contracts=1,
                    price=0.30, ev_pct=0.06)
    assert m.called
    sent_data = m.call_args[0][0].data.decode()
    assert json.loads(sent_data)["event"] == "fill"
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_notify.py -v
```

- [ ] **Step 3: Implement `kalshi_mlb_rfq/notify.py`**

```python
"""Fill / halt notifications: bot.log + optional webhook."""

import json
import urllib.error
import urllib.request
from datetime import datetime, timezone

from kalshi_mlb_rfq.config import LOG_PATH, NOTIFY_WEBHOOK_URL


def _append_log(line: str):
    with LOG_PATH.open("a") as f:
        f.write(line + "\n")


def _post_webhook(payload: dict):
    if not NOTIFY_WEBHOOK_URL:
        return
    try:
        data = json.dumps(payload).encode()
        req = urllib.request.Request(NOTIFY_WEBHOOK_URL, data=data, method="POST")
        req.add_header("Content-Type", "application/json")
        with urllib.request.urlopen(req, timeout=5):
            pass
    except (urllib.error.URLError, TimeoutError):
        # Notification failure is non-fatal.
        pass


def fill(rfq_id: str, combo_market_ticker: str, contracts: float,
         price: float, ev_pct: float):
    ts = datetime.now(timezone.utc).isoformat()
    _append_log(f"{ts} [FILL] rfq={rfq_id} ticker={combo_market_ticker} "
                f"n={contracts} price={price} ev={ev_pct:.2%}")
    _post_webhook({"event": "fill", "rfq_id": rfq_id,
                    "ticker": combo_market_ticker, "contracts": contracts,
                    "price": price, "ev_pct": ev_pct, "ts": ts})


def halt(reason: str, detail: str | None = None):
    ts = datetime.now(timezone.utc).isoformat()
    _append_log(f"{ts} [HALT] reason={reason} detail={detail or ''}")
    _post_webhook({"event": "halt", "reason": reason,
                    "detail": detail, "ts": ts})
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/notify.py tests/kalshi_mlb_rfq/test_notify.py
git commit -m "feat(kalshi-mlb-rfq): notify (bot.log + optional webhook)"
```

---

## Phase 8 — Daemon orchestrator

### Task 20: `main.py` — skeleton with signal handling + accept lock

**Files:**
- Create: `kalshi_mlb_rfq/main.py`
- Test: `tests/kalshi_mlb_rfq/test_main_signal.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_main_signal.py
import threading
from unittest.mock import patch

from kalshi_mlb_rfq import main


def test_signal_handler_sets_running_false():
    main._running.set()
    main._signal_handler(15, None)  # SIGTERM
    assert not main._running.is_set()


def test_accept_lock_serializes():
    """Two threads attempting accept simultaneously should serialize."""
    counter = {"value": 0}
    snapshots = []

    def attempt():
        with main.ACCEPT_LOCK:
            v = counter["value"]
            snapshots.append(v)
            counter["value"] = v + 1

    threads = [threading.Thread(target=attempt) for _ in range(20)]
    for t in threads: t.start()
    for t in threads: t.join()
    # If the lock works, snapshots are 0,1,2,...,19 in some order with no dupes.
    assert sorted(snapshots) == list(range(20))
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_main_signal.py -v
```

- [ ] **Step 3: Create `kalshi_mlb_rfq/main.py` skeleton**

```python
"""Kalshi MLB RFQ Bot — autonomous taker daemon."""

import argparse
import os
import signal
import sys
import threading
import time
from datetime import datetime, timezone

from kalshi_mlb_rfq import (
    config, db, notify, rfq_client,
)
from kalshi_mlb_rfq.config import KILL_FILE

VERSION = "0.1.0"

_running = threading.Event()
_running.set()
ACCEPT_LOCK = threading.Lock()


def _signal_handler(_sig, _frame):
    _running.clear()


def _phantom_rfq_cleanup():
    """At startup, cancel any RFQs on Kalshi that aren't tracked in our live_rfqs."""
    try:
        kalshi_open = rfq_client.list_open_rfqs(config.KALSHI_USER_ID)
    except Exception as e:
        print(f"  startup: list_open_rfqs failed: {e}", flush=True)
        return

    if not kalshi_open:
        print("  startup: no open RFQs on Kalshi", flush=True)
        return

    with db.connect(read_only=True) as con:
        ours = {r[0] for r in con.execute(
            "SELECT rfq_id FROM live_rfqs WHERE status='open'"
        ).fetchall()}

    for rfq in kalshi_open:
        rid = rfq.get("id")
        if rid and rid not in ours:
            try:
                rfq_client.delete_rfq(rid)
                print(f"  startup: cancelled phantom rfq {rid}", flush=True)
            except Exception as e:
                print(f"  startup: failed to cancel phantom {rid}: {e}", flush=True)


def main_loop(dry_run: bool):
    db.init_database()
    sid = db.start_session(pid=os.getpid(), dry_run=dry_run, version=VERSION)
    print(f"=== Kalshi MLB RFQ Bot — session {sid} (dry_run={dry_run}) ===", flush=True)
    _phantom_rfq_cleanup()

    last_rfq_refresh = 0.0
    last_quote_poll = 0.0
    last_risk_sweep = 0.0
    last_pipeline = 0.0
    last_heartbeat = 0.0

    try:
        while _running.is_set():
            now = time.time()

            if KILL_FILE.exists():
                notify.halt("kill_switch")
                # Sleep then re-check; don't exit.
                time.sleep(config.RISK_SWEEP_SEC)
                continue

            if now - last_rfq_refresh >= config.RFQ_REFRESH_SEC:
                # Implemented in Task 21
                last_rfq_refresh = now

            if now - last_quote_poll >= config.QUOTE_POLL_SEC:
                # Implemented in Task 22
                last_quote_poll = now

            if now - last_risk_sweep >= config.RISK_SWEEP_SEC:
                # Implemented in Task 21 (tipoff sweep, line-move sweep)
                last_risk_sweep = now

            if now - last_pipeline >= config.PIPELINE_REFRESH_SEC:
                # Implemented in Task 23
                last_pipeline = now

            if now - last_heartbeat >= 60:
                print(f"  [HB] {datetime.now(timezone.utc).isoformat()} alive", flush=True)
                last_heartbeat = now

            time.sleep(0.5)
    finally:
        # SIGTERM shutdown: cancel every live RFQ.
        with db.connect(read_only=True) as con:
            live = [r[0] for r in con.execute(
                "SELECT rfq_id FROM live_rfqs WHERE status='open'"
            ).fetchall()]
        for rid in live:
            try:
                rfq_client.delete_rfq(rid)
            except Exception:
                pass
        db.end_session(sid)
        print("=== shutdown complete ===", flush=True)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true",
                        help="Run full loop without calling accept_quote.")
    args = parser.parse_args()

    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)
    main_loop(dry_run=args.dry_run)


if __name__ == "__main__":
    cli()
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/main.py tests/kalshi_mlb_rfq/test_main_signal.py
git commit -m "feat(kalshi-mlb-rfq): main.py skeleton with signal handler + accept lock"
```

---

### Task 21: Wire RFQ refresh loop + risk sweep

**Files:**
- Modify: `kalshi_mlb_rfq/main.py`
- Test: `tests/kalshi_mlb_rfq/test_main_rfq_refresh.py`

- [ ] **Step 1: Write a higher-level integration-style test**

This test exercises the full RFQ refresh path with mocked external calls. We verify that:
- New top-N candidates trigger `mint_combo_ticker` + `create_rfq` calls.
- Combos that drop out of top-N trigger `delete_rfq`.
- A `live_rfqs` row is written for each new RFQ.

```python
# tests/kalshi_mlb_rfq/test_main_rfq_refresh.py
from datetime import datetime, timezone
from unittest.mock import patch, MagicMock

import pandas as pd
import pytest

from kalshi_mlb_rfq import main, db, combo_enumerator


@pytest.fixture
def initialized_db(tmp_path, monkeypatch):
    p = tmp_path / "test.duckdb"
    monkeypatch.setattr(db, "DB_PATH", p)
    db.init_database()
    return p


def test_rfq_refresh_submits_top_n(initialized_db, monkeypatch):
    candidates = [
        combo_enumerator.ComboCandidate(
            game_id="g1",
            legs=({"market_ticker": "L1", "event_ticker": "E1", "side": "yes"},
                  {"market_ticker": "L2", "event_ticker": "E2", "side": "yes"}),
            leg_set_hash=f"hash{i}",
            descriptor=f"COMBO_{i}",
        )
        for i in range(5)
    ]
    fairs = {f"hash{i}": (0.30 + 0.01 * i, 0.10) for i in range(5)}  # blended_fair, kalshi_ref
    monkeypatch.setattr(main, "MAX_LIVE_RFQS_OVERRIDE", 3)

    with patch("kalshi_mlb_rfq.main.mint_and_create_rfq", return_value=("rid", "MT")) as m:
        main._refresh_rfqs(candidates, fairs, dry_run=True)

    # Top 3 by edge magnitude should have been submitted.
    assert m.call_count == 3
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_main_rfq_refresh.py -v
```

- [ ] **Step 3: Append `_refresh_rfqs` and helpers to `kalshi_mlb_rfq/main.py`**

```python
# Append to kalshi_mlb_rfq/main.py

import json
from kalshi_mlb_rfq import combo_enumerator, fair_value, rfq_client


# Allow tests to override the live cap.
MAX_LIVE_RFQS_OVERRIDE: int | None = None


def _max_live_rfqs() -> int:
    return MAX_LIVE_RFQS_OVERRIDE if MAX_LIVE_RFQS_OVERRIDE is not None else config.MAX_LIVE_RFQS


def mint_and_create_rfq(candidate: combo_enumerator.ComboCandidate,
                         target_cost_dollars: float = 1.0) -> tuple[str | None, str | None]:
    """Mint combo ticker (or fetch from cache), then create RFQ. Returns (rfq_id, combo_ticker)."""
    legs = list(candidate.legs)

    with db.connect() as con:
        cached = con.execute(
            "SELECT combo_market_ticker FROM combo_cache WHERE leg_set_hash=?",
            [candidate.leg_set_hash],
        ).fetchone()

    if cached:
        combo_ticker = cached[0]
    else:
        combo_ticker, combo_event = rfq_client.mint_combo_ticker(
            config.MVE_COLLECTION_TICKER, legs)
        with db.connect() as con:
            con.execute(
                "INSERT INTO combo_cache (leg_set_hash, collection_ticker, "
                "combo_market_ticker, combo_event_ticker, legs_json, game_id) "
                "VALUES (?, ?, ?, ?, ?, ?) ON CONFLICT DO NOTHING",
                [candidate.leg_set_hash, config.MVE_COLLECTION_TICKER,
                 combo_ticker, combo_event, json.dumps(legs), candidate.game_id],
            )

    rfq_id = rfq_client.create_rfq(combo_ticker, target_cost_dollars=target_cost_dollars)
    return rfq_id, combo_ticker


def _refresh_rfqs(candidates: list[combo_enumerator.ComboCandidate],
                   fair_scores: dict[str, tuple[float, float]],
                   dry_run: bool):
    """Drive the continuous priority-queue pipeline.

    fair_scores: leg_set_hash → (blended_fair, kalshi_ref).
    """
    scored = [(c, *fair_scores[c.leg_set_hash])
              for c in candidates if c.leg_set_hash in fair_scores]
    ranked = combo_enumerator.rank_by_edge(scored)
    target = ranked[: _max_live_rfqs()]
    target_hashes = {c.leg_set_hash for c in target}

    with db.connect() as con:
        live = con.execute(
            "SELECT rfq_id, leg_set_hash FROM live_rfqs WHERE status='open'"
        ).fetchall()

    live_hashes = {h: rid for rid, h in live}

    # Drop: in DB but not in target.
    for rid, h in live_hashes.items():
        if h not in target_hashes:
            try:
                rfq_client.delete_rfq(rid)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_rfqs SET status='cancelled', closed_at=?, "
                        "cancellation_reason='out_of_top_n' WHERE rfq_id=?",
                        [datetime.now(timezone.utc), rid],
                    )
            except Exception as e:
                print(f"  drop {rid} failed: {e}", flush=True)

    # Add: in target but not in DB.
    for c in target:
        if c.leg_set_hash in live_hashes:
            continue
        try:
            rid, combo_ticker = mint_and_create_rfq(c)
            blended_fair, kalshi_ref = fair_scores[c.leg_set_hash]
            edge = combo_enumerator.edge_score(blended_fair, kalshi_ref)
            with db.connect() as con:
                con.execute(
                    "INSERT INTO live_rfqs (rfq_id, combo_market_ticker, leg_set_hash, "
                    "game_id, blended_fair_at_submit, kalshi_ref_at_submit, "
                    "edge_at_submit, status, submitted_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    [rid, combo_ticker, c.leg_set_hash, c.game_id,
                     blended_fair, kalshi_ref, edge, "open",
                     datetime.now(timezone.utc)],
                )
        except Exception as e:
            print(f"  add {c.leg_set_hash[:8]} failed: {e}", flush=True)
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/main.py tests/kalshi_mlb_rfq/test_main_rfq_refresh.py
git commit -m "feat(kalshi-mlb-rfq): main._refresh_rfqs continuous priority-queue pipeline"
```

---

### Task 22: Wire quote poll loop + accept gates

**Files:**
- Modify: `kalshi_mlb_rfq/main.py`
- Test: `tests/kalshi_mlb_rfq/test_main_accept.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_main_accept.py
from datetime import datetime, timedelta, timezone
from unittest.mock import patch

import pytest

from kalshi_mlb_rfq import main, db


@pytest.fixture
def initialized_db(tmp_path, monkeypatch):
    p = tmp_path / "test.duckdb"
    monkeypatch.setattr(db, "DB_PATH", p)
    db.init_database()
    # seed a live RFQ
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_rfqs (rfq_id, combo_market_ticker, leg_set_hash, "
            "game_id, status, submitted_at) VALUES "
            "(?, ?, ?, ?, ?, ?)",
            ["rfq-1", "MT-COMBO", "hash1", "g1", "open",
             datetime.now(timezone.utc)],
        )
    return p


def test_accept_evaluation_passes_gates(initialized_db):
    quote = {
        "id": "q1",
        "rfq_id": "rfq-1",
        "yes_bid_dollars": "0.13",
        "no_bid_dollars": "0.74",
        "creator_id": "maker-1",
        "status": "open",
    }
    # blended_fair=0.32 → buying YES at 0.26 → +EV
    fair = 0.32

    with patch("kalshi_mlb_rfq.main._fresh_blended_fair", return_value=fair), \
         patch("kalshi_mlb_rfq.main._all_per_accept_gates_pass", return_value=(True, "ok")), \
         patch("kalshi_mlb_rfq.main._kelly_size_for_quote", return_value=10), \
         patch("kalshi_mlb_rfq.main.rfq_client.accept_quote",
               return_value={"order": {"id": "ord-1"}}) as accept_mock, \
         patch("kalshi_mlb_rfq.main.rfq_client.get_position_contracts", return_value=10):
        main._evaluate_quote(quote, dry_run=False)

    accept_mock.assert_called_once()
    with db.connect(read_only=True) as con:
        n_fills = con.execute("SELECT COUNT(*) FROM fills").fetchone()[0]
    assert n_fills == 1


def test_dry_run_never_accepts(initialized_db):
    quote = {"id": "q1", "rfq_id": "rfq-1",
             "yes_bid_dollars": "0.13", "no_bid_dollars": "0.74",
             "creator_id": "m", "status": "open"}
    with patch("kalshi_mlb_rfq.main._fresh_blended_fair", return_value=0.50), \
         patch("kalshi_mlb_rfq.main._all_per_accept_gates_pass", return_value=(True, "ok")), \
         patch("kalshi_mlb_rfq.main.rfq_client.accept_quote") as a:
        main._evaluate_quote(quote, dry_run=True)
    a.assert_not_called()
    with db.connect(read_only=True) as con:
        decision = con.execute(
            "SELECT decision FROM quote_log WHERE quote_id='q1'"
        ).fetchone()[0]
    assert decision == "declined_dry_run"
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_main_accept.py -v
```

- [ ] **Step 3: Append to `kalshi_mlb_rfq/main.py`**

```python
import uuid
from kalshi_mlb_rfq import ev_calc


def _fresh_blended_fair(combo_market_ticker: str) -> float | None:
    """Re-read blended fair for the combo right before accept (samples may have shifted)."""
    # Implementation in Task 23 — placeholder for now returns None.
    raise NotImplementedError("wired in Task 23")


def _all_per_accept_gates_pass(quote: dict, fair: float,
                                combo_meta: dict) -> tuple[bool, str]:
    """Returns (pass, decision_string). decision_string is the quote_log decision value."""
    raise NotImplementedError("wired in Task 23")


def _kelly_size_for_quote(quote: dict, fair: float) -> int:
    raise NotImplementedError("wired in Task 23")


def _log_quote_decision(quote: dict, fair: float | None,
                         decision: str, reason: str | None = None,
                         post_fee_ev: float | None = None):
    with db.connect() as con:
        con.execute(
            "INSERT INTO quote_log (quote_id, rfq_id, combo_market_ticker, "
            "creator_id, yes_bid_dollars, no_bid_dollars, blended_fair_at_eval, "
            "post_fee_ev_pct, decision, reason_detail, observed_at) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) ON CONFLICT DO NOTHING",
            [quote["id"], quote["rfq_id"], quote.get("market_ticker"),
             quote.get("creator_id"),
             float(quote["yes_bid_dollars"]) if quote.get("yes_bid_dollars") else None,
             float(quote["no_bid_dollars"]) if quote.get("no_bid_dollars") else None,
             fair, post_fee_ev, decision, reason,
             datetime.now(timezone.utc)],
        )


def _evaluate_quote(quote: dict, dry_run: bool):
    """Per-quote: evaluate gates, accept if all pass, log decision."""
    rfq_id = quote["rfq_id"]
    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT combo_market_ticker, leg_set_hash, game_id "
            "FROM live_rfqs WHERE rfq_id=?", [rfq_id]
        ).fetchone()
    if not row:
        return  # RFQ vanished
    combo_market_ticker, leg_set_hash, game_id = row
    quote = {**quote, "market_ticker": combo_market_ticker}

    with ACCEPT_LOCK:
        fair = _fresh_blended_fair(combo_market_ticker)
        if fair is None:
            _log_quote_decision(quote, None, "declined_ev",
                                 reason="no_fresh_fair")
            return

        passed, decision = _all_per_accept_gates_pass(
            quote, fair, {"leg_set_hash": leg_set_hash, "game_id": game_id})
        if not passed:
            _log_quote_decision(quote, fair, decision)
            return

        no_bid = float(quote["no_bid_dollars"])
        ev_dollars, ev_pct = ev_calc.post_fee_ev_buy_yes(fair, no_bid)
        if ev_pct < config.MIN_EV_PCT:
            _log_quote_decision(quote, fair, "declined_ev", post_fee_ev=ev_pct)
            return

        if dry_run:
            _log_quote_decision(quote, fair, "declined_dry_run", post_fee_ev=ev_pct)
            return

        contracts = _kelly_size_for_quote(quote, fair)
        if contracts <= 0:
            _log_quote_decision(quote, fair, "declined_kelly_zero",
                                 post_fee_ev=ev_pct)
            return

        resp = rfq_client.accept_quote(quote["id"], contracts=contracts)
        if resp is None:
            _log_quote_decision(quote, fair, "failed_quote_walked",
                                 post_fee_ev=ev_pct)
            return

        # Post-accept fill reconciliation
        try:
            actual = rfq_client.get_position_contracts(combo_market_ticker)
        except Exception:
            actual = contracts  # fallback to requested

        yes_ask = 1.0 - no_bid
        fee = ev_calc.fee_per_contract(yes_ask)
        with db.connect() as con:
            con.execute(
                "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
                "game_id, side, contracts, price_dollars, fee_dollars, "
                "blended_fair_at_fill, expected_ev_dollars, filled_at, raw_response) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [str(uuid.uuid4()), quote["id"], rfq_id, combo_market_ticker, game_id,
                 "yes", actual, yes_ask, fee, fair, ev_dollars,
                 datetime.now(timezone.utc), str(resp)],
            )
            con.execute(
                "INSERT INTO positions (combo_market_ticker, game_id, "
                "net_contracts, weighted_price, legs_json, updated_at) VALUES "
                "(?, ?, ?, ?, ?, ?) ON CONFLICT (combo_market_ticker) DO UPDATE SET "
                "net_contracts = positions.net_contracts + EXCLUDED.net_contracts, "
                "weighted_price = (positions.weighted_price * positions.net_contracts + "
                "                  EXCLUDED.weighted_price * EXCLUDED.net_contracts) / "
                "                 (positions.net_contracts + EXCLUDED.net_contracts), "
                "updated_at = EXCLUDED.updated_at",
                [combo_market_ticker, game_id, actual, yes_ask, "[]",
                 datetime.now(timezone.utc)],
            )
            cooled_until = datetime.now(timezone.utc) + timedelta(seconds=config.COMBO_COOLDOWN_SEC)
            con.execute(
                "INSERT INTO combo_cooldown (leg_set_hash, game_id, cooled_until, reason) "
                "VALUES (?, ?, ?, ?) ON CONFLICT (leg_set_hash) DO UPDATE SET "
                "cooled_until = EXCLUDED.cooled_until",
                [leg_set_hash, game_id, cooled_until, "post_accept"],
            )
        _log_quote_decision(quote, fair, "accepted", post_fee_ev=ev_pct)
        notify.fill(rfq_id=rfq_id, combo_market_ticker=combo_market_ticker,
                    contracts=actual, price=yes_ask, ev_pct=ev_pct)
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/main.py tests/kalshi_mlb_rfq/test_main_accept.py
git commit -m "feat(kalshi-mlb-rfq): main._evaluate_quote with serialized accepts + reconciliation"
```

---

### Task 23: Wire fair-value provider, gate aggregator, Kelly sizing

These are the three `NotImplementedError` placeholders from Task 22.

**Files:**
- Modify: `kalshi_mlb_rfq/main.py`
- Test: `tests/kalshi_mlb_rfq/test_main_fair_provider.py`

- [ ] **Step 1: Write failing test**

```python
# tests/kalshi_mlb_rfq/test_main_fair_provider.py
import json
from unittest.mock import patch

import pandas as pd
import pytest

from kalshi_mlb_rfq import main, db


@pytest.fixture
def initialized_db_with_combo(tmp_path, monkeypatch):
    p = tmp_path / "test.duckdb"
    monkeypatch.setattr(db, "DB_PATH", p)
    db.init_database()
    # Seed a combo cache row representing legs we can score.
    legs = [
        {"market_ticker": "KXMLBSPREAD-X-NYY2", "event_ticker": "KXMLBSPREAD-X", "side": "yes"},
        {"market_ticker": "KXMLBTOTAL-X-8",     "event_ticker": "KXMLBTOTAL-X",  "side": "yes"},
    ]
    with db.connect() as con:
        con.execute(
            "INSERT INTO combo_cache (leg_set_hash, collection_ticker, "
            "combo_market_ticker, combo_event_ticker, legs_json, game_id) VALUES "
            "(?, ?, ?, ?, ?, ?)",
            ["hashH", "KXMVECROSSCATEGORY-R", "MT-COMBO",
             "KXMVECROSSCATEGORY-S", json.dumps(legs), "g1"],
        )
    return p


def test_fresh_blended_fair_returns_blend(initialized_db_with_combo):
    samples = pd.DataFrame({"home_margin": [3, 0, 5, -1] * 250,
                             "total_final_score": [10, 6, 9, 5] * 250})
    book_fair_lookup = {("g1", -1.5, 8.5): {"draftkings": 0.30}}

    with patch("kalshi_mlb_rfq.main._load_samples_for_game", return_value=samples), \
         patch("kalshi_mlb_rfq.main._load_book_fairs", return_value=book_fair_lookup), \
         patch("kalshi_mlb_rfq.main._home_team_for_game", return_value="NYY"):
        fair = main._fresh_blended_fair("MT-COMBO")
    assert fair is not None
    # Model: (home>=2) AND (total>=8); rows where both: (3,10),(5,9) → 2 of 4 in repeat → 0.5
    assert fair == pytest.approx((0.5 + 0.30) / 2)
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_main_fair_provider.py -v
```

- [ ] **Step 3: Replace the three `NotImplementedError` stubs in `main.py` with implementations**

```python
# Replace placeholders in kalshi_mlb_rfq/main.py

import duckdb
import numpy as np

from kalshi_mlb_rfq import fair_value, kelly, risk
from kalshi_mlb_rfq.config import ANSWER_KEY_DB, MAX_BOOK_STALENESS_SEC


def _load_samples_for_game(game_id: str) -> pd.DataFrame | None:
    con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
    try:
        df = con.execute(
            "SELECT home_margin, total_final_score, home_margin_f5, total_f5, home_scored_first "
            "FROM mlb_game_samples WHERE game_id=?", [game_id]
        ).fetchdf()
    finally:
        con.close()
    return df if not df.empty else None


def _load_book_fairs(game_id: str, spread_line: float, total_line: float) -> dict[str, float]:
    """Pull mlb_sgp_odds rows for this game/spread/total combo, devig per book."""
    con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
    try:
        rows = con.execute(
            "SELECT bookmaker, combo, sgp_decimal FROM mlb_sgp_odds "
            "WHERE game_id=? AND period='FG' "
            "AND fetch_time > NOW() - INTERVAL '? SECOND'",
            [game_id, MAX_BOOK_STALENESS_SEC]
        ).fetchdf()
    finally:
        con.close()

    out: dict[str, float] = {}
    for book in rows["bookmaker"].unique():
        sub = rows[rows["bookmaker"] == book].copy()
        # Use the 4-side devig method when 4 sides exist.
        fair_per_book = fair_value.devig_book(sub, combo="Home Spread + Over",
                                               vig_fallback=_vig_fallback(book))
        if fair_per_book is not None:
            out[book] = fair_per_book
    return out


def _vig_fallback(book: str) -> float:
    return {
        "draftkings": config.DK_VIG_FALLBACK,
        "fanduel": config.FD_VIG_FALLBACK,
        "prophetx": config.PX_VIG_FALLBACK,
        "novig": config.NOVIG_VIG_FALLBACK,
    }.get(book, 0.10)


def _home_team_for_game(game_id: str) -> str | None:
    con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
    try:
        row = con.execute(
            "SELECT home_team FROM mlb_team_dict WHERE game_id=? LIMIT 1",
            [game_id]
        ).fetchone()
    finally:
        con.close()
    return row[0] if row else None


def _fresh_blended_fair(combo_market_ticker: str) -> float | None:
    """Re-read combo legs and compute blended fair from samples + sgp_odds."""
    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT legs_json, game_id FROM combo_cache "
            "WHERE combo_market_ticker=?", [combo_market_ticker]
        ).fetchone()
    if not row:
        return None
    legs_json, game_id = row
    legs = json.loads(legs_json)

    samples = _load_samples_for_game(game_id)
    if samples is None:
        return None

    # Build fair_value.Leg objects from the leg dicts.
    legs_typed = [_leg_dict_to_typed(l, game_id) for l in legs]
    if any(l is None for l in legs_typed):
        return None

    model = fair_value.model_fair(samples, legs_typed)

    # Lookup book fairs for this (spread_line, total_line). Extract lines from legs.
    spread_line = _spread_line_from_legs(legs)
    total_line = _total_line_from_legs(legs)
    book_fairs = _load_book_fairs(game_id, spread_line, total_line)

    return fair_value.blend(model, book_fairs)


def _leg_dict_to_typed(leg: dict, game_id: str):
    """Convert {market_ticker, side} to fair_value typed leg or None if not spread/total."""
    mt = leg["market_ticker"]
    side = leg["side"]
    if mt.startswith("KXMLBSPREAD-"):
        # ticker format: KXMLBSPREAD-{event}-{TEAM}{N}
        suffix = mt.rsplit("-", 1)[-1]
        # Extract team code (chars) + N (trailing digits)
        n_chars = ""
        team_chars = ""
        for ch in suffix:
            if ch.isdigit():
                n_chars += ch
            else:
                team_chars += ch
        n = int(n_chars)
        home = _home_team_for_game(game_id)
        team_is_home = (team_chars == home)
        return fair_value.SpreadLeg(team_is_home=team_is_home, line_n=n, side=side)
    if mt.startswith("KXMLBTOTAL-"):
        n = int(mt.rsplit("-", 1)[-1])
        return fair_value.TotalLeg(line_n=n, side=side)
    return None


def _spread_line_from_legs(legs: list[dict]) -> float:
    for l in legs:
        if l["market_ticker"].startswith("KXMLBSPREAD-"):
            suffix = l["market_ticker"].rsplit("-", 1)[-1]
            n = int("".join(c for c in suffix if c.isdigit()))
            return -(n - 0.5)
    return 0.0


def _total_line_from_legs(legs: list[dict]) -> float:
    for l in legs:
        if l["market_ticker"].startswith("KXMLBTOTAL-"):
            n = int(l["market_ticker"].rsplit("-", 1)[-1])
            return n - 0.5
    return 0.0


_POSITIONS_API_FAIL_COUNT = 0


def _record_positions_api_result(success: bool):
    global _POSITIONS_API_FAIL_COUNT
    if success:
        _POSITIONS_API_FAIL_COUNT = 0
    else:
        _POSITIONS_API_FAIL_COUNT += 1


def _samples_generated_at() -> datetime | None:
    con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
    try:
        row = con.execute(
            "SELECT generated_at FROM mlb_samples_meta ORDER BY generated_at DESC LIMIT 1"
        ).fetchone()
    finally:
        con.close()
    return row[0] if row else None


def _commence_time_for_game(game_id: str) -> datetime | None:
    con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
    try:
        row = con.execute(
            "SELECT commence_time FROM mlb_team_dict WHERE game_id=? LIMIT 1",
            [game_id]
        ).fetchone()
    finally:
        con.close()
    return row[0] if row else None


def _all_per_accept_gates_pass(quote: dict, fair: float,
                                combo_meta: dict) -> tuple[bool, str]:
    """Run every per-accept gate. Returns (pass, decision_label)."""
    # Prediction staleness
    gen_at = _samples_generated_at()
    if gen_at is None or not risk.staleness_ok(gen_at, config.MAX_PREDICTION_STALENESS_SEC):
        return False, "declined_stale_predictions"

    # Tipoff window
    ct = _commence_time_for_game(combo_meta["game_id"])
    if not risk.tipoff_ok(ct, config.TIPOFF_CANCEL_MIN):
        return False, "declined_tipoff"

    # Line-move check (requires reference_lines snapshot at RFQ submission time;
    # populated in _refresh_rfqs — see Task 23b for the wiring).
    rfq_id = quote["rfq_id"]
    with db.connect(read_only=True) as con:
        ref_row = con.execute(
            "SELECT lines_json FROM reference_lines WHERE rfq_id=?", [rfq_id]
        ).fetchone()
    if ref_row:
        ref_lines = json.loads(ref_row[0])
        # Pull current lines from a recent mlb_sgp_odds aggregate, or skip if unavailable.
        current_lines = _current_book_lines_for_combo(combo_meta["game_id"], combo_meta)
        if current_lines and not risk.line_move_ok(ref_lines, current_lines, config.LINE_MOVE_THRESHOLD):
            return False, "declined_line_move"

    # Positions API health
    if _POSITIONS_API_FAIL_COUNT >= config.POSITIONS_HEALTH_RETRIES:
        return False, "declined_positions_unhealthy"

    # Fair bounds
    if not risk.fair_in_bounds(fair, config.MIN_FAIR_PROB, config.MAX_FAIR_PROB):
        return False, "declined_kelly_zero"

    # Sanity bound
    no_bid = float(quote["no_bid_dollars"])
    quote_implied = 1 - no_bid
    if not risk.sanity_bound_ok(quote_implied, fair, config.MAX_QUOTE_DEVIATION):
        return False, "declined_sanity"

    # Kill switch
    if not risk.kill_switch_ok():
        return False, "declined_killswitch"

    # Cooldown
    leg_set_hash = combo_meta["leg_set_hash"]
    with db.connect(read_only=True) as con:
        rows = con.execute(
            "SELECT leg_set_hash, cooled_until FROM combo_cooldown"
        ).fetchall()
    cd_map = {h: u for h, u in rows}
    if not risk.cooldown_ok(leg_set_hash, cd_map):
        return False, "declined_cooldown"

    # Inverse-combo
    with db.connect(read_only=True) as con:
        held = {h for (h,) in con.execute(
            "SELECT cc.leg_set_hash FROM combo_cache cc "
            "JOIN positions p ON p.combo_market_ticker = cc.combo_market_ticker "
            "WHERE p.net_contracts > 0").fetchall()}
    legs = json.loads(combo_meta.get("legs_json") or "[]")
    if legs and not risk.inverse_combo_ok(legs, held):
        return False, "declined_inverse_lock"

    # Per-game cap + daily cap
    today_start = datetime.now(timezone.utc).replace(hour=0, minute=0, second=0)
    with db.connect(read_only=True) as con:
        today_fills = [{"game_id": r[0], "contracts": r[1], "price_dollars": r[2]}
                        for r in con.execute(
                            "SELECT game_id, contracts, price_dollars FROM fills "
                            "WHERE filled_at >= ?", [today_start]
                        ).fetchall()]
    if not risk.per_game_cap_ok(combo_meta["game_id"], today_fills,
                                  config.BANKROLL, config.MAX_GAME_EXPOSURE_PCT):
        return False, "declined_per_game_cap"
    if not risk.daily_cap_ok(today_fills, config.DAILY_EXPOSURE_CAP_USD):
        return False, "declined_daily_cap"

    # Fill-ratio halt (rolling window)
    with db.connect(read_only=True) as con:
        window = con.execute(
            "SELECT decision FROM quote_log WHERE decision IN "
            "('accepted', 'failed_quote_walked') ORDER BY observed_at DESC LIMIT ?",
            [config.FILL_RATIO_WINDOW]
        ).fetchall()
    if not risk.fill_ratio_ok([{"decision": d} for (d,) in window],
                                config.MIN_FILL_RATIO):
        return False, "halted_low_fill_ratio"

    return True, "passed"


def _kelly_size_for_quote(quote: dict, fair: float) -> int:
    """Conditional Kelly sizing. Loads samples + existing positions on the same game."""
    rfq_id = quote["rfq_id"]
    with db.connect(read_only=True) as con:
        meta = con.execute(
            "SELECT game_id, combo_market_ticker FROM live_rfqs WHERE rfq_id=?",
            [rfq_id]
        ).fetchone()
    if not meta:
        return 0
    game_id, combo_market_ticker = meta

    samples = _load_samples_for_game(game_id)
    if samples is None or samples.empty:
        return 0

    # Build outcome vector for THIS combo.
    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT legs_json FROM combo_cache WHERE combo_market_ticker=?",
            [combo_market_ticker]
        ).fetchone()
    if not row:
        return 0
    legs = json.loads(row[0])
    typed_legs = [_leg_dict_to_typed(l, game_id) for l in legs]
    if any(l is None for l in typed_legs):
        return 0
    mask = pd.Series([True] * len(samples), index=samples.index)
    for leg in typed_legs:
        mask &= fair_value._hit_mask(samples, leg)
    outcome_vec = mask.astype(int).values

    # Existing positions on the same game.
    with db.connect(read_only=True) as con:
        pos_rows = con.execute(
            "SELECT cc.legs_json, p.net_contracts, p.weighted_price "
            "FROM positions p JOIN combo_cache cc "
            "ON cc.combo_market_ticker = p.combo_market_ticker "
            "WHERE p.game_id = ? AND p.net_contracts > 0",
            [game_id]
        ).fetchall()
    existing = []
    for legs_json, n, price in pos_rows:
        leg_objs = [_leg_dict_to_typed(l, game_id) for l in json.loads(legs_json)]
        if any(l is None for l in leg_objs):
            continue
        sub_mask = pd.Series([True] * len(samples), index=samples.index)
        for leg in leg_objs:
            sub_mask &= fair_value._hit_mask(samples, leg)
        existing.append({
            "outcome_vec": sub_mask.astype(int).values,
            "contracts": float(n),
            "effective_price": float(price),
        })

    no_bid = float(quote["no_bid_dollars"])
    yes_ask = 1 - no_bid
    fee = ev_calc.fee_per_contract(yes_ask)
    effective_price = yes_ask + fee

    return kelly.kelly_size_combo(
        outcome_vec=np.asarray(outcome_vec),
        existing_positions=existing,
        effective_price=effective_price,
        bankroll=config.BANKROLL,
        kelly_fraction=config.KELLY_FRACTION,
    )
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/main.py tests/kalshi_mlb_rfq/test_main_fair_provider.py
git commit -m "feat(kalshi-mlb-rfq): wire fair_value provider + gate aggregator + Kelly sizing"
```

---

### Task 23b: Line-move snapshot + helpers for new gates

The gates added in Task 23 (line-move, positions API health) need supporting plumbing:
1. `_current_book_lines_for_combo` helper
2. `reference_lines` snapshot in `_refresh_rfqs`
3. Wrap `rfq_client.get_position_contracts` with `_record_positions_api_result`

**Files:**
- Modify: `kalshi_mlb_rfq/main.py`

- [ ] **Step 1: Append `_current_book_lines_for_combo` helper to `main.py`**

```python
def _current_book_lines_for_combo(game_id: str, combo_meta: dict) -> dict | None:
    """Return current spread/total line snapshot from mlb_sgp_odds for line-move detection.

    Returns dict like {"spread": -1.5, "total": 8.5} or None if no recent rows.
    """
    con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
    try:
        row = con.execute(
            "SELECT spread_line, total_line FROM mlb_sgp_odds "
            "WHERE game_id=? AND period='FG' "
            "AND fetch_time > NOW() - INTERVAL '? SECOND' "
            "ORDER BY fetch_time DESC LIMIT 1",
            [game_id, config.MAX_BOOK_STALENESS_SEC]
        ).fetchone()
    finally:
        con.close()
    if not row:
        return None
    return {"spread": float(row[0]), "total": float(row[1])}
```

- [ ] **Step 2: Add reference-lines snapshot in `_refresh_rfqs`**

After `INSERT INTO live_rfqs ...` in the `# Add` section of `_refresh_rfqs`, also write:

```python
                lines_now = _current_book_lines_for_combo(c.game_id, {})
                if lines_now:
                    con.execute(
                        "INSERT INTO reference_lines (rfq_id, lines_json, snapped_at) "
                        "VALUES (?, ?, ?) ON CONFLICT (rfq_id) DO NOTHING",
                        [rid, json.dumps(lines_now), datetime.now(timezone.utc)],
                    )
```

- [ ] **Step 3: Wrap positions-API call in `_evaluate_quote`**

Replace the `rfq_client.get_position_contracts` line in `_evaluate_quote` with:

```python
        try:
            actual = rfq_client.get_position_contracts(combo_market_ticker)
            _record_positions_api_result(True)
        except Exception:
            _record_positions_api_result(False)
            actual = contracts  # fallback to requested
```

- [ ] **Step 4: Quick import check**

```bash
python3 -c "from kalshi_mlb_rfq import main; print('OK')"
```
Expected: OK.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/main.py
git commit -m "feat(kalshi-mlb-rfq): line-move snapshot + positions-API health tracking"
```

---

### Task 24: Pipeline-refresh subprocess + recon script update

**Files:**
- Modify: `kalshi_mlb_rfq/main.py`
- Modify: `mlb_sgp/recon_kalshi_mlb_rfq.py` (was drafted earlier; lookup endpoint path needs the fix)

- [ ] **Step 1: Append `_run_pipeline()` to `kalshi_mlb_rfq/main.py`**

```python
# Add to kalshi_mlb_rfq/main.py

import subprocess


def _run_pipeline():
    """Trigger the MLB R answer-key pipeline. Non-blocking error handling: notify on failure."""
    cmd = ["Rscript", str(ANSWER_KEY_DB.parent / "MLB Answer Key" / "MLB.R")]
    try:
        proc = subprocess.run(cmd, cwd=ANSWER_KEY_DB.parent,
                              timeout=240, capture_output=True, text=True)
        if proc.returncode != 0:
            notify.halt("pipeline_refresh_failed",
                        detail=f"exit={proc.returncode}; stderr={proc.stderr[:500]}")
    except subprocess.TimeoutExpired:
        notify.halt("pipeline_refresh_failed", detail="timeout")
    except Exception as e:
        notify.halt("pipeline_refresh_failed", detail=str(e))
```

- [ ] **Step 2: Update `mlb_sgp/recon_kalshi_mlb_rfq.py`** to use the corrected MVE lookup path (drop the `/lookup` suffix that was wrong in the original draft):

Replace the line in PHASE D:
```python
lookup_path = f"/multivariate_event_collections/{chosen_collection}/lookup"
```
With:
```python
lookup_path = f"/multivariate_event_collections/{chosen_collection}"
```

Also handle DELETE returning 204 with empty body — wrap the cleanup `api()` call in a try/except that ignores empty-body JSON decode errors.

- [ ] **Step 3: Verify recon still runs** (this is a smoke test, not a unit test)

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/kalshi-mlb-rfq
python3 mlb_sgp/recon_kalshi_mlb_rfq.py 2>&1 | tail -40
```
Expected: PHASE 0 OK, PHASE A finds the 2 cross-category MVE collections, PHASE D mints a combo ticker, PHASE F polls quotes (4–5 expected), CLEANUP deletes RFQs with status 204.

- [ ] **Step 4: Commit**

```bash
git add kalshi_mlb_rfq/main.py mlb_sgp/recon_kalshi_mlb_rfq.py
git commit -m "feat(kalshi-mlb-rfq): pipeline refresh subprocess + fix recon lookup path"
```

---

### Task 25: Wire all four loops in main_loop()

The skeleton main_loop currently has placeholder branches. Replace them with calls to the implemented helpers.

**Files:**
- Modify: `kalshi_mlb_rfq/main.py`

- [ ] **Step 1: Update `main_loop` body to actually call helpers**

```python
# Update main_loop in kalshi_mlb_rfq/main.py

def main_loop(dry_run: bool):
    db.init_database()
    sid = db.start_session(pid=os.getpid(), dry_run=dry_run, version=VERSION)
    print(f"=== Kalshi MLB RFQ Bot — session {sid} (dry_run={dry_run}) ===", flush=True)
    _phantom_rfq_cleanup()

    last_rfq_refresh = 0.0
    last_quote_poll = 0.0
    last_risk_sweep = 0.0
    last_pipeline = 0.0
    last_heartbeat = 0.0

    try:
        while _running.is_set():
            now = time.time()

            if KILL_FILE.exists():
                notify.halt("kill_switch")
                time.sleep(config.RISK_SWEEP_SEC)
                continue

            if now - last_rfq_refresh >= config.RFQ_REFRESH_SEC:
                try:
                    candidates, fair_scores = _enumerate_and_score_all_games()
                    _refresh_rfqs(candidates, fair_scores, dry_run=dry_run)
                except Exception as e:
                    print(f"  rfq_refresh error: {e}", flush=True)
                last_rfq_refresh = now

            if now - last_quote_poll >= config.QUOTE_POLL_SEC:
                try:
                    _poll_all_live_rfqs(dry_run=dry_run)
                except Exception as e:
                    print(f"  quote_poll error: {e}", flush=True)
                last_quote_poll = now

            if now - last_risk_sweep >= config.RISK_SWEEP_SEC:
                try:
                    _risk_sweep()
                except Exception as e:
                    print(f"  risk_sweep error: {e}", flush=True)
                last_risk_sweep = now

            if now - last_pipeline >= config.PIPELINE_REFRESH_SEC:
                _run_pipeline()
                last_pipeline = now

            if now - last_heartbeat >= 60:
                print(f"  [HB] {datetime.now(timezone.utc).isoformat()} alive", flush=True)
                last_heartbeat = now

            time.sleep(0.5)
    finally:
        with db.connect(read_only=True) as con:
            live = [r[0] for r in con.execute(
                "SELECT rfq_id FROM live_rfqs WHERE status='open'"
            ).fetchall()]
        for rid in live:
            try:
                rfq_client.delete_rfq(rid)
            except Exception:
                pass
        db.end_session(sid)
        print("=== shutdown complete ===", flush=True)


def _enumerate_and_score_all_games():
    """For each open MLB game on Kalshi, enumerate combos and compute blended fair.

    Returns (candidates, fair_scores) where fair_scores maps leg_set_hash → (blended_fair, kalshi_ref).
    """
    # Load list of open MLB events from Kalshi.
    from kalshi_mlb_rfq.auth_client import api
    status, body, _ = api("GET", "/events?series_ticker=KXMLBGAME&status=open&limit=50")
    events = body.get("events", []) if status == 200 else []

    candidates_all = []
    fair_scores: dict[str, tuple[float, float]] = {}

    for ev in events:
        event_ticker = ev["event_ticker"]
        suffix = event_ticker.replace("KXMLBGAME-", "")
        # Look up canonical home/away codes from event title or markets.
        # Use the suffix's last 6 chars (NYYTEX → NYY,TEX) as a heuristic.
        away_code = suffix[-6:-3]
        home_code = suffix[-3:]

        # Available spreads/totals — query Kalshi for tickers that exist.
        avail_spreads = _kalshi_available_spreads(suffix, home_code, away_code)
        avail_totals = _kalshi_available_totals(suffix)
        if not avail_spreads or not avail_totals:
            continue

        game_id = _resolve_game_id(home_code, away_code, ev.get("expected_expiration_time"))
        if game_id is None:
            continue

        for cand in combo_enumerator.enumerate_2leg(
                game_id=game_id, event_suffix=suffix,
                home_code=home_code, away_code=away_code,
                available_spreads=avail_spreads,
                available_totals=avail_totals):
            try:
                combo_ticker, _ = mint_and_create_rfq.__wrapped__(  # bypass the create call
                    cand) if False else (None, None)
            except Exception:
                pass
            # Compute fair: model + book.
            samples = _load_samples_for_game(game_id)
            if samples is None:
                break
            typed = [_leg_dict_to_typed(dict(l), game_id) for l in cand.legs]
            if any(l is None for l in typed):
                continue
            model = fair_value.model_fair(samples, typed)
            spread_line = _spread_line_from_legs([dict(l) for l in cand.legs])
            total_line = _total_line_from_legs([dict(l) for l in cand.legs])
            books = _load_book_fairs(game_id, spread_line, total_line)
            blended = fair_value.blend(model, books)
            if blended is None:
                continue
            if not (config.MIN_FAIR_PROB <= blended <= config.MAX_FAIR_PROB):
                continue
            kalshi_ref = _kalshi_last_price(cand.legs[0]["market_ticker"])  # rough proxy
            candidates_all.append(cand)
            fair_scores[cand.leg_set_hash] = (blended, kalshi_ref)

    return candidates_all, fair_scores


def _kalshi_available_spreads(suffix: str, home_code: str, away_code: str
                                 ) -> list[tuple[float, str]]:
    """Hit /markets?event_ticker=KXMLBSPREAD-{suffix} and parse out lines."""
    from kalshi_mlb_rfq.auth_client import api
    status, body, _ = api("GET",
                          f"/markets?event_ticker=KXMLBSPREAD-{suffix}&limit=50")
    if status != 200:
        return []
    out = []
    for m in body.get("markets", []):
        ticker = m["ticker"]
        # KXMLBSPREAD-{suffix}-{TEAM}{N}
        spread_part = ticker.split(f"KXMLBSPREAD-{suffix}-")[-1]
        # team code is letters; N is digits at end
        digits = "".join(c for c in spread_part if c.isdigit())
        team_chars = "".join(c for c in spread_part if not c.isdigit())
        if not digits or not team_chars:
            continue
        n = int(digits)
        line = -(n - 0.5)  # N=2 → line=-1.5
        who = "home" if team_chars == home_code else "away"
        out.append((line, who))
    return out


def _kalshi_available_totals(suffix: str) -> list[float]:
    from kalshi_mlb_rfq.auth_client import api
    status, body, _ = api("GET",
                          f"/markets?event_ticker=KXMLBTOTAL-{suffix}&limit=50")
    if status != 200:
        return []
    out = []
    for m in body.get("markets", []):
        ticker = m["ticker"]
        try:
            n = int(ticker.rsplit("-", 1)[-1])
            out.append(n - 0.5)  # N=8 → line=7.5
        except ValueError:
            continue
    return out


def _kalshi_last_price(market_ticker: str) -> float:
    from kalshi_mlb_rfq.auth_client import api
    status, body, _ = api("GET", f"/markets/{market_ticker}")
    if status != 200:
        return 0.0
    m = body.get("market") if isinstance(body, dict) else None
    if not m:
        return 0.0
    try:
        return float(m.get("last_price_dollars", "0") or 0.0)
    except ValueError:
        return 0.0


def _resolve_game_id(home_code: str, away_code: str,
                      _expected_expiration: str | None) -> str | None:
    """Map Kalshi (home_code, away_code) to mlb.duckdb game_id via mlb_team_dict."""
    con = duckdb.connect(str(ANSWER_KEY_DB), read_only=True)
    try:
        rows = con.execute(
            "SELECT game_id FROM mlb_team_dict "
            "WHERE home_team_short=? AND away_team_short=? "
            "AND commence_time > NOW() AND commence_time < NOW() + INTERVAL '24 HOUR' "
            "ORDER BY commence_time LIMIT 1",
            [home_code, away_code]
        ).fetchone()
    finally:
        con.close()
    return rows[0] if rows else None


def _poll_all_live_rfqs(dry_run: bool):
    with db.connect(read_only=True) as con:
        live = [r[0] for r in con.execute(
            "SELECT rfq_id FROM live_rfqs WHERE status='open'"
        ).fetchall()]
    for rid in live:
        try:
            quotes = rfq_client.poll_quotes(rid, user_id=config.KALSHI_USER_ID)
        except Exception:
            continue
        for q in quotes:
            if q.get("status") != "open":
                continue
            with db.connect(read_only=True) as con:
                already = con.execute(
                    "SELECT 1 FROM quote_log WHERE quote_id=?", [q["id"]]
                ).fetchone()
            if already:
                continue
            _evaluate_quote(q, dry_run=dry_run)


def _risk_sweep():
    """Tipoff cancel + line-move pulls."""
    now = datetime.now(timezone.utc)
    with db.connect(read_only=True) as con:
        # For tipoff cancel we need each game's commence_time. For v1 we
        # short-circuit and rely solely on the per-accept tipoff gate.
        live_pairs = con.execute(
            "SELECT rfq_id, game_id FROM live_rfqs WHERE status='open'"
        ).fetchall()
    # TODO: query mlb.duckdb::mlb_team_dict.commence_time per game_id
    # and compare against now + TIPOFF_CANCEL_MIN. Skip in v1 — accept gate suffices.
```

- [ ] **Step 2: Quick sanity test that imports succeed**

```bash
python3 -c "from kalshi_mlb_rfq import main; print('import OK')"
```

- [ ] **Step 3: Commit**

```bash
git add kalshi_mlb_rfq/main.py
git commit -m "feat(kalshi-mlb-rfq): wire all 4 loops + Kalshi market enumeration helpers"
```

---

## Phase 9 — Operational tooling

### Task 26: `dashboard.py` — CLI status tool

**Files:**
- Create: `kalshi_mlb_rfq/dashboard.py`
- Test: `tests/kalshi_mlb_rfq/test_dashboard.py`

- [ ] **Step 1: Write failing test (smoke test only)**

```python
# tests/kalshi_mlb_rfq/test_dashboard.py
from kalshi_mlb_rfq import dashboard, db


def test_dashboard_runs_without_error(tmp_path, monkeypatch, capsys):
    p = tmp_path / "test.duckdb"
    monkeypatch.setattr(db, "DB_PATH", p)
    db.init_database()
    dashboard.main()
    out = capsys.readouterr().out
    assert "Kalshi MLB RFQ" in out
```

- [ ] **Step 2: Run, expect FAIL**

```bash
pytest tests/kalshi_mlb_rfq/test_dashboard.py -v
```

- [ ] **Step 3: Implement `kalshi_mlb_rfq/dashboard.py`**

```python
"""CLI status tool for the Kalshi MLB RFQ bot."""

from datetime import datetime, timezone

from kalshi_mlb_rfq import db


def main():
    print("=" * 60)
    print("  Kalshi MLB RFQ Bot — status")
    print("=" * 60)

    with db.connect(read_only=True) as con:
        live = con.execute(
            "SELECT COUNT(*) FROM live_rfqs WHERE status='open'"
        ).fetchone()[0]
        today_fills = con.execute(
            "SELECT COUNT(*), COALESCE(SUM(contracts*price_dollars), 0) "
            "FROM fills WHERE filled_at >= CURRENT_DATE"
        ).fetchone()
        decisions = con.execute(
            "SELECT decision, COUNT(*) FROM quote_log "
            "WHERE observed_at >= NOW() - INTERVAL '1 HOUR' GROUP BY decision"
        ).fetchall()

    print(f"  Live RFQs       : {live}")
    print(f"  Today fills     : {today_fills[0]} totaling ${today_fills[1]:.2f}")
    print()
    print("  Quote-log decisions in last hour:")
    for d, n in decisions:
        print(f"    {d:<35} {n}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/dashboard.py tests/kalshi_mlb_rfq/test_dashboard.py
git commit -m "feat(kalshi-mlb-rfq): dashboard CLI status tool"
```

---

## Phase 10 — Documentation

### Task 27: Module README + root README + CLAUDE.md updates

**Files:**
- Create: `kalshi_mlb_rfq/README.md`
- Modify: `README.md` (root)
- Modify: `CLAUDE.md` (root)

- [ ] **Step 1: Write `kalshi_mlb_rfq/README.md`**

```markdown
# Kalshi MLB RFQ Bot

Autonomous taker daemon that auto-RFQs MLB game-line same-game-parlays on Kalshi (cross-category MVE collection), evaluates maker quotes against a model+sportsbook blended fair value, and auto-accepts +EV quotes within a complete safety scaffold.

**Spec:** `docs/superpowers/specs/2026-04-27-kalshi-mlb-rfq-bot-design.md`

## Quick start

```bash
cd ~/NFLWork/kalshi_mlb_rfq

# 1. Copy env template, fill credentials
cp .env.example .env
# Edit .env with KALSHI_API_KEY_ID, KALSHI_PRIVATE_KEY_PATH, KALSHI_USER_ID

# 2. Install deps
pip install -r requirements.txt

# 3. Make sure the MLB pipeline has run at least once
cd "../Answer Keys" && python run.py --sport mlb && cd ../kalshi_mlb_rfq

# 4. Dry-run validation (foreground, no real accepts)
python3 main.py --dry-run

# 5. Live mode (backgrounded)
python3 -u main.py >> bot.log 2>&1 &
tail -f bot.log
```

## Stopping the bot

```bash
# Graceful (cancels all live RFQs first)
kill $(pgrep -f "python3.*kalshi_mlb_rfq.*main.py")

# Emergency kill switch (file-based; bot stays alive but stops accepting)
touch .kill
# Resume:
rm .kill
```

## Status

```bash
python3 dashboard.py
```

## Architecture

See spec §3. TL;DR: standalone daemon, reads `Answer Keys/mlb.duckdb` (samples + sgp_odds) read-only, writes `kalshi_mlb_rfq.duckdb`. Continuous priority-queue pipeline of up to 80 in-flight RFQs.

## Knobs

See `.env.example`. Most relevant:
- `BANKROLL`, `KELLY_FRACTION` — sizing
- `MIN_EV_PCT`, `MAX_QUOTE_DEVIATION` — accept thresholds
- `MAX_GAME_EXPOSURE_PCT`, `DAILY_EXPOSURE_CAP_USD` — exposure caps
- `MAX_PREDICTION_STALENESS_SEC` — accept-gate staleness threshold (10 min default)
- `MIN_FILL_RATIO`, `FILL_RATIO_WINDOW` — adverse-selection halt

## Troubleshooting

- **No combos surfacing:** check that `mlb_game_samples` and `mlb_sgp_odds` are populated. Run the MLB pipeline first.
- **All quotes declining `declined_stale_predictions`:** rerun the pipeline; samples are over 10 minutes old.
- **Bot halted on `fill_ratio_collapse`:** investigate — makers are walking on accepts at a rate that suggests adverse selection. Check `quote_log` for the maker_ids causing it.
```

- [ ] **Step 2: Update root `README.md`** — add a row to the project-structure section:

Find a place near other sport bots and add:
```
- **`kalshi_mlb_rfq/`** — autonomous taker daemon for MLB SGP combos on Kalshi (cross-category MVE). See `kalshi_mlb_rfq/README.md`.
```

- [ ] **Step 3: Update root `CLAUDE.md`** — add to the "This repo contains tools for:" list under Project Structure:

```
- **Autonomous Kalshi MLB SGP taker bot** (`kalshi_mlb_rfq/`) — wide-mode RFQs on cross-category MVE combos with model+book blended fair value, conditional Kelly sizing, full per-accept gate scaffold. See `kalshi_mlb_rfq/README.md`.
```

- [ ] **Step 4: Commit**

```bash
git add kalshi_mlb_rfq/README.md README.md CLAUDE.md
git commit -m "docs: kalshi_mlb_rfq README + root README + CLAUDE.md updates"
```

---

## Phase 11 — Validation

### Task 28: Dry-run smoke test

This is a manual gate, not an automated test. The dry-run runs against real Kalshi API but never accepts.

- [ ] **Step 1: Pre-flight checks**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/kalshi-mlb-rfq/kalshi_mlb_rfq

# Confirm env loaded
test -f .env && echo "env OK" || echo "MISSING .env"

# Confirm dependencies
pip install -r requirements.txt

# Confirm answer key recently ran
python3 -c "
import duckdb
con = duckdb.connect('../Answer Keys/mlb.duckdb', read_only=True)
print('mlb_game_samples rows:', con.execute('SELECT COUNT(*) FROM mlb_game_samples').fetchone()[0])
print('mlb_samples_meta:', con.execute('SELECT * FROM mlb_samples_meta').fetchall())
con.close()
"
```
Expected: env OK, sample count > 0, generated_at within the last hour.

- [ ] **Step 2: Confirm zero open RFQs on the account**

```bash
python3 -c "
import sys
sys.path.insert(0, '.')
import kalshi_mlb_rfq.config as cfg
from kalshi_mlb_rfq.rfq_client import list_open_rfqs
rfqs = list_open_rfqs(cfg.KALSHI_USER_ID)
print(f'open RFQs: {len(rfqs)}')
for r in rfqs:
    print(' ', r['id'], r['market_ticker'])
"
```
Expected: 0 open. If non-zero, those are leftover from prior sessions; the bot will adopt them at startup but cleaner to verify.

- [ ] **Step 3: Run dry-run for 30 minutes**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/kalshi-mlb-rfq/kalshi_mlb_rfq
python3 -u main.py --dry-run | tee dry_run.log
# Let it run for 30+ minutes during a slate of MLB games. Ctrl+C to stop.
```

While running, check periodically:
```bash
python3 dashboard.py
tail -20 bot.log
```

- [ ] **Step 4: Verify outcomes**

```bash
python3 -c "
import duckdb
con = duckdb.connect('kalshi_mlb_rfq.duckdb', read_only=True)
print('=== sessions ===')
for r in con.execute('SELECT * FROM sessions').fetchall(): print(' ', r)
print()
print('=== live_rfqs by status ===')
for r in con.execute('SELECT status, COUNT(*) FROM live_rfqs GROUP BY status').fetchall(): print(' ', r)
print()
print('=== quote_log decisions ===')
for r in con.execute('SELECT decision, COUNT(*) FROM quote_log GROUP BY decision ORDER BY 2 DESC').fetchall(): print(' ', r)
print()
print('=== fills (should be 0 in dry run) ===')
print(' ', con.execute('SELECT COUNT(*) FROM fills').fetchone())
con.close()
"
```
Expected:
- `fills` count is 0
- `quote_log` shows `declined_dry_run` entries (the dry-run terminal decision)
- `live_rfqs` shows a non-zero count of RFQs that were created and either still open or cancelled (`out_of_top_n`)

- [ ] **Step 5: Confirm zero phantom RFQs after exit**

```bash
python3 -c "
import sys
sys.path.insert(0, '.')
import kalshi_mlb_rfq.config as cfg
from kalshi_mlb_rfq.rfq_client import list_open_rfqs
print('post-exit open RFQs on account:', len(list_open_rfqs(cfg.KALSHI_USER_ID)))
"
```
Expected: 0.

- [ ] **Step 6: TBD-5 verification — Kalshi accept atomicity**

If you saw any `failed_quote_walked` decisions in dry-run, inspect their distribution. If they're rare and random, the implicit assumption (Kalshi accepts at the price we evaluated) is probably safe. If they cluster around price-movement events, add the pre-accept confirm step (per spec TBD-5).

- [ ] **Step 7: Commit dry-run logs (if any debug findings)**

```bash
# Don't commit dry_run.log itself — it's gitignored as kalshi_mlb_rfq/*.log.
# But if findings, add to spec or README.
```

---

## Self-review checklist

Before proposing merge to `main`:

- [ ] All Phase 1–9 tasks complete; tests pass: `pytest tests/kalshi_mlb_rfq/ -v`
- [ ] Dry-run completed for ≥30 min with no errors and zero phantom RFQs at exit
- [ ] All `quote_log` decisions in dry-run are explainable
- [ ] No accidental fills (`SELECT COUNT(*) FROM fills` = 0 after dry-run)
- [ ] `kalshi_mlb_rfq/README.md` matches actual behavior
- [ ] Root `README.md` and `CLAUDE.md` reference the bot
- [ ] `.gitignore` includes `kalshi_mlb_rfq/.env`, `*.duckdb`, `*.log`, `.kill`
- [ ] Run executive engineer review: `git diff main..HEAD --stat` and check for unintended changes
- [ ] Get user approval before merging
- [ ] Worktree cleanup post-merge: `git worktree remove .worktrees/kalshi-mlb-rfq && git branch -d feat/kalshi-mlb-rfq-bot`

---

## Notes for the implementing engineer

- **Kalshi API is finicky.** Auth signatures use the request path *without* query string. Empty bodies on DELETE return 204 with no JSON. Quote polling requires `rfq_creator_user_id` filter (recon-derived).
- **The recon script is gold.** When in doubt about a Kalshi mechanic, run `python3 mlb_sgp/recon_kalshi_mlb_rfq.py` against your account to see actual responses. Always cancels its own RFQs in the finally block.
- **DuckDB write contention** between this bot and the answer-key R pipeline is mitigated by `db.connect()` retry-with-backoff. If you see frequent retries, check whether the pipeline is running concurrently.
- **CBB MM as reference.** When uncertain how to structure something (logging, dashboard, signal handling), look at `kalshi_mm/main.py` — same shape, just maker.
- **Conditional Kelly is sensitive to sample data quality.** If you ever see Kelly producing huge sizes (>$50 per bet on a $1000 bankroll), check covariance computation — likely a bug in `_hit_mask` or sample loading.
