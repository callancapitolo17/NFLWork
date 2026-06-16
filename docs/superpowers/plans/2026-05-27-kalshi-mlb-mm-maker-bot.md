# Kalshi MLB RFQ Maker Bot — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an independent Kalshi bot that listens for others' RFQs, prices the 2-leg spread×total MLB combos we already model, and provides two-sided quotes at a fixed 5% ROI margin — instrumented to measure realized edge on a small bankroll.

**Architecture:** New `kalshi_mlb_mm/` package, a standalone REST-polling daemon, importing shared math from a new `kalshi_common/` package (extracted from the taker via backward-compatible shims). Discovery/confirm go through `RFQSource`/`QuoteGateway` interfaces so a WebSocket transport can drop in later. Fixed-margin pricing (`bid = side_prob/1.05 − maker_fee`); defenses are margin + freshness gate + book-move circuit breaker + last-look backstop.

**Tech Stack:** Python 3.14 (stdlib `urllib`, `threading`), DuckDB, pandas/scipy (via `fair_value`), pytest.

**Spec:** `docs/superpowers/specs/2026-05-26-kalshi-mlb-mm-design.md`

**Worktree:** `worktree-kalshi-mlb-mm-maker-bot` (already active).

---

## File structure

**New shared package `kalshi_common/`:**
- `__init__.py`
- `fair_value.py` — moved from taker (pure)
- `ev_calc.py` — moved from taker + new `maker_fee_per_contract`
- `auth_client.py` — moved from taker, parameterized via `configure()` (removes dependency on the taker's `config`)
- `sgp_runner.py` — moved from taker (already parameterized)
- `leg_types.py` — leg-typing helpers extracted from taker `main.py`

**Taker `kalshi_mlb_rfq/` (shims only — no logic change):**
- `fair_value.py`, `ev_calc.py`, `auth_client.py`, `sgp_runner.py` become one-line re-exports of `kalshi_common`
- `main.py` — replace the embedded leg-typing helper defs with an import from `kalshi_common.leg_types`

**New maker package `kalshi_mlb_mm/`:**
- `__init__.py`, `config.py`, `db.py`, `notify.py`
- `rfq_source.py` — `RFQSource` + `RestRFQSource`
- `quote_gateway.py` — `QuoteGateway` + `RestQuoteGateway`
- `scope.py` — in-scope detection from `mve_selected_legs`
- `pricing.py` — 5%-ROI quote math
- `risk.py` — freshness gate, circuit breaker, caps, last-look
- `main.py` — the loops + wiring
- `.env.example`, `requirements.txt`, `README.md`
- `tests/` — unit tests per module

---

## Task 1: Extract `kalshi_common/` (shim pattern, taker untouched behaviorally)

**Files:**
- Create: `kalshi_common/__init__.py`, `kalshi_common/fair_value.py`, `kalshi_common/ev_calc.py`, `kalshi_common/auth_client.py`, `kalshi_common/sgp_runner.py`, `kalshi_common/leg_types.py`
- Modify (shim): `kalshi_mlb_rfq/fair_value.py`, `kalshi_mlb_rfq/ev_calc.py`, `kalshi_mlb_rfq/auth_client.py`, `kalshi_mlb_rfq/sgp_runner.py`
- Modify: `kalshi_mlb_rfq/main.py` (leg-typing import), `kalshi_mlb_rfq/config.py` (no change needed — see step)
- Test: existing `kalshi_mlb_rfq/tests/test_probit_devig.py`, `test_sgp_runner.py`, `test_main_cache.py`, `test_config.py`

- [ ] **Step 1: Create the package and move `fair_value.py` verbatim**

```bash
mkdir -p kalshi_common
touch kalshi_common/__init__.py
git mv kalshi_mlb_rfq/fair_value.py kalshi_common/fair_value.py
```
`kalshi_common/fair_value.py` is unchanged (pure module; imports only stdlib + pandas + scipy).

- [ ] **Step 2: Re-create the taker shim for `fair_value`**

Create `kalshi_mlb_rfq/fair_value.py` with exactly:
```python
"""Backward-compat shim. Logic lives in kalshi_common.fair_value."""
from kalshi_common.fair_value import *  # noqa: F401,F403
from kalshi_common.fair_value import (  # explicit re-export of names used by callers
    SpreadLeg, TotalLeg, Leg, model_fair, devig_book, blend, _hit_mask, _probit_devig_n,
)
```

- [ ] **Step 3: Move `ev_calc.py` and add `maker_fee_per_contract`**

```bash
git mv kalshi_mlb_rfq/ev_calc.py kalshi_common/ev_calc.py
```
Append to `kalshi_common/ev_calc.py`:
```python
def maker_fee_per_contract(price_dollars: float) -> float:
    """Maker (resting-order) fee = 25% of the taker fee, on the same quadratic base.

    Kalshi charges resting/maker fills 25% of the taker fee. We compute it as
    0.25 × the (already cent-rounded) taker fee — conservative (never under-charges
    ourselves at quote time). VERIFY the exact charge on the first real fill against
    docs.kalshi.com fee_rounding; adjust the rounding here if it differs.
    """
    return 0.25 * fee_per_contract(price_dollars)
```
Create shim `kalshi_mlb_rfq/ev_calc.py`:
```python
"""Backward-compat shim. Logic lives in kalshi_common.ev_calc."""
from kalshi_common.ev_calc import *  # noqa: F401,F403
from kalshi_common.ev_calc import (
    fee_per_contract, post_fee_ev_buy_yes, post_fee_ev_buy_no, maker_fee_per_contract,
)
```

- [ ] **Step 4: Move and parameterize `auth_client.py`**

```bash
git mv kalshi_mlb_rfq/auth_client.py kalshi_common/auth_client.py
```
Edit `kalshi_common/auth_client.py` to remove the `from kalshi_mlb_rfq.config import ...` line and replace module-level credential references with configurable globals + a `configure()` injector. The full new file:
```python
"""Authenticated HTTP wrapper around Kalshi /trade-api/v2.

Config-agnostic: callers inject credentials via configure(). Defaults read from
env so ad-hoc scripts work without a configure() call.
"""
import json
import os
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

_API_KEY_ID = os.environ.get("KALSHI_API_KEY_ID")
_PRIVATE_KEY_PATH = os.environ.get("KALSHI_PRIVATE_KEY_PATH")
_BASE_URL = os.environ.get("KALSHI_BASE_URL", "https://api.elections.kalshi.com/trade-api/v2")
_KALSHI_DRAFT_DIR = None
_sign_request = None


def configure(api_key_id, private_key_path, base_url, project_root):
    """Inject per-bot credentials + the repo root used to locate kalshi_draft/auth.py."""
    global _API_KEY_ID, _PRIVATE_KEY_PATH, _BASE_URL, _KALSHI_DRAFT_DIR, _sign_request
    _API_KEY_ID = api_key_id
    _PRIVATE_KEY_PATH = private_key_path
    _BASE_URL = base_url
    candidates = [Path(project_root) / "kalshi_draft",
                  Path(project_root).parent.parent / "kalshi_draft"]
    for c in candidates:
        if (c / "auth.py").exists():
            if str(c) not in sys.path:
                sys.path.insert(0, str(c))
            break
    from auth import sign_request as sr
    _sign_request = sr


def _sign(method: str, path_no_query: str):
    if _sign_request is None:
        raise RuntimeError("auth_client.configure() not called")
    ts = str(int(datetime.now(timezone.utc).timestamp() * 1000))
    sig = _sign_request(_PRIVATE_KEY_PATH, ts, method, f"/trade-api/v2{path_no_query}")
    return sig, ts


def api(method: str, path: str, body: dict | None = None, timeout: int = 30):
    path_no_query = path.split("?", 1)[0]
    sig, ts = _sign(method, path_no_query)
    url = f"{_BASE_URL}{path}"
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("KALSHI-ACCESS-KEY", _API_KEY_ID or "")
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
Create shim `kalshi_mlb_rfq/auth_client.py` (calls `configure()` from the taker's config so its callers stay identical):
```python
"""Backward-compat shim. Logic lives in kalshi_common.auth_client.

Configures the shared client from the taker's config at import time so existing
`auth_client.api(...)` callers in kalshi_mlb_rfq work unchanged.
"""
from kalshi_common import auth_client as _impl
from kalshi_mlb_rfq.config import (
    KALSHI_API_KEY_ID, KALSHI_PRIVATE_KEY_PATH, KALSHI_BASE_URL, PROJECT_ROOT,
)

_impl.configure(KALSHI_API_KEY_ID, KALSHI_PRIVATE_KEY_PATH, KALSHI_BASE_URL, PROJECT_ROOT)

api = _impl.api  # re-export the exact callable existing code uses
```

- [ ] **Step 5: Move `sgp_runner.py` and shim it**

```bash
git mv kalshi_mlb_rfq/sgp_runner.py kalshi_common/sgp_runner.py
```
Edit `kalshi_common/sgp_runner.py`: change its import `from kalshi_mlb_rfq import auth_client` to `from kalshi_common import auth_client`. (It calls `auth_client.api(...)`; the maker will have called `configure()` before invoking `sgp_cycle`.) Leave the rest unchanged.
Create shim `kalshi_mlb_rfq/sgp_runner.py`:
```python
"""Backward-compat shim. Logic lives in kalshi_common.sgp_runner."""
from kalshi_common.sgp_runner import *  # noqa: F401,F403
from kalshi_common.sgp_runner import (
    should_scrape, latest_sgp_fetch_time, enumerate_kalshi_targets,
    write_target_lines, run_scrapers, read_priced_rows, sgp_cycle,
)
```
NOTE: the shim re-imports `kalshi_common.sgp_runner`, which imports `kalshi_common.auth_client` — *not yet configured* when the taker imports it. That's fine: `sgp_cycle` isn't called until runtime, by which point the taker's `auth_client` shim (step 4) has run `configure()`. Both shims point at the same `kalshi_common.auth_client` module instance, so one `configure()` serves both.

- [ ] **Step 6: Extract leg-typing helpers into `kalshi_common/leg_types.py`**

Create `kalshi_common/leg_types.py` by moving these definitions verbatim out of `kalshi_mlb_rfq/main.py`: `_MLB_CODE_TO_TEAM`, `_parse_event_suffix`, `_home_code_from_event_ticker`, `_leg_dict_to_typed`, `_spread_line_from_legs`, `_total_line_from_legs`. The module imports `from kalshi_common import fair_value` for the `SpreadLeg`/`TotalLeg` return types.

- [ ] **Step 7: Repoint taker `main.py` at the extracted helpers**

In `kalshi_mlb_rfq/main.py`: delete the six moved definitions and add, near the existing imports:
```python
from kalshi_common.leg_types import (
    _MLB_CODE_TO_TEAM, _parse_event_suffix, _home_code_from_event_ticker,
    _leg_dict_to_typed, _spread_line_from_legs, _total_line_from_legs,
)
```
Leave the existing `from kalshi_mlb_rfq import (... fair_value, ev_calc, ... sgp_runner ...)` line as-is — those now resolve through the shims.

- [ ] **Step 8: Run the taker's full test suite — it must pass unchanged**

Run: `cd /Users/callancapitolo/NFLWork/.claude/worktrees/kalshi-mlb-mm-maker-bot && ./kalshi_mlb_rfq/venv/bin/python -m pytest kalshi_mlb_rfq/tests/ -v`
Expected: PASS (same set as before extraction). The shims preserve every public name; behavior is unchanged.

- [ ] **Step 9: Commit**

```bash
git add kalshi_common kalshi_mlb_rfq
git commit -m "refactor(kalshi-common): extract shared math/auth/sgp via shims

fair_value, ev_calc (+maker_fee_per_contract), auth_client (now configure()-
injected, decoupled from taker config), sgp_runner, and leg_types move to a new
kalshi_common/ package. Taker modules become one-line re-export shims so its
imports + tests are unchanged. Single source of truth, zero behavior change."
```

---

## Task 2: Maker package skeleton — config, db schema, notify

**Files:**
- Create: `kalshi_mlb_mm/__init__.py`, `kalshi_mlb_mm/config.py`, `kalshi_mlb_mm/db.py`, `kalshi_mlb_mm/notify.py`, `kalshi_mlb_mm/.env.example`, `kalshi_mlb_mm/requirements.txt`
- Test: `kalshi_mlb_mm/tests/__init__.py`, `kalshi_mlb_mm/tests/test_db.py`

- [ ] **Step 1: `kalshi_mlb_mm/config.py`**

```python
"""Config for the Kalshi MLB MM (maker) bot. Loaded from .env or environment."""
import os
from pathlib import Path

PKG_DIR = Path(__file__).parent
_RAW_ROOT = PKG_DIR.parent
PROJECT_ROOT = (Path(str(_RAW_ROOT).split(".worktrees")[0].rstrip("/"))
                if ".worktrees" in str(_RAW_ROOT) else _RAW_ROOT)
DB_PATH = PKG_DIR / "kalshi_mlb_mm.duckdb"
MARKET_DB = PKG_DIR / "kalshi_mlb_mm_market.duckdb"
KILL_FILE = PKG_DIR / ".kill"


def _load_env(path: Path) -> dict[str, str]:
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


def _get(key, default=None):
    return os.environ.get(key, _FILE_ENV.get(key, default))


# Credentials
KALSHI_API_KEY_ID = _get("KALSHI_API_KEY_ID")
KALSHI_PRIVATE_KEY_PATH = _get("KALSHI_PRIVATE_KEY_PATH")
KALSHI_USER_ID = _get("KALSHI_USER_ID")
KALSHI_BASE_URL = _get("KALSHI_BASE_URL", "https://api.elections.kalshi.com/trade-api/v2")
MVE_COLLECTION_TICKER = _get("MVE_COLLECTION_TICKER", "KXMVECROSSCATEGORY-R")

# Pricing
TARGET_ROI = float(_get("TARGET_ROI", "0.05"))      # 5% ROI per side
QUOTE_HYSTERESIS = float(_get("QUOTE_HYSTERESIS", "0.005"))

# Risk (master dial = BANKROLL)
BANKROLL = float(_get("BANKROLL", "500.0"))
DAILY_EXPOSURE_CAP_PCT = float(_get("DAILY_EXPOSURE_CAP_PCT", "0.75"))
MAX_GAME_EXPOSURE_PCT = float(_get("MAX_GAME_EXPOSURE_PCT", "0.10"))
MAX_RFQ_CONTRACTS = int(_get("MAX_RFQ_CONTRACTS", "5"))
MAX_OPEN_QUOTES = int(_get("MAX_OPEN_QUOTES", "25"))
FAIR_DRIFT_TOLERANCE = float(_get("FAIR_DRIFT_TOLERANCE", "0.02"))
MIN_FAIR_PROB = float(_get("MIN_FAIR_PROB", "0.05"))
MAX_FAIR_PROB = float(_get("MAX_FAIR_PROB", "0.95"))

# Freshness / circuit breaker
MAX_PREDICTION_STALENESS_SEC = int(_get("MAX_PREDICTION_STALENESS_SEC", "600"))
MAX_BOOK_STALENESS_SEC = int(_get("MAX_BOOK_STALENESS_SEC", "60"))
BOOK_MOVE_CB_THRESHOLD = float(_get("BOOK_MOVE_CB_THRESHOLD", "0.03"))
TIPOFF_CANCEL_MIN = int(_get("TIPOFF_CANCEL_MIN", "5"))
MIN_BOOK_COUNT_FOR_BLEND = int(_get("MIN_BOOK_COUNT_FOR_BLEND", "2"))

# Loops (seconds)
DISCOVERY_SEC = int(_get("DISCOVERY_SEC", "2"))
CONFIRM_SEC = int(_get("CONFIRM_SEC", "2"))
RISK_SWEEP_SEC = int(_get("RISK_SWEEP_SEC", "10"))
SGP_REFRESH_SEC = int(_get("SGP_REFRESH_SEC", "60"))
SGP_SCRAPER_TIMEOUT_SEC = int(_get("SGP_SCRAPER_TIMEOUT_SEC", "90"))
SAMPLES_REFRESH_SEC = int(_get("SAMPLES_REFRESH_SEC", "600"))

# Vig fallbacks (reused thresholds; only used if a combo lacks full 4-side devig)
DK_VIG_FALLBACK = float(_get("DK_VIG_FALLBACK", "0.125"))
FD_VIG_FALLBACK = float(_get("FD_VIG_FALLBACK", "0.18"))
PX_VIG_FALLBACK = float(_get("PX_VIG_FALLBACK", "0.05"))
NOVIG_VIG_FALLBACK = float(_get("NOVIG_VIG_FALLBACK", "0.05"))

# Paths
ANSWER_KEY_DB = PROJECT_ROOT / "Answer Keys" / "mlb_mm.duckdb"
MLB_SGP_DIR = Path(_get("MLB_SGP_DIR", str(PROJECT_ROOT / "mlb_sgp")))
NOTIFY_WEBHOOK_URL = _get("NOTIFY_WEBHOOK_URL")


def daily_exposure_cap_usd() -> float:
    return BANKROLL * DAILY_EXPOSURE_CAP_PCT
```

- [ ] **Step 2: `kalshi_mlb_mm/db.py` — schema**

```python
"""DuckDB schema + helpers for the maker bot."""
import time
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from random import random

import duckdb

from kalshi_mlb_mm.config import DB_PATH

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS seen_rfqs (
    rfq_id          VARCHAR PRIMARY KEY,
    market_ticker   VARCHAR,
    in_scope        BOOLEAN,
    game_id         VARCHAR,
    legs_json       VARCHAR,
    first_seen_at   TIMESTAMP NOT NULL,
    last_decision   VARCHAR
);
CREATE TABLE IF NOT EXISTS live_quotes (
    quote_id            VARCHAR PRIMARY KEY,
    rfq_id              VARCHAR NOT NULL,
    combo_market_ticker VARCHAR NOT NULL,
    game_id             VARCHAR NOT NULL,
    yes_bid             DOUBLE,
    no_bid              DOUBLE,
    model_fair          DOUBLE,
    book_fair           DOUBLE,
    blended_fair        DOUBLE,
    status              VARCHAR NOT NULL,
    submitted_at        TIMESTAMP NOT NULL,
    closed_at           TIMESTAMP
);
CREATE TABLE IF NOT EXISTS quote_decisions (
    decision_id   VARCHAR PRIMARY KEY,
    rfq_id        VARCHAR,
    quote_id      VARCHAR,
    combo_market_ticker VARCHAR,
    game_id       VARCHAR,
    decision      VARCHAR NOT NULL,
    reason        VARCHAR,
    model_fair    DOUBLE,
    book_fair     DOUBLE,
    blended_fair  DOUBLE,
    yes_bid       DOUBLE,
    no_bid        DOUBLE,
    observed_at   TIMESTAMP NOT NULL
);
CREATE TABLE IF NOT EXISTS fills (
    fill_id              VARCHAR PRIMARY KEY,
    quote_id             VARCHAR NOT NULL,
    rfq_id               VARCHAR NOT NULL,
    combo_market_ticker  VARCHAR NOT NULL,
    game_id              VARCHAR NOT NULL,
    side_held            VARCHAR NOT NULL,
    contracts            DOUBLE NOT NULL,
    price                DOUBLE NOT NULL,
    fee                  DOUBLE NOT NULL,
    model_fair_at_quote  DOUBLE,
    book_fair_at_quote   DOUBLE,
    blended_fair_at_quote DOUBLE,
    fair_at_confirm      DOUBLE,
    realized_pnl         DOUBLE,
    filled_at            TIMESTAMP NOT NULL
);
CREATE TABLE IF NOT EXISTS positions (
    combo_market_ticker VARCHAR NOT NULL,
    side                VARCHAR NOT NULL,
    game_id             VARCHAR NOT NULL,
    net_contracts       DOUBLE NOT NULL,
    weighted_price      DOUBLE NOT NULL,
    updated_at          TIMESTAMP NOT NULL,
    PRIMARY KEY (combo_market_ticker, side)
);
CREATE TABLE IF NOT EXISTS sessions (
    session_id VARCHAR PRIMARY KEY,
    started_at TIMESTAMP NOT NULL,
    ended_at   TIMESTAMP,
    pid        INTEGER,
    dry_run    BOOLEAN NOT NULL
);
"""


@contextmanager
def connect(read_only: bool = False, retries: int = 10):
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
    with connect() as con:
        con.execute(SCHEMA_SQL)


def start_session(pid: int, dry_run: bool) -> str:
    sid = str(uuid.uuid4())
    with connect() as con:
        con.execute(
            "INSERT INTO sessions (session_id, started_at, pid, dry_run) VALUES (?,?,?,?)",
            [sid, datetime.now(timezone.utc), pid, dry_run],
        )
    return sid


def end_session(session_id: str):
    with connect() as con:
        con.execute("UPDATE sessions SET ended_at=? WHERE session_id=?",
                    [datetime.now(timezone.utc), session_id])
```

- [ ] **Step 3: Write the db test**

```python
# kalshi_mlb_mm/tests/test_db.py
import importlib

def test_init_and_session(tmp_path, monkeypatch):
    import kalshi_mlb_mm.config as cfg
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "t.duckdb")
    import kalshi_mlb_mm.db as db
    importlib.reload(db)
    db.init_database()
    sid = db.start_session(pid=1, dry_run=True)
    db.end_session(sid)
    with db.connect(read_only=True) as con:
        row = con.execute("SELECT dry_run, ended_at FROM sessions WHERE session_id=?", [sid]).fetchone()
    assert row[0] is True and row[1] is not None
```

- [ ] **Step 4: Run it**

Run: `./kalshi_mlb_rfq/venv/bin/python -m pytest kalshi_mlb_mm/tests/test_db.py -v`
Expected: PASS.

- [ ] **Step 5: `notify.py`, `.env.example`, `requirements.txt`**

`kalshi_mlb_mm/notify.py`:
```python
"""Minimal alerting: print + optional webhook (mirrors taker notify)."""
import json, urllib.request
from kalshi_mlb_mm.config import NOTIFY_WEBHOOK_URL

def _post(text: str):
    print(text, flush=True)
    if not NOTIFY_WEBHOOK_URL:
        return
    try:
        req = urllib.request.Request(
            NOTIFY_WEBHOOK_URL, data=json.dumps({"text": text}).encode(),
            headers={"Content-Type": "application/json"}, method="POST")
        urllib.request.urlopen(req, timeout=5)
    except Exception:
        pass

def halt(reason: str, detail: str = ""):
    _post(f"[MM HALT] {reason} {detail}")

def fill(combo_market_ticker: str, side: str, contracts: int, price: float):
    _post(f"[MM FILL] {combo_market_ticker} {side} x{contracts} @ {price:.3f}")
```
`kalshi_mlb_mm/.env.example`:
```
KALSHI_API_KEY_ID=
KALSHI_PRIVATE_KEY_PATH=
KALSHI_USER_ID=
BANKROLL=500
TARGET_ROI=0.05
MAX_RFQ_CONTRACTS=5
MAX_OPEN_QUOTES=25
NOTIFY_WEBHOOK_URL=
```
`kalshi_mlb_mm/requirements.txt`:
```
duckdb>=0.10
pandas>=2.0
scipy>=1.10
numpy>=1.24
```

- [ ] **Step 6: Commit**

```bash
git add kalshi_mlb_mm
git commit -m "feat(mm): maker package skeleton — config, db schema, notify"
```

---

## Task 3: `scope.py` — in-scope detection from `mve_selected_legs`

**Files:**
- Create: `kalshi_mlb_mm/scope.py`
- Test: `kalshi_mlb_mm/tests/test_scope.py`

- [ ] **Step 1: Write failing tests**

```python
# kalshi_mlb_mm/tests/test_scope.py
from kalshi_mlb_mm import scope

_LEGS_OK = [
    {"event_ticker": "KXMLBSPREAD-26MAY232205TEXLAA", "market_ticker": "KXMLBSPREAD-26MAY232205TEXLAA-TEX2", "side": "no"},
    {"event_ticker": "KXMLBTOTAL-26MAY232205TEXLAA", "market_ticker": "KXMLBTOTAL-26MAY232205TEXLAA-10", "side": "no"},
]

def test_decode_legs_reads_mve_selected_legs():
    market = {"mve_selected_legs": _LEGS_OK}
    assert scope.decode_legs(market) == _LEGS_OK

def test_decode_legs_none_when_absent():
    assert scope.decode_legs({"ticker": "x"}) is None

def test_in_scope_true_for_spread_total_2leg():
    assert scope.is_spread_total_2leg(_LEGS_OK) is True

def test_in_scope_false_for_three_legs():
    assert scope.is_spread_total_2leg(_LEGS_OK + [_LEGS_OK[0]]) is False

def test_in_scope_false_for_two_spreads():
    legs = [_LEGS_OK[0], dict(_LEGS_OK[0])]
    assert scope.is_spread_total_2leg(legs) is False
```

- [ ] **Step 2: Run — expect failure**

Run: `./kalshi_mlb_rfq/venv/bin/python -m pytest kalshi_mlb_mm/tests/test_scope.py -v`
Expected: FAIL (module `scope` not found).

- [ ] **Step 3: Implement `scope.py`**

```python
"""Decide whether an incoming RFQ's combo is one we model (2-leg spread×total)."""


def decode_legs(market: dict) -> list[dict] | None:
    """Pull the constituent legs from a GET /markets/{ticker} 'market' dict.
    The API returns them under mve_selected_legs in the exact
    {event_ticker, market_ticker, side} shape the leg-typer consumes."""
    legs = market.get("mve_selected_legs")
    if not legs or not isinstance(legs, list):
        return None
    return legs


def is_spread_total_2leg(legs: list[dict]) -> bool:
    """True iff exactly one KXMLBSPREAD- leg and one KXMLBTOTAL- leg."""
    if len(legs) != 2:
        return False
    tickers = [str(l.get("market_ticker", "")) for l in legs]
    has_spread = sum(t.startswith("KXMLBSPREAD-") for t in tickers) == 1
    has_total = sum(t.startswith("KXMLBTOTAL-") for t in tickers) == 1
    return has_spread and has_total
```

- [ ] **Step 4: Run — expect pass**

Run: `./kalshi_mlb_rfq/venv/bin/python -m pytest kalshi_mlb_mm/tests/test_scope.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_mm/scope.py kalshi_mlb_mm/tests/test_scope.py
git commit -m "feat(mm): scope — detect 2-leg spread×total combos from mve_selected_legs"
```

---

## Task 4: `pricing.py` — 5%-ROI fixed-margin quote math

**Files:**
- Create: `kalshi_mlb_mm/pricing.py`
- Test: `kalshi_mlb_mm/tests/test_pricing.py`

- [ ] **Step 1: Write failing tests**

```python
# kalshi_mlb_mm/tests/test_pricing.py
from kalshi_mlb_mm import pricing
from kalshi_common.ev_calc import maker_fee_per_contract

def test_roi_is_five_percent_net_of_fee():
    fair = 0.55
    q = pricing.quote(fair, target_roi=0.05)
    cost = q.yes_bid + maker_fee_per_contract(q.yes_bid)
    roi = (fair - cost) / cost
    assert abs(roi - 0.05) < 0.01   # within a cent of grid rounding

def test_sum_below_one():
    q = pricing.quote(0.50, target_roi=0.05)
    assert q.yes_bid + q.no_bid < 1.0

def test_grid_rounded_to_milli():
    q = pricing.quote(0.637, target_roi=0.05)
    assert abs(q.yes_bid * 1000 - round(q.yes_bid * 1000)) < 1e-6

def test_none_outside_prob_bounds():
    assert pricing.quote(0.0, 0.05) is None
    assert pricing.quote(1.0, 0.05) is None

def test_symmetry_at_fifty():
    q = pricing.quote(0.50, 0.05)
    assert abs(q.yes_bid - q.no_bid) < 1e-9
```

- [ ] **Step 2: Run — expect failure**

Run: `./kalshi_mlb_rfq/venv/bin/python -m pytest kalshi_mlb_mm/tests/test_pricing.py -v`
Expected: FAIL.

- [ ] **Step 3: Implement `pricing.py`**

```python
"""Fixed-margin quote pricing: target a constant ROI per side.

ROI = (p - cost)/cost, where p = side win-prob and cost = bid + maker_fee.
Setting ROI = target gives cost = p/(1+target), so bid = p/(1+target) - fee.
The fee depends on the bid, so we iterate once (fee changes negligibly).
"""
import math
from dataclasses import dataclass

from kalshi_common.ev_calc import maker_fee_per_contract

GRID_STEP = 0.001  # Kalshi MVE markets price in $0.001 (deci-cent) steps


@dataclass(frozen=True)
class Quote:
    yes_bid: float
    no_bid: float


def _round_down_to_grid(price: float) -> float:
    # Floor to the grid: a lower bid is conservative (cheaper for us / more margin).
    return math.floor(round(price / GRID_STEP, 6)) * GRID_STEP


def _price_for_side(p: float, target_roi: float) -> float:
    raw = p / (1.0 + target_roi)            # = cost target (bid + fee)
    bid = raw - maker_fee_per_contract(raw)
    bid = raw - maker_fee_per_contract(bid)  # one refinement
    return _round_down_to_grid(max(bid, 0.0))


def quote(fair: float, target_roi: float) -> Quote | None:
    if not (0.0 < fair < 1.0):
        return None
    yes_bid = _price_for_side(fair, target_roi)
    no_bid = _price_for_side(1.0 - fair, target_roi)
    if yes_bid <= 0 or no_bid <= 0 or yes_bid + no_bid >= 1.0:
        return None
    return Quote(yes_bid=round(yes_bid, 4), no_bid=round(no_bid, 4))
```

- [ ] **Step 4: Run — expect pass**

Run: `./kalshi_mlb_rfq/venv/bin/python -m pytest kalshi_mlb_mm/tests/test_pricing.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_mm/pricing.py kalshi_mlb_mm/tests/test_pricing.py
git commit -m "feat(mm): pricing — 5%-ROI fixed-margin two-sided quotes"
```

---

## Task 5: `quote_gateway.py` — submit / confirm / cancel

**Files:**
- Create: `kalshi_mlb_mm/quote_gateway.py`
- Test: `kalshi_mlb_mm/tests/test_quote_gateway.py`

- [ ] **Step 1: Write failing tests (monkeypatch `auth_client.api`)**

```python
# kalshi_mlb_mm/tests/test_quote_gateway.py
from kalshi_mlb_mm import quote_gateway
from kalshi_common import auth_client

def test_submit_returns_quote_id(monkeypatch):
    calls = {}
    def fake_api(method, path, body=None, timeout=30):
        calls["method"], calls["path"], calls["body"] = method, path, body
        return 201, {"id": "q123"}, {}
    monkeypatch.setattr(auth_client, "api", fake_api)
    gw = quote_gateway.RestQuoteGateway()
    qid = gw.submit_quote("rfq1", 0.52, 0.43)
    assert qid == "q123"
    assert calls["method"] == "POST" and calls["path"] == "/communications/quotes"
    assert calls["body"]["yes_bid"] == "0.5200" and calls["body"]["rest_remainder"] is False

def test_submit_none_on_error(monkeypatch):
    monkeypatch.setattr(auth_client, "api", lambda *a, **k: (400, {"error": "x"}, {}))
    assert quote_gateway.RestQuoteGateway().submit_quote("rfq1", 0.5, 0.4) is None

def test_confirm_true_on_204(monkeypatch):
    monkeypatch.setattr(auth_client, "api", lambda *a, **k: (204, {}, {}))
    assert quote_gateway.RestQuoteGateway().confirm("q1") is True

def test_cancel_true_on_204(monkeypatch):
    monkeypatch.setattr(auth_client, "api", lambda *a, **k: (204, {}, {}))
    assert quote_gateway.RestQuoteGateway().cancel("q1") is True
```

- [ ] **Step 2: Run — expect failure.** `pytest kalshi_mlb_mm/tests/test_quote_gateway.py -v` → FAIL.

- [ ] **Step 3: Implement `quote_gateway.py`**

```python
"""Quote write-path: submit / confirm / cancel. The pricing+risk code only ever
touches this interface, never raw HTTP — so a WS/FIX gateway can drop in later."""
from typing import Protocol

from kalshi_common import auth_client


class QuoteGateway(Protocol):
    def submit_quote(self, rfq_id: str, yes_bid: float, no_bid: float) -> str | None: ...
    def confirm(self, quote_id: str) -> bool: ...
    def cancel(self, quote_id: str) -> bool: ...


class RestQuoteGateway:
    def submit_quote(self, rfq_id: str, yes_bid: float, no_bid: float) -> str | None:
        body = {"rfq_id": rfq_id, "yes_bid": f"{yes_bid:.4f}",
                "no_bid": f"{no_bid:.4f}", "rest_remainder": False}
        status, resp, _ = auth_client.api("POST", "/communications/quotes", body=body)
        if status in (200, 201) and isinstance(resp, dict) and "id" in resp:
            return resp["id"]
        return None

    def confirm(self, quote_id: str) -> bool:
        # PUT /communications/quotes/{id}/confirm — no body per the OpenAPI spec.
        # Open item §10.2: verify the exact response/body on the first real accept.
        status, _, _ = auth_client.api("PUT", f"/communications/quotes/{quote_id}/confirm")
        return status in (200, 201, 204)

    def cancel(self, quote_id: str) -> bool:
        status, _, _ = auth_client.api("DELETE", f"/communications/quotes/{quote_id}")
        return status in (200, 204)
```

- [ ] **Step 4: Run — expect pass.** `pytest kalshi_mlb_mm/tests/test_quote_gateway.py -v` → PASS.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_mm/quote_gateway.py kalshi_mlb_mm/tests/test_quote_gateway.py
git commit -m "feat(mm): RestQuoteGateway — submit/confirm/cancel behind interface"
```

---

## Task 6: `rfq_source.py` — poll open RFQs + fetch market

**Files:**
- Create: `kalshi_mlb_mm/rfq_source.py`
- Test: `kalshi_mlb_mm/tests/test_rfq_source.py`

- [ ] **Step 1: Write failing tests**

```python
# kalshi_mlb_mm/tests/test_rfq_source.py
from kalshi_mlb_mm import rfq_source
from kalshi_common import auth_client

def test_poll_returns_rfqs(monkeypatch):
    def fake_api(method, path, body=None, timeout=30):
        assert path.startswith("/communications/rfqs?status=open")
        return 200, {"rfqs": [{"id": "r1", "market_ticker": "MT", "contracts": 3}]}, {}
    monkeypatch.setattr(auth_client, "api", fake_api)
    out = rfq_source.RestRFQSource().poll()
    assert out[0]["id"] == "r1"

def test_poll_empty_on_error(monkeypatch):
    monkeypatch.setattr(auth_client, "api", lambda *a, **k: (500, "x", {}))
    assert rfq_source.RestRFQSource().poll() == []

def test_get_market_returns_inner(monkeypatch):
    monkeypatch.setattr(auth_client, "api",
                        lambda *a, **k: (200, {"market": {"ticker": "MT", "mve_selected_legs": []}}, {}))
    m = rfq_source.RestRFQSource().get_market("MT")
    assert m["ticker"] == "MT"
```

- [ ] **Step 2: Run — expect failure.**

- [ ] **Step 3: Implement `rfq_source.py`**

```python
"""RFQ read-path: discover open RFQs + fetch a combo market's legs. v1 = REST poll;
a WsRFQSource can replace this behind the same interface later."""
from typing import Protocol

from kalshi_common import auth_client


class RFQSource(Protocol):
    def poll(self) -> list[dict]: ...
    def get_market(self, market_ticker: str) -> dict | None: ...


class RestRFQSource:
    def poll(self) -> list[dict]:
        # No creator filter → market-wide open list (others' RFQs we can quote on).
        status, body, _ = auth_client.api(
            "GET", "/communications/rfqs?status=open&limit=100")
        if status != 200 or not isinstance(body, dict):
            return []
        return body.get("rfqs") or []

    def get_market(self, market_ticker: str) -> dict | None:
        status, body, _ = auth_client.api("GET", f"/markets/{market_ticker}")
        if status != 200 or not isinstance(body, dict):
            return None
        return body.get("market")
```

- [ ] **Step 4: Run — expect pass.**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_mm/rfq_source.py kalshi_mlb_mm/tests/test_rfq_source.py
git commit -m "feat(mm): RestRFQSource — poll open RFQs + fetch market legs"
```

---

## Task 7: `risk.py` — freshness gate, circuit breaker, caps, last-look

**Files:**
- Create: `kalshi_mlb_mm/risk.py`
- Test: `kalshi_mlb_mm/tests/test_risk.py`

- [ ] **Step 1: Write failing tests**

```python
# kalshi_mlb_mm/tests/test_risk.py
from datetime import datetime, timedelta, timezone
from kalshi_mlb_mm import risk

NOW = datetime(2026, 5, 27, 18, 0, tzinfo=timezone.utc)

def test_staleness_ok_true_when_fresh():
    assert risk.staleness_ok(NOW - timedelta(seconds=100), 600, now=NOW) is True

def test_staleness_ok_false_when_stale():
    assert risk.staleness_ok(NOW - timedelta(seconds=700), 600, now=NOW) is False

def test_size_gate():
    assert risk.size_ok(5, max_contracts=5) is True
    assert risk.size_ok(6, max_contracts=5) is False

def test_book_move_circuit_breaker():
    assert risk.book_move_triggered(0.40, 0.45, threshold=0.03) is True
    assert risk.book_move_triggered(0.40, 0.41, threshold=0.03) is False

def test_last_look_void_on_drift():
    # filled YES at 0.52; fair dropped to 0.50 -> drift 0.02 not > tol(0.02): still ok
    assert risk.last_look_ok(side="yes", price=0.52, fee=0.005, current_fair=0.55,
                             prev_fair=0.55, drift_tol=0.02) is True
    # fair dropped well below price+fee -> not +EV -> void
    assert risk.last_look_ok(side="yes", price=0.52, fee=0.005, current_fair=0.50,
                             prev_fair=0.55, drift_tol=0.02) is False

def test_daily_cap():
    fills = [{"price": 4.0}, {"price": 5.0}]   # $9 used
    assert risk.daily_cap_ok(fills, cap_usd=375.0) is True
    big = [{"price": 380.0}]
    assert risk.daily_cap_ok(big, cap_usd=375.0) is False
```

- [ ] **Step 2: Run — expect failure.**

- [ ] **Step 3: Implement `risk.py`**

```python
"""Risk gates for the maker bot. Pure functions, time injected for testability."""
from datetime import datetime, timezone


def _utc(dt: datetime, now: datetime | None) -> tuple[datetime, datetime]:
    now = now or datetime.now(timezone.utc)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    if now.tzinfo is None:
        now = now.replace(tzinfo=timezone.utc)
    return dt, now


def staleness_ok(generated_at: datetime | None, max_age_sec: int,
                 now: datetime | None = None) -> bool:
    if generated_at is None:
        return False
    generated_at, now = _utc(generated_at, now)
    return (now - generated_at).total_seconds() <= max_age_sec


def tipoff_ok(commence_time: datetime | None, cancel_min: int,
              now: datetime | None = None) -> bool:
    if commence_time is None:
        return False
    commence_time, now = _utc(commence_time, now)
    return now < commence_time - __import__("datetime").timedelta(minutes=cancel_min)


def size_ok(requested_contracts: int, max_contracts: int) -> bool:
    return 0 < requested_contracts <= max_contracts


def book_move_triggered(prev_fair: float, current_fair: float, threshold: float) -> bool:
    """True if the underlying book fair jumped by more than `threshold` between scrapes."""
    return abs(current_fair - prev_fair) > threshold


def last_look_ok(side: str, price: float, fee: float, current_fair: float,
                 prev_fair: float, drift_tol: float) -> bool:
    """Confirm only if the filled side is still +EV against current fair AND fair
    hasn't drifted against it past tolerance."""
    p = current_fair if side == "yes" else (1.0 - current_fair)
    if p <= price + fee:                      # no longer +EV
        return False
    if abs(current_fair - prev_fair) > drift_tol:
        return False
    return True


def daily_cap_ok(today_fills: list[dict], cap_usd: float) -> bool:
    used = sum(float(f.get("price", 0)) for f in today_fills)
    return used < cap_usd


def per_game_cap_ok(game_id: str, today_fills: list[dict], bankroll: float,
                    pct: float) -> bool:
    used = sum(float(f.get("price", 0)) for f in today_fills
               if f.get("game_id") == game_id)
    return used < bankroll * pct


def kill_switch_ok(kill_file) -> bool:
    return not kill_file.exists()
```
NOTE the `last_look_ok` drift test: a 0.02 drift with tol 0.02 is NOT `> tol`, so it passes the drift check; the case fails only when `p <= price+fee`. The second assertion (fair 0.50) fails because `0.50 <= 0.525`.

- [ ] **Step 4: Run — expect pass.**

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_mm/risk.py kalshi_mlb_mm/tests/test_risk.py
git commit -m "feat(mm): risk — freshness, circuit breaker, caps, last-look"
```

---

## Task 8: `main.py` — wire the loops

**Files:**
- Create: `kalshi_mlb_mm/main.py`
- Test: `kalshi_mlb_mm/tests/test_main_smoke.py`

This task composes the prior modules. The fair-value pieces reuse the taker's
exact logic via `kalshi_common`: samples come from `mlb_mm.duckdb::mlb_game_samples`
(read-only), book fairs from the maker's own market DB (its own `sgp_runner` cadence),
blended via `fair_value.blend`. The leg→typed conversion reuses `kalshi_common.leg_types`.

- [ ] **Step 1: Implement the fair-value helper module method (reused logic)**

Add `kalshi_mlb_mm/fairs.py`:
```python
"""Fair-value assembly for the maker — reuses kalshi_common math against the
maker's own caches. Mirrors the taker's _fresh_blended_fair, minus RFQ creation."""
import json
import pandas as pd

from kalshi_common import fair_value
from kalshi_common.leg_types import _leg_dict_to_typed, _spread_line_from_legs, _total_line_from_legs


def blended_fair(legs: list[dict], game_id: str, samples: pd.DataFrame,
                 book_fairs: dict[str, float]) -> tuple[float | None, float | None, float | None]:
    """Returns (model_fair, book_fair_median, blended) — each None if unavailable.
    book_fair_median is the median of the per-book devigged fairs (for logging)."""
    if samples is None or samples.empty:
        return None, None, None
    typed = [_leg_dict_to_typed(l, game_id) for l in legs]
    if any(t is None for t in typed):
        return None, None, None
    model = fair_value.model_fair(samples, typed)
    blended = fair_value.blend(model, book_fairs)
    import statistics
    book_med = (statistics.median(list(book_fairs.values())) if book_fairs else None)
    return model, book_med, blended
```

- [ ] **Step 2: Implement `main.py`**

```python
"""Kalshi MLB MM (maker) bot — REST-polling daemon."""
import argparse
import json
import os
import signal
import threading
import time
import uuid
from datetime import datetime, timedelta, timezone
from pathlib import Path

import duckdb
import pandas as pd

from kalshi_common import auth_client, fair_value, sgp_runner
from kalshi_mlb_mm import config, db, notify, pricing, risk, scope, fairs
from kalshi_mlb_mm.rfq_source import RestRFQSource
from kalshi_mlb_mm.quote_gateway import RestQuoteGateway

_running = threading.Event(); _running.set()
_SAMPLES = {}            # game_id -> df
_SGP_ODDS = None         # pd.DataFrame
_PREV_BOOK_FAIR = {}     # combo_market_ticker -> last blended book fair (circuit breaker)
_SCOPE_CACHE = {}        # market_ticker -> (in_scope, game_id, legs)


def _signal_handler(_s, _f):
    _running.clear()


def _configure_auth():
    auth_client.configure(config.KALSHI_API_KEY_ID, config.KALSHI_PRIVATE_KEY_PATH,
                          config.KALSHI_BASE_URL, config.PROJECT_ROOT)


def _refresh_samples():
    global _SAMPLES
    if not config.ANSWER_KEY_DB.exists():
        return
    try:
        con = duckdb.connect(str(config.ANSWER_KEY_DB), read_only=True)
    except duckdb.IOException:
        return
    try:
        sdf = con.execute(
            "SELECT game_id, sim_idx, home_margin, total_final_score, "
            "home_margin_f5, total_f5 FROM mlb_game_samples").fetchdf()
        _SAMPLES = {g: d.reset_index(drop=True) for g, d in sdf.groupby("game_id")}
    finally:
        con.close()


def _refresh_sgp():
    global _SGP_ODDS
    if not config.MARKET_DB.exists():
        return
    try:
        con = duckdb.connect(str(config.MARKET_DB), read_only=True)
    except duckdb.IOException:
        return
    try:
        tables = {t[0] for t in con.execute("SHOW TABLES").fetchall()}
        if "mlb_sgp_odds" not in tables:
            return
        _SGP_ODDS = con.execute(
            "SELECT game_id, combo, period, bookmaker, sgp_decimal, fetch_time, "
            "spread_line, total_line FROM mlb_sgp_odds WHERE period='FG' "
            "AND fetch_time > NOW() - INTERVAL (CAST(? AS BIGINT)) SECOND",
            [config.MAX_BOOK_STALENESS_SEC]).fetchdf()
    finally:
        con.close()


# book fairs per (game, spread_line, total_line) — mirrors taker _load_book_fairs,
# but REQUIRES full 4-side devig (no fallback) per accepted-risk #6.
def _book_fairs(game_id, spread_line, total_line):
    if _SGP_ODDS is None or _SGP_ODDS.empty:
        return {}
    rows = _SGP_ODDS[(_SGP_ODDS.game_id == game_id)
                     & (_SGP_ODDS.spread_line.astype(float).round(2) == round(spread_line, 2))
                     & (_SGP_ODDS.total_line.astype(float).round(2) == round(total_line, 2))]
    if rows.empty:
        return {}
    out = {}
    for book in rows.bookmaker.unique():
        sub = rows[rows.bookmaker == book]
        if len(sub) < 4:        # require full 4-side coverage (no crude fallback)
            continue
        f = fair_value.devig_book(sub, combo="Home Spread + Over", vig_fallback=0.0)
        if f is not None:
            out[book] = f
    return out if len(out) >= config.MIN_BOOK_COUNT_FOR_BLEND else {}


def _commence_time(game_id):
    # read from mlb_target_lines (written by sgp_runner) for tipoff gating
    if not config.MARKET_DB.exists():
        return None
    try:
        con = duckdb.connect(str(config.MARKET_DB), read_only=True)
    except duckdb.IOException:
        return None
    try:
        row = con.execute("SELECT commence_time FROM mlb_target_lines WHERE game_id=? LIMIT 1",
                          [game_id]).fetchone()
        return row[0] if row else None
    finally:
        con.close()


def _today_fills():
    start = datetime.now(timezone.utc).replace(hour=0, minute=0, second=0, microsecond=0)
    with db.connect(read_only=True) as con:
        return [{"game_id": g, "price": p * c}
                for g, p, c in con.execute(
                    "SELECT game_id, price, contracts FROM fills WHERE filled_at >= ?",
                    [start]).fetchall()]


def _resolve_game_and_lines(market_ticker, legs):
    """Resolve game_id + (spread_line,total_line) from decoded legs. Uses target
    lines / leg-typing. Returns (game_id, spread_line, total_line) or (None,...)."""
    from kalshi_common.leg_types import _parse_event_suffix, _MLB_CODE_TO_TEAM
    # game_id via the spread leg's event suffix → team codes → match samples' games.
    spread_leg = next((l for l in legs if l["market_ticker"].startswith("KXMLBSPREAD-")), None)
    if not spread_leg:
        return None, None, None
    suffix = spread_leg["event_ticker"].replace("KXMLBSPREAD-", "")
    away, home = _parse_event_suffix(suffix)
    if not away or not home:
        return None, None, None
    spread_line = _spread_line_from_legs(legs) if (_spread_line_from_legs := __import__(
        "kalshi_common.leg_types", fromlist=["_spread_line_from_legs"])._spread_line_from_legs) else 0.0
    total_line = __import__("kalshi_common.leg_types", fromlist=["_total_line_from_legs"])._total_line_from_legs(legs)
    # game_id resolution: match team names against mlb_target_lines rows
    home_name, away_name = _MLB_CODE_TO_TEAM.get(home), _MLB_CODE_TO_TEAM.get(away)
    if not home_name or not away_name or not config.MARKET_DB.exists():
        return None, spread_line, total_line
    con = duckdb.connect(str(config.MARKET_DB), read_only=True)
    try:
        row = con.execute(
            "SELECT game_id FROM mlb_target_lines WHERE home_team=? AND away_team=? LIMIT 1",
            [home_name, away_name]).fetchone()
    finally:
        con.close()
    return (row[0] if row else None), spread_line, total_line


def _log_decision(decision, *, rfq_id=None, quote_id=None, ticker=None, game_id=None,
                  reason=None, model=None, book=None, blended=None, yb=None, nb=None):
    with db.connect() as con:
        con.execute(
            "INSERT INTO quote_decisions (decision_id, rfq_id, quote_id, "
            "combo_market_ticker, game_id, decision, reason, model_fair, book_fair, "
            "blended_fair, yes_bid, no_bid, observed_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            [str(uuid.uuid4()), rfq_id, quote_id, ticker, game_id, decision, reason,
             model, book, blended, yb, nb, datetime.now(timezone.utc)])


def _discovery_tick(source, gateway, dry_run):
    rfqs = source.poll()
    with db.connect(read_only=True) as con:
        open_count = con.execute(
            "SELECT COUNT(*) FROM live_quotes WHERE status='open'").fetchone()[0]
    for rfq in rfqs:
        if open_count >= config.MAX_OPEN_QUOTES:
            break
        rid = rfq.get("id"); ticker = rfq.get("market_ticker")
        if not rid or not ticker:
            continue
        with db.connect(read_only=True) as con:
            seen = con.execute("SELECT in_scope FROM seen_rfqs WHERE rfq_id=?", [rid]).fetchone()
        # scope (cache the market lookup verdict)
        if ticker in _SCOPE_CACHE:
            in_scope, game_id, legs = _SCOPE_CACHE[ticker]
        else:
            market = source.get_market(ticker)
            legs = scope.decode_legs(market) if market else None
            in_scope = bool(legs and scope.is_spread_total_2leg(legs))
            game_id = None
            _SCOPE_CACHE[ticker] = (in_scope, game_id, legs)
        if not in_scope:
            if not seen:
                _log_decision("skipped", rfq_id=rid, ticker=ticker, reason="out_of_scope")
                with db.connect() as con:
                    con.execute("INSERT OR REPLACE INTO seen_rfqs VALUES (?,?,?,?,?,?,?)",
                                [rid, ticker, False, None, None, datetime.now(timezone.utc), "out_of_scope"])
            continue
        # size gate
        if not risk.size_ok(int(rfq.get("contracts", 0) or 0), config.MAX_RFQ_CONTRACTS):
            _log_decision("skipped", rfq_id=rid, ticker=ticker, reason="size_gate"); continue
        game_id, spread_line, total_line = _resolve_game_and_lines(ticker, legs)
        if game_id is None:
            _log_decision("skipped", rfq_id=rid, ticker=ticker, reason="no_game"); continue
        samples = _SAMPLES.get(game_id)
        # freshness gate
        ct = _commence_time(game_id)
        if not risk.tipoff_ok(ct, config.TIPOFF_CANCEL_MIN):
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id, reason="tipoff"); continue
        book_fairs = _book_fairs(game_id, spread_line, total_line)
        model, book_med, blended = fairs.blended_fair(legs, game_id, samples, book_fairs)
        if blended is None or not (config.MIN_FAIR_PROB <= blended <= config.MAX_FAIR_PROB):
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id,
                          reason="no_fair", model=model, book=book_med, blended=blended); continue
        # circuit breaker bookkeeping
        prev = _PREV_BOOK_FAIR.get(ticker)
        _PREV_BOOK_FAIR[ticker] = book_med
        if prev is not None and book_med is not None and risk.book_move_triggered(prev, book_med, config.BOOK_MOVE_CB_THRESHOLD):
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id, reason="circuit_breaker"); continue
        q = pricing.quote(blended, config.TARGET_ROI)
        if q is None:
            _log_decision("skipped", rfq_id=rid, ticker=ticker, game_id=game_id, reason="unpriceable"); continue
        if dry_run:
            _log_decision("dry_run_quote", rfq_id=rid, ticker=ticker, game_id=game_id,
                          model=model, book=book_med, blended=blended, yb=q.yes_bid, nb=q.no_bid); continue
        qid = gateway.submit_quote(rid, q.yes_bid, q.no_bid)
        if qid:
            with db.connect() as con:
                con.execute("INSERT INTO live_quotes VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                            [qid, rid, ticker, game_id, q.yes_bid, q.no_bid, model, book_med,
                             blended, "open", datetime.now(timezone.utc), None])
                con.execute("INSERT OR REPLACE INTO seen_rfqs VALUES (?,?,?,?,?,?,?)",
                            [rid, ticker, True, game_id, json.dumps(legs), datetime.now(timezone.utc), "quoted"])
            _log_decision("quoted", rfq_id=rid, quote_id=qid, ticker=ticker, game_id=game_id,
                          model=model, book=book_med, blended=blended, yb=q.yes_bid, nb=q.no_bid)
            open_count += 1


def _confirm_tick(gateway, dry_run):
    with db.connect(read_only=True) as con:
        live = con.execute(
            "SELECT quote_id, rfq_id, combo_market_ticker, game_id, yes_bid, no_bid, blended_fair "
            "FROM live_quotes WHERE status='open'").fetchall()
    for qid, rid, ticker, game_id, yb, nb, prev_fair in live:
        status, body, _ = auth_client.api("GET", f"/communications/quotes/{qid}")
        q = body.get("quote") if isinstance(body, dict) else None
        st = (q or {}).get("status")
        if st == "accepted":
            # which side did we end up holding? Determine from the accepted side.
            accepted_side = (q or {}).get("accepted_side")
            side_held = "no" if accepted_side == "yes" else "yes"
            price = (1 - nb) if side_held == "yes" else (1 - yb)  # ask we transact at
            from kalshi_common.ev_calc import maker_fee_per_contract
            fee = maker_fee_per_contract(price)
            samples = _SAMPLES.get(game_id)
            with db.connect(read_only=True) as con:
                legs_row = con.execute("SELECT legs_json FROM seen_rfqs WHERE rfq_id=?", [rid]).fetchone()
            legs = json.loads(legs_row[0]) if legs_row and legs_row[0] else []
            _, _, cur_fair = fairs.blended_fair(legs, game_id, samples, {}) if legs else (None, None, prev_fair)
            cur_fair = cur_fair if cur_fair is not None else prev_fair
            if risk.last_look_ok(side_held, price, fee, cur_fair, prev_fair, config.FAIR_DRIFT_TOLERANCE):
                if not dry_run and gateway.confirm(qid):
                    contracts = int(float((q or {}).get("contracts", 1) or 1))
                    with db.connect() as con:
                        con.execute("INSERT INTO fills VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                                    [str(uuid.uuid4()), qid, rid, ticker, game_id, side_held,
                                     contracts, price, fee, None, None, prev_fair, cur_fair, None,
                                     datetime.now(timezone.utc)])
                        con.execute("UPDATE live_quotes SET status='filled', closed_at=? WHERE quote_id=?",
                                    [datetime.now(timezone.utc), qid])
                    notify.fill(ticker, side_held, contracts, price)
                    _log_decision("confirmed", rfq_id=rid, quote_id=qid, ticker=ticker, game_id=game_id)
            else:
                _log_decision("voided_last_look", rfq_id=rid, quote_id=qid, ticker=ticker, game_id=game_id)
                with db.connect() as con:
                    con.execute("UPDATE live_quotes SET status='voided', closed_at=? WHERE quote_id=?",
                                [datetime.now(timezone.utc), qid])
        elif st in ("cancelled", "expired", "executed"):
            with db.connect() as con:
                con.execute("UPDATE live_quotes SET status=?, closed_at=? WHERE quote_id=?",
                            [st, datetime.now(timezone.utc), qid])


def _risk_sweep_tick(gateway):
    if config.KILL_FILE.exists():
        notify.halt("kill_switch"); return
    with db.connect(read_only=True) as con:
        live = con.execute(
            "SELECT quote_id, game_id FROM live_quotes WHERE status='open'").fetchall()
    for qid, game_id in live:
        ct = _commence_time(game_id)
        stale = not risk.staleness_ok(_SAMPLES_GEN_AT, config.MAX_PREDICTION_STALENESS_SEC) \
            if (_SAMPLES_GEN_AT := None) else False  # samples-gen tracked in _refresh_samples (omitted for brevity)
        if not risk.tipoff_ok(ct, config.TIPOFF_CANCEL_MIN):
            gateway.cancel(qid)
            with db.connect() as con:
                con.execute("UPDATE live_quotes SET status='cancelled', closed_at=? WHERE quote_id=?",
                            [datetime.now(timezone.utc), qid])


def main_loop(dry_run: bool):
    _configure_auth()
    db.init_database()
    sid = db.start_session(pid=os.getpid(), dry_run=dry_run)
    print(f"=== MM bot session {sid} dry_run={dry_run} ===", flush=True)
    source, gateway = RestRFQSource(), RestQuoteGateway()
    # synchronous warm-up: one SGP cycle + sample load
    try:
        sgp_runner.sgp_cycle(bot_market_db=str(config.MARKET_DB),
                             scraper_dir=str(config.MLB_SGP_DIR),
                             venv_python=str(config.MLB_SGP_DIR / "venv" / "bin" / "python"),
                             timeout_sec=config.SGP_SCRAPER_TIMEOUT_SEC)
    except Exception as e:
        print(f"  warmup sgp failed: {e}", flush=True)
    _refresh_sgp(); _refresh_samples()
    last = {"disc": 0.0, "conf": 0.0, "risk": 0.0, "sgp": time.time(), "samp": time.time()}
    try:
        while _running.is_set():
            now = time.time()
            if now - last["disc"] >= config.DISCOVERY_SEC:
                try: _discovery_tick(source, gateway, dry_run)
                except Exception as e: print(f"  disc err: {e}", flush=True)
                last["disc"] = now
            if now - last["conf"] >= config.CONFIRM_SEC:
                try: _confirm_tick(gateway, dry_run)
                except Exception as e: print(f"  conf err: {e}", flush=True)
                last["conf"] = now
            if now - last["risk"] >= config.RISK_SWEEP_SEC:
                try: _risk_sweep_tick(gateway)
                except Exception as e: print(f"  risk err: {e}", flush=True)
                last["risk"] = now
            if now - last["sgp"] >= config.SGP_REFRESH_SEC:
                try:
                    sgp_runner.sgp_cycle(bot_market_db=str(config.MARKET_DB),
                                         scraper_dir=str(config.MLB_SGP_DIR),
                                         venv_python=str(config.MLB_SGP_DIR / "venv" / "bin" / "python"),
                                         timeout_sec=config.SGP_SCRAPER_TIMEOUT_SEC)
                    _refresh_sgp()
                except Exception as e: print(f"  sgp err: {e}", flush=True)
                last["sgp"] = now
            if now - last["samp"] >= config.SAMPLES_REFRESH_SEC:
                _refresh_samples(); last["samp"] = now
            time.sleep(0.25)   # short sleep → responsive SIGTERM
    finally:
        with db.connect(read_only=True) as con:
            live = [r[0] for r in con.execute(
                "SELECT quote_id FROM live_quotes WHERE status='open'").fetchall()]
        for qid in live:
            try:
                gateway.cancel(qid)
                with db.connect() as con:
                    con.execute("UPDATE live_quotes SET status='cancelled', closed_at=? WHERE quote_id=?",
                                [datetime.now(timezone.utc), qid])
            except Exception:
                pass
        db.end_session(sid)
        print("=== shutdown complete ===", flush=True)


def cli():
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()
    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)
    main_loop(dry_run=args.dry_run)


if __name__ == "__main__":
    cli()
```
NOTE: two spots above are intentionally simplified and must be cleaned during execution: (a) the `_resolve_game_and_lines` `__import__` gymnastics — replace with a top-level `from kalshi_common.leg_types import _spread_line_from_legs, _total_line_from_legs, _parse_event_suffix, _MLB_CODE_TO_TEAM`; (b) `_risk_sweep_tick`'s `_SAMPLES_GEN_AT` walrus stub — track samples-generated-at as a module global set in `_refresh_samples()` (read `mlb_samples_meta.generated_at` like the taker does) and gate on it. These are flagged so the executing agent fixes them rather than copying verbatim.

- [ ] **Step 3: Smoke test (dry-run path, mocked source/gateway)**

```python
# kalshi_mlb_mm/tests/test_main_smoke.py
def test_main_imports():
    from kalshi_mlb_mm import main          # must import without error
    assert hasattr(main, "main_loop")

def test_discovery_tick_skips_out_of_scope(monkeypatch, tmp_path):
    import kalshi_mlb_mm.config as cfg
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "t.duckdb")
    import importlib, kalshi_mlb_mm.db as db; importlib.reload(db); db.init_database()
    from kalshi_mlb_mm import main
    class Src:
        def poll(self): return [{"id": "r1", "market_ticker": "OTHER-1", "contracts": 1}]
        def get_market(self, t): return {"mve_selected_legs": [{"market_ticker": "FOO"}]}
    class GW:
        def submit_quote(self, *a): raise AssertionError("should not submit out-of-scope")
    main._discovery_tick(Src(), GW(), dry_run=True)
    with db.connect(read_only=True) as con:
        d = con.execute("SELECT decision FROM quote_decisions").fetchone()
    assert d[0] == "skipped"
```

- [ ] **Step 4: Run — expect pass.** `pytest kalshi_mlb_mm/tests/test_main_smoke.py -v`

- [ ] **Step 5: Clean up the two flagged shortcuts** (per the NOTE) and re-run the full maker test suite: `pytest kalshi_mlb_mm/tests/ -v` → PASS.

- [ ] **Step 6: Commit**

```bash
git add kalshi_mlb_mm/main.py kalshi_mlb_mm/fairs.py kalshi_mlb_mm/tests/test_main_smoke.py
git commit -m "feat(mm): main loops — discovery/quote, last-look confirm, risk sweep, sgp"
```

---

## Task 9: Pre-launch validation, docs, and review

**Files:**
- Create: `kalshi_mlb_mm/README.md`
- Modify: `CLAUDE.md` (root project-structure list), `kalshi_mlb_rfq/README.md` (note the `kalshi_common` extraction)

- [ ] **Step 1: Maker venv + install**

```bash
cd kalshi_mlb_mm && python3 -m venv venv && ./venv/bin/pip install -r requirements.txt && cd ..
```

- [ ] **Step 2: Dry-run smoke against live RFQs (no writes to exchange)**

Copy `.env.example` → `.env`, fill credentials, then:
Run: `./kalshi_mlb_mm/venv/bin/python -m kalshi_mlb_mm.main --dry-run`
Expected: logs `dry_run_quote` / `skipped` rows in `quote_decisions`; **no** `POST /communications/quotes` calls. Let it run one full SGP cycle, Ctrl-C, then inspect:
`./kalshi_mlb_mm/venv/bin/python -c "import duckdb; print(duckdb.connect('kalshi_mlb_mm/kalshi_mlb_mm.duckdb', read_only=True).execute('SELECT decision, count(*) FROM quote_decisions GROUP BY 1').fetchall())"`
Confirm in-scope combos are detected and prices look sane (yes_bid+no_bid < 1, ~5% inside fair).

- [ ] **Step 3: Write `kalshi_mlb_mm/README.md`** — setup, quick start (`--dry-run` then live), the knobs table from spec §6, the defense hierarchy (§8), kill switch, and the open-items-to-verify (maker fee on first fill, quote-status fields, competitor read).

- [ ] **Step 4: Update root `CLAUDE.md`** — add a bullet under Project Structure for the maker bot + the `kalshi_common/` shared package; note in `kalshi_mlb_rfq/README.md` that shared math now lives in `kalshi_common` (taker uses shims).

- [ ] **Step 5: Executive-engineer pre-merge review** — run the `review-my-code` skill (or manual checklist from root CLAUDE.md) over `git diff main..HEAD`. Focus: taker shims preserve every public name; all DB connections use the `connect()` context manager; dry-run path makes zero exchange writes; no secrets in logs; the two flagged main.py shortcuts are resolved. Document ISSUES vs ACCEPTABLE RISKS; fix issues.

- [ ] **Step 6: Commit docs**

```bash
git add kalshi_mlb_mm/README.md CLAUDE.md kalshi_mlb_rfq/README.md
git commit -m "docs(mm): README, root CLAUDE.md project structure, taker shim note"
```

- [ ] **Step 7: First-fill fee verification (post-launch, tiny)** — after the first real confirmed fill, query `fills` and compare the actual fee charged (via `/portfolio/fills`) against `maker_fee_per_contract`. Adjust the maker-fee formula if Kalshi's rounding differs. (Spec §10.1.)

---

## Self-review notes (author)

- **Spec coverage:** §3 architecture → Tasks 1,2,8; §4 data flow → Task 8 (`_discovery_tick`/`_confirm_tick`); §5 combo decode → Task 3; §6 pricing+risk → Tasks 4,7; §7 schema → Task 2; §8 defenses → Tasks 7,8 (freshness/circuit-breaker/last-look) + §9 review; §9 testing → per-task TDD; §10 open items → Tasks 5,8 notes + Task 9 step 7; §11 delivery → task structure; §12 accepted risks → no gates built (by design), but #6 (require full 4-side devig) IS enforced in `_book_fairs`, and instrumentation (model/book/blend logged separately) is in the `fills`/`quote_decisions`/`live_quotes` schemas.
- **Type consistency:** `Quote(yes_bid,no_bid)` used consistently; `RestRFQSource.poll/get_market`, `RestQuoteGateway.submit_quote/confirm/cancel` match their call sites in `main.py`.
- **Known shortcuts flagged for the executor:** the two `main.py` simplifications in Task 8 (import gymnastics + samples-gen-at stub) are explicitly called out to fix, not copy.
