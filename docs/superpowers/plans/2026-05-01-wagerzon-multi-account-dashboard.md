# Wagerzon Multi-Account Dashboard Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Display available balances for all configured Wagerzon accounts on the MLB correlated parlay dashboard, and let the user choose which account each placement uses via a single global selector.

**Architecture:** A new account registry (`wagerzon_accounts.py`) discovers WZ accounts from env vars by suffix pattern. A consolidated auth module (`wagerzon_auth.py`) owns logged-in session caching per account. A balance fetcher (`wagerzon_balance.py`) returns `BalanceSnapshot` objects (parallel fan-out across accounts). `parlay_placer.place_parlays()` is refactored to take an explicit `account` parameter. The Flask server adds three new endpoints (balances list, last-used GET/POST) and threads `account` through `/api/place-parlay`. The R dashboard adds a header strip with per-account balance pills + a global selector.

**Tech Stack:** Python 3.10+ (stdlib `concurrent.futures`, `dataclasses`), `requests`, `python-dotenv`, DuckDB, Flask, R Shiny, vanilla JS, `pytest`.

**Worktree:** This plan is executed inside `.worktrees/wagerzon-multi-account-dashboard` on branch `feature/wagerzon-multi-account-dashboard`. The worktree is already created and the spec is already committed (commits `996233a`, `c874187`).

**Spec reference:** `docs/superpowers/specs/2026-05-01-wagerzon-multi-account-dashboard-design.md`

---

## File Plan

### New files

| Path | Responsibility | Approx LOC |
|---|---|---|
| `wagerzon_odds/wagerzon_accounts.py` | Account registry; env discovery; `WagerzonAccount` dataclass + `list_accounts()` + `get_account(label)` | ~80 |
| `wagerzon_odds/wagerzon_auth.py` | Logged-in session cache keyed by account label; ASP.NET form login; transparent re-login on 401 | ~120 |
| `wagerzon_odds/wagerzon_balance.py` | `BalanceSnapshot` dataclass; `fetch_available_balance(account)`; `fetch_all(accounts)` parallel fan-out | ~100 |
| `wagerzon_odds/test_wagerzon_accounts.py` | pytest unit tests for the registry (env discovery, ordering, error handling) | ~120 |
| `wagerzon_odds/test_wagerzon_balance.py` | pytest unit tests for the balance module (using requests-mock) | ~100 |
| `wagerzon_odds/CLAUDE.md` | Module-level orientation: registry pattern, auth helper, balance, parlay_placer's account parameter | ~60 |

### Modified files

| Path | Reason | Approx delta |
|---|---|---|
| `wagerzon_odds/parlay_placer.py` | Accept explicit `account: WagerzonAccount`; remove module-level env reads + `_get_session()`; call `wagerzon_auth.get_session(account)` | ~30 lines changed |
| `wagerzon_odds/parlay_pricer.py` | Replace `get_wz_session()` body to delegate to `wagerzon_auth.get_session(default_account)` | ~10 lines changed |
| `wagerzon_odds/test_parlay_placer.py` | Update existing tests to pass an account into placement calls | ~40 lines changed |
| `wagerzon_odds/README.md` | Add "Multi-account support" section | ~50 lines added |
| `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` | New endpoints; modified `/api/place-parlay`; schema migrations at startup; balance/session caches | ~250 lines added |
| `Answer Keys/MLB Dashboard/mlb_dashboard.R` | Header strip with pills + selector; updated `placeParlay()` JS to send `account` and refresh pill; insufficient-balance warning | ~200 lines added |
| `Answer Keys/CLAUDE.md` | Note new endpoints + `placed_parlays.account` column + `dashboard_settings` table + how to add a WZ account | ~30 lines added |

### Out of scope (per spec)

- `bet_logger/` — no changes (`--account c` already shipped via parallel worktree).
- `nfl_draft/` portal account selection.
- `/api/place-bet` (single-bet flow is tracking-only).

---

## Phase 0 — Prerequisites

### Task 0.1: Verify worktree state

**Files:** none (read-only check)

- [ ] **Step 1: Confirm worktree branch and clean state**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git branch --show-current
git status --short
```

Expected:
- `feature/wagerzon-multi-account-dashboard`
- Empty status (or only the in-progress plan file).

- [ ] **Step 2: Confirm Python version and pytest availability**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
python3 --version
python3 -c "import pytest, requests, duckdb, dotenv; print('ok')"
```

Expected: Python ≥ 3.10 and `ok`.

If `dotenv` or `pytest` is missing, install into the same venv used by `wagerzon_odds`:

```bash
pip install pytest requests-mock python-dotenv
```

- [ ] **Step 3: Confirm WagerzonC env vars exist**

```bash
grep -E '^WAGERZON[A-Z]?_USERNAME' /Users/callancapitolo/NFLWork/bet_logger/.env | cut -d= -f1
```

Expected: three lines — `WAGERZON_USERNAME`, `WAGERZONJ_USERNAME`, `WAGERZONC_USERNAME`.

If `WAGERZONC_USERNAME` is missing, **stop**. The parallel WagerzonC worktree must merge first.

---

## Phase 1 — Account registry (no integration yet)

### Task 1.1: Write the registry tests

**Files:**
- Create: `wagerzon_odds/test_wagerzon_accounts.py`

- [ ] **Step 1: Write the test file**

Create `wagerzon_odds/test_wagerzon_accounts.py` with:

```python
"""Unit tests for wagerzon_accounts registry."""
from __future__ import annotations

import pytest

from wagerzon_accounts import (
    WagerzonAccount,
    list_accounts,
    get_account,
    AccountNotFoundError,
)


@pytest.fixture
def env_two_accounts(monkeypatch):
    monkeypatch.setenv("WAGERZON_USERNAME", "primary_user")
    monkeypatch.setenv("WAGERZON_PASSWORD", "primary_pw")
    monkeypatch.setenv("WAGERZONJ_USERNAME", "j_user")
    monkeypatch.setenv("WAGERZONJ_PASSWORD", "j_pw")
    yield


@pytest.fixture
def env_three_accounts(env_two_accounts, monkeypatch):
    monkeypatch.setenv("WAGERZONC_USERNAME", "c_user")
    monkeypatch.setenv("WAGERZONC_PASSWORD", "c_pw")
    yield


@pytest.fixture
def env_clean(monkeypatch):
    """Strip every WAGERZON* env var so tests start from a clean slate."""
    for k in list(__import__("os").environ.keys()):
        if k.startswith("WAGERZON"):
            monkeypatch.delenv(k, raising=False)
    yield


def test_primary_only(env_clean, monkeypatch):
    monkeypatch.setenv("WAGERZON_USERNAME", "u")
    monkeypatch.setenv("WAGERZON_PASSWORD", "p")
    accounts = list_accounts()
    assert len(accounts) == 1
    assert accounts[0] == WagerzonAccount(label="Wagerzon", suffix="", username="u", password="p")


def test_two_accounts_primary_first(env_clean, env_two_accounts):
    accounts = list_accounts()
    labels = [a.label for a in accounts]
    assert labels == ["Wagerzon", "WagerzonJ"]


def test_three_accounts_alphabetical_after_primary(env_clean, env_three_accounts):
    accounts = list_accounts()
    labels = [a.label for a in accounts]
    assert labels == ["Wagerzon", "WagerzonC", "WagerzonJ"]


def test_skip_pair_missing_password(env_clean, monkeypatch, caplog):
    monkeypatch.setenv("WAGERZON_USERNAME", "u")
    monkeypatch.setenv("WAGERZON_PASSWORD", "p")
    monkeypatch.setenv("WAGERZONJ_USERNAME", "j_user")
    # No WAGERZONJ_PASSWORD on purpose
    accounts = list_accounts()
    assert [a.label for a in accounts] == ["Wagerzon"]
    assert "WAGERZONJ" in caplog.text


def test_skip_lowercase_suffix(env_clean, monkeypatch):
    monkeypatch.setenv("WAGERZON_USERNAME", "u")
    monkeypatch.setenv("WAGERZON_PASSWORD", "p")
    monkeypatch.setenv("WAGERZONj_USERNAME", "junk")
    monkeypatch.setenv("WAGERZONj_PASSWORD", "junk")
    accounts = list_accounts()
    assert [a.label for a in accounts] == ["Wagerzon"]


def test_skip_multichar_suffix(env_clean, monkeypatch):
    monkeypatch.setenv("WAGERZON_USERNAME", "u")
    monkeypatch.setenv("WAGERZON_PASSWORD", "p")
    monkeypatch.setenv("WAGERZONJX_USERNAME", "junk")
    monkeypatch.setenv("WAGERZONJX_PASSWORD", "junk")
    accounts = list_accounts()
    assert [a.label for a in accounts] == ["Wagerzon"]


def test_get_account_by_label(env_clean, env_two_accounts):
    acct = get_account("WagerzonJ")
    assert acct.label == "WagerzonJ"
    assert acct.username == "j_user"


def test_get_account_unknown_raises(env_clean, env_two_accounts):
    with pytest.raises(AccountNotFoundError):
        get_account("WagerzonZ")


def test_no_accounts_returns_empty_list(env_clean):
    assert list_accounts() == []


def test_password_not_in_repr(env_clean, monkeypatch):
    monkeypatch.setenv("WAGERZON_USERNAME", "u")
    monkeypatch.setenv("WAGERZON_PASSWORD", "supersecret123")
    acct = list_accounts()[0]
    assert "supersecret123" not in repr(acct)
    assert "supersecret123" not in str(acct)


def test_skip_pair_missing_username(env_clean, monkeypatch, caplog):
    import logging
    caplog.set_level(logging.WARNING)
    monkeypatch.setenv("WAGERZON_USERNAME", "u")
    monkeypatch.setenv("WAGERZON_PASSWORD", "p")
    monkeypatch.setenv("WAGERZONJ_PASSWORD", "j_pw")
    # No WAGERZONJ_USERNAME on purpose
    accounts = list_accounts()
    assert [a.label for a in accounts] == ["Wagerzon"]
    assert "WAGERZONJ" in caplog.text
    assert "WAGERZONJ_USERNAME" in caplog.text
```

- [ ] **Step 2: Run tests — they should FAIL (module doesn't exist yet)**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard/wagerzon_odds
python3 -m pytest test_wagerzon_accounts.py -v
```

Expected: `ERRORS / collection failure: cannot import name 'WagerzonAccount' from 'wagerzon_accounts'` or `ModuleNotFoundError`.

### Task 1.2: Implement the registry

**Files:**
- Create: `wagerzon_odds/wagerzon_accounts.py`

- [ ] **Step 1: Write the implementation**

Create `wagerzon_odds/wagerzon_accounts.py`:

```python
"""Registry of configured Wagerzon accounts.

Discovers accounts from environment variables matching the pattern
`WAGERZON{SUFFIX}_USERNAME` / `WAGERZON{SUFFIX}_PASSWORD`, where SUFFIX is
either empty (primary account) or a single uppercase ASCII letter (e.g. J, C).

Usage:
    from wagerzon_accounts import list_accounts, get_account
    for acct in list_accounts():
        print(acct.label, acct.username)
    j = get_account("WagerzonJ")
"""
from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

try:
    from dotenv import load_dotenv
    _ENV_PATH = Path(__file__).resolve().parent.parent / "bet_logger" / ".env"
    if _ENV_PATH.exists():
        load_dotenv(_ENV_PATH)
except ImportError:
    pass

logger = logging.getLogger(__name__)

_SUFFIX_RE = re.compile(r"^WAGERZON([A-Z]?)_USERNAME$")
_PW_RE = re.compile(r"^WAGERZON([A-Z]?)_PASSWORD$")


class AccountNotFoundError(KeyError):
    """Raised by get_account when the label is not in the current registry."""


@dataclass(frozen=True)
class WagerzonAccount:
    """A Wagerzon login config discovered from env vars.

    suffix: "" for the primary account, otherwise a single uppercase letter
            (e.g. "J", "C").
    label:  "Wagerzon" + suffix; stable identifier used by /api/place-parlay etc.
    """
    label: str
    suffix: str
    username: str
    password: str = field(repr=False)


def list_accounts() -> list[WagerzonAccount]:
    """Discover WZ accounts from environment.

    Returns the primary account (suffix '') first if present, then alphabetical
    by suffix. Pairs missing one half (USERNAME or PASSWORD) are skipped with a
    warning.
    """
    found: dict[str, WagerzonAccount] = {}
    for key, val in os.environ.items():
        m = _SUFFIX_RE.match(key)
        if not m:
            continue
        suffix = m.group(1)
        pw_key = f"WAGERZON{suffix}_PASSWORD"
        password = os.environ.get(pw_key)
        if not val:
            logger.warning(
                "Wagerzon account WAGERZON%s skipped: %s is empty",
                suffix or "(primary)", key,
            )
            continue
        if not password:
            logger.warning(
                "Wagerzon account WAGERZON%s skipped: %s set but %s missing/empty",
                suffix or "(primary)", key, pw_key,
            )
            continue
        label = "Wagerzon" + suffix
        found[suffix] = WagerzonAccount(
            label=label, suffix=suffix, username=val, password=password,
        )

    # Second pass: warn on orphan PASSWORD vars (USERNAME missing for that suffix).
    for key in os.environ:
        m = _PW_RE.match(key)
        if not m:
            continue
        suffix = m.group(1)
        if suffix in found:
            continue
        un_key = f"WAGERZON{suffix}_USERNAME"
        if not os.environ.get(un_key):
            logger.warning(
                "Wagerzon account WAGERZON%s skipped: %s set but %s missing/empty",
                suffix or "(primary)", key, un_key,
            )

    ordered: list[WagerzonAccount] = []
    if "" in found:
        ordered.append(found.pop(""))
    for suffix in sorted(found.keys()):
        ordered.append(found[suffix])
    return ordered


def get_account(label: str) -> WagerzonAccount:
    """Look up an account by label (e.g. 'WagerzonJ'). Raises AccountNotFoundError."""
    for acct in list_accounts():
        if acct.label == label:
            return acct
    raise AccountNotFoundError(label)
```

- [ ] **Step 2: Run tests — they should PASS**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard/wagerzon_odds
python3 -m pytest test_wagerzon_accounts.py -v
```

Expected: 11 tests pass, 0 failures.

If a test fails, do NOT proceed. Read the failure, fix the implementation, rerun.

- [ ] **Step 3: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add wagerzon_odds/wagerzon_accounts.py wagerzon_odds/test_wagerzon_accounts.py
git commit -m "$(cat <<'EOF'
feat(wagerzon): add account registry module

Discovers Wagerzon accounts from env vars by suffix pattern. Returns
primary first, then alphabetical. Pairs with one half missing are
skipped with a warning. Includes pytest unit tests covering primary-only,
multi-account ordering, malformed suffix filtering, and lookup-by-label.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 2 — Auth helper

### Task 2.1: Inspect the existing auth code paths

**Files:** none (read-only)

- [ ] **Step 1: Read the two existing auth implementations and write notes**

Open both for reference (you'll consolidate them):

- `wagerzon_odds/parlay_placer.py` lines 21–33 (env load), lines 131–165 (`_get_session()`)
- `wagerzon_odds/parlay_pricer.py` lines 34–50 (`get_wz_session()`), and `scraper_v2.login()` it delegates to

The new `wagerzon_auth.py` must:
1. Take a `WagerzonAccount` (not env vars).
2. Login via the same ASP.NET-form mechanism (`__VIEWSTATE`, `__EVENTVALIDATION`, etc.) used by `_get_session()`.
3. Cache one `requests.Session` per account label.
4. Detect 401 / login-page-bounce and re-login transparently.

### Task 2.2: Write the auth helper tests

**Files:**
- Create: `wagerzon_odds/test_wagerzon_auth.py`

- [ ] **Step 1: Write the test file (mocking HTTP)**

Create `wagerzon_odds/test_wagerzon_auth.py`:

```python
"""Unit tests for wagerzon_auth (login + session caching)."""
from __future__ import annotations

import re
import pytest
import requests_mock

from wagerzon_accounts import WagerzonAccount
import wagerzon_auth


@pytest.fixture
def acct():
    return WagerzonAccount(label="Wagerzon", suffix="", username="u", password="p")


@pytest.fixture(autouse=True)
def reset_cache():
    wagerzon_auth.clear_session_cache()
    yield
    wagerzon_auth.clear_session_cache()


LOGIN_PAGE_HTML = (
    '<html><form>'
    '<input name="__VIEWSTATE" value="vs1"/>'
    '<input name="__VIEWSTATEGENERATOR" value="vsg1"/>'
    '<input name="__EVENTVALIDATION" value="ev1"/>'
    '</form></html>'
)
LOGGED_IN_URL = wagerzon_auth.WAGERZON_BASE_URL.rstrip("/") + "/NewSchedule.aspx"


def test_first_call_logs_in(acct):
    with requests_mock.Mocker() as m:
        m.get(wagerzon_auth.WAGERZON_BASE_URL, text=LOGIN_PAGE_HTML)
        m.post(wagerzon_auth.WAGERZON_BASE_URL, text="ok",
               headers={"Location": LOGGED_IN_URL}, status_code=302)
        sess = wagerzon_auth.get_session(acct)
        assert isinstance(sess, requests_mock.Adapter) is False  # got a Session
        assert m.call_count == 2  # GET login page + POST login form


def test_second_call_returns_cached(acct):
    with requests_mock.Mocker() as m:
        m.get(wagerzon_auth.WAGERZON_BASE_URL, text=LOGIN_PAGE_HTML)
        m.post(wagerzon_auth.WAGERZON_BASE_URL, text="ok",
               headers={"Location": LOGGED_IN_URL}, status_code=302)
        s1 = wagerzon_auth.get_session(acct)
        s2 = wagerzon_auth.get_session(acct)
        assert s1 is s2
        assert m.call_count == 2  # no second login


def test_already_logged_in_skips_form_post(acct):
    """If GET base URL bounces straight to NewSchedule, no POST required."""
    with requests_mock.Mocker() as m:
        # The cookie jar already has a valid session — GET redirects.
        m.get(wagerzon_auth.WAGERZON_BASE_URL,
              status_code=302, headers={"Location": LOGGED_IN_URL})
        m.get(LOGGED_IN_URL, text="hello")
        sess = wagerzon_auth.get_session(acct)
        # Only the redirect chain happened; no form POST.
        post_calls = [r for r in m.request_history if r.method == "POST"]
        assert len(post_calls) == 0


def test_two_accounts_get_separate_sessions():
    a = WagerzonAccount(label="Wagerzon",  suffix="",  username="u1", password="p1")
    b = WagerzonAccount(label="WagerzonJ", suffix="J", username="u2", password="p2")
    with requests_mock.Mocker() as m:
        m.get(wagerzon_auth.WAGERZON_BASE_URL, text=LOGIN_PAGE_HTML)
        m.post(wagerzon_auth.WAGERZON_BASE_URL, text="ok",
               headers={"Location": LOGGED_IN_URL}, status_code=302)
        s_a = wagerzon_auth.get_session(a)
        s_b = wagerzon_auth.get_session(b)
        assert s_a is not s_b


def test_clear_session_cache_forces_relogin(acct):
    with requests_mock.Mocker() as m:
        m.get(wagerzon_auth.WAGERZON_BASE_URL, text=LOGIN_PAGE_HTML)
        m.post(wagerzon_auth.WAGERZON_BASE_URL, text="ok",
               headers={"Location": LOGGED_IN_URL}, status_code=302)
        wagerzon_auth.get_session(acct)
        wagerzon_auth.clear_session_cache(label=acct.label)
        wagerzon_auth.get_session(acct)
        # Two full login sequences (GET + POST) = 4 calls
        assert m.call_count == 4
```

- [ ] **Step 2: Run tests — they should FAIL (module not yet implemented)**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard/wagerzon_odds
python3 -m pytest test_wagerzon_auth.py -v
```

Expected: ImportError or `AttributeError` on `wagerzon_auth.get_session`.

### Task 2.3: Implement the auth helper

**Files:**
- Create: `wagerzon_odds/wagerzon_auth.py`

- [ ] **Step 1: Write the implementation**

Create `wagerzon_odds/wagerzon_auth.py`:

```python
"""Wagerzon session/auth helper, keyed by WagerzonAccount.

Consolidates the ASP.NET form-login mechanism previously duplicated in
parlay_placer._get_session() and parlay_pricer.get_wz_session().

Public API:
    get_session(account)         -> logged-in requests.Session
    clear_session_cache(label?)  -> drop one or all cached sessions
"""
from __future__ import annotations

import logging
import re
import threading
from typing import Optional

import requests
from requests.adapters import HTTPAdapter

from wagerzon_accounts import WagerzonAccount

logger = logging.getLogger(__name__)

WAGERZON_BASE_URL = "https://backend.wagerzon.com"
DEFAULT_TIMEOUT = 15
POOL_SIZE = 8

# label -> Session
_SESSION_CACHE: dict[str, requests.Session] = {}
_LOCK = threading.RLock()


def get_session(account: WagerzonAccount) -> requests.Session:
    """Return a logged-in Session for `account`, logging in on first call."""
    with _LOCK:
        cached = _SESSION_CACHE.get(account.label)
        if cached is not None:
            return cached
        session = _build_session()
        _login(session, account)
        _SESSION_CACHE[account.label] = session
        return session


def clear_session_cache(label: Optional[str] = None) -> None:
    """Drop cached sessions. Pass a label to drop one; omit to drop all."""
    with _LOCK:
        if label is None:
            _SESSION_CACHE.clear()
        else:
            _SESSION_CACHE.pop(label, None)


def _build_session() -> requests.Session:
    s = requests.Session()
    adapter = HTTPAdapter(pool_connections=POOL_SIZE, pool_maxsize=POOL_SIZE)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    return s


def _login(session: requests.Session, account: WagerzonAccount) -> None:
    """Run the ASP.NET form-login flow for account."""
    resp = session.get(WAGERZON_BASE_URL, timeout=DEFAULT_TIMEOUT,
                       allow_redirects=True)
    resp.raise_for_status()
    if "NewSchedule" in resp.url or "Welcome" in resp.url:
        # Server-side cookie reuse already logged us in — nothing to do.
        return

    html = resp.text
    fields: dict[str, str] = {}
    for name in ("__VIEWSTATE", "__VIEWSTATEGENERATOR", "__EVENTVALIDATION",
                 "__EVENTTARGET", "__EVENTARGUMENT"):
        m = re.search(rf'(?:name|id)="{name}"[^>]*value="([^"]*)"', html)
        if m:
            fields[name] = m.group(1)
    fields["Account"] = account.username
    fields["Password"] = account.password
    fields["BtnSubmit"] = ""

    resp = session.post(WAGERZON_BASE_URL, data=fields,
                        timeout=DEFAULT_TIMEOUT, allow_redirects=True)
    resp.raise_for_status()
    if "NewSchedule" not in resp.url and "Welcome" not in resp.url:
        # Login appeared to succeed at the HTTP level but the response
        # didn't land us on a logged-in page. Log + continue; callers
        # will get a 401-equivalent on the next real request.
        logger.warning(
            "Wagerzon login for %s did not redirect to a logged-in page (url=%s)",
            account.label, resp.url,
        )
```

- [ ] **Step 2: Run tests — should PASS**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard/wagerzon_odds
python3 -m pytest test_wagerzon_auth.py -v
```

Expected: 6 tests pass (5 from the original list + the new password-not-logged test).

- [ ] **Step 3: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add wagerzon_odds/wagerzon_auth.py wagerzon_odds/test_wagerzon_auth.py
git commit -m "$(cat <<'EOF'
feat(wagerzon): consolidate auth into wagerzon_auth module

New module owns the ASP.NET form-login flow and a per-account session
cache. Replaces the two ad-hoc auth implementations in parlay_placer
(_get_session) and parlay_pricer (get_wz_session). Sessions are keyed
by account label and protected by an RLock for thread safety.

Includes pytest tests using requests-mock covering first-call login,
cache hit, redirect-when-already-logged-in, per-account separation, and
cache invalidation.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 3 — Balance fetcher

### Task 3.1: Wagerzon balance endpoint (already discovered)

**Files:** none (reference for Tasks 3.2 + 3.3)

No browser/devtools work required — the endpoint was found in `wagerzon_odds/recon_place_parlay.json` (an existing recon dump from the parlay-placer build).

| Field | Value |
|---|---|
| HTTP method | `GET` |
| URL | `https://backend.wagerzon.com/wager/PlayerInfoHelper.aspx` |
| Required query params | none |
| Required POST body | none (it's a GET) |
| Required headers beyond cookies | none — session cookies from `wagerzon_auth.get_session()` carry identity |
| Response body | JSON: `{"result": {"AvailBalance": "1,245.32 ", "CurrentBalance": "1,300 ", "AmountAtRisk": "715 ", "RealAvailBalance": "...", "CreditLimit": "...", "BonusPoints": 0.0, "FreePlayAmount": "0 ", "Player": "<acct>", "Password": "<plaintext>", "ErrorCode": {}, "ErrorMsg": "", ...}}` |
| Path to available balance | `result.AvailBalance` (string with thousand-separator comma + trailing space; parse with `float(s.strip().replace(",", ""))`) |
| Path to cash balance | `result.CurrentBalance` (same string format) |

**Why `AvailBalance` and not `RealAvailBalance`:** `AvailBalance` = `CurrentBalance - AmountAtRisk` (the wagerable amount from your own cash). `RealAvailBalance` adds `CreditLimit` on top — not what we want for the gating decision.

**Security gotcha (CRITICAL):** the response body also includes `result.Player` and `result.Password` (the user's WZ password in plaintext). The balance parser MUST extract only the two numeric fields above and never log the raw response body. Tests in Task 3.2 include an explicit assertion that the password value never appears in `caplog`.

### Task 3.2: Write the balance module tests

**Files:**
- Create: `wagerzon_odds/test_wagerzon_balance.py`

> **Note to executor:** the test cases below assume the balance endpoint is JSON. If Task 3.1 found HTML scraping is required instead, adapt the parser tests but keep the same `BalanceSnapshot` shape and the same public functions.

- [ ] **Step 1: Write the test file (replace `BALANCE_URL` and `MOCK_RESPONSE` with values from Task 3.1)**

Create `wagerzon_odds/test_wagerzon_balance.py`:

```python
"""Unit tests for wagerzon_balance."""
from __future__ import annotations

from datetime import datetime, timezone
import pytest
import requests_mock

from wagerzon_accounts import WagerzonAccount
import wagerzon_auth
import wagerzon_balance


@pytest.fixture(autouse=True)
def clean_caches():
    wagerzon_auth.clear_session_cache()
    yield
    wagerzon_auth.clear_session_cache()


@pytest.fixture
def acct():
    return WagerzonAccount(label="Wagerzon", suffix="", username="u", password="p")


_LOGGED_IN_URL = wagerzon_auth.WAGERZON_BASE_URL.rstrip("/") + "/NewSchedule.aspx"


def _mock_login(m):
    """Helper: mock the login endpoints so get_session() succeeds."""
    LOGIN_HTML = (
        '<html><form>'
        '<input name="__VIEWSTATE" value="x"/>'
        '<input name="__VIEWSTATEGENERATOR" value="x"/>'
        '<input name="__EVENTVALIDATION" value="x"/>'
        '</form></html>'
    )
    m.get(wagerzon_auth.WAGERZON_BASE_URL, text=LOGIN_HTML)
    m.post(wagerzon_auth.WAGERZON_BASE_URL, text="",
           headers={"Location": _LOGGED_IN_URL},
           status_code=302)
    m.get(_LOGGED_IN_URL, text="logged in")


# Real WZ balance endpoint and a sample of its response shape.
BALANCE_URL = "https://backend.wagerzon.com/wager/PlayerInfoHelper.aspx"
MOCK_BALANCE_RESPONSE = {
    "result": {
        "AmountAtRisk": "715 ",
        "AvailBalance": "1,245.32 ",       # the gating number we display
        "BonusPoints": 0.0000,
        "CreditLimit": "2,000 ",
        "CurrentBalance": "1,300.00 ",     # the "cash" we expose for tooltip/debug
        "FreePlayAmount": "0 ",
        "RealAvailBalance": "3,245.32 ",
        # NOTE: real responses also include "Player" and "Password" fields
        # (yes, the password in plaintext). The parser MUST extract only the
        # numeric balance fields and never log the raw response body.
        "ErrorCode": {},
        "ErrorMsg": "",
    }
}


def test_fetch_returns_snapshot(acct):
    with requests_mock.Mocker() as m:
        _mock_login(m)
        m.get(BALANCE_URL, json=MOCK_BALANCE_RESPONSE)
        snap = wagerzon_balance.fetch_available_balance(acct)
        assert snap.label == "Wagerzon"
        assert snap.available == 1245.32
        assert snap.cash == 1300.00
        assert snap.error is None
        assert snap.fetched_at.tzinfo == timezone.utc


def test_fetch_does_not_log_password_from_response(acct, caplog):
    """The WZ balance response includes the user's password in plaintext.
    The fetcher must never log the raw response body."""
    import logging
    response_with_password = {
        "result": {
            "AvailBalance": "100 ", "CurrentBalance": "100 ",
            "Player": "MYACCT", "Password": "supersecret-do-not-log",
            "ErrorCode": {}, "ErrorMsg": "",
        }
    }
    caplog.set_level(logging.DEBUG)
    with requests_mock.Mocker() as m:
        _mock_login(m)
        m.get(BALANCE_URL, json=response_with_password)
        wagerzon_balance.fetch_available_balance(acct)
    assert "supersecret-do-not-log" not in caplog.text


def test_fetch_timeout_returns_error_snapshot(acct):
    import requests
    with requests_mock.Mocker() as m:
        _mock_login(m)
        m.get(BALANCE_URL, exc=requests.exceptions.ConnectTimeout)
        snap = wagerzon_balance.fetch_available_balance(acct)
        assert snap.error == "timeout"
        assert snap.available is None


def test_fetch_5xx_returns_wz_error_snapshot(acct):
    with requests_mock.Mocker() as m:
        _mock_login(m)
        m.get(BALANCE_URL, status_code=503, text="oops")
        snap = wagerzon_balance.fetch_available_balance(acct)
        assert snap.error == "wz_error"
        assert snap.available is None


def test_fetch_401_triggers_relogin_then_retries_once(acct):
    with requests_mock.Mocker() as m:
        _mock_login(m)
        # First balance call: 401. Second: success.
        m.get(BALANCE_URL, [
            {"status_code": 401, "text": "unauthorized"},
            {"json": MOCK_BALANCE_RESPONSE, "status_code": 200},
        ])
        snap = wagerzon_balance.fetch_available_balance(acct)
        assert snap.error is None
        assert snap.available == 1245.32


def test_fetch_all_runs_in_parallel():
    a = WagerzonAccount(label="Wagerzon",  suffix="",  username="u1", password="p1")
    b = WagerzonAccount(label="WagerzonJ", suffix="J", username="u2", password="p2")
    with requests_mock.Mocker() as m:
        _mock_login(m)
        m.get(BALANCE_URL, json=MOCK_BALANCE_RESPONSE)
        snaps = wagerzon_balance.fetch_all([a, b])
        labels = sorted(s.label for s in snaps)
        assert labels == ["Wagerzon", "WagerzonJ"]
```

- [ ] **Step 2: Run tests — they should FAIL (module not yet written)**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard/wagerzon_odds
python3 -m pytest test_wagerzon_balance.py -v
```

Expected: ImportError on `wagerzon_balance`.

### Task 3.3: Implement the balance module

**Files:**
- Create: `wagerzon_odds/wagerzon_balance.py`

> **Endpoint:** `GET https://backend.wagerzon.com/wager/PlayerInfoHelper.aspx`. Auth = the same session cookies established by `wagerzon_auth.get_session()`. Response is JSON shaped `{"result": {"AvailBalance": "1,245.32 ", "CurrentBalance": "1,300 ", "AmountAtRisk": "715 ", "RealAvailBalance": "...", "CreditLimit": "...", "Player": ..., "Password": ..., ...}}`. Numeric fields arrive as strings with thousand-separator commas and a trailing space — parse with `float(s.strip().replace(",", ""))`. **Security: the response body includes the user's password in plaintext — the parser must extract only `AvailBalance` + `CurrentBalance` and never log the raw response.**

- [ ] **Step 1: Write the implementation**

Create `wagerzon_odds/wagerzon_balance.py`:

```python
"""Fetch available balance for one or more Wagerzon accounts.

Public API:
    fetch_available_balance(account)  -> BalanceSnapshot
    fetch_all(accounts)               -> list[BalanceSnapshot]   (parallel)
"""
from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Optional

import requests

import wagerzon_auth
from wagerzon_accounts import WagerzonAccount

logger = logging.getLogger(__name__)

_BALANCE_URL = "https://backend.wagerzon.com/wager/PlayerInfoHelper.aspx"

_TIMEOUT_SECONDS = 5
_MAX_PARALLEL = 8


@dataclass
class BalanceSnapshot:
    label: str
    available: Optional[float]
    cash: Optional[float]
    fetched_at: datetime
    error: Optional[str]  # one of: timeout, auth_failed, wz_error, parse_error, None


def fetch_available_balance(account: WagerzonAccount) -> BalanceSnapshot:
    """Fetch the available balance. Network/auth failures return a snapshot
    with `error` set rather than raising."""
    try:
        snap = _fetch_once(account)
        return snap
    except _UnauthorizedError:
        # Force re-login and try once more.
        wagerzon_auth.clear_session_cache(label=account.label)
        try:
            return _fetch_once(account)
        except _UnauthorizedError:
            return _error_snapshot(account, "auth_failed")
        except requests.exceptions.Timeout:
            return _error_snapshot(account, "timeout")
        except requests.exceptions.RequestException as e:
            logger.warning("balance fetch retry failed for %s: %s", account.label, e)
            return _error_snapshot(account, "wz_error")


def fetch_all(accounts: list[WagerzonAccount]) -> list[BalanceSnapshot]:
    """Fan out balance fetches in parallel. Order follows the input list."""
    if not accounts:
        return []
    workers = min(_MAX_PARALLEL, len(accounts))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(fetch_available_balance, accounts))
    return results


# ----- internals -----

class _UnauthorizedError(Exception):
    pass


def _fetch_once(account: WagerzonAccount) -> BalanceSnapshot:
    session = wagerzon_auth.get_session(account)
    try:
        # Adapt method (GET/POST), URL, and any required headers/body
        # to whatever Task 3.1 discovered.
        resp = session.get(_BALANCE_URL, timeout=_TIMEOUT_SECONDS,
                           allow_redirects=False)
    except requests.exceptions.Timeout:
        return _error_snapshot(account, "timeout")
    except requests.exceptions.RequestException as e:
        logger.warning("balance fetch transport error for %s: %s", account.label, e)
        return _error_snapshot(account, "wz_error")

    if resp.status_code == 401 or _looks_like_login_redirect(resp):
        raise _UnauthorizedError()
    if resp.status_code >= 500:
        return _error_snapshot(account, "wz_error")
    if resp.status_code != 200:
        logger.warning("balance fetch unexpected status for %s: %s",
                       account.label, resp.status_code)
        return _error_snapshot(account, "wz_error")

    try:
        available, cash = _parse_balance_response(resp)
    except Exception as e:
        logger.warning("balance parse error for %s: %s", account.label, e)
        return _error_snapshot(account, "parse_error")

    return BalanceSnapshot(
        label=account.label,
        available=available,
        cash=cash,
        fetched_at=datetime.now(timezone.utc),
        error=None,
    )


def _parse_balance_response(resp: requests.Response) -> tuple[float, Optional[float]]:
    """Parse the WZ PlayerInfoHelper response.

    Response shape:
        {"result": {"AvailBalance": "1,245.32 ", "CurrentBalance": "1,300 ", ...,
                    "Player": "ACCT", "Password": "<plaintext>", ...}}

    Numbers come as strings with comma thousand-separators and a trailing
    space. We extract ONLY AvailBalance (gating value) and CurrentBalance
    (cash, optional). The response body also contains the user's password
    in plaintext — never log the full body and never return any field
    other than the two below.
    """
    data = resp.json()
    result = data.get("result") or {}

    raw_avail = result.get("AvailBalance")
    if raw_avail is None:
        raise ValueError("balance response missing 'AvailBalance'")
    available = _parse_money_string(raw_avail)

    raw_cash = result.get("CurrentBalance")
    cash = _parse_money_string(raw_cash) if raw_cash is not None else None

    return available, cash


def _parse_money_string(s) -> float:
    """Parse a Wagerzon-formatted money string like '1,245.32 ' or 0.0 (number)."""
    if isinstance(s, (int, float)):
        return float(s)
    return float(str(s).strip().replace(",", ""))


def _looks_like_login_redirect(resp: requests.Response) -> bool:
    loc = resp.headers.get("Location", "")
    return "Default.aspx" in loc or "Login" in loc


def _error_snapshot(account: WagerzonAccount, error_code: str) -> BalanceSnapshot:
    return BalanceSnapshot(
        label=account.label,
        available=None,
        cash=None,
        fetched_at=datetime.now(timezone.utc),
        error=error_code,
    )
```

- [ ] **Step 2: Run tests — they should PASS**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard/wagerzon_odds
python3 -m pytest test_wagerzon_balance.py -v
```

Expected: 6 tests pass (5 from the original list + the new password-not-logged test).

- [ ] **Step 3: Smoke test against the real WZ endpoint**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard/wagerzon_odds
python3 -c "
from wagerzon_accounts import list_accounts
from wagerzon_balance import fetch_all
for snap in fetch_all(list_accounts()):
    print(snap)
"
```

Expected: one `BalanceSnapshot` per configured account, each with a real `available` number and `error=None`.

If any account returns `error="parse_error"`: the response shape didn't match the parser. Compare the live response (from a one-off curl with the session cookies) to the shape documented in Task 3.1, fix `_parse_balance_response`, rerun.

If any returns `error="auth_failed"`: the auth flow is broken for that account — likely a credential typo in env, or WZ requires a CAPTCHA. Investigate before proceeding.

- [ ] **Step 4: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add wagerzon_odds/wagerzon_balance.py wagerzon_odds/test_wagerzon_balance.py
git commit -m "$(cat <<'EOF'
feat(wagerzon): add balance fetcher

Adds fetch_available_balance(account) and fetch_all(accounts) returning
BalanceSnapshot objects (label, available, cash, fetched_at, error).
Wraps the Wagerzon balance endpoint discovered in Task 3.1. Network and
auth failures return a snapshot with error set, never raise. fetch_all
runs in parallel via ThreadPoolExecutor (5s per-account timeout).

Includes pytest tests using requests-mock for success, timeout, 5xx,
401-then-retry, and parallel fan-out.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 4 — Refactor parlay_placer to take an explicit account

### Task 4.1: Read existing tests and understand the surface

**Files:** none (read-only)

- [ ] **Step 1: Read `wagerzon_odds/test_parlay_placer.py` cover-to-cover**

Identify which tests call `place_parlays()` directly. These will need an account argument.

### Task 4.2: Refactor parlay_placer to require an account

**Files:**
- Modify: `wagerzon_odds/parlay_placer.py`
- Modify: `wagerzon_odds/test_parlay_placer.py`

- [ ] **Step 1: Replace the module-level env load + `_get_session()`**

In `wagerzon_odds/parlay_placer.py`:

Find lines 21–33 (the dotenv load + try/except around it) and **leave them in place** — `wagerzon_auth` reuses the same env pattern via `wagerzon_accounts`, but other parts of `parlay_placer` may still read env vars (logging, debugging). Remove only the parts that won't be needed.

Find `_CACHED_SESSION = None` (the module-level cache variable, search for it). Delete it.

Find `_get_session()` (lines ~131–165). Delete the entire function.

At the top of the file, after existing imports, add:

```python
import wagerzon_auth
from wagerzon_accounts import WagerzonAccount
```

- [ ] **Step 2: Update `place_parlays()` signature**

Find the current signature (line ~379):

```python
def place_parlays(specs: list[ParlaySpec]) -> list[PlacementResult]:
```

Change to:

```python
def place_parlays(
    specs: list[ParlaySpec],
    account: WagerzonAccount,
) -> list[PlacementResult]:
```

In the docstring, add: `account: which Wagerzon login to use; required.`

- [ ] **Step 3: Replace every `_get_session()` call inside `place_parlays()` (and its helpers)**

Search the file:

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
grep -n "_get_session" wagerzon_odds/parlay_placer.py
```

For each call, replace with `wagerzon_auth.get_session(account)`. Helpers that previously took no args may need an explicit `account` parameter — pass it through from `place_parlays`.

After replacement, `_get_session` should appear nowhere in the file.

- [ ] **Step 4: Update existing tests**

In `wagerzon_odds/test_parlay_placer.py`, every test that calls `place_parlays(specs)` must be updated to pass an account.

Add a test fixture at the top of the file:

```python
import pytest
from wagerzon_accounts import WagerzonAccount

@pytest.fixture
def primary_acct():
    return WagerzonAccount(
        label="Wagerzon", suffix="",
        username="test_user", password="test_pw",
    )
```

For each test calling `place_parlays(specs)`, change to:

```python
def test_xxx(primary_acct, ...):
    ...
    result = place_parlays(specs, primary_acct)
```

If a test patches `parlay_placer._get_session`, change to patch `wagerzon_auth.get_session` instead.

- [ ] **Step 5: Run the full parlay_placer test suite**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard/wagerzon_odds
python3 -m pytest test_parlay_placer.py -v
```

Expected: all 22 tests pass.

If any fail: read the failure carefully. The most common cause is a helper function that still calls `_get_session` (now undefined) — search and fix.

- [ ] **Step 6: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add wagerzon_odds/parlay_placer.py wagerzon_odds/test_parlay_placer.py
git commit -m "$(cat <<'EOF'
refactor(wagerzon): parlay_placer takes explicit account parameter

place_parlays() now requires an account: WagerzonAccount argument and
uses wagerzon_auth.get_session(account) instead of the deleted internal
_get_session() / module-level cache. Tests updated to pass a primary
account fixture. No behavior change for the primary account.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task 4.3: Update parlay_pricer to use shared auth (cleanup)

**Files:**
- Modify: `wagerzon_odds/parlay_pricer.py`

- [ ] **Step 1: Replace the body of `get_wz_session()`**

In `wagerzon_odds/parlay_pricer.py` (lines ~34–50), replace the function with:

```python
import wagerzon_auth
from wagerzon_accounts import list_accounts


def get_wz_session() -> requests.Session:
    """Return an authenticated Wagerzon session for the primary account.

    Pricing is account-agnostic (the `RiskWin=2` trick skips balance
    validation), so we always use the primary login. Callers that need
    a specific account should call wagerzon_auth.get_session(acct) directly.
    """
    accounts = list_accounts()
    if not accounts:
        raise RuntimeError("no Wagerzon accounts configured (env vars missing)")
    return wagerzon_auth.get_session(accounts[0])
```

Remove the now-unused `from scraper_v2 import login` import if it's used nowhere else in the file.

- [ ] **Step 2: Smoke-test the pricer still works**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard/wagerzon_odds
python3 -c "
from parlay_pricer import get_wz_session
s = get_wz_session()
print('session created:', type(s).__name__)
"
```

Expected: `session created: Session`.

- [ ] **Step 3: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add wagerzon_odds/parlay_pricer.py
git commit -m "$(cat <<'EOF'
refactor(wagerzon): parlay_pricer.get_wz_session delegates to wagerzon_auth

Replaces the bespoke login flow (via scraper_v2.login) with a thin
wrapper around wagerzon_auth.get_session(primary_account). Pricing
remains account-agnostic — RiskWin=2 skips balance validation server-
side, so we always use the primary account.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 5 — Backend / Flask endpoints

### Task 5.1: Add schema migrations at server startup

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`

- [ ] **Step 1: Locate the existing `placed_parlays` CREATE TABLE block**

Search for it:

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
grep -n "CREATE TABLE.*placed_parlays\|CREATE TABLE IF NOT EXISTS placed_parlays" "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
```

Expected: a hit around line 187. Read lines 187–220 for context.

- [ ] **Step 2: Add `account` to the CREATE TABLE column list**

Inside the `CREATE TABLE IF NOT EXISTS placed_parlays (...)` block, add (alphabetically near other text columns, but in practice append before the closing `)` since DuckDB doesn't care about order for new tables):

```sql
account TEXT,
```

- [ ] **Step 3: Add an idempotent ALTER for existing databases**

Immediately after the `CREATE TABLE` execute, add:

```python
# Multi-account support: add `account` column to existing placed_parlays.
try:
    con.execute("ALTER TABLE placed_parlays ADD COLUMN account TEXT")
except duckdb.CatalogException:
    pass  # Column already exists.
```

Make sure `duckdb` is imported at the top of the file (it should be).

- [ ] **Step 4: Add the `dashboard_settings` table**

In the same startup-migration block (where `placed_parlays` is created), add:

```python
con.execute("""
    CREATE TABLE IF NOT EXISTS dashboard_settings (
        key   TEXT PRIMARY KEY,
        value TEXT NOT NULL,
        updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
""")
```

- [ ] **Step 5: Restart the server and verify schema**

```bash
cd /Users/callancapitolo/NFLWork
# Kill any running dashboard server first.
pkill -f mlb_dashboard_server.py || true
# Start it (background; redirect logs).
"Answer Keys/MLB Dashboard"/venv/bin/python "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" > /tmp/mlb_dashboard_server.log 2>&1 &
sleep 4
# Verify the schema.
"Answer Keys/MLB Dashboard"/venv/bin/python -c "
import duckdb
con = duckdb.connect('Answer Keys/MLB Dashboard/mlb_dashboard.duckdb')
print('placed_parlays columns:')
for row in con.execute(\"PRAGMA table_info('placed_parlays')\").fetchall():
    print(' ', row[1], row[2])
print('dashboard_settings columns:')
for row in con.execute(\"PRAGMA table_info('dashboard_settings')\").fetchall():
    print(' ', row[1], row[2])
"
```

Expected:
- `placed_parlays` shows an `account` column of type `VARCHAR` (DuckDB's `TEXT`).
- `dashboard_settings` exists with `key`, `value`, `updated_at` columns.

- [ ] **Step 6: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): schema migrations for multi-account support

Adds placed_parlays.account column (idempotent ALTER) and a generic
dashboard_settings(key,value,updated_at) table for last-used-account
persistence. Both run at server startup.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task 5.2: Add `/api/wagerzon/balances` endpoint with cache

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`

- [ ] **Step 1: Add the imports near the top of the file**

Find the existing `from parlay_pricer import (...)` block (around line 31–37). After it, add:

```python
from wagerzon_accounts import list_accounts as wz_list_accounts, get_account as wz_get_account, AccountNotFoundError
import wagerzon_balance
```

- [ ] **Step 2: Add the in-memory balance cache near other module-level state**

Search for an existing module-level cache or constants section. Add:

```python
import threading
from datetime import datetime, timezone

# Last balance snapshot per account label. Powers the "stale" tag.
_BALANCE_CACHE: dict[str, wagerzon_balance.BalanceSnapshot] = {}
_BALANCE_CACHE_LOCK = threading.RLock()


def _balance_snapshot_to_json(snap: wagerzon_balance.BalanceSnapshot) -> dict:
    now = datetime.now(timezone.utc)
    stale_seconds = int((now - snap.fetched_at).total_seconds())
    return {
        "label": snap.label,
        "available": snap.available,
        "cash": snap.cash,
        "fetched_at": snap.fetched_at.isoformat().replace("+00:00", "Z"),
        "error": snap.error,
        "stale_seconds": stale_seconds,
    }


def _refresh_one_balance(account) -> wagerzon_balance.BalanceSnapshot:
    snap = wagerzon_balance.fetch_available_balance(account)
    with _BALANCE_CACHE_LOCK:
        # On error, prefer keeping the prior successful snapshot (with its
        # original `available` value) but bump fetched_at to "now-but-failed"
        # — actually: keep the OLD snapshot's available + cash, but flag it.
        prior = _BALANCE_CACHE.get(snap.label)
        if snap.error and prior is not None and prior.error is None:
            snap = wagerzon_balance.BalanceSnapshot(
                label=prior.label,
                available=prior.available,
                cash=prior.cash,
                fetched_at=prior.fetched_at,  # keep the OLD timestamp so stale_seconds is honest
                error=snap.error,
            )
        _BALANCE_CACHE[snap.label] = snap
    return snap


def _refresh_all_balances() -> list[wagerzon_balance.BalanceSnapshot]:
    accounts = wz_list_accounts()
    if not accounts:
        return []
    # fetch_all already parallelises.
    snaps = wagerzon_balance.fetch_all(accounts)
    with _BALANCE_CACHE_LOCK:
        for snap in snaps:
            prior = _BALANCE_CACHE.get(snap.label)
            if snap.error and prior is not None and prior.error is None:
                snap = wagerzon_balance.BalanceSnapshot(
                    label=prior.label,
                    available=prior.available,
                    cash=prior.cash,
                    fetched_at=prior.fetched_at,
                    error=snap.error,
                )
            _BALANCE_CACHE[snap.label] = snap
    return list(_BALANCE_CACHE.values())
```

- [ ] **Step 3: Add the endpoint**

Find an existing `@app.route("/api/...")` for context. Add a new route (anywhere among the other routes):

```python
@app.route("/api/wagerzon/balances", methods=["GET"])
def api_wagerzon_balances():
    snaps = _refresh_all_balances()
    return jsonify({"balances": [_balance_snapshot_to_json(s) for s in snaps]})
```

If `jsonify` isn't already imported at the top, add it: `from flask import jsonify`.

- [ ] **Step 4: Test the endpoint**

```bash
# Restart the server first.
pkill -f mlb_dashboard_server.py || true
cd /Users/callancapitolo/NFLWork
"Answer Keys/MLB Dashboard"/venv/bin/python "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" > /tmp/mlb_dashboard_server.log 2>&1 &
sleep 4
curl -s http://localhost:8083/api/wagerzon/balances | python3 -m json.tool
```

Expected: JSON like

```json
{
  "balances": [
    {"label": "Wagerzon", "available": 1245.32, "cash": 1300.00,
     "fetched_at": "...Z", "error": null, "stale_seconds": 0},
    {"label": "WagerzonC", "available": 2010.55, "cash": 2010.55, ...},
    {"label": "WagerzonJ", "available": 890.10, "cash": 890.10, ...}
  ]
}
```

If `error` is populated for any account: check `/tmp/mlb_dashboard_server.log` for the underlying message; most likely cause is a credential issue or the balance endpoint URL still being a placeholder.

- [ ] **Step 5: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): GET /api/wagerzon/balances + balance cache

Returns one BalanceSnapshot per discovered Wagerzon account, fetched in
parallel. Errors are surfaced per-account in the response (HTTP always
200). Cache preserves the prior successful snapshot's value when a fetch
fails — UI then shows it with stale_seconds indicating freshness.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task 5.3: Add `/api/wagerzon/last-used` GET + POST

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`

- [ ] **Step 1: Add helper functions for the settings table**

Near the other helper functions (above the routes), add:

```python
DB_PATH = ...  # already defined elsewhere in the file; reuse the existing constant


def _get_setting(key: str) -> str | None:
    con = duckdb.connect(str(DB_PATH))
    try:
        row = con.execute(
            "SELECT value FROM dashboard_settings WHERE key = ?", [key]
        ).fetchone()
        return row[0] if row else None
    finally:
        con.close()


def _set_setting(key: str, value: str) -> None:
    con = duckdb.connect(str(DB_PATH))
    try:
        con.execute(
            """
            INSERT INTO dashboard_settings (key, value, updated_at)
            VALUES (?, ?, CURRENT_TIMESTAMP)
            ON CONFLICT (key)
            DO UPDATE SET value = excluded.value, updated_at = CURRENT_TIMESTAMP
            """,
            [key, value],
        )
    finally:
        con.close()


def _get_default_account_label() -> str | None:
    """Return last-used label if it's still in the registry, else fall back
    to the registry primary, else None."""
    accounts = wz_list_accounts()
    if not accounts:
        return None
    valid = {a.label for a in accounts}
    saved = _get_setting("wagerzon_last_used")
    if saved in valid:
        return saved
    return accounts[0].label
```

> **Note:** if `DB_PATH` isn't already a module-level constant in the file, search for how the file currently opens the DB (likely `duckdb.connect("...mlb_dashboard.duckdb")`) and reuse that pattern.

- [ ] **Step 2: Add the endpoints**

```python
@app.route("/api/wagerzon/last-used", methods=["GET"])
def api_wagerzon_last_used_get():
    return jsonify({"label": _get_default_account_label()})


@app.route("/api/wagerzon/last-used", methods=["POST"])
def api_wagerzon_last_used_post():
    data = request.get_json(silent=True) or {}
    label = data.get("label")
    if not isinstance(label, str) or not label:
        return jsonify({"error": "label required"}), 400
    valid = {a.label for a in wz_list_accounts()}
    if label not in valid:
        return jsonify({"error": f"unknown label: {label}"}), 400
    _set_setting("wagerzon_last_used", label)
    return jsonify({"ok": True})
```

If `request` isn't imported at top, add: `from flask import request`.

- [ ] **Step 3: Test both endpoints**

```bash
pkill -f mlb_dashboard_server.py || true
cd /Users/callancapitolo/NFLWork
"Answer Keys/MLB Dashboard"/venv/bin/python "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" > /tmp/mlb_dashboard_server.log 2>&1 &
sleep 4
# Should return {"label": "Wagerzon"} (the registry primary, since nothing's been saved yet).
curl -s http://localhost:8083/api/wagerzon/last-used
echo
# Persist a choice.
curl -s -X POST -H "Content-Type: application/json" \
  -d '{"label": "WagerzonJ"}' http://localhost:8083/api/wagerzon/last-used
echo
# Should now return {"label": "WagerzonJ"}.
curl -s http://localhost:8083/api/wagerzon/last-used
echo
# Should 400.
curl -s -X POST -H "Content-Type: application/json" \
  -d '{"label": "WagerzonZ"}' http://localhost:8083/api/wagerzon/last-used
echo
```

Expected outputs (in order): `{"label": "Wagerzon"}`, `{"ok": true}`, `{"label": "WagerzonJ"}`, `{"error": "unknown label: WagerzonZ"}`.

- [ ] **Step 4: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): GET/POST /api/wagerzon/last-used persistence

Stores the user's last selected account in the dashboard_settings table
so the dashboard re-loads onto the same account across sessions. Helper
_get_default_account_label() falls back to the registry primary when no
saved value exists or when the saved label is no longer configured.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task 5.4: Modify `/api/place-parlay` to accept `account`

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`

- [ ] **Step 1: Re-read the handler at lines 1258–1392**

Identify:
- Where `parlay_placer.place_parlays([spec])` is called (line ~1345).
- Where `_upsert_placed_parlay(result, row)` is called (line ~1375).
- Where the request body is parsed (the `request.get_json()` line near the top of the handler).

- [ ] **Step 2: Parse `account` from the request body**

In the handler, where the body is parsed, add:

```python
account_label = data.get("account")
if not account_label:
    return jsonify({"error": "account required"}), 400
try:
    account = wz_get_account(account_label)
except AccountNotFoundError:
    return jsonify({"error": f"unknown account: {account_label}"}), 400
```

`data` here is whatever variable holds the parsed JSON body — match the existing name (likely `data` or `body`).

- [ ] **Step 3: Pass account into `place_parlays`**

Change:

```python
results = parlay_placer.place_parlays([spec])
```

to:

```python
results = parlay_placer.place_parlays([spec], account)
```

- [ ] **Step 4: Persist `account` on the placed_parlays row**

Find `_upsert_placed_parlay` (the function it delegates to). Modify its INSERT/UPDATE to write `account = ?`.

If `_upsert_placed_parlay` doesn't currently take `account`, change its signature:

```python
def _upsert_placed_parlay(result, row, account_label: str) -> None:
    ...
```

Inside, ensure the INSERT statement (around line 1314-1340 originally) includes `account` in the column list and `?` in the values, and the parameter list passes `account_label`.

Update the caller in the route handler:

```python
_upsert_placed_parlay(result, row, account.label)
```

> **Important:** the existing INSERT uses `'placing'` as the literal status. The same row may later be UPDATEd by another path. Search for all `UPDATE placed_parlays` statements and add `account = ?` to those that should preserve it (most likely none — `account` is set on insert and never changes).

- [ ] **Step 5: Trigger a balance refresh after successful placement**

In the success branch of the handler (where `status == "placed"`), after persisting, add:

```python
# Refresh the balance for the placed-on account so the UI shows the new value.
try:
    fresh = _refresh_one_balance(account)
    balance_after = _balance_snapshot_to_json(fresh)
except Exception as e:
    app.logger.warning("post-placement balance refresh failed: %s", e)
    balance_after = None
```

In the JSON response (whatever the handler returns), add the `balance_after` field:

```python
response_payload["balance_after"] = balance_after
```

(Match the existing response variable name — adapt accordingly.)

- [ ] **Step 6: Smoke-test with a tiny placement**

> **DO NOT actually place a $100 bet for testing.** Test the auth/route plumbing only — either against a stub or by sending a request that you expect to fail at WZ price-validation (which is a clean failure, not a placed bet).

```bash
pkill -f mlb_dashboard_server.py || true
cd /Users/callancapitolo/NFLWork
"Answer Keys/MLB Dashboard"/venv/bin/python "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" > /tmp/mlb_dashboard_server.log 2>&1 &
sleep 4
# Missing account → 400
curl -s -X POST -H "Content-Type: application/json" \
  -d '{"parlay_hash": "fakehash"}' http://localhost:8083/api/place-parlay
echo
# Unknown account → 400
curl -s -X POST -H "Content-Type: application/json" \
  -d '{"parlay_hash": "fakehash", "account": "WagerzonZ"}' http://localhost:8083/api/place-parlay
echo
```

Expected: both return 400 with the appropriate error string. Server logs should NOT show a stack trace.

- [ ] **Step 7: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): /api/place-parlay accepts account parameter

The endpoint now requires an `account` field in the body, resolves it
through the wagerzon_accounts registry, threads it into
parlay_placer.place_parlays(), and writes the resolved label onto the
placed_parlays row. After a successful placement, the response includes
a fresh balance_after snapshot for the placed-on account so the UI can
update its pill without a separate fetch.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 6 — Dashboard UI (R + JS)

### Task 6.1: Locate the right insertion point for the header strip

**Files:** none (read-only)

- [ ] **Step 1: Read `Answer Keys/MLB Dashboard/mlb_dashboard.R` for layout structure**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
grep -n "tabsetPanel\|navbarPage\|fluidPage\|titlePanel\|headerPanel\|tags\$h1\|tags\$header" "Answer Keys/MLB Dashboard/mlb_dashboard.R" | head -30
```

Identify where the top-of-page UI is constructed. The new header strip needs to render above the existing tab panel so it's visible on every tab.

### Task 6.2: Render balance pills + selector in the header

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R`

- [ ] **Step 1: Add a new `div` near the top of the UI**

In the UI block (NOT the server block), insert above the `tabsetPanel(...)` (or whatever the top-level layout container is):

```r
tags$div(
  id = "wz-account-bar",
  style = paste(
    "display:flex; align-items:center; gap:14px;",
    "padding:8px 16px; background:#f7f7f7;",
    "border-bottom:1px solid #ddd; font-family:sans-serif;"
  ),
  tags$div(id = "wz-account-pills", style = "display:flex; gap:8px;"),
  tags$button(
    id = "wz-refresh-btn", type = "button",
    style = "border:none; background:transparent; cursor:pointer; font-size:18px;",
    title = "Refresh balances",
    HTML("&#x21bb;")  # ↻
  ),
  tags$span("Placing on:", style = "margin-left:auto; font-weight:bold;"),
  tags$select(id = "wz-account-select", style = "padding:4px 8px; font-size:14px;")
)
```

- [ ] **Step 2: Add the JS that populates the bar on page load**

Find an existing `tags$script(HTML(...))` block. Append (or add a new one near the bottom of the UI, before the closing of the layout):

```r
tags$script(HTML("
(function() {
  var BAR_ID    = 'wz-account-pills';
  var SELECT_ID = 'wz-account-select';
  var REFRESH_BTN_ID = 'wz-refresh-btn';

  // Currently-selected account label, source of truth for the dashboard.
  window.WZ_SELECTED_ACCOUNT = null;
  window.WZ_BALANCES = {};   // label -> snapshot

  function fmtMoney(n) {
    if (n === null || n === undefined) return '—';
    return '$' + Number(n).toLocaleString('en-US', {minimumFractionDigits:2, maximumFractionDigits:2});
  }

  function renderPills() {
    var bar = document.getElementById(BAR_ID);
    if (!bar) return;
    bar.innerHTML = '';
    Object.keys(WZ_BALANCES).forEach(function(label) {
      var snap = WZ_BALANCES[label];
      var pill = document.createElement('span');
      pill.className = 'wz-pill';
      pill.dataset.label = label;
      var stale = snap.error && snap.stale_seconds > 600;
      pill.style.cssText = 'padding:4px 10px; border-radius:14px; ' +
        'background:' + (stale ? '#fdd' : '#fff') + '; ' +
        'border:' + (label === window.WZ_SELECTED_ACCOUNT ? '2px solid #333' : '1px solid #bbb') + ';' +
        'font-size:13px;';
      var text = label + ': ' + fmtMoney(snap.available);
      if (snap.error) {
        text += ' ⚠';  // ⚠
        var sub = snap.stale_seconds < 60
          ? snap.stale_seconds + 's'
          : Math.floor(snap.stale_seconds/60) + 'm';
        text += ' (stale ' + sub + ' ago)';
      }
      pill.textContent = text;
      bar.appendChild(pill);
    });
  }

  function renderSelect(orderedLabels) {
    var sel = document.getElementById(SELECT_ID);
    if (!sel) return;
    sel.innerHTML = '';
    orderedLabels.forEach(function(label) {
      var opt = document.createElement('option');
      opt.value = label;
      opt.textContent = label;
      if (label === window.WZ_SELECTED_ACCOUNT) opt.selected = true;
      sel.appendChild(opt);
    });
    if (orderedLabels.length === 0) {
      var opt = document.createElement('option');
      opt.value = '';
      opt.textContent = 'No Wagerzon accounts configured';
      opt.disabled = true; opt.selected = true;
      sel.appendChild(opt);
      sel.disabled = true;
    } else if (orderedLabels.length === 1) {
      sel.disabled = true;
    }
  }

  function refreshBalances() {
    return fetch('/api/wagerzon/balances')
      .then(function(r) { return r.json(); })
      .then(function(payload) {
        var orderedLabels = [];
        WZ_BALANCES = {};
        (payload.balances || []).forEach(function(snap) {
          WZ_BALANCES[snap.label] = snap;
          orderedLabels.push(snap.label);
        });
        renderPills();
        renderSelect(orderedLabels);
        recomputeAllInsufficiencyWarnings();
      });
  }

  function loadLastUsed() {
    return fetch('/api/wagerzon/last-used')
      .then(function(r) { return r.json(); })
      .then(function(payload) {
        window.WZ_SELECTED_ACCOUNT = payload.label || null;
      });
  }

  function persistSelection(label) {
    return fetch('/api/wagerzon/last-used', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({label: label})
    });
  }

  function recomputeAllInsufficiencyWarnings() {
    // Implemented in Task 6.4 once placeParlay learns about risk.
    if (typeof window._wzRecomputeWarnings === 'function') {
      window._wzRecomputeWarnings();
    }
  }

  // Public hook used by placeParlay() after a successful placement
  // to update the pill without a full re-fetch.
  window.wzApplyBalanceAfter = function(label, snap) {
    if (!snap) return;
    WZ_BALANCES[label] = snap;
    renderPills();
    recomputeAllInsufficiencyWarnings();
  };

  document.addEventListener('DOMContentLoaded', function() {
    var sel = document.getElementById(SELECT_ID);
    sel.addEventListener('change', function() {
      window.WZ_SELECTED_ACCOUNT = sel.value;
      persistSelection(sel.value);
      renderPills();  // re-highlight the selected pill
      recomputeAllInsufficiencyWarnings();
    });
    document.getElementById(REFRESH_BTN_ID).addEventListener('click', refreshBalances);

    loadLastUsed().then(refreshBalances);
  });
})();
"))
```

- [ ] **Step 2: Restart and visually verify**

```bash
pkill -f mlb_dashboard_server.py || true
pkill -f mlb_dashboard.R || true
cd /Users/callancapitolo/NFLWork
"Answer Keys/MLB Dashboard"/venv/bin/python "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" > /tmp/mlb_dashboard_server.log 2>&1 &
sleep 3
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R" > /tmp/mlb_dashboard_R.log 2>&1 &
sleep 6
echo "Open http://localhost:8083 in a browser."
```

Visually confirm:
- Pills render with three balances at the top of the page.
- Dropdown is selectable and shows three options.
- Clicking a pill or refresh ↻ re-fetches.
- Switching the dropdown highlights the new pill.

- [ ] **Step 3: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): header bar with balance pills + account selector

Adds a sticky strip above the tab panel showing one pill per configured
Wagerzon account (label + available balance), plus a global dropdown
that drives every subsequent placement. Stale fetches show a ⚠ tag and
flip red after 10 minutes. Selection persists to dashboard_settings via
the /api/wagerzon/last-used endpoint.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task 6.3: Update `placeParlay()` JS to send the selected account

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R`

- [ ] **Step 1: Locate `placeParlay()` (lines 3689–3744)**

```bash
grep -n "function placeParlay" "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

- [ ] **Step 2: Add account to the request body, handle balance_after in the response**

Find the `body: JSON.stringify({ parlay_hash: hash })` line and replace it with:

```javascript
body: JSON.stringify({
  parlay_hash: hash,
  account: window.WZ_SELECTED_ACCOUNT
})
```

Just above that, add a guard (so we don't POST without an account):

```javascript
if (!window.WZ_SELECTED_ACCOUNT) {
  showToast("No Wagerzon account selected", "error");
  btn.disabled = false;
  btn.textContent = originalText;
  return;
}
```

In the success branch (`if (result.status === "placed") {`), after the existing pill update, add:

```javascript
if (result.balance_after && window.wzApplyBalanceAfter) {
  window.wzApplyBalanceAfter(result.balance_after.label, result.balance_after);
}
```

In the success toast, append the account label:

```javascript
showToast(ticket
  ? ("Placed on " + window.WZ_SELECTED_ACCOUNT + " (#" + ticket + ")")
  : ("Placed on " + window.WZ_SELECTED_ACCOUNT),
  "success");
```

- [ ] **Step 3: Manual smoke test (no real bet)**

Restart dashboard. In the browser console:

```javascript
window.WZ_SELECTED_ACCOUNT = "WagerzonZ";  // intentionally invalid
// Click any Place button; should toast "Place failed: unknown account: WagerzonZ"
```

Reset:

```javascript
window.WZ_SELECTED_ACCOUNT = "Wagerzon";
```

- [ ] **Step 4: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): placeParlay sends selected account + updates pill

placeParlay() now sends the selected account in the request body, blocks
the placement if no account is selected, and applies the balance_after
snapshot returned by the server to the corresponding pill. Success toast
includes the account label so the user can confirm at a glance which
account a ticket landed on.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task 6.4: Add transient insufficient-balance warning per parlay

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R`

- [ ] **Step 1: Identify how parlay risk is rendered today**

```bash
grep -n "data-risk\|riskAmount\|riskValue\|recommended_size\|kelly_bet" "Answer Keys/MLB Dashboard/mlb_dashboard.R" | head -30
```

Find where the Place button is rendered for each parlay row. The `risk` value (the dollars at stake) needs to be on the button as a `data-risk` attribute so JS can read it. If it isn't already, add `` `data-risk="' + risk + '"' `` to the rendered HTML.

- [ ] **Step 2: Add the warning element next to each Place button**

Where the Place button is built, add immediately after it (in the same parent cell):

```r
# In the R/HTML that builds each row's actions:
tags$span(class = "wz-insufficient-warning", style = "margin-left:8px; color:#c0392b; font-size:12px;")
```

For dynamic JS-rendered rows, do the equivalent in JS — find where row HTML is appended.

- [ ] **Step 3: Implement the recompute function in the existing JS block**

In the same `tags$script(HTML(...))` block from Task 6.2, replace the placeholder hook with a real implementation. Add (within the IIFE, before `document.addEventListener`):

```javascript
window._wzRecomputeWarnings = function() {
  var label = window.WZ_SELECTED_ACCOUNT;
  var snap  = label ? WZ_BALANCES[label] : null;
  var available = snap ? snap.available : null;

  document.querySelectorAll('button[data-hash][data-risk]').forEach(function(btn) {
    var risk = parseFloat(btn.dataset.risk);
    var warn = btn.parentElement.querySelector('.wz-insufficient-warning');
    if (!warn) return;

    if (available === null || isNaN(risk)) {
      // Suppress the warning when we don't know the balance OR risk.
      warn.textContent = '';
      return;
    }
    if (risk > available) {
      warn.textContent = '⚠ insufficient on ' + label +
        ' ($' + Number(available).toFixed(2) + ' < $' +
        risk.toFixed(2) + ' risk)';
    } else {
      warn.textContent = '';
    }
  });
};
```

- [ ] **Step 4: Trigger recompute when parlay tables re-render**

If the dashboard re-renders parlay tables on a Shiny event (e.g., new parlays appearing), call `_wzRecomputeWarnings()` afterwards. Find the existing post-render hook (search for `Shiny.addCustomMessageHandler` or DOM observer code) and add a call. If no hook exists, install a `MutationObserver` watching the parlay-table container that calls `_wzRecomputeWarnings()` on changes.

- [ ] **Step 5: Visual verification**

Restart the dashboard. Pick an account whose balance is small (or temporarily edit `WZ_BALANCES.Wagerzon.available = 0` in the browser console). Confirm the warning text appears next to Place buttons whose risk exceeds the balance, and disappears when you switch to an account with sufficient balance.

- [ ] **Step 6: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
feat(mlb-dashboard): transient insufficient-balance warning per parlay

Renders a small ⚠ message beside each Place button when the selected
account's available balance is less than the parlay's risk. Recomputes
on selector change, balance refresh, post-placement balance update, and
parlay-table re-render. Warning is suppressed when balance is unknown
(no successful fetch ever) so the UI doesn't show misleading text.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 7 — Documentation

### Task 7.1: Update `wagerzon_odds/README.md`

**Files:**
- Modify: `wagerzon_odds/README.md`

- [ ] **Step 1: Add a new "Multi-account support" section**

After the existing "Auth" section, append:

```markdown
## Multi-account support

The Wagerzon modules now support multiple accounts simultaneously.

### Account discovery

Accounts are discovered from environment variables in `bet_logger/.env`
matching the pattern `WAGERZON{SUFFIX}_USERNAME` / `WAGERZON{SUFFIX}_PASSWORD`,
where `{SUFFIX}` is empty (primary) or a single uppercase letter:

```
WAGERZON_USERNAME=...     WAGERZON_PASSWORD=...      # primary, label = "Wagerzon"
WAGERZONJ_USERNAME=...    WAGERZONJ_PASSWORD=...     # label = "WagerzonJ"
WAGERZONC_USERNAME=...    WAGERZONC_PASSWORD=...     # label = "WagerzonC"
```

Add a new account by adding two env vars and restarting any consumer
process (the MLB dashboard server). No code changes required.

### Modules

- `wagerzon_accounts.py` — registry. `list_accounts()`, `get_account(label)`.
- `wagerzon_auth.py` — logged-in `requests.Session` cache keyed by account label. `get_session(account)`.
- `wagerzon_balance.py` — `fetch_available_balance(account)` and `fetch_all(accounts)` returning `BalanceSnapshot`.

### parlay_placer

`place_parlays(specs, account)` now requires an explicit account argument.
The previous module-level cache and `_get_session()` helper are removed.

### Dashboard integration

The MLB correlated parlay dashboard (`Answer Keys/MLB Dashboard/`) exposes:
- `GET /api/wagerzon/balances` — all accounts' current balances.
- `GET/POST /api/wagerzon/last-used` — persisted selector value.
- `POST /api/place-parlay` — accepts an `account` field that resolves through the registry.
```

- [ ] **Step 2: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add wagerzon_odds/README.md
git commit -m "$(cat <<'EOF'
docs(wagerzon): document multi-account discovery + new modules

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task 7.2: Create `wagerzon_odds/CLAUDE.md`

**Files:**
- Create: `wagerzon_odds/CLAUDE.md`

- [ ] **Step 1: Write the file**

Create `wagerzon_odds/CLAUDE.md`:

```markdown
# wagerzon_odds module

This directory holds the Wagerzon-specific scrapers, pricer, placer, and
shared auth/account/balance helpers used by the MLB correlated parlay
dashboard.

## Quick map

- **Scraping odds:** `scraper_v2.py`, `scraper_specials.py` write to
  `wagerzon_odds/wagerzon.duckdb`.
- **Pricing parlays:** `parlay_pricer.py::get_parlay_price()` calls the
  WZ `ConfirmWagerHelper` endpoint with `RiskWin=2` (skips balance
  validation server-side).
- **Placing parlays:** `parlay_placer.py::place_parlays(specs, account)`
  calls `PostWagerMultipleHelper`. **Requires** an explicit account.
- **Balance:** `wagerzon_balance.py::fetch_available_balance(account)`
  returns a `BalanceSnapshot`.
- **Account registry:** `wagerzon_accounts.py` discovers accounts by env
  var suffix (see `README.md` for the pattern).
- **Auth:** `wagerzon_auth.py::get_session(account)` caches one
  logged-in `requests.Session` per account label, in-memory only.

## Adding a new Wagerzon account

1. Pick a single uppercase suffix letter (e.g. `D`).
2. Add `WAGERZOND_USERNAME` and `WAGERZOND_PASSWORD` to `bet_logger/.env`.
3. Add the same pair to `bet_logger/.env.example` (placeholder values).
4. Restart any process that uses the registry (MLB dashboard server).
5. Optional: add `--account d` support to `bet_logger/scraper_wagerzon.py`
   if you want the bet-logging side to also pick up tickets from the
   account (separate concern; see `bet_logger/CLAUDE.md`).

## Pitfalls

- The auth flow scrapes `__VIEWSTATE` etc. from the Default.aspx form.
  If WZ changes the login page, both `wagerzon_auth.py` and any
  dependent test fixtures need updating.
- The balance endpoint URL is hardcoded in `wagerzon_balance.py` — if
  WZ moves it, update the constant + retest.
- Pricing always uses the primary account regardless of which account
  will eventually place the bet. This is intentional (`RiskWin=2`) but
  means a primary-account login failure breaks pricing for all accounts.
```

- [ ] **Step 2: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add wagerzon_odds/CLAUDE.md
git commit -m "$(cat <<'EOF'
docs(wagerzon): add module-level CLAUDE.md

Quick map of the new account registry / auth / balance helpers, the
multi-account-aware place_parlays signature, and a step-by-step for
adding a new Wagerzon account.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

### Task 7.3: Update `Answer Keys/CLAUDE.md`

**Files:**
- Modify: `Answer Keys/CLAUDE.md`

- [ ] **Step 1: Add a Wagerzon multi-account section**

Append to the file:

```markdown
## MLB Dashboard — Wagerzon multi-account

The MLB correlated parlay dashboard (port 8083) supports multiple
Wagerzon accounts. Configuration lives in `bet_logger/.env` and is
discovered by `wagerzon_odds/wagerzon_accounts.py`. See
`wagerzon_odds/CLAUDE.md` for the discovery pattern and how to add an
account.

### New endpoints (mlb_dashboard_server.py)

- `GET /api/wagerzon/balances` — list all account snapshots.
- `GET/POST /api/wagerzon/last-used` — persist the selector value.
- `POST /api/place-parlay` — now requires `{"account": "<label>"}` in
  the body in addition to `parlay_hash`.

### Schema

- `placed_parlays.account` — TEXT, label of the WZ account the parlay
  was placed on. NULL for pre-feature rows (treat as primary).
- `dashboard_settings(key, value, updated_at)` — generic key/value
  preferences. Currently only `wagerzon_last_used`.
```

- [ ] **Step 2: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git add "Answer Keys/CLAUDE.md"
git commit -m "$(cat <<'EOF'
docs(answer-keys): document Wagerzon multi-account dashboard endpoints + schema

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 8 — Pre-merge review and cleanup

### Task 8.1: Executive engineer review

**Files:** none (read-only diff review)

- [ ] **Step 1: Run the full diff against main**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
git fetch origin main
git diff origin/main..HEAD --stat
git diff origin/main..HEAD | wc -l
```

Capture the file list and total line delta — useful for the review writeup.

- [ ] **Step 2: Run the review checklist**

For each item below, check the diff and write FINDING (or ACCEPTABLE-AS-IS).

| Check | Method |
|---|---|
| **Data integrity** | Does any new INSERT/UPDATE create duplicate rows? Does the `placed_parlays.account` write happen exactly once per placement (not on retries)? |
| **Resource safety** | Every `duckdb.connect(...)` paired with `try / finally con.close()`? Any `requests.Session` created without ever being closed? |
| **Edge cases** | Off-season (no games, no parlays) — does the dashboard still load with the header bar? Empty registry — does the server start without crashing? Auth failure on dashboard load — does the UI degrade gracefully? |
| **Dead code** | Any `_get_session` references left? Any unused imports? Any test fixtures unused after the refactor? |
| **Log/disk hygiene** | Any new `print()` calls that should be `logger.info/warning`? Any log line that prints credentials? |
| **Security** | Search the diff for the literal string `password` — every match should be either a config key name or the `WagerzonAccount.password` field. No print statements logging the value. |

```bash
git diff origin/main..HEAD | grep -i "password" | grep -v "WagerzonAccount\|env_key\|getenv\|os\.environ"
```

Expected: no hits.

- [ ] **Step 3: Write the review findings to the conversation**

Document each FINDING with:
- File + line
- What's wrong
- Suggested fix

Then split into ISSUES TO FIX vs ACCEPTABLE RISKS.

- [ ] **Step 4: Fix all ISSUES TO FIX**

For each issue, edit, re-run tests, commit:

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-multi-account-dashboard
python3 -m pytest wagerzon_odds/ -v
git add ...
git commit -m "fix(...): description"
```

### Task 8.2: End-to-end smoke test

**Files:** none (manual test)

- [ ] **Step 1: Restart everything from a clean state**

```bash
cd /Users/callancapitolo/NFLWork
pkill -f mlb_dashboard_server.py || true
pkill -f mlb_dashboard.R || true
sleep 1
"Answer Keys/MLB Dashboard"/venv/bin/python "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" > /tmp/mlb_dashboard_server.log 2>&1 &
sleep 4
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R" > /tmp/mlb_dashboard_R.log 2>&1 &
sleep 6
```

- [ ] **Step 2: Open the dashboard and walk through this checklist**

Open `http://localhost:8083`. Verify:

- [ ] All three account pills render with real numbers.
- [ ] Selector dropdown shows all three labels; defaults to whatever was last used (or primary on first run).
- [ ] Refresh button (↻) re-fetches.
- [ ] Switching the dropdown highlights the new pill and persists (refresh the page — the selection should stick).
- [ ] At least one parlay is showing; its Place button is enabled.
- [ ] Switch to an account with low balance; if any parlay's risk exceeds it, the warning appears beside Place.

- [ ] **Step 3: Place ONE small real parlay (e.g. $1 risk)**

> **Talk to the user first** — placing a real bet is a real-money action.

After placement:
- [ ] Toast says `"Placed on <label> (#<ticket>)"`.
- [ ] The pill for the placed-on account updates immediately to show the lower balance.
- [ ] `placed_parlays` table has a new row with `account = '<label>'`:

```bash
"Answer Keys/MLB Dashboard"/venv/bin/python -c "
import duckdb
con = duckdb.connect('Answer Keys/MLB Dashboard/mlb_dashboard.duckdb')
print(con.execute(\"SELECT parlay_hash, status, account, ticket_number FROM placed_parlays ORDER BY updated_at DESC LIMIT 3\").fetchall())
"
```

- [ ] **Step 4: Repeat the placement on a different account**

Switch the dropdown to a second account. Place another small parlay. Verify the new row in `placed_parlays` has the correct `account` value.

- [ ] **Step 5: Stop the dashboard processes for now (don't merge yet)**

```bash
pkill -f mlb_dashboard_server.py || true
pkill -f mlb_dashboard.R || true
```

### Task 8.3: User approval to merge

**Files:** none

- [ ] **Step 1: Summarize for the user**

Write up: what shipped, what was tested (the smoke-test results above), what's NOT tested (e.g., 401 auth-failure path may not have been exercised), and any ACCEPTABLE-AS-IS items from the review.

- [ ] **Step 2: Wait for explicit user approval**

Per CLAUDE.md: never merge to `main` without explicit user approval. Do not proceed to the next step until the user says "merge it" or equivalent.

### Task 8.4: Merge to main and clean up worktree

**Files:** none

- [ ] **Step 1: Merge to main with `--no-ff`**

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git pull origin main
git merge --no-ff feature/wagerzon-multi-account-dashboard -m "$(cat <<'EOF'
Merge feature/wagerzon-multi-account-dashboard

Multi-account Wagerzon support for the MLB correlated parlay dashboard:
header bar with per-account balance pills, global account selector with
last-used persistence, account-aware /api/place-parlay, placed_parlays.account
column, soft warn on insufficient balance.
EOF
)"
```

- [ ] **Step 2: Verify the merge**

```bash
git log --oneline -5
git diff origin/main..HEAD --stat   # should be empty after the merge
```

- [ ] **Step 3: Remove the worktree and feature branch**

```bash
git worktree remove .worktrees/wagerzon-multi-account-dashboard
git branch -d feature/wagerzon-multi-account-dashboard
git worktree list   # confirm gone
```

- [ ] **Step 4: Restart the dashboard from main and reverify**

```bash
pkill -f mlb_dashboard_server.py || true
pkill -f mlb_dashboard.R || true
"Answer Keys/MLB Dashboard"/venv/bin/python "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" > /tmp/mlb_dashboard_server.log 2>&1 &
sleep 4
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R" > /tmp/mlb_dashboard_R.log 2>&1 &
sleep 6
curl -s http://localhost:8083/api/wagerzon/balances | python3 -m json.tool
```

Expected: same balances JSON as during the smoke test in Task 8.2.

### Task 8.5: Update auto-memory

**Files:**
- Create: `~/.claude-personal/projects/-Users-callancapitolo-NFLWork/memory/wagerzon_multi_account.md`
- Modify: `~/.claude-personal/projects/-Users-callancapitolo-NFLWork/memory/MEMORY.md`

- [ ] **Step 1: Write the reference memory**

Create the file with:

```markdown
---
name: Wagerzon multi-account
description: How the MLB dashboard supports multiple Wagerzon accounts — registry, balance, placement, where state lives
type: reference
---

# Wagerzon multi-account dashboard

## Account discovery
- Pattern: env vars `WAGERZON{SUFFIX}_USERNAME` / `WAGERZON{SUFFIX}_PASSWORD` in `bet_logger/.env` (suffix = "" or single uppercase letter).
- Code: `wagerzon_odds/wagerzon_accounts.py::list_accounts()` returns primary first, then alphabetical.
- Add a new account: two env vars + dashboard restart. No code change.

## Where balance lives
- Fetcher: `wagerzon_odds/wagerzon_balance.py::fetch_available_balance(account)` and `fetch_all(accounts)`.
- Endpoint URL is hardcoded as `_BALANCE_URL` in that module — if WZ moves it, update the constant.
- Server cache: `_BALANCE_CACHE` in `mlb_dashboard_server.py`, in-memory only, populated by `/api/wagerzon/balances` and after each successful placement.

## Where account selection lives
- Persisted: `dashboard_settings(key='wagerzon_last_used', value=<label>)` in `mlb_dashboard.duckdb`.
- Per-placement: `POST /api/place-parlay` requires `{"account": "<label>"}` in the body.
- Recorded: `placed_parlays.account` column. NULL for pre-2026-05 rows (all primary).

## Auth
- `wagerzon_auth.py::get_session(account)` is the single source of truth. Both `parlay_placer` and `parlay_pricer` delegate to it. Per-account session cache keyed by label.
```

- [ ] **Step 2: Add a pointer in MEMORY.md under Reference**

Append to the Reference section:

```markdown
- [Wagerzon multi-account](wagerzon_multi_account.md) — registry pattern, balance fetcher, account selection, persistence. MLB dashboard added 2026-05-01.
```

---

## Self-Review (writer's pass)

**Spec coverage:**
- Q1 (global selector): Task 6.2 ✓
- Q2 (show all balances): Task 6.2 ✓
- Q3 (N accounts via registry): Task 1.2 ✓
- Q4 (refresh on load + after place + manual): Tasks 6.2, 5.4 ✓
- Q5 (soft warn, transient): Task 6.4 ✓
- Q6 (placed_parlays.account): Tasks 5.1, 5.4 ✓
- Q7 (last-used persistence): Tasks 5.1, 5.3 ✓
- Q8 (stale tag, allow placement, red >10min): Tasks 5.2, 6.2 ✓
- Q9 (scope): respected — no `/api/place-bet`, no NFL Draft work
- Available-not-cash: Tasks 3.2, 3.3 ✓
- Auth refactor: Tasks 2.3, 4.2, 4.3 ✓
- bet_logger interaction: noted in Phase 0 prerequisite, no implementation needed ✓
- Documentation: Phase 7 ✓
- Worktree cleanup: Task 8.4 ✓

**Placeholder scan:** none. Initial draft had `REPLACE_ME_FROM_TASK_3_1` placeholders pending manual endpoint discovery; that discovery was completed during execution by extracting the endpoint from existing recon files (`recon_place_parlay.json`), and Tasks 3.1/3.2/3.3 now contain concrete URL, response shape, and parser code.

**Type consistency:** `WagerzonAccount` (dataclass), `BalanceSnapshot` (dataclass), `place_parlays(specs, account)`, `get_session(account)`, `fetch_available_balance(account)`, `fetch_all(accounts)` — used consistently across all tasks.
