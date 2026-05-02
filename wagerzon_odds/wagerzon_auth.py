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
# Guards mutations of _SESSION_CACHE itself and of _LABEL_LOCKS.
_CACHE_LOCK = threading.Lock()
# Per-label lock so concurrent logins for *different* accounts can run in
# parallel. A cache miss for label X holds only the lock for X — other
# labels are unaffected.
_LABEL_LOCKS: dict[str, threading.Lock] = {}


def _lock_for_label(label: str) -> threading.Lock:
    with _CACHE_LOCK:
        return _LABEL_LOCKS.setdefault(label, threading.Lock())


def get_session(account: WagerzonAccount) -> requests.Session:
    """Return a logged-in Session for `account`, logging in on first call.

    Concurrent calls for *different* accounts run in parallel (each label
    has its own lock). Concurrent calls for the *same* account serialize
    on that label's lock; the loser sees the winner's cached session.
    """
    # Fast path: cache hit, no lock acquisition.
    cached = _SESSION_CACHE.get(account.label)
    if cached is not None:
        return cached
    label_lock = _lock_for_label(account.label)
    with label_lock:
        # Re-check under the per-label lock — another thread may have
        # logged in for this same label while we were waiting.
        cached = _SESSION_CACHE.get(account.label)
        if cached is not None:
            return cached
        session = _build_session()
        _login(session, account)  # network I/O — held only by THIS label's lock
        with _CACHE_LOCK:
            _SESSION_CACHE[account.label] = session
        return session


def clear_session_cache(label: Optional[str] = None) -> None:
    """Drop cached sessions. Pass a label to drop one; omit to drop all.

    Closes the underlying connection pool of any session dropped, so
    Phase 3's "401 -> clear -> retry" path doesn't leak idle connections.
    """
    with _CACHE_LOCK:
        if label is None:
            dropped = list(_SESSION_CACHE.values())
            _SESSION_CACHE.clear()
        else:
            dropped_one = _SESSION_CACHE.pop(label, None)
            dropped = [dropped_one] if dropped_one is not None else []
    for sess in dropped:
        try:
            sess.close()
        except Exception:
            pass


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
    if "__VIEWSTATE" not in fields:
        raise RuntimeError(
            "Could not find __VIEWSTATE on Wagerzon login page — "
            "page structure may have changed (account=%s)" % account.label
        )
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
