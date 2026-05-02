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
