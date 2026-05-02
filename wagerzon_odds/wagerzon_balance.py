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
        {"result": {"RealAvailBalance": "256 ", "AvailBalance": "-1,744 ",
                    "CurrentBalance": "-959 ", "CreditLimit": "2,000 ",
                    "AmountAtRisk": "785 ", "Player": "ACCT",
                    "Password": "<plaintext>", ...}}

    Returns:
        (available, cash)
        - available: RealAvailBalance — the actual wagerable amount
          (= AvailBalance + CreditLimit, precomputed by WZ). This is
          what determines whether WZ will accept a new bet, so it's
          what the dashboard's insufficient-balance warning gates on.
          AvailBalance alone (cash minus open exposure) goes negative
          whenever the user is using their credit line, which would
          cause the warning to fire on almost every bet.
        - cash: CurrentBalance — raw cash on deposit, exposed for
          tooltip / debugging. Optional.

    Numbers come as strings with comma thousand-separators and a
    trailing space. We extract ONLY the two numeric fields above.
    The response body also contains the user's password in plaintext
    (result.Password); never log the raw body.
    """
    data = resp.json()
    result = data.get("result") or {}

    raw_avail = result.get("RealAvailBalance")
    if raw_avail is None:
        raise ValueError("balance response missing 'RealAvailBalance'")
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
