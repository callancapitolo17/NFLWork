# wagerzon_odds/parlay_placer.py
"""Wagerzon parlay placer — pure REST, no browser.

See docs/superpowers/specs/2026-04-26-mlb-parlay-auto-placement-design.md
for the full design. This module provides:
    - ParlaySpec / Leg / PlacementResult dataclasses
    - encode_sel / encode_detail_data leg-encoding helpers
    - _get_session / _clear_session session management
    - _confirm_preflight / _drift_ok price safety checks
    - place_parlays(specs, dry_run=False) — top-level entry point
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
import os
import re
import requests
import json as _json

WAGERZON_BASE_URL = "https://backend.wagerzon.com"

_CACHED_SESSION: Optional[requests.Session] = None


@dataclass
class Leg:
    """One leg of a parlay.

    play codes: 0=away spread, 1=home spread, 2=over, 3=under,
                4=away ML, 5=home ML
    points: signed line value (negative for fav spread / over)
    odds: integer American odds (e.g. 117 = +117, -130 = -130)
    pitcher: 0 = action / 3 = listed (mirrors what ConfirmWagerHelper expects)
    """
    idgm: int
    play: int
    points: float
    odds: int
    pitcher: int = 0


@dataclass
class ParlaySpec:
    parlay_hash: str
    legs: list[Leg]
    amount: float
    expected_win: float
    expected_risk: float


@dataclass
class PlacementResult:
    parlay_hash: str
    status: str
    ticket_number: Optional[str] = None
    idwt: Optional[int] = None
    actual_win: Optional[float] = None
    actual_risk: Optional[float] = None
    error_msg: str = ""
    error_msg_key: str = ""
    raw_response: Optional[str] = None


def encode_sel(legs: list[Leg]) -> str:
    """Encode legs as Wagerzon's `sel` parameter: play_idgm_points_odds, ..."""
    parts = []
    for leg in legs:
        # points printed without trailing .0 if integer
        pts = int(leg.points) if leg.points == int(leg.points) else leg.points
        parts.append(f"{leg.play}_{leg.idgm}_{pts}_{leg.odds}")
    return ",".join(parts)


def encode_detail_data(legs: list[Leg], amount: float) -> list[dict]:
    """Build the `detailData` JSON array for ConfirmWagerHelper / PostWager."""
    return [
        {
            "Amount": str(int(amount)) if amount == int(amount) else str(amount),
            "RiskWin": 0,
            "TeaserPointsPurchased": 0,
            "IdGame": leg.idgm,
            "Play": leg.play,
            "Pitcher": leg.pitcher,
            "Points": {
                "BuyPoints": 0,
                "BuyPointsDesc": "",
                "LineDesc": "",
                "selected": True,
            },
        }
        for leg in legs
    ]


def _get_session() -> requests.Session:
    """Return a logged-in Wagerzon session, reusing a cached one if present.

    Mirrors the login pattern in wagerzon_odds/parlay_pricer.py: GET base URL,
    if not already authenticated, parse __VIEWSTATE etc. and form-POST
    Account/Password.

    Resets via _clear_session() (called from auth_error retry path).
    """
    global _CACHED_SESSION
    if _CACHED_SESSION is not None:
        return _CACHED_SESSION

    session = requests.Session()
    resp = session.get(WAGERZON_BASE_URL, timeout=15)
    resp.raise_for_status()
    if "NewSchedule" in resp.url or "Welcome" in resp.url:
        _CACHED_SESSION = session
        return session

    html = resp.text
    fields = {}
    for name in ("__VIEWSTATE", "__VIEWSTATEGENERATOR", "__EVENTVALIDATION",
                 "__EVENTTARGET", "__EVENTARGUMENT"):
        m = re.search(rf'(?:name|id)="{name}"[^>]*value="([^"]*)"', html)
        if m:
            fields[name] = m.group(1)
    fields["Account"] = os.environ["WAGERZON_USERNAME"]
    fields["Password"] = os.environ["WAGERZON_PASSWORD"]
    fields["BtnSubmit"] = ""

    resp = session.post(WAGERZON_BASE_URL, data=fields, timeout=15)
    resp.raise_for_status()
    _CACHED_SESSION = session
    return session


def _clear_session() -> None:
    """Force the next _get_session() call to re-login."""
    global _CACHED_SESSION
    _CACHED_SESSION = None


CONFIRM_URL = f"{WAGERZON_BASE_URL}/wager/ConfirmWagerHelper.aspx"
DRIFT_TOLERANCE_USD = 0.01


def _drift_ok(expected: float, actual: float) -> bool:
    """True if returned Win matches expected within $0.01."""
    return abs(actual - expected) <= DRIFT_TOLERANCE_USD


class AuthExpired(Exception):
    """Wagerzon returned an HTML login page instead of JSON."""


def _raise_if_html(resp) -> None:
    ct = resp.headers.get("content-type", "")
    if "json" not in ct:
        raise AuthExpired(f"Non-JSON response: content-type={ct!r}")


def _confirm_preflight(legs: list[Leg], amount: float) -> tuple[float, float]:
    """Call ConfirmWagerHelper, return (Win, Risk).

    Raises ValueError on malformed response. Auth/HTML responses raise
    AuthExpired (handled by caller).
    """
    session = _get_session()
    detail_data = encode_detail_data(legs, amount)
    params = {
        "IDWT": "0",
        "WT": "1",
        "amountType": "0",
        "open": "0",
        "sameAmount": "false",
        "sameAmountNumber": str(int(amount)),
        "useFreePlayAmount": "false",
        "sel": encode_sel(legs),
        "detailData": _json.dumps(detail_data),
    }
    resp = session.post(CONFIRM_URL, data=params, timeout=15,
                        headers={"Accept": "application/json"})
    _raise_if_html(resp)
    data = resp.json()
    details = data["result"]["details"][0]
    return float(details["Win"]), float(details["Risk"])


POST_URL = f"{WAGERZON_BASE_URL}/wager/PostWagerMultipleHelper.aspx"

ERROR_KEY_MAP = {
    "insufficient_funds": "insufficient balance",
    "bet_too_large":      "exceeds limit",
    "line_unavailable":   "line pulled",
}


def _build_post_request(spec: ParlaySpec, password: str) -> dict:
    """Build a single parlay POST request for PostWagerMultipleHelper."""
    detail_data = encode_detail_data(spec.legs, spec.amount)
    return {
        "WT": "1",
        "open": 0,
        "IDWT": "0",
        "sel": encode_sel(spec.legs),
        "sameAmount": False,
        "amountType": "0",
        "detailData": _json.dumps(detail_data),
        "confirmPassword": password,
        "sameAmountNumber": str(int(spec.amount)) if spec.amount == int(spec.amount) else str(spec.amount),
        "useFreePlayAmount": False,
        "roundRobinCombinations": "",
    }


def _post_wagers(specs: list[ParlaySpec]) -> list[PlacementResult]:
    """POST one bulk request to PostWagerMultipleHelper. Returns one
    PlacementResult per input spec, in order.

    Raises AuthExpired on HTML response.
    """
    session = _get_session()
    password = os.environ["WAGERZON_PASSWORD"]
    payload = [_build_post_request(s, password) for s in specs]
    body = {"postWagerRequests": _json.dumps(payload)}

    resp = session.post(POST_URL, data=body, timeout=30,
                        headers={"Accept": "application/json"})
    _raise_if_html(resp)
    data = resp.json()
    raw = resp.text

    results = []
    for spec, item in zip(specs, data["result"]):
        wpr = item["WagerPostResult"]
        if wpr.get("Confirm") and not wpr.get("ErrorMsgKey"):
            results.append(PlacementResult(
                parlay_hash=spec.parlay_hash,
                status="placed",
                ticket_number=str(wpr.get("TicketNumber")) if wpr.get("TicketNumber") else None,
                idwt=int(wpr.get("IDWT")) if wpr.get("IDWT") else None,
                actual_win=float(wpr.get("Win", 0)),
                actual_risk=float(wpr.get("Risk", 0)),
                raw_response=raw,
            ))
        else:
            key = wpr.get("ErrorMsgKey", "") or ""
            friendly = ERROR_KEY_MAP.get(key, key) or "unknown error"
            results.append(PlacementResult(
                parlay_hash=spec.parlay_hash,
                status="rejected",
                error_msg=f"rejected: {friendly}" if key else (wpr.get("ErrorMsg") or "rejected"),
                error_msg_key=key,
                raw_response=raw,
            ))
    return results
