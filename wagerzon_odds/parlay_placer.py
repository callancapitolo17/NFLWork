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
from pathlib import Path
from typing import Optional
import os
import re
import requests
import json as _json

# Load WAGERZON_USERNAME / WAGERZON_PASSWORD from bet_logger/.env at import
# time, matching the pattern in wagerzon_odds/scraper_v2.py. This way the
# placer works whether it's imported by the dashboard server, run as a
# script, or called from a Python REPL — no shell-export dance required.
# Falls through silently if dotenv isn't installed; in that case the caller
# is expected to have set the env vars some other way.
try:
    from dotenv import load_dotenv
    _ENV_PATH = Path(__file__).resolve().parent.parent / "bet_logger" / ".env"
    if _ENV_PATH.exists():
        load_dotenv(_ENV_PATH)
except ImportError:
    pass

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
    """Build the `detailData` JSON array for ConfirmWagerHelper / PostWager.

    On RiskWin:
        We send `RiskWin: 0` for both ConfirmWagerHelper (preflight) and
        PostWagerMultipleHelper (placement). This matches the Wagerzon UI's
        captured behavior — the browser sends RiskWin=0 to both endpoints
        when a real user clicks Place.

        wagerzon_odds/parlay_pricer.py uses RiskWin="2" instead, but for a
        different use case: it queries the parlay-payout curve at many
        stakes including stakes the user couldn't fund, and "2" suppresses
        the balance check so it gets a price at any amount. The placer
        only ever queries at the user's actual chosen wager, so the
        balance check is appropriate.

        Edge case (TODO): if account balance < wager at preflight time,
        Wagerzon returns an error response without a `details` array and
        _confirm_preflight raises an unhandled IndexError/KeyError. This
        is rare in practice (user has $3k+ balance, bets $15-$100) but
        worth handling cleanly in a future hardening pass — surface as
        a `rejected: insufficient_balance` status instead of a 500.
    """
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
            # Defensive parse: IDWT can come back as int, str, or None depending
            # on Wagerzon-side serialization quirks. Use the truthy-check
            # pattern but log when something unexpected happens so we can
            # diagnose post-hoc orphans without losing the audit trail.
            raw_idwt = wpr.get("IDWT")
            raw_ticket = wpr.get("TicketNumber")
            try:
                idwt_val = int(raw_idwt) if raw_idwt not in (None, "", 0) else None
            except (TypeError, ValueError):
                idwt_val = None
            ticket_val = str(raw_ticket) if raw_ticket not in (None, "") else None
            if idwt_val is None or ticket_val is None:
                # Not fatal — the bet is placed, but the local audit will
                # be incomplete. Caller (api_place_parlay) writes an orphan
                # row with raw_response so we can reconstruct after the fact.
                import sys as _sys
                print(
                    f"!! parlay_placer: confirmed placement missing "
                    f"identifiers — parlay_hash={spec.parlay_hash} "
                    f"raw_idwt={raw_idwt!r} raw_ticket={raw_ticket!r} "
                    f"wpr_keys={list(wpr.keys())}",
                    file=_sys.stderr, flush=True,
                )
            results.append(PlacementResult(
                parlay_hash=spec.parlay_hash,
                status="placed",
                ticket_number=ticket_val,
                idwt=idwt_val,
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


# ---------------------------------------------------------------------------
# Top-level place_parlays with retry, drift, dry-run, network handling
# ---------------------------------------------------------------------------

def _drift_error_msg(expected: float, actual: float) -> str:
    return f"Expected ${expected:.2f}, Wagerzon offered ${actual:.2f}"


def _preflight_with_retry(spec: ParlaySpec) -> tuple[float, float]:
    """ConfirmWagerHelper with one re-login retry on AuthExpired."""
    try:
        return _confirm_preflight(spec.legs, spec.amount)
    except AuthExpired:
        _clear_session()
        return _confirm_preflight(spec.legs, spec.amount)


def _post_with_retry(specs: list[ParlaySpec]) -> list[PlacementResult]:
    """PostWagerMultipleHelper with one re-login + re-preflight retry.

    On AuthExpired: re-login, re-preflight EACH spec independently, mark
    drifted specs as price_moved (per-spec, not batch-wide), then place
    only the survivors. Bare _confirm_preflight here is intentional — we
    have already consumed the one allowed auth retry.
    """
    try:
        return _post_wagers(specs)
    except AuthExpired:
        _clear_session()
        retry_results: dict[str, PlacementResult] = {}
        retry_survivors: list[ParlaySpec] = []
        for spec in specs:
            win, risk = _confirm_preflight(spec.legs, spec.amount)
            if not _drift_ok(spec.expected_win, win):
                retry_results[spec.parlay_hash] = PlacementResult(
                    parlay_hash=spec.parlay_hash, status="price_moved",
                    error_msg=_drift_error_msg(spec.expected_win, win),
                    actual_win=win, actual_risk=risk,
                )
            else:
                retry_survivors.append(spec)
        if retry_survivors:
            for r in _post_wagers(retry_survivors):
                retry_results[r.parlay_hash] = r
        return [retry_results[spec.parlay_hash] for spec in specs]


def place_parlays(specs: list[ParlaySpec], dry_run: bool = False) -> list[PlacementResult]:
    """Place a batch of parlays at Wagerzon.

    For each spec:
      1) ConfirmWagerHelper preflight (with one auth-retry)
      2) Drift check vs expected_win — abort that spec on drift > $0.01
    Then for surviving specs:
      3) If dry_run: return would_place results
      4) Else: PostWagerMultipleHelper (with one auth-retry that re-runs
         preflight defensively before placing).
    """
    seen = set()
    for spec in specs:
        if spec.parlay_hash in seen:
            raise ValueError(f"duplicate parlay_hash in batch: {spec.parlay_hash!r}")
        seen.add(spec.parlay_hash)

    results_by_hash: dict[str, PlacementResult] = {}
    survivors: list[ParlaySpec] = []

    # Preflight + drift check per spec
    for spec in specs:
        try:
            win, risk = _preflight_with_retry(spec)
        except AuthExpired:
            results_by_hash[spec.parlay_hash] = PlacementResult(
                parlay_hash=spec.parlay_hash, status="auth_error",
                error_msg="auth_error: session expired",
            )
            continue
        except requests.exceptions.RequestException as e:
            results_by_hash[spec.parlay_hash] = PlacementResult(
                parlay_hash=spec.parlay_hash, status="network_error",
                error_msg=f"network_error: {type(e).__name__}",
            )
            continue
        if not _drift_ok(spec.expected_win, win):
            results_by_hash[spec.parlay_hash] = PlacementResult(
                parlay_hash=spec.parlay_hash, status="price_moved",
                error_msg=_drift_error_msg(spec.expected_win, win),
                actual_win=win, actual_risk=risk,
            )
            continue
        if dry_run:
            results_by_hash[spec.parlay_hash] = PlacementResult(
                parlay_hash=spec.parlay_hash, status="would_place",
                actual_win=win, actual_risk=risk,
            )
            continue
        survivors.append(spec)

    # Live placement for specs that passed preflight
    if survivors:
        try:
            live_results = _post_with_retry(survivors)
        except AuthExpired:
            live_results = [
                PlacementResult(parlay_hash=s.parlay_hash, status="auth_error",
                                error_msg="auth_error: re-login failed")
                for s in survivors
            ]
        except requests.exceptions.RequestException as e:
            live_results = [
                PlacementResult(parlay_hash=s.parlay_hash, status="network_error",
                                error_msg=f"network_error: {type(e).__name__}")
                for s in survivors
            ]
        for r in live_results:
            results_by_hash[r.parlay_hash] = r

    # Preserve input order
    return [results_by_hash[s.parlay_hash] for s in specs]
