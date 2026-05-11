# wagerzon_odds/single_placer.py
"""Wagerzon direct-API placer for SINGLE (straight) bets.

Mirrors parlay_placer.py's session/preflight/drift/error-classification
patterns but for 1-leg payloads.

Public API:
    build_sel_for_single(bet: dict) -> str
    place_single(account: str, bet: dict, session=None) -> dict

Returns a dict with keys:
    status:         'placed' | 'price_moved' | 'rejected' | 'auth_error'
                    | 'network_error' | 'orphaned'
    ticket_number:  WZ ticket string when status == 'placed', else None
    balance_after:  Available balance snapshot from WZ when status == 'placed'
    error_msg:      Human-readable error text (None on success)
    error_msg_key:  Short key for pill rendering (None on success)

Session management:
    In production, pass account=<label string>; the function resolves it to a
    WagerzonAccount via wagerzon_accounts.get_account() and loads a cached
    requests.Session via wagerzon_auth.get_session().
    In tests, pass session=<MagicMock> to bypass the live auth layer entirely.

Endpoints:
    Preflight:   /wager/ConfirmWagerHelper.aspx  (same as parlay placer)
    Submission:  /wager/MakeWagerHelper.aspx      (single-bet endpoint;
                 parlay placer uses PostWagerMultipleHelper.aspx instead)

Tech debt note:
    _is_html_response() and _network_error() duplicate logic from
    parlay_placer.py (_raise_if_html / no equivalent helper). If the
    pattern stabilises, these could be extracted to a shared wagerzon_helpers.py.
"""
from __future__ import annotations
from pathlib import Path
from typing import Optional
import json as _json

import requests

# Load Wagerzon env vars from bet_logger/.env at import time so callers
# importing only single_placer still get .env loaded without depending on
# import order. Falls through silently if dotenv isn't installed.
try:
    from dotenv import load_dotenv
    _ENV_PATH = Path(__file__).resolve().parent.parent / "bet_logger" / ".env"
    if _ENV_PATH.exists():
        load_dotenv(_ENV_PATH)
except ImportError:
    pass

import wagerzon_auth
from wagerzon_accounts import WagerzonAccount, get_account

WAGERZON_BASE_URL = "https://backend.wagerzon.com"
CONFIRM_URL = f"{WAGERZON_BASE_URL}/wager/ConfirmWagerHelper.aspx"
MAKE_URL    = f"{WAGERZON_BASE_URL}/wager/MakeWagerHelper.aspx"

# Tolerance for American-odds drift check. If WZ's preflight returns an Odds
# value that differs from what the dashboard showed the user by more than this,
# we abort and return price_moved.
DRIFT_TOLERANCE_AMERICAN = 1


# ---------------------------------------------------------------------------
# Sel encoding — mirrors parlay_placer.encode_sel for a single leg
# ---------------------------------------------------------------------------

def build_sel_for_single(bet: dict) -> str:
    """Encode a single-leg sel string matching parlay_placer's per-leg format.

    Format: "<play>_<idgm>_<pts>_<odds>"
    Points are printed without trailing .0 for integer values (e.g. -3 not -3.0),
    and with the decimal for non-integers (e.g. -1.5).

    This matches encode_sel() in parlay_placer.py exactly for 1-leg use.
    """
    play = bet["play"]
    idgm = bet["idgm"]
    pts  = bet["line"]
    odds = bet["american_odds"]
    # Drop trailing .0 for integer lines (e.g. -3.0 → -3); keep decimal for fractions
    pts_fmt = int(pts) if pts == int(pts) else pts
    return f"{play}_{idgm}_{pts_fmt}_{odds}"


# ---------------------------------------------------------------------------
# Auth detection — mirrors parlay_placer._raise_if_html logic
# ---------------------------------------------------------------------------

def _is_html_response(resp) -> bool:
    """Return True if the response is HTML (session expired), not JSON."""
    ct = resp.headers.get("content-type", "")
    return "json" not in ct


# ---------------------------------------------------------------------------
# Error helpers
# ---------------------------------------------------------------------------

def _network_error(msg: str) -> dict:
    return {
        "status": "network_error",
        "ticket_number": None,
        "balance_after": None,
        "error_msg": msg,
        "error_msg_key": "network_error",
    }


def _auth_error(msg: str = "session expired") -> dict:
    return {
        "status": "auth_error",
        "ticket_number": None,
        "balance_after": None,
        "error_msg": msg,
        "error_msg_key": "session_expired",
    }


def _rejected(msg: str, key: str = "") -> dict:
    return {
        "status": "rejected",
        "ticket_number": None,
        "balance_after": None,
        "error_msg": msg,
        "error_msg_key": key,
    }


def _price_moved(expected: int, actual: int) -> dict:
    return {
        "status": "price_moved",
        "ticket_number": None,
        "balance_after": None,
        "error_msg": f"Odds drifted from {expected} to {actual}",
        "error_msg_key": "drift",
    }


# ---------------------------------------------------------------------------
# Preflight payload builder — mirrors parlay_placer.encode_detail_data
# ---------------------------------------------------------------------------

def _build_confirm_payload(bet: dict) -> dict:
    """Build the ConfirmWagerHelper POST parameters for a single leg."""
    amount = bet["actual_size"]
    amount_str = str(int(amount)) if amount == int(amount) else str(amount)
    detail_data = [
        {
            "Amount": amount_str,
            "RiskWin": 0,
            "TeaserPointsPurchased": 0,
            "IdGame": bet["idgm"],
            "Play": bet["play"],
            "Pitcher": bet.get("pitcher", 0),
            "Points": {
                "BuyPoints": 0,
                "BuyPointsDesc": "",
                "LineDesc": "",
                "selected": True,
            },
        }
    ]
    return {
        "IDWT": "0",
        "WT": "0",          # WT=0 = straight/single (WT=1 = parlay)
        "amountType": "0",
        "open": "0",
        "sameAmount": "false",
        "sameAmountNumber": amount_str,
        "useFreePlayAmount": "false",
        "sel": build_sel_for_single(bet),
        "detailData": _json.dumps(detail_data),
    }


# ---------------------------------------------------------------------------
# Top-level place_single
# ---------------------------------------------------------------------------

def place_single(account: str, bet: dict, session=None) -> dict:
    """Submit a single straight bet at Wagerzon via REST.

    Args:
        account:  Wagerzon account label (e.g. "primary"). Ignored if
                  `session` is provided (test injection path).
        bet:      Dict with keys: idgm, play, line, american_odds,
                  actual_size, wz_odds_at_place, and optionally pitcher.
        session:  For testing only — pass a MagicMock requests.Session.
                  Production code should leave this as None.

    Returns:
        dict with status / ticket_number / balance_after / error_msg /
        error_msg_key as documented in the module docstring.
    """
    # Session resolution: injected (tests) or from live auth cache (production)
    if session is not None:
        sess = session
    else:
        wz_account: WagerzonAccount = get_account(account)
        sess = wagerzon_auth.get_session(wz_account)

    confirm_payload = _build_confirm_payload(bet)

    # -----------------------------------------------------------------------
    # Step 1: ConfirmWagerHelper (preflight + drift check)
    # -----------------------------------------------------------------------
    try:
        confirm_resp = sess.post(
            CONFIRM_URL,
            data=confirm_payload,
            timeout=15,
            headers={"Accept": "application/json"},
        )
    except requests.RequestException as e:
        return _network_error(f"{type(e).__name__}: {e}")

    if _is_html_response(confirm_resp):
        return _auth_error("Wagerzon returned HTML at preflight (session expired)")

    confirm_json = confirm_resp.json()
    result_block = confirm_json.get("result") or {}

    # Check for an error code in the preflight response itself
    if "ErrorCode" in result_block:
        key = result_block["ErrorCode"]
        msg = result_block.get("ErrorMessage") or key
        return _rejected(msg, key)

    details = result_block.get("details") or []
    if not details:
        # Line pulled, bet too large, or other server-side refusal
        return _rejected(
            msg="Wagerzon returned empty details (line pulled?)",
            key="empty_details",
        )

    first_detail = details[0]
    wz_odds_now = first_detail.get("Odds")
    if wz_odds_now is None:
        return _rejected(
            msg="Wagerzon preflight details missing Odds field",
            key="missing_odds",
        )

    expected_odds = int(bet["wz_odds_at_place"])
    if abs(int(wz_odds_now) - expected_odds) > DRIFT_TOLERANCE_AMERICAN:
        return _price_moved(expected_odds, int(wz_odds_now))

    # -----------------------------------------------------------------------
    # Step 2: MakeWagerHelper (real submission)
    # Network failure here is an orphan candidate — we log forensics and
    # return network_error so the dashboard can surface it.
    # -----------------------------------------------------------------------
    make_payload = dict(confirm_payload)
    try:
        make_resp = sess.post(
            MAKE_URL,
            data=make_payload,
            timeout=20,
            headers={"Accept": "application/json"},
        )
    except requests.RequestException as e:
        return _network_error(f"{type(e).__name__}: {e}")

    if _is_html_response(make_resp):
        return _auth_error("Wagerzon returned HTML at submission (session expired)")

    make_json = make_resp.json()
    make_result = make_json.get("result") or {}

    # Post-submission rejection (e.g. INSUFFICIENT_BALANCE at placement time)
    if "ErrorCode" in make_result:
        key = make_result["ErrorCode"]
        msg = make_result.get("ErrorMessage") or key
        return _rejected(msg, key)

    ticket = make_result.get("WagerNumber")
    if not ticket:
        # WZ confirmed but returned no ticket; treat as orphan candidate.
        return {
            "status": "orphaned",
            "ticket_number": None,
            "balance_after": None,
            "error_msg": "Submitted to WZ but no WagerNumber returned",
            "error_msg_key": "no_ticket",
        }

    return {
        "status": "placed",
        "ticket_number": str(ticket),
        "balance_after": make_result.get("AvailBalance"),
        "error_msg": None,
        "error_msg_key": None,
    }
