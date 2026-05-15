"""Wagerzon single-bet preview / pricer.

Mirror of `wagerzon_odds.parlay_pricer.get_parlay_price` for single
straight bets. Calls ConfirmWagerHelper with `RiskWin="2"` so WZ
returns a price quote without balance validation. Used by the MLB
Dashboard's editable-Risk feature to:

  1. Verify the actual `Win` (to-win) WZ would credit at any user-typed
     amount, instead of relying on local American-odds math.
  2. Detect line drift before the user clicks Place — the response
     includes WZ's *current* odds for the leg.
  3. Surface other WZ rejections (integer-only amounts, MAXRISK, line
     pulled, etc.) directly via `ErrorMsg` / `ErrorMsgKey`.

Does NOT place the bet. For actual placement, see `single_placer.py`.
"""

import json as _json

import requests

from config import WAGERZON_BASE_URL
from single_placer import build_sel_for_single

CONFIRM_URL = f"{WAGERZON_BASE_URL}/wager/ConfirmWagerHelper.aspx"


def _network_error(msg: str) -> dict:
    return {
        "win": None,
        "current_wz_odds": None,
        "error_msg": msg,
        "error_msg_key": "network_error",
    }


def _auth_error(msg: str = "Wagerzon returned HTML at preflight (session expired)") -> dict:
    return {
        "win": None,
        "current_wz_odds": None,
        "error_msg": msg,
        "error_msg_key": "session_expired",
    }


def get_single_price(session: requests.Session, bet: dict, amount: float) -> dict:
    """Get WZ's current price for a single bet at the given amount.

    Args:
        session: An authenticated requests.Session for WZ. In tests, pass a
            MagicMock — the function only exercises `.post`, `.json()`,
            `.headers`, `.status_code`.
        bet: Dict with keys `idgm, play, line, american_odds, pitcher`
            (pitcher optional, defaults to 0).
        amount: Risk amount to query at. Passed verbatim to WZ — if WZ
            rejects (e.g. integer-only rule), the rejection surfaces in
            the returned `error_msg_key`.

    Returns:
        Dict with `{win, current_wz_odds, error_msg, error_msg_key}`.
        `win` is the to-win amount WZ would credit; `current_wz_odds` is
        the leg's current American odds at WZ. Both None on failure.
    """
    # ConfirmWagerHelper payload mirrors single_placer._build_confirm_payload
    # but with RiskWin="2" (preview, no balance check) instead of "0".
    amount_str = str(int(amount)) if amount == int(amount) else str(amount)
    detail_data = [
        {
            "Amount": amount_str,
            "RiskWin": "2",
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
    payload = {
        "IDWT": "0",
        "WT": "0",
        "amountType": "0",
        "open": "0",
        "sameAmount": "false",
        "sameAmountNumber": amount_str,
        "useFreePlayAmount": "false",
        "sel": build_sel_for_single(bet),
        "detailData": _json.dumps(detail_data),
    }

    try:
        resp = session.post(
            CONFIRM_URL,
            data=payload,
            timeout=15,
            headers={"Accept": "application/json"},
        )
    except requests.RequestException as e:
        return _network_error(f"{type(e).__name__}: {e}")

    if "json" not in resp.headers.get("content-type", ""):
        return _auth_error()

    try:
        body = resp.json()
    except ValueError as e:
        return _network_error(f"json decode failed: {e}")

    result = body.get("result") or {}
    err_key = result.get("ErrorMsgKey") or result.get("ErrorCode") or ""
    err_msg = result.get("ErrorMsg") or result.get("ErrorMessage") or ""
    # Unlike parlay_pricer (which still surfaces Win/Risk when err_key is
    # something like MINWAGERONLINE), single_pricer treats any error key as
    # fatal — the caller's UI shows the key directly, so we don't need the
    # valid-but-flagged data.
    if err_key:
        return {
            "win": None,
            "current_wz_odds": None,
            "error_msg": err_msg or err_key,
            "error_msg_key": err_key,
        }

    details = result.get("details") or []
    if not details:
        return {
            "win": None,
            "current_wz_odds": None,
            "error_msg": "Wagerzon returned empty details (line pulled?)",
            "error_msg_key": "empty_details",
        }
    outer = details[0]
    win = outer.get("Win")
    odds_now = outer.get("Odds")
    return {
        "win": win,
        "current_wz_odds": int(odds_now) if odds_now is not None else None,
        "error_msg": "",
        "error_msg_key": "",
    }
