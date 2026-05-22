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
    Submission:  /wager/PostWagerMultipleHelper.aspx  (same endpoint parlays
                 use; WT=0 distinguishes singles from parlays which use WT=1)

Tech debt note:
    _is_html_response() and _network_error() duplicate logic from
    parlay_placer.py (_raise_if_html / no equivalent helper). If the
    pattern stabilises, these could be extracted to a shared wagerzon_helpers.py.

DEPLOYMENT NOTE (2026-05-11): The submission flow uses PostWagerMultipleHelper.aspx
(same endpoint parlays use) with WT=0 to distinguish singles. This payload shape
is inferred from recon — it has NOT been verified end-to-end against a live
Wagerzon account for singles. Before enabling this module in /api/place-bet
for production traffic, run wagerzon_odds/recon_place_parlay.py (or analogous
recon for singles) and confirm:
  1. PostWagerMultipleHelper accepts a single-leg payload with WT=0
  2. CreateWagerHelper.aspx is NOT called in the parlay placement flow
     (grep confirms it is absent from parlay_placer.py) — so it is likely
     a Wagerzon UI-only step and is not required in our API path
  3. The actual ticket-extraction key (WagerNumber vs different name)

Password / confirmPassword:
  PostWagerMultipleHelper requires `confirmPassword` in the submission payload
  to authorize the wager (mirrors parlay_placer._build_post_request). In the
  production path the field is populated from wz_account.password resolved via
  get_account(). In the test path (session= injected), the field is omitted
  because mock sessions do not validate credentials.
"""
from __future__ import annotations
from datetime import datetime, timezone
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
# Single source of truth for the drift tolerance — same value parlay placement
# uses ($0.01). Importing instead of re-declaring keeps the two placers in
# lockstep if we ever revisit it.
from parlay_placer import DRIFT_TOLERANCE_USD

WAGERZON_BASE_URL = "https://backend.wagerzon.com"
CONFIRM_URL = f"{WAGERZON_BASE_URL}/wager/ConfirmWagerHelper.aspx"
MAKE_URL    = f"{WAGERZON_BASE_URL}/wager/PostWagerMultipleHelper.aspx"


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


def _price_moved_by_win(expected: float, actual: float) -> dict:
    """Drift detected via Win-on-Win comparison. Message names both dollar
    values so a real drift is debuggable from the error toast alone."""
    return {
        "status": "price_moved",
        "ticket_number": None,
        "balance_after": None,
        "error_msg": f"Expected ${expected:.2f}, Wagerzon offered ${actual:.2f}",
        "error_msg_key": "drift",
    }


def _compute_expected_win(odds: int, risk: float) -> float:
    """American-odds → expected Win at the given stake. Used as fallback
    when the dashboard didn't supply a verified-from-WZ `expected_win`
    (e.g. user clicked Place without ever editing the Risk field).

    Mirrors the math in mlb_dashboard_server.py::_resolve_amount_and_win's
    fallback branch so the placer and dashboard agree by construction."""
    decimal = (odds / 100 + 1) if odds > 0 else (100 / -odds + 1)
    return round(risk * (decimal - 1), 2)


def _first_present(*candidates):
    """Return the first candidate that isn't None/""/0/0.0. Mirrors the
    helper in parlay_placer._post_wagers so single + parlay handle WZ's
    variable response shapes identically."""
    for c in candidates:
        if c not in (None, "", 0, 0.0):
            return c
    return None


# Where placement-response forensics get written. Gitignored so we can
# capture full WZ responses for orphan / rejection paths without leaking
# them to source control. .gitignore pattern: .placement_debug_*.json
_DEBUG_DIR = Path(__file__).resolve().parent
_DEBUG_SENSITIVE_KEYS = {"confirmPassword", "password", "Password"}


def _sanitize_for_log(obj):
    """Recursively redact known-sensitive keys before writing to disk.
    Defensive — WZ's balance endpoint is known to echo plaintext passwords,
    so we never trust an arbitrary WZ response with raw on-disk persistence."""
    if isinstance(obj, dict):
        return {
            k: ("<redacted>" if k in _DEBUG_SENSITIVE_KEYS else _sanitize_for_log(v))
            for k, v in obj.items()
        }
    if isinstance(obj, list):
        return [_sanitize_for_log(x) for x in obj]
    return obj


def _log_placement_response(
    bet: dict,
    request_payload: dict,
    response_body: dict,
    classification: str,
    reason: str,
    endpoint: str = "PostWagerMultipleHelper",
) -> None:
    """Dump WZ's request/response pair to a gitignored file. Best-effort:
    never raises, never blocks the placement flow.

    Called whenever placement did NOT end in 'placed' — at EITHER the
    preflight step (ConfirmWagerHelper: price_moved, empty_details,
    missing_win, rejected) OR the submission step (PostWagerMultipleHelper:
    rejected, orphaned). The captured file is the forensic trail for any
    surprise; `endpoint` and `classification` together identify which
    branch fired. Filename includes classification so multiple captures
    for the same bet are triageable at a glance.
    """
    try:
        ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        bet_hash = (bet.get("bet_hash") or "no-hash")
        # Strip path separators defensively; cap length to keep filenames sane.
        safe_hash = "".join(c for c in bet_hash if c.isalnum() or c in "-_")[:40]
        # Include classification in the filename so multiple captures for the
        # same bet (e.g. retry attempts that each hit price_moved) are easy to
        # triage at a glance without opening each file.
        safe_class = "".join(c for c in (classification or "unknown") if c.isalnum() or c in "-_")[:24]
        path = _DEBUG_DIR / f".placement_debug_{ts}_{safe_class}_{safe_hash}.json"
        payload = {
            "classification": classification,
            "reason":         reason,
            "endpoint":       endpoint,
            "bet_summary": {
                "bet_hash":      bet.get("bet_hash"),
                "idgm":          bet.get("idgm"),
                "play":          bet.get("play"),
                "line":          bet.get("line"),
                "american_odds": bet.get("american_odds"),
                "actual_size":   bet.get("actual_size"),
                "expected_win":  bet.get("expected_win"),
                "market":        bet.get("market"),
                "bet_on":        bet.get("bet_on"),
            },
            "request_payload":  _sanitize_for_log(request_payload),
            "response_body":    _sanitize_for_log(response_body),
        }
        path.write_text(_json.dumps(payload, indent=2, default=str))
    except Exception:
        # NEVER raise from logging — placement flow must not break on disk I/O.
        pass


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
        wz_account = None  # test path — no live account object
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
        # No auth_error retry: single placer is stateless; caller is expected to
        # retry on auth_error. The parlay placer retries internally because parlay
        # state is more complex to reconstruct. Singles are cheap to retry.
        return _auth_error("Wagerzon returned HTML at preflight (session expired)")

    confirm_json = confirm_resp.json()
    result_block = confirm_json.get("result") or {}

    # Check for an error code in the preflight response itself.
    # WZ endpoints are inconsistent: ConfirmWagerHelper uses ErrorMsgKey/ErrorMsg
    # (same convention as parlay_placer) while other endpoints use ErrorCode/ErrorMessage.
    # Check both so we handle either shape defensively.
    err_key = result_block.get("ErrorMsgKey") or result_block.get("ErrorCode")
    if err_key:
        err_msg = result_block.get("ErrorMsg") or result_block.get("ErrorMessage") or err_key
        _log_placement_response(bet, confirm_payload, confirm_json,
                                classification="rejected",
                                reason=str(err_key),
                                endpoint="ConfirmWagerHelper")
        return _rejected(err_msg, err_key)

    details = result_block.get("details") or []
    if not details:
        # Line pulled, bet too large, or other server-side refusal
        _log_placement_response(bet, confirm_payload, confirm_json,
                                classification="empty_details",
                                reason="line_pulled_or_refused",
                                endpoint="ConfirmWagerHelper")
        return _rejected(
            msg="Wagerzon returned empty details (line pulled?)",
            key="empty_details",
        )

    # Price-match check: Win-on-Win, identical protocol to parlay_placer.
    # WZ reliably returns `Win` in details[0]; `Odds` is per-market-conditional
    # (omitted for some single shapes like alt-spreads and F-period totals).
    # When the dashboard verified the bet via /api/wz-quote-single, it round-
    # trips the verified Win as `expected_win` so the comparison is WZ-said-X
    # against WZ-says-X. If absent, we fall back to math from the user's
    # stake + dashboard-quoted American odds.
    first_detail = details[0]
    wz_win = first_detail.get("Win")
    if wz_win is None:
        _log_placement_response(bet, confirm_payload, confirm_json,
                                classification="missing_win",
                                reason="no_win_in_details",
                                endpoint="ConfirmWagerHelper")
        return _rejected(
            msg="Wagerzon preflight details missing Win field",
            key="missing_win",
        )

    expected_win = bet.get("expected_win")
    if expected_win is None:
        expected_win = _compute_expected_win(
            odds=int(bet["wz_odds_at_place"]),
            risk=float(bet["actual_size"]),
        )

    if abs(float(wz_win) - float(expected_win)) > DRIFT_TOLERANCE_USD:
        # Log the preflight body so we can see WZ's exact Win + Odds + line
        # info. This is the diagnostic gap we're closing — without this, a
        # price_moved rejection gives us no visibility into what WZ actually
        # returned, making it impossible to distinguish genuine line drift
        # from a precision/rounding issue from a protocol mismatch.
        _log_placement_response(bet, confirm_payload, confirm_json,
                                classification="price_moved",
                                reason=f"win_drift_{abs(float(wz_win) - float(expected_win)):.4f}",
                                endpoint="ConfirmWagerHelper")
        return _price_moved_by_win(
            expected=float(expected_win),
            actual=float(wz_win),
        )

    # -----------------------------------------------------------------------
    # Step 2: PostWagerMultipleHelper (real submission; same endpoint parlays
    # use, distinguished by WT=0 in the confirm_payload built above)
    # Network failure here is an orphan candidate — we log forensics and
    # return network_error so the dashboard can surface it.
    # -----------------------------------------------------------------------
    wager_request = dict(confirm_payload)
    # Normalize types + required keys to match parlay_placer._build_post_request
    # exactly. confirm_payload uses string-typed booleans ("false") and string
    # zeros ("0") because ConfirmWagerHelper (preflight) tolerates them when
    # the payload is form-encoded. But PostWagerMultipleHelper JSON-decodes
    # each wager inside `postWagerRequests`, where JSON "false" (string) is
    # not the same as JSON false (bool). Parlay's working request uses bool
    # types and includes `roundRobinCombinations: ""` — mirror it exactly.
    wager_request["open"]                   = 0
    wager_request["sameAmount"]             = False
    wager_request["useFreePlayAmount"]      = False
    wager_request["roundRobinCombinations"] = ""
    # confirmPassword is required by WZ to authorize the placement (mirrors
    # parlay_placer._build_post_request). The test path (wz_account is None)
    # skips it because mock sessions do not validate credentials.
    if wz_account is not None:
        wager_request["confirmPassword"] = wz_account.password

    # WZ's PostWagerMultipleHelper expects the wager(s) wrapped in a
    # `postWagerRequests` JSON-serialized array — exactly the same envelope
    # parlay_placer._post_wagers uses. Sending the wager dict flat returns
    # `result: {ErrorMessage: "No Payload"}` (verified 2026-05-22 against
    # a live account). For a single bet, the array has exactly one element.
    make_body = {"postWagerRequests": _json.dumps([wager_request])}
    try:
        make_resp = sess.post(
            MAKE_URL,
            data=make_body,
            timeout=20,
            headers={"Accept": "application/json"},
        )
    except requests.RequestException as e:
        return _network_error(f"{type(e).__name__}: {e}")

    if _is_html_response(make_resp):
        # No auth_error retry: single placer is stateless; caller is expected to
        # retry on auth_error. The parlay placer retries internally because parlay
        # state is more complex to reconstruct. Singles are cheap to retry.
        return _auth_error("Wagerzon returned HTML at submission (session expired)")

    make_json = make_resp.json()

    # Navigate WZ's PostWagerMultipleHelper response defensively, the same
    # way parlay_placer._post_wagers does. The shape we've observed for
    # parlays is:
    #     {"result": [ {"WagerPostResult": {Confirm, ErrorMsgKey,
    #                                        TicketNumber, IDWT, Win, Risk,
    #                                        details: [{...}]}} ]}
    # but for singles (WT=0) we don't have an end-to-end verified sample,
    # and the original prototype assumed a flat dict with `WagerNumber`.
    # Handle all of:
    #   - result as list (parlay-style)        → take result[0]
    #   - result as dict (legacy single)       → use directly
    #   - WagerPostResult wrapper present      → unwrap it
    #   - fields at wpr top-level OR inside    → check both via _first_present
    #     details[0]
    #   - legacy "WagerNumber" key alongside   → accept either as ticket
    #     the parlay-style "TicketNumber"
    raw_result = make_json.get("result")
    if isinstance(raw_result, list):
        item = raw_result[0] if raw_result else {}
    elif isinstance(raw_result, dict):
        item = raw_result
    else:
        item = {}

    wpr = item.get("WagerPostResult") if isinstance(item, dict) else None
    if not isinstance(wpr, dict):
        wpr = item if isinstance(item, dict) else {}

    inner_list = wpr.get("details") if isinstance(wpr, dict) else None
    inner = inner_list[0] if isinstance(inner_list, list) and inner_list and isinstance(inner_list[0], dict) else {}

    # Explicit rejection takes priority — surface WZ's real reason
    # (LINE_PULLED, MINWAGERONLINE, INSUFFICIENT_BALANCE, "No Payload",
    # etc.) so the user sees what's blocking the placement instead of a
    # generic orphan. WZ may surface a rejection with either a machine
    # key (ErrorMsgKey/ErrorCode) or a human message (ErrorMsg/
    # ErrorMessage), or both — accept any combination as a rejection
    # signal. When only the message is present, use a generic "wz_error"
    # key so the dashboard knows this is a rejection.
    err_key = _first_present(
        wpr.get("ErrorMsgKey"),   wpr.get("ErrorCode"),
        inner.get("ErrorMsgKey"), inner.get("ErrorCode"),
    )
    err_msg = _first_present(
        wpr.get("ErrorMsg"),   wpr.get("ErrorMessage"),
        inner.get("ErrorMsg"), inner.get("ErrorMessage"),
    )
    if err_key or err_msg:
        final_msg = str(err_msg) if err_msg else str(err_key)
        final_key = str(err_key) if err_key else "wz_error"
        _log_placement_response(bet, wager_request, make_json,
                                classification="rejected",
                                reason=final_key,
                                endpoint="PostWagerMultipleHelper")
        return _rejected(final_msg, final_key)

    # Ticket lookup at every plausible location. Includes the legacy
    # `WagerNumber` key as a fallback so a partially-correct WZ response
    # variant still recovers the ticket.
    ticket = _first_present(
        wpr.get("TicketNumber"),   inner.get("TicketNumber"),
        wpr.get("WagerNumber"),    inner.get("WagerNumber"),
    )
    balance_after = _first_present(
        wpr.get("AvailBalance"),   inner.get("AvailBalance"),
    )

    if ticket is None:
        # No error key, no ticket — WZ returned SOMETHING but neither a
        # success-confirmation nor a rejection reason we know how to read.
        # Could be a third response shape variant we haven't sampled yet.
        # Log the body so we can finally see what came back.
        _log_placement_response(bet, wager_request, make_json,
                                classification="orphaned",
                                reason="no_ticket_no_error_key",
                                endpoint="PostWagerMultipleHelper")
        return {
            "status": "orphaned",
            "ticket_number": None,
            "balance_after": balance_after,
            "error_msg": "Submitted to WZ but response had no ticket or error key",
            "error_msg_key": "no_ticket",
        }

    return {
        "status": "placed",
        "ticket_number": str(ticket),
        "balance_after": balance_after,
        "error_msg": None,
        "error_msg_key": None,
    }
