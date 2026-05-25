"""Tests for single_placer.place_single's Step 2 response navigation.

Background: After the Win-on-Win drift-check fix (2026-05-19), single
placements got PAST preflight but the first end-to-end attempt
(ARI -2.5 alt-spread, 2026-05-20) came back as `orphaned: Submitted to
WZ but no WagerNumber returned`. The user confirmed the bet was NOT in
WZ's ticket history — meaning WZ took the submission and silently
rejected it.

The root cause was the same family of bug we just fixed at preflight:
the placer was guessing at WZ's response field names. Specifically:
  - parlay_placer reads `TicketNumber`, NOT `WagerNumber`.
  - parlay_placer navigates a `result: [{WagerPostResult: {...}}]`
    wrapper; the original single placer assumed a flat `result: {...}`.
  - Real rejections (e.g. LINE_PULLED) carry `ErrorMsgKey` inside
    `WagerPostResult`, not at the top level — so they were being
    masked into the orphan path.

These tests verify the defensive navigation handles every shape we know
about: parlay-style list-with-wrapper, dict-with-wrapper, and the
legacy flat shape with `WagerNumber`. Plus the orphan and rejection
paths now trigger response-body logging to gitignored files for any
future surprise.
"""
import json as _json
import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock

import pytest

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3] / "wagerzon_odds"))

import single_placer


# ---------------------------------------------------------------------------
# Mock-session helpers
# ---------------------------------------------------------------------------

def _json_resp(body: dict):
    resp = MagicMock()
    resp.headers = {"content-type": "application/json"}
    resp.json.return_value = body
    return resp


def _make_session(confirm_body: dict, make_body: dict):
    """First POST returns preflight, second returns placement."""
    sess = MagicMock()
    sess.post.side_effect = [_json_resp(confirm_body), _json_resp(make_body)]
    return sess


def _wz_bet(**overrides) -> dict:
    payload = dict(
        bet_hash="abc123",
        idgm=5685842, play=1, line=-2.5,
        american_odds=215, wz_odds_at_place=215,
        actual_size=205.0, expected_win=440.75,
        market="alternate_spreads_fg", bet_on="Cleveland Guardians",
        pitcher=0,
    )
    payload.update(overrides)
    return payload


# Standard preflight response that passes the Win-on-Win drift check.
# All Step-2 tests start from a clean preflight so they only exercise
# the placement-response navigation.
PASS_PREFLIGHT = {
    "result": {"details": [{"Win": 440.75, "Risk": 205.0}]},
}


# ---------------------------------------------------------------------------
# Response-shape navigation
# ---------------------------------------------------------------------------

def test_parlay_style_response_list_with_wrapper_recovers_ticket():
    """The shape parlay_placer relies on: result is a list of items,
    each wrapping a `WagerPostResult` that holds TicketNumber + IDWT."""
    make_body = {
        "result": [{
            "WagerPostResult": {
                "Confirm":      True,
                "TicketNumber": "T-1001",
                "IDWT":         42,
                "AvailBalance": 1000.0,
                "details":      [{"Win": 440.75, "Risk": 205.0}],
            }
        }]
    }
    sess = _make_session(PASS_PREFLIGHT, make_body)
    result = single_placer.place_single(account=None, bet=_wz_bet(), session=sess)
    assert result["status"] == "placed"
    assert result["ticket_number"] == "T-1001"
    assert result["balance_after"] == 1000.0


def test_dict_with_wrapper_recovers_ticket():
    """Same WagerPostResult wrapper but result is a dict, not a list.
    Defensive: WZ may serve singles with a non-list `result`."""
    make_body = {
        "result": {
            "WagerPostResult": {
                "Confirm":      True,
                "TicketNumber": "T-2002",
                "AvailBalance": 750.0,
            }
        }
    }
    sess = _make_session(PASS_PREFLIGHT, make_body)
    result = single_placer.place_single(account=None, bet=_wz_bet(), session=sess)
    assert result["status"] == "placed"
    assert result["ticket_number"] == "T-2002"


def test_legacy_flat_response_with_wagernumber_still_recovers():
    """The shape the original prototype guessed: flat result dict with
    `WagerNumber` (not `TicketNumber`). Kept as a fallback so we don't
    regress if WZ does happen to use this shape for some single shapes."""
    make_body = {
        "result": {
            "WagerNumber":  "T-LEGACY",
            "AvailBalance": 500.0,
        }
    }
    sess = _make_session(PASS_PREFLIGHT, make_body)
    result = single_placer.place_single(account=None, bet=_wz_bet(), session=sess)
    assert result["status"] == "placed"
    assert result["ticket_number"] == "T-LEGACY"


def test_ticket_inside_details_array_recovers():
    """parlay_placer sometimes finds TicketNumber inside details[0]
    rather than at the wpr top level — handle that case for singles too."""
    make_body = {
        "result": [{
            "WagerPostResult": {
                "Confirm": True,
                "details": [{
                    "TicketNumber": "T-3003",
                    "IDWT":         99,
                    "AvailBalance": 600.0,
                    "Win":          440.75,
                }],
            }
        }]
    }
    sess = _make_session(PASS_PREFLIGHT, make_body)
    result = single_placer.place_single(account=None, bet=_wz_bet(), session=sess)
    assert result["status"] == "placed"
    assert result["ticket_number"] == "T-3003"
    assert result["balance_after"] == 600.0


# ---------------------------------------------------------------------------
# Rejection unmasking (the main user-facing benefit)
# ---------------------------------------------------------------------------

def test_rejection_via_errormessage_only_returns_rejected(tmp_path, monkeypatch):
    """REGRESSION for 2026-05-22: WZ rejected our flat (un-wrapped)
    submission with `{result: {ErrorMessage: "No Payload"}}` — note no
    ErrorMsgKey. Old logic required a machine error key, so this was
    masked as 'orphaned' for two days. Now an ErrorMessage alone is a
    sufficient rejection signal; the placer surfaces the message and
    tags it with a generic 'wz_error' key so the dashboard renders it
    as a rejection in the user-facing toast."""
    monkeypatch.setattr(single_placer, "_DEBUG_DIR", tmp_path)
    make_body = {"result": {"ErrorMessage": "No Payload"}}
    sess = _make_session(PASS_PREFLIGHT, make_body)
    result = single_placer.place_single(account=None, bet=_wz_bet(), session=sess)
    assert result["status"] == "rejected"
    assert "No Payload" in result["error_msg"]
    # No machine key in the response → we fall back to a generic key so
    # the dashboard still routes this as rejected, not orphan.
    assert result["error_msg_key"] == "wz_error"


def test_request_body_uses_postwagerrequests_envelope():
    """REGRESSION for 2026-05-22 root cause: WZ's PostWagerMultipleHelper
    only accepts the wager(s) wrapped in `postWagerRequests=[...]`.
    Sending the wager dict flat returned `ErrorMessage: "No Payload"`.
    Mirror parlay_placer._post_wagers' envelope exactly.

    This test inspects what the placer actually POSTs to the MAKE_URL,
    asserts the envelope structure, and verifies the wager dict (one
    element) carries the same per-leg fields the preflight used."""
    preflight_body = PASS_PREFLIGHT
    make_body = {
        "result": [{
            "WagerPostResult": {
                "Confirm":      True,
                "TicketNumber": "T-9999",
            }
        }]
    }
    sess = _make_session(preflight_body, make_body)
    single_placer.place_single(account=None, bet=_wz_bet(), session=sess)

    # Two POSTs happened: [confirm, make]. Inspect the make POST's body.
    assert sess.post.call_count == 2
    make_call = sess.post.call_args_list[1]
    # Sent under the `data=` kwarg.
    data = make_call.kwargs.get("data") or (make_call.args[1] if len(make_call.args) > 1 else None)
    assert isinstance(data, dict), f"expected dict body, got {type(data).__name__}"
    assert "postWagerRequests" in data, \
        f"placement body missing postWagerRequests envelope: keys={list(data.keys())}"

    # The envelope value is a JSON-encoded array of wager dicts.
    decoded = _json.loads(data["postWagerRequests"])
    assert isinstance(decoded, list) and len(decoded) == 1
    wager = decoded[0]
    # Sanity: the single wager carries the WT=0 marker (singles vs parlays)
    # and the same `sel` shape parlay_placer's encode_sel produces per leg.
    assert wager["WT"] == "0"
    assert wager["IDWT"] == "0"
    assert wager["sel"] == "1_5685842_-2.5_215"  # matches _wz_bet defaults
    # Types and required keys must exactly match parlay_placer._build_post_request.
    # PostWagerMultipleHelper JSON-decodes each wager strictly, so string-typed
    # booleans ("false") would not parse as the bool `false` and could be the
    # difference between accept and reject.
    assert wager["open"] == 0 and isinstance(wager["open"], int)
    assert wager["sameAmount"] is False
    assert wager["useFreePlayAmount"] is False
    assert wager["roundRobinCombinations"] == ""


def test_rejection_inside_wagerpostresult_surfaces_real_error(tmp_path, monkeypatch):
    """The bug-of-the-day: WZ returns a real rejection inside
    WagerPostResult.ErrorMsgKey, but the old code only checked the top
    level → silently became an orphan. Now we navigate into it and
    surface the real reason (e.g., LINE_PULLED, MINWAGERONLINE)."""
    monkeypatch.setattr(single_placer, "_DEBUG_DIR", tmp_path)
    make_body = {
        "result": [{
            "WagerPostResult": {
                "Confirm":      False,
                "ErrorMsgKey":  "LINE_PULLED",
                "ErrorMsg":     "Line is no longer available",
            }
        }]
    }
    sess = _make_session(PASS_PREFLIGHT, make_body)
    result = single_placer.place_single(account=None, bet=_wz_bet(), session=sess)
    assert result["status"] == "rejected"
    assert result["error_msg_key"] == "LINE_PULLED"
    assert "no longer available" in result["error_msg"]
    # A debug-file should also have been written so we have the body.
    debug_files = list(tmp_path.glob(".placement_debug_*.json"))
    assert len(debug_files) == 1


def test_rejection_at_top_level_still_works(tmp_path, monkeypatch):
    """Older code path: ErrorMsgKey sits at the result top level (no
    WagerPostResult wrapper). The defensive read picks it up via the
    'item itself as wpr' fallback. monkeypatch _DEBUG_DIR so the
    rejection-path body-log doesn't litter the real wagerzon_odds dir."""
    monkeypatch.setattr(single_placer, "_DEBUG_DIR", tmp_path)
    make_body = {
        "result": {
            "ErrorMsgKey":  "INSUFFICIENT_BALANCE",
            "ErrorMessage": "Not enough funds",
        }
    }
    sess = _make_session(PASS_PREFLIGHT, make_body)
    result = single_placer.place_single(account=None, bet=_wz_bet(), session=sess)
    assert result["status"] == "rejected"
    assert result["error_msg_key"] == "INSUFFICIENT_BALANCE"


# ---------------------------------------------------------------------------
# Orphan path + body logging
# ---------------------------------------------------------------------------

def test_no_ticket_no_error_key_returns_orphaned_and_logs_body(tmp_path, monkeypatch):
    """The exact case the user hit on 2026-05-20: submission succeeded,
    but the response had neither a ticket nor an error key in any
    location we know to check. Return orphaned and dump the response
    body to a debug file so we can finally see what came back."""
    monkeypatch.setattr(single_placer, "_DEBUG_DIR", tmp_path)
    make_body = {"result": {"some_unexpected_key": "some_value"}}
    sess = _make_session(PASS_PREFLIGHT, make_body)
    result = single_placer.place_single(account=None, bet=_wz_bet(), session=sess)
    assert result["status"] == "orphaned"
    assert result["error_msg_key"] == "no_ticket"

    debug_files = list(tmp_path.glob(".placement_debug_*.json"))
    assert len(debug_files) == 1
    body = _json.loads(debug_files[0].read_text())
    assert body["classification"] == "orphaned"
    assert body["response_body"]["result"]["some_unexpected_key"] == "some_value"
    # bet_hash from the bet payload is recorded for cross-referencing.
    assert body["bet_summary"]["bet_hash"] == "abc123"


def test_log_redacts_sensitive_keys(tmp_path, monkeypatch):
    """If WZ ever echoes a password (the balance endpoint is known to do
    this), our debug log MUST NOT persist it. confirmPassword on the
    request side must also be stripped — we send the user's password
    every placement."""
    monkeypatch.setattr(single_placer, "_DEBUG_DIR", tmp_path)
    # Simulate WZ echoing a password back (balance-endpoint pattern).
    make_body = {
        "result": {
            "WagerPostResult": {
                "Confirm":  False,
                "ErrorMsgKey": "AUTH",
                "Password":    "super-secret",
                "Player":      "username",
            }
        }
    }
    sess = _make_session(PASS_PREFLIGHT, make_body)
    single_placer.place_single(account=None, bet=_wz_bet(), session=sess)
    debug_file = next(tmp_path.glob(".placement_debug_*.json"))
    persisted = _json.loads(debug_file.read_text())
    persisted_str = _json.dumps(persisted)
    assert "super-secret" not in persisted_str
    assert "<redacted>" in persisted_str


def test_log_failure_does_not_break_placement(tmp_path, monkeypatch):
    """If the debug-file write fails for any reason (disk full, perms,
    path nonexistent), the placement flow must still complete cleanly.
    Set _DEBUG_DIR to an unwriteable path and confirm we still get a
    proper orphaned status."""
    bad_dir = tmp_path / "nonexistent" / "does_not_exist"
    monkeypatch.setattr(single_placer, "_DEBUG_DIR", bad_dir)
    make_body = {"result": {"empty": "no_ticket_no_error"}}
    sess = _make_session(PASS_PREFLIGHT, make_body)
    result = single_placer.place_single(account=None, bet=_wz_bet(), session=sess)
    # Function still returns the orphaned dict cleanly even if logging blew up.
    assert result["status"] == "orphaned"


# ---------------------------------------------------------------------------
# _first_present unit checks (the helper backing the defensive read)
# ---------------------------------------------------------------------------

def test_first_present_returns_first_non_empty():
    assert single_placer._first_present(None, "", 0, "x", "y") == "x"


def test_first_present_returns_none_when_all_empty():
    assert single_placer._first_present(None, "", 0, 0.0) is None


def test_first_present_treats_string_zero_as_present():
    """Edge case: a ticket number serialized as the string '0' is still
    a real value, not the int sentinel. Cover the boundary explicitly."""
    assert single_placer._first_present(None, "0") == "0"


# ---------------------------------------------------------------------------
# Preflight (ConfirmWagerHelper) failure logging — added 2026-05-22 to close
# the diagnostic gap on price_moved / empty_details / missing_win paths.
# Without these, when a placement fails at preflight we have no visibility
# into what WZ actually returned (vs Step 2 failures, which were already
# logged).
# ---------------------------------------------------------------------------

def test_price_moved_writes_debug_file_capturing_preflight_response(tmp_path, monkeypatch):
    """REGRESSION for 2026-05-22: user reported repeated 'Price moved' on
    bets where the line had not visibly moved at WZ. Previously the
    price_moved path wrote nothing to disk, so we had no way to see what
    Win WZ actually returned vs what expected_win the dashboard sent.

    This test verifies that price_moved now writes the full preflight
    response body (request + response) to a gitignored debug file with
    classification='price_moved' so a future investigation can read it.
    """
    monkeypatch.setattr(single_placer, "_DEBUG_DIR", tmp_path)
    # Preflight returns Win=$50 but expected_win is $440.75 -> drift fires.
    preflight = {"result": {"details": [{"Win": 50.0, "Risk": 205.0}]}}
    # Use a session that only returns the preflight response — the placer
    # should refuse to call Step 2 because the drift check trips first.
    sess = MagicMock()
    sess.post.side_effect = [_json_resp(preflight),
                              AssertionError("Step 2 must not fire on price_moved")]
    result = single_placer.place_single(
        account=None, bet=_wz_bet(expected_win=440.75), session=sess)
    assert result["status"] == "price_moved"

    # Find the debug file written for this rejection.
    debug_files = list(tmp_path.glob(".placement_debug_*_price_moved_*.json"))
    assert len(debug_files) == 1, \
        f"expected exactly one price_moved debug file, got {debug_files}"
    body = _json.loads(debug_files[0].read_text())
    assert body["classification"] == "price_moved"
    assert body["endpoint"] == "ConfirmWagerHelper"
    # The captured response_body should be exactly the preflight JSON,
    # giving the investigator everything WZ returned at preflight.
    assert body["response_body"]["result"]["details"][0]["Win"] == 50.0
    # The reason field encodes the drift amount so we can see at a glance
    # how far off we were without opening the file.
    assert body["reason"].startswith("win_drift_")


def test_empty_details_writes_debug_file(tmp_path, monkeypatch):
    """Symmetric coverage: empty_details at preflight also logs."""
    monkeypatch.setattr(single_placer, "_DEBUG_DIR", tmp_path)
    preflight = {"result": {"details": []}}  # WZ returned empty details
    sess = MagicMock()
    sess.post.side_effect = [_json_resp(preflight),
                              AssertionError("Step 2 must not fire")]
    result = single_placer.place_single(
        account=None, bet=_wz_bet(), session=sess)
    assert result["status"] == "rejected"
    assert result["error_msg_key"] == "empty_details"
    debug_files = list(tmp_path.glob(".placement_debug_*_empty_details_*.json"))
    assert len(debug_files) == 1
    body = _json.loads(debug_files[0].read_text())
    assert body["endpoint"] == "ConfirmWagerHelper"


def test_missing_win_at_preflight_writes_debug_file(tmp_path, monkeypatch):
    """If WZ's preflight returns details[0] without a Win field at all,
    we now log the body so we can finally see what was returned."""
    monkeypatch.setattr(single_placer, "_DEBUG_DIR", tmp_path)
    preflight = {"result": {"details": [{"Risk": 100.0}]}}  # no Win
    sess = MagicMock()
    sess.post.side_effect = [_json_resp(preflight),
                              AssertionError("Step 2 must not fire")]
    result = single_placer.place_single(
        account=None, bet=_wz_bet(), session=sess)
    assert result["status"] == "rejected"
    assert result["error_msg_key"] == "missing_win"
    debug_files = list(tmp_path.glob(".placement_debug_*_missing_win_*.json"))
    assert len(debug_files) == 1
    body = _json.loads(debug_files[0].read_text())
    assert body["endpoint"] == "ConfirmWagerHelper"


def test_preflight_error_key_writes_debug_file(tmp_path, monkeypatch):
    """A preflight that returns an ErrorMsgKey (e.g. INSUFFICIENT_BALANCE)
    also logs now — previously this path returned _rejected with no
    forensics. Useful when WZ refuses for a reason we haven't seen before."""
    monkeypatch.setattr(single_placer, "_DEBUG_DIR", tmp_path)
    preflight = {"result": {"ErrorMsgKey": "MINWAGERONLINE",
                              "ErrorMsg": "Below minimum wager"}}
    sess = MagicMock()
    sess.post.side_effect = [_json_resp(preflight),
                              AssertionError("Step 2 must not fire")]
    result = single_placer.place_single(
        account=None, bet=_wz_bet(), session=sess)
    assert result["status"] == "rejected"
    assert result["error_msg_key"] == "MINWAGERONLINE"
    debug_files = list(tmp_path.glob(".placement_debug_*_rejected_*.json"))
    assert len(debug_files) == 1
    body = _json.loads(debug_files[0].read_text())
    assert body["endpoint"] == "ConfirmWagerHelper"
    assert body["reason"] == "MINWAGERONLINE"


# ---------------------------------------------------------------------------
# Amount-field semantics — added 2026-05-22 after capturing live evidence
# that WZ at RiskWin=0 interprets `Amount` as the SMALLER of (Risk, Win).
# At negative odds (favorites) the smaller value is Win, not Risk; sending
# the user's Risk as Amount caused WZ to compute the "wrong" Risk and Win
# back, which our drift check then false-flagged as price_moved.
# ---------------------------------------------------------------------------

def test_compute_wz_amount_negative_odds_returns_win():
    """At -110 odds and $22 risk: Win = 22 × 100/110 = $20.00. WZ wants
    that as Amount when RiskWin=0 because Win < Risk."""
    assert single_placer._compute_wz_amount(22.0, -110) == 20.0


def test_compute_wz_amount_positive_odds_returns_risk_unchanged():
    """At +125 odds and $22 risk: Win would be $27.50; Risk ($22) is the
    smaller value, so Amount = Risk = $22 unchanged."""
    assert single_placer._compute_wz_amount(22.0, 125) == 22.0


def test_compute_wz_amount_even_money_returns_risk():
    """At ±100 odds, Risk == Win exactly; either works as Amount. Use
    Risk for consistency with the positive-odds branch."""
    assert single_placer._compute_wz_amount(22.0, 100) == 22.0
    assert single_placer._compute_wz_amount(22.0, -100) == 22.0


def test_compute_wz_amount_heavy_favorite():
    """At -300 odds, $30 risk wins $10. Amount = $10."""
    assert single_placer._compute_wz_amount(30.0, -300) == 10.0


def test_build_confirm_payload_sends_win_as_amount_at_negative_odds():
    """REGRESSION for the 2026-05-22 price_moved storm: at negative odds
    the placer must send Win (not Risk) in Amount and sameAmountNumber,
    so WZ's preflight returns the matching (Risk, Win) pair our verified
    expected_win was computed against."""
    bet = dict(
        idgm=5697583, play=3, line=4.5,
        american_odds=-110, wz_odds_at_place=-110,
        actual_size=22.0,
        pitcher=0,
    )
    payload = single_placer._build_confirm_payload(bet)
    detail = _json.loads(payload["detailData"])[0]
    # At -110 with risk=$22, WZ wants Amount=Win=$20, not Amount=Risk=$22.
    assert detail["Amount"] == "20"
    assert payload["sameAmountNumber"] == "20"


def test_build_confirm_payload_sends_risk_as_amount_at_positive_odds():
    """At positive odds (underdog), Risk is the smaller value — keep the
    behavior that already worked (CLE -2.5 at +215, ARI -2.5 at +190
    both reached Step 2 with the old code)."""
    bet = dict(
        idgm=5685842, play=1, line=-2.5,
        american_odds=215, wz_odds_at_place=215,
        actual_size=205.0,
        pitcher=0,
    )
    payload = single_placer._build_confirm_payload(bet)
    detail = _json.loads(payload["detailData"])[0]
    assert detail["Amount"] == "205"
    assert payload["sameAmountNumber"] == "205"


def test_negative_odds_placement_no_longer_false_fires_price_moved():
    """End-to-end: simulate the exact 2026-05-22 case. User risks $22 at
    -110, expected_win=$20. With the fix, the placer sends Amount=$20 to
    WZ; WZ's preflight returns Risk=$22, Win=$20; drift = $0 → proceeds
    to Step 2 cleanly. Without the fix, WZ would have returned Win=$22
    (the old bug), tripping the drift check."""
    # Mirror what the WZ recon-style happy path returns when Amount=Win=20.
    preflight = {"result": {"details": [{"Win": 20.0, "Risk": 22.0,
                                          "Odds": -110}]}}
    make_body = {"result": [{"WagerPostResult":
                              {"Confirm": True, "TicketNumber": "T-FAV-OK"}}]}
    sess = _make_session(preflight, make_body)
    bet = _wz_bet(
        idgm=5697583, play=3, line=4.5,
        american_odds=-110, wz_odds_at_place=-110,
        actual_size=22.0, expected_win=20.0,
    )
    result = single_placer.place_single(account=None, bet=bet, session=sess)
    assert result["status"] == "placed", \
        f"favorite-odds placement should now succeed, got {result}"
    assert result["ticket_number"] == "T-FAV-OK"

    # And confirm the request body actually sent Amount=20 (not 22).
    confirm_call = sess.post.call_args_list[0]
    confirm_body = confirm_call.kwargs.get("data") or confirm_call.args[1]
    detail = _json.loads(confirm_body["detailData"])[0]
    assert detail["Amount"] == "20", \
        f"placer should have sent Amount=20 (= Win) for -110, sent {detail['Amount']}"
