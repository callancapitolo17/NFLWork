"""Unit tests for single_placer.place_single's Win-on-Win drift check.

Background: WZ's ConfirmWagerHelper preflight response for SOME single-bet
shapes (e.g. alt-spreads, F-period totals) returns `details[0]` without an
`Odds` key. The placer used to require `Odds` and rejected all such bets.
It now drift-checks on `Win` (dollars), identical to the parlay placer.

These tests mock the requests.Session that single_placer would normally
get from wagerzon_auth, so we can control exactly what WZ "returns" at
each step. The placer's `session=` injection path makes this clean.

Two-step protocol the tests exercise:
  Step 1: POST /wager/ConfirmWagerHelper.aspx (preflight, price-match)
  Step 2: POST /wager/PostWagerMultipleHelper.aspx (actual placement)

Tests:
  1. expected_win supplied + matches preflight -> proceeds to Step 2
  2. expected_win absent + math fallback matches preflight -> proceeds
  3. expected_win supplied + Win drifted -> price_moved with both values
  4. Win absent from details[0] -> missing_win
  5. REGRESSION: Win present, Odds absent -> proceeds (the CLE -2.5 case)
  6-8. _compute_expected_win unit cases for +/-/even American odds.
"""
import json
import sys
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
    """Build a MagicMock that quacks like a requests.Response carrying JSON."""
    resp = MagicMock()
    resp.headers = {"content-type": "application/json"}
    resp.json.return_value = body
    return resp


def _make_session(confirm_body: dict, make_body: dict | None = None):
    """Build a MagicMock requests.Session that returns `confirm_body` on the
    first POST (preflight) and `make_body` on the second POST (placement).
    If `make_body` is None, the second POST raises — the test should not get
    that far."""
    sess = MagicMock()
    side_effects = [_json_resp(confirm_body)]
    if make_body is not None:
        side_effects.append(_json_resp(make_body))
    else:
        side_effects.append(AssertionError("placement POST should not be reached"))
    sess.post.side_effect = side_effects
    return sess


def _wz_bet(**overrides) -> dict:
    """A baseline WZ single-bet payload matching the CLE -2.5 card shape."""
    payload = dict(
        idgm=5685842,
        play=1,
        line=-2.5,
        american_odds=215,
        wz_odds_at_place=215,
        actual_size=205.0,
        pitcher=0,
    )
    payload.update(overrides)
    return payload


PLACED_BODY_OK = {
    "result": {
        "WagerNumber": "T-1001",
        "AvailBalance": 1000.0,
    }
}


# ---------------------------------------------------------------------------
# Drift-check behavior
# ---------------------------------------------------------------------------

def test_expected_win_supplied_and_matches_proceeds_to_placement():
    """Happy path: dashboard supplied WZ-verified expected_win=440.75 and
    preflight returns Win=440.75 (same number, same endpoint). Drift check
    passes; placer proceeds to Step 2 and returns 'placed'."""
    preflight = {"result": {"details": [{"Win": 440.75, "Risk": 205.0}]}}
    sess = _make_session(preflight, PLACED_BODY_OK)
    result = single_placer.place_single(
        account=None, bet=_wz_bet(expected_win=440.75), session=sess)
    assert result["status"] == "placed"
    assert result["ticket_number"] == "T-1001"


def test_expected_win_absent_falls_back_to_math_and_matches():
    """When the dashboard didn't supply expected_win (user clicked Place
    without verifying), the placer computes from odds + risk. For 215 odds
    at $205 risk: 205 * 2.15 = $440.75. Preflight returns the same → place."""
    preflight = {"result": {"details": [{"Win": 440.75, "Risk": 205.0}]}}
    sess = _make_session(preflight, PLACED_BODY_OK)
    bet = _wz_bet()  # no expected_win
    assert "expected_win" not in bet
    result = single_placer.place_single(account=None, bet=bet, session=sess)
    assert result["status"] == "placed"


def test_expected_win_supplied_with_real_drift_returns_price_moved():
    """Dashboard quoted user $440.75 to win; preflight now says $410. Drift
    is $30.75 — well above the $0.01 tolerance. Refuse to place; surface
    both dollar values in the error so the user sees what changed."""
    preflight = {"result": {"details": [{"Win": 410.00, "Risk": 205.0}]}}
    sess = _make_session(preflight, make_body=None)  # Step 2 must not fire
    result = single_placer.place_single(
        account=None, bet=_wz_bet(expected_win=440.75), session=sess)
    assert result["status"] == "price_moved"
    assert result["error_msg_key"] == "drift"
    assert "440.75" in result["error_msg"]
    assert "410.00" in result["error_msg"]


def test_win_absent_returns_missing_win():
    """If WZ's preflight returns details[0] with no Win field at all, we
    can't price-match — refuse to place. (Parlay placer relies on Win
    existing too; if this ever fires it's a genuinely new WZ behavior.)"""
    preflight = {"result": {"details": [{"Risk": 205.0}]}}  # no Win
    sess = _make_session(preflight, make_body=None)
    result = single_placer.place_single(
        account=None, bet=_wz_bet(expected_win=440.75), session=sess)
    assert result["status"] == "rejected"
    assert result["error_msg_key"] == "missing_win"


def test_odds_absent_but_win_present_proceeds_regression():
    """REGRESSION TEST for the CLE -2.5 bug: WZ's preflight returned
    details[0] with Win but no Odds. The old placer rejected with
    missing_odds. The new placer ignores the absent Odds field and uses
    Win for drift-check."""
    # Same shape as the actual user-reported failure: Win present, Odds absent.
    preflight = {"result": {"details": [{"Win": 440.75, "Risk": 205.0}]}}
    assert "Odds" not in preflight["result"]["details"][0]
    sess = _make_session(preflight, PLACED_BODY_OK)
    result = single_placer.place_single(
        account=None, bet=_wz_bet(expected_win=440.75), session=sess)
    assert result["status"] == "placed", \
        f"expected 'placed', got {result['status']}: {result.get('error_msg')}"


def test_drift_within_tolerance_still_proceeds():
    """Exactly at $0.01 difference should be considered matching (the
    parlay tolerance is inclusive: abs(diff) <= 0.01)."""
    preflight = {"result": {"details": [{"Win": 440.76, "Risk": 205.0}]}}
    sess = _make_session(preflight, PLACED_BODY_OK)
    result = single_placer.place_single(
        account=None, bet=_wz_bet(expected_win=440.75), session=sess)
    assert result["status"] == "placed"


def test_supplied_expected_win_takes_priority_over_math_fallback():
    """Tightening test: prove the dashboard-supplied expected_win is the
    value actually used, not just coincidentally equal to the math.

    Setup: bet has wz_odds_at_place=215 + actual_size=205, so the math
    fallback would compute $440.75. We supply expected_win=300 (a value
    the math would NEVER produce for this bet) and have WZ return Win=300.

      - If supplied is used (correct): expected=300 vs WZ=300 -> place.
      - If math is silently used (broken): expected=440.75 vs WZ=300
        -> $140.75 drift -> price_moved.

    Asserting 'placed' proves the supplied value is what the drift check
    consumed."""
    preflight = {"result": {"details": [{"Win": 300.00, "Risk": 205.0}]}}
    sess = _make_session(preflight, PLACED_BODY_OK)
    result = single_placer.place_single(
        account=None,
        bet=_wz_bet(expected_win=300.00),
        session=sess,
    )
    assert result["status"] == "placed", \
        f"supplied expected_win was ignored; got {result}"


# ---------------------------------------------------------------------------
# _compute_expected_win — pure math
# ---------------------------------------------------------------------------

def test_compute_expected_win_positive_odds():
    """+215 at $205 risk → $440.75 (matches the screenshot card)."""
    assert single_placer._compute_expected_win(odds=215, risk=205.0) == 440.75


def test_compute_expected_win_negative_odds():
    """-110 at $100 risk → $90.91 (standard juice example)."""
    assert single_placer._compute_expected_win(odds=-110, risk=100.0) == 90.91


def test_compute_expected_win_even_money():
    """+100 or -100 at $100 risk → $100.00 either way."""
    assert single_placer._compute_expected_win(odds=100, risk=100.0) == 100.0
    assert single_placer._compute_expected_win(odds=-100, risk=100.0) == 100.0
