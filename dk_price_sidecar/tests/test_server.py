"""dk_price_sidecar — pure-function tests. No browser, no network."""
from __future__ import annotations

import json

from dk_price_sidecar.server import build_calculate_bets_body, parse_calculate_bets

YOURBET_200 = json.dumps({
    "bets": [
        {"type": "Single", "selectionsMapped": [{"id": "a"}], "trueOdds": 3.95},
        {"type": "Single", "selectionsMapped": [{"id": "b"}], "trueOdds": 2.42},
        {"type": "YourBet", "selectionsMapped": [{"id": "a"}, {"id": "b"}],
         "trueOdds": 7.5, "displayOdds": "+650"},
    ],
    "combinabilityRestrictions": [],
})
SINGLES_ONLY_200 = json.dumps({
    "bets": [
        {"type": "Single", "selectionsMapped": [{"id": "a"}], "trueOdds": 2.22},
        {"type": "Single", "selectionsMapped": [{"id": "b"}], "trueOdds": 12.4},
    ],
    "combinabilityRestrictions": [
        {"selections": ["a", "b"], "restrictionType": "NonCombinableGroup"}],
})


def test_body_is_the_production_shape():
    body = build_calculate_bets_body(["a", "b"])
    assert body["selectionsForYourBet"] == [{"id": "a", "yourBetGroup": 0},
                                            {"id": "b", "yourBetGroup": 0}]
    assert body["selections"] == []
    assert body["oddsStyle"] == "american"


def test_yourbet_price_is_extracted_for_the_full_leg_set():
    out = parse_calculate_bets(200, YOURBET_200, n_legs=2)
    assert out["true_odds"] == 7.5
    assert out["display"] == "+650"
    assert out["restrictions"] is None


def test_singles_only_with_restriction_is_a_decline_not_a_price():
    """Captured live 2026-09-02: a non-combinable pair returns two Singles
    and a NonCombinableGroup. That must never be reported as a price."""
    out = parse_calculate_bets(200, SINGLES_ONLY_200, n_legs=2)
    assert out["true_odds"] is None
    assert out["restrictions"][0]["restrictionType"] == "NonCombinableGroup"


def test_single_leg_bet_does_not_satisfy_a_two_leg_request():
    body = json.dumps({"bets": [{"type": "Single", "selectionsMapped": [{"id": "a"}],
                                 "trueOdds": 2.0}]})
    assert parse_calculate_bets(200, body, n_legs=2)["true_odds"] is None


def test_non_200_carries_status_and_error_snippet():
    out = parse_calculate_bets(403, "<HTML>Access Denied", n_legs=2)
    assert out["status"] == 403 and out["true_odds"] is None
    assert "Access Denied" in out["error"]


def test_thrown_fetch_is_status_zero():
    out = parse_calculate_bets(0, "TypeError: Failed to fetch", n_legs=2)
    assert out["status"] == 0 and out["true_odds"] is None


def test_non_json_200_is_not_a_price():
    out = parse_calculate_bets(200, "<html>challenge</html>", n_legs=2)
    assert out["true_odds"] is None and out["error"] == "non-JSON body"
