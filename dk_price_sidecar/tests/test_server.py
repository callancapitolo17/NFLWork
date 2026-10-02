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


class _ScriptedBrowser:
    """DkBrowser with the page replaced by a script of fetch statuses."""

    def __new__(cls, statuses):
        from dk_price_sidecar.server import DkBrowser

        class Scripted(DkBrowser):
            def __init__(self):
                super().__init__()
                self._page = object()          # ready() -> True
                self.statuses = list(statuses)
                self.events = []

            def _reload(self):
                self.events.append("reload")
                self.reloads += 1
                self.calls_since_reload = 0

            def _call(self, selection_ids):
                self.events.append("call")
                self.calls_since_reload += 1
                status = self.statuses.pop(0)
                text = YOURBET_200 if status == 200 else "TypeError: Failed to fetch"
                return parse_calculate_bets(status, text, len(selection_ids))

        browser = Scripted()
        browser._page_loaded_at = float("inf")  # never reload on the timer
        return browser


def test_reloads_before_the_page_budget_runs_out():
    from dk_price_sidecar.server import CALLS_PER_PAGE_LOAD
    browser = _ScriptedBrowser([200] * (CALLS_PER_PAGE_LOAD + 1))
    for _ in range(CALLS_PER_PAGE_LOAD + 1):
        assert browser.price(["a", "b"])["true_odds"] == 7.5
    assert browser.events == ["call"] * CALLS_PER_PAGE_LOAD + ["reload", "call"]


def test_a_denied_call_reloads_and_retries_once():
    browser = _ScriptedBrowser([0, 200])
    assert browser.price(["a", "b"])["true_odds"] == 7.5
    assert browser.events == ["call", "reload", "call"]
    assert browser.blocks == 1


def test_a_second_denial_is_returned_not_retried_forever():
    browser = _ScriptedBrowser([0, 0])
    assert browser.price(["a", "b"])["status"] == 0
    assert browser.events == ["call", "reload", "call"]
