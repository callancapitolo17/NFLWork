"""FanDuel on-demand pricing (kalshi_mlb_mm Phase 2).

Covers three additions:
  1. ``fetch_event_runners`` captures each runner's single decimal odds —
     buckets now store ``(marketId, selectionId, decimal_or_None)``.
     Odds path (verified against fixtures/fd_runners.json and the SGP
     odds shape at scraper_fanduel_sgp.py:482):
     ``winRunnerOdds.trueOdds.decimalOdds.decimalOdds``.
  2. ``resolve_legs(structure, legs, home_team, away_team)`` — canonical
     legs -> ResolvedLeg refs against the fg structure. Spread keys are
     ``(side, signed_line)`` with the away key at the NEGATED
     home-perspective line (mirrors price_sgps: home_line=t.spread,
     away_line=-t.spread).
  3. ``price_selection_set(client, refs)`` — N-leg implyBets, decimal or
     None on any failure.

No network: fake sessions/clients throughout (style of test_fd_moneyline.py).
"""
import json
import sys
from pathlib import Path
from unittest.mock import MagicMock

sys.path.insert(0, str(Path(__file__).parent.parent))

from kalshi_common.legset import CanonicalLeg  # noqa: E402
from mlb_sgp._shared import ResolvedLeg  # noqa: E402


# --------------------------------------------------------------------------
# 1. fetch_event_runners odds capture
# --------------------------------------------------------------------------

_EVENT_PAGE = {
    "attachments": {
        "markets": {
            "1": {"marketId": "734.rl", "marketName": "Run Line", "runners": [
                {"selectionId": 101, "runnerName": "Home Team", "handicap": -1.5,
                 "winRunnerOdds": {"trueOdds": {"decimalOdds": {"decimalOdds": 2.05}}}},
                {"selectionId": 102, "runnerName": "Away Team", "handicap": 1.5,
                 "winRunnerOdds": {"trueOdds": {"decimalOdds": {"decimalOdds": 1.78}}}},
            ]},
            "2": {"marketId": "734.tot", "marketName": "Total Runs", "runners": [
                {"selectionId": 201, "runnerName": "Over", "handicap": 8.5,
                 "winRunnerOdds": {"trueOdds": {"decimalOdds": {"decimalOdds": 1.91}}}},
                # No winRunnerOdds at all -> odds slot must be None.
                {"selectionId": 202, "runnerName": "Under", "handicap": 8.5},
            ]},
            "3": {"marketId": "734.ml", "marketName": "Moneyline", "runners": [
                {"selectionId": 301, "runnerName": "Home Team",
                 "winRunnerOdds": {"trueOdds": {"decimalOdds": {"decimalOdds": 1.65}}}},
                # winRunnerOdds present but trueOdds nesting absent -> None.
                {"selectionId": 302, "runnerName": "Away Team",
                 "winRunnerOdds": {"americanDisplayOdds": {"americanOdds": 195}}},
            ]},
        }
    }
}


def _fake_session(payload, status=200):
    resp = MagicMock()
    resp.status_code = status
    resp.json.return_value = payload
    session = MagicMock()
    session.get.return_value = resp
    session.post.return_value = resp
    return session


def test_fetch_event_runners_captures_single_odds():
    from scraper_fanduel_sgp import fetch_event_runners
    out = fetch_event_runners(_fake_session(_EVENT_PAGE), "ev1",
                              "Home Team", "Away Team")
    fg = out["fg"]
    assert fg["spreads"][("home", -1.5)] == ("734.rl", 101, 2.05)
    assert fg["spreads"][("away", 1.5)] == ("734.rl", 102, 1.78)
    assert fg["totals"][("O", 8.5)] == ("734.tot", 201, 1.91)
    # Runner without winRunnerOdds: ref still captured, odds slot None.
    assert fg["totals"][("U", 8.5)] == ("734.tot", 202, None)
    assert fg["moneyline"]["home"] == ("734.ml", 301, 1.65)
    # winRunnerOdds present but missing the trueOdds nesting -> None.
    assert fg["moneyline"]["away"] == ("734.ml", 302, None)


def test_fetch_event_runners_malformed_odds_is_none():
    """Defensive: garbage winRunnerOdds types must not raise."""
    from scraper_fanduel_sgp import fetch_event_runners
    page = {"markets": [
        {"marketId": "m1", "marketName": "Moneyline", "runners": [
            {"selectionId": 1, "runnerName": "Home Team",
             "winRunnerOdds": "not-a-dict"},
            {"selectionId": 2, "runnerName": "Away Team",
             "winRunnerOdds": {"trueOdds": {"decimalOdds": "nope"}}},
        ]},
    ]}
    out = fetch_event_runners(_fake_session(page), "ev1",
                              "Home Team", "Away Team")
    assert out["fg"]["moneyline"]["home"] == ("m1", 1, None)
    assert out["fg"]["moneyline"]["away"] == ("m1", 2, None)


# --------------------------------------------------------------------------
# 2. resolve_legs
# --------------------------------------------------------------------------

# fetch_event_runners(...)["fg"] shape, post-odds-capture. Spreads keyed
# (side, SIGNED line): home -1.5 favored pairs with away +1.5.
STRUCT = {
    "spreads": {
        ("home", -1.5): ("734.rl", 101, 2.05),
        ("away", 1.5): ("734.rl", 102, 1.78),
        # One-sided alt: home +2.5 offered, away -2.5 missing.
        ("home", 2.5): ("734.arl", 103, 1.30),
    },
    "totals": {
        ("O", 8.5): ("734.tot", 201, 1.91),
        ("U", 8.5): ("734.tot", 202, 1.91),
        # One-sided alt total: Over 10.5 offered, Under missing.
        ("O", 10.5): ("734.atot", 203, 3.1),
    },
    "moneyline": {
        "home": ("734.ml", 301, 1.65),
        "away": ("734.ml", 302, 2.30),
    },
}

GID = "KXMLBGAME-26JUL10NYYBOS"
HOME, AWAY = "Home Team", "Away Team"


def _wrap(bucket):
    """The two-period structure resolve_legs consumes since #85 (F5
    coverage lives in test_f5_on_demand_resolvers.py)."""
    return {"FG": bucket, "F5": None}


def _resolve(legs, struct=STRUCT):
    from mlb_sgp.fanduel import resolve_legs
    return resolve_legs(_wrap(struct), legs, HOME, AWAY)


def test_resolve_spread_home_side():
    """Home spread -1.5: chosen key ('home', -1.5); opposite ('away', +1.5)."""
    out = _resolve([CanonicalLeg(GID, "spread", -1.5, "home")])
    assert out == [ResolvedLeg(ref=("734.rl", 101),
                               opposite_ref=("734.rl", 102),
                               single_decimal=2.05,
                               opposite_decimal=1.78)]


def test_resolve_spread_away_side_negates_line():
    """Away side of the SAME home-perspective line -1.5 maps to the
    ('away', +1.5) key — mirrors price_sgps's away_line = -t.spread."""
    out = _resolve([CanonicalLeg(GID, "spread", -1.5, "away")])
    assert out == [ResolvedLeg(ref=("734.rl", 102),
                               opposite_ref=("734.rl", 101),
                               single_decimal=1.78,
                               opposite_decimal=2.05)]


def test_resolve_total_over_and_under():
    out = _resolve([CanonicalLeg(GID, "total", 8.5, "over"),
                    CanonicalLeg(GID, "total", 8.5, "under")])
    assert out is not None and len(out) == 2
    assert out[0].ref == ("734.tot", 201)
    assert out[0].opposite_ref == ("734.tot", 202)
    assert out[0].single_decimal == 1.91
    assert out[1].ref == ("734.tot", 202)
    assert out[1].opposite_ref == ("734.tot", 201)


def test_resolve_moneyline_both_sides():
    out = _resolve([CanonicalLeg(GID, "ml", None, "home"),
                    CanonicalLeg(GID, "ml", None, "away")])
    assert out is not None
    assert out[0].ref == ("734.ml", 301)
    assert out[0].opposite_ref == ("734.ml", 302)
    assert out[0].single_decimal == 1.65
    assert out[0].opposite_decimal == 2.30
    assert out[1].ref == ("734.ml", 302)
    assert out[1].opposite_ref == ("734.ml", 301)


def test_resolve_chosen_side_missing_returns_none():
    """Any leg whose CHOSEN side is unresolvable fails the whole set."""
    legs = [CanonicalLeg(GID, "total", 8.5, "over"),      # resolvable
            CanonicalLeg(GID, "spread", -3.5, "home")]    # not offered
    assert _resolve(legs) is None


def test_resolve_opposite_missing_yields_opposite_ref_none():
    """home +2.5 offered but away -2.5 absent -> opposite_ref=None."""
    out = _resolve([CanonicalLeg(GID, "spread", 2.5, "home")])
    assert out == [ResolvedLeg(ref=("734.arl", 103),
                               opposite_ref=None,
                               single_decimal=1.30,
                               opposite_decimal=None)]


def test_resolve_opposite_missing_total():
    out = _resolve([CanonicalLeg(GID, "total", 10.5, "over")])
    assert out is not None
    assert out[0].ref == ("734.atot", 203)
    assert out[0].opposite_ref is None
    assert out[0].opposite_decimal is None


def test_resolve_none_odds_slot_is_none_safe():
    """A tuple whose odds slot is None resolves with single_decimal=None."""
    struct = {"spreads": {}, "totals": {},
              "moneyline": {"home": ("m", 1, None), "away": ("m", 2, 2.0)}}
    out = _resolve([CanonicalLeg(GID, "ml", None, "home")], struct)
    assert out == [ResolvedLeg(ref=("m", 1), opposite_ref=("m", 2),
                               single_decimal=None, opposite_decimal=2.0)]


def test_resolve_never_raises_on_garbage_structure():
    from mlb_sgp.fanduel import resolve_legs
    legs = [CanonicalLeg(GID, "spread", -1.5, "home")]
    assert resolve_legs(None, legs, HOME, AWAY) is None
    assert resolve_legs({}, legs, HOME, AWAY) is None
    assert resolve_legs({"spreads": "junk"}, legs, HOME, AWAY) is None
    assert resolve_legs(_wrap(STRUCT), [], HOME, AWAY) is None
    assert resolve_legs(
        _wrap(STRUCT), [CanonicalLeg(GID, "anytime_hr", None, "home")],
        HOME, AWAY) is None


def test_resolve_multi_leg_order_preserved():
    legs = [CanonicalLeg(GID, "spread", -1.5, "home"),
            CanonicalLeg(GID, "total", 8.5, "under"),
            CanonicalLeg(GID, "ml", None, "away")]
    out = _resolve(legs)
    assert [r.ref for r in out] == [("734.rl", 101), ("734.tot", 202),
                                    ("734.ml", 302)]


# --------------------------------------------------------------------------
# 3. price_selection_set
# --------------------------------------------------------------------------

def _imply_bets_response(decimal=5.6, is_sgm=True, n_legs=None):
    # Live FD responses carry one legCombinations entry per covered leg
    # (verified 2026-07-10); the leg-coverage guard requires it to equal
    # len(refs). Callers pass n_legs to match their refs.
    return {"betCombinations": [
        {"isSGM": is_sgm,
         "legCombinations": [{"legType": "SIMPLE_SELECTION"}] * (n_legs or 1),
         "winAvgOdds": {"trueOdds": {"decimalOdds": {"decimalOdds": decimal}},
                        "americanDisplayOdds": {"americanOddsInt": 460}}},
    ]}


def _fake_client(payload, status=200):
    client = MagicMock()
    client.session = _fake_session(payload, status=status)
    return client


def test_price_selection_set_three_leg_body_shape():
    from mlb_sgp.fanduel import price_selection_set
    refs = [("734.rl", 101), ("734.tot", 201), ("734.ml", 302)]
    client = _fake_client(_imply_bets_response(decimal=7.25, n_legs=3))

    out = price_selection_set(client, refs)
    assert out == 7.25

    # Body: one SIMPLE_SELECTION betLeg per ref, order preserved —
    # generalizes the 2-leg body of scraper_fanduel_sgp.price_combo.
    (_, kwargs) = client.session.post.call_args
    body = json.loads(kwargs["data"])
    assert len(body["betLegs"]) == 3
    for leg, (mid, sid) in zip(body["betLegs"], refs):
        assert leg["legType"] == "SIMPLE_SELECTION"
        assert leg["betRunners"] == [
            {"runner": {"marketId": mid, "selectionId": sid}}]


def test_price_selection_set_http_error_returns_none():
    from mlb_sgp.fanduel import price_selection_set
    client = _fake_client(_imply_bets_response(), status=403)
    assert price_selection_set(client, [("m", 1), ("m", 2)]) is None


def test_price_selection_set_no_sgm_entry_returns_none():
    from mlb_sgp.fanduel import price_selection_set
    client = _fake_client(_imply_bets_response(is_sgm=False))
    assert price_selection_set(client, [("m", 1), ("m", 2)]) is None


def test_price_selection_set_single_leg_accepts_non_sgm():
    """1-element refs prices a single: FD tags nothing isSGM for one leg,
    so the single combination's odds are accepted."""
    from mlb_sgp.fanduel import price_selection_set
    client = _fake_client(_imply_bets_response(decimal=1.91, is_sgm=False))
    assert price_selection_set(client, [("m", 1)]) == 1.91


def test_price_selection_set_rejects_partial_leg_coverage():
    """Adversarial review #10: an SGM combination that covers fewer legs
    than requested (FD silently dropped a suspended leg) must not be taken
    as the N-leg price."""
    from mlb_sgp.fanduel import price_selection_set
    client = _fake_client(_imply_bets_response(decimal=7.25, n_legs=2))
    refs = [("734.rl", 101), ("734.tot", 201), ("734.ml", 302)]
    assert price_selection_set(client, refs) is None


def test_price_selection_set_never_raises():
    from mlb_sgp.fanduel import price_selection_set
    client = MagicMock()
    client.session.post.side_effect = RuntimeError("boom")
    assert price_selection_set(client, [("m", 1), ("m", 2)]) is None

    bad_json = MagicMock()
    bad_json.session.post.return_value.status_code = 200
    bad_json.session.post.return_value.json.side_effect = ValueError("bad json")
    assert price_selection_set(bad_json, [("m", 1), ("m", 2)]) is None

    assert price_selection_set(_fake_client(_imply_bets_response()), []) is None
