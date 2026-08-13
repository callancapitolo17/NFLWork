"""Issue #87: SGPService.price_on_demand admits I1 (1st inning) legs.

RFI legs re-encode at parse as I1 totals at 0.5 (legset), so the period
whitelist must include "I1" for them to reach the book hooks. Unknown
periods (and F5 moneylines) stay CLOSED with zero wire calls — same
contract as test_sgp_service_f5_guard.
"""
from kalshi_common.legset import CanonicalLeg
from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import GameRef, ResolvedLeg

EVT = "26AUG132210MILLAD"
GAME = GameRef(game_id="g1", home_team="Los Angeles Dodgers",
               away_team="Milwaukee Brewers", commence_time=None)

I1_OVER = CanonicalLeg(EVT, "total", 0.5, "over", "I1")


def _tracking_hooks(calls):
    def price(client, refs, event):
        calls.append(list(refs))
        return 2.0

    return {
        "match_event": lambda client, game: calls.append("match") or "EV",
        "build_structure": lambda client, ev, game: {"fake": "structure"},
        "resolve": lambda structure, legs, home, away: [
            ResolvedLeg(ref=f"r{j}", opposite_ref=f"o{j}",
                        single_decimal=1.9, opposite_decimal=1.9)
            for j in range(len(legs))],
        "price": price,
    }


def _service(calls):
    return SGPService(books=("novig",),
                      on_demand_hooks={"novig": _tracking_hooks(calls)})


def test_i1_total_plus_fg_ml_prices_through_hooks():
    calls = []
    svc = _service(calls)
    legs = [I1_OVER, CanonicalLeg(EVT, "ml", None, "home")]
    res = svc.price_on_demand("novig", GAME, legs)
    assert res is not None
    assert len(calls) > 0


def test_lone_i1_leg_prices_through_hooks():
    """A cross-game combo's lone RFI partition leg prices as a 1-leg set."""
    calls = []
    svc = _service(calls)
    res = svc.price_on_demand("novig", GAME, [I1_OVER])
    assert res is not None
    assert len(calls) > 0


def test_unknown_period_still_declines_without_touching_book_hooks():
    calls = []
    svc = _service(calls)
    legs = [CanonicalLeg(EVT, "total", 1.5, "over", "I2"),
            I1_OVER]
    assert svc.price_on_demand("novig", GAME, legs) is None
    assert calls == []


def test_f5_moneyline_still_declines():
    calls = []
    svc = _service(calls)
    legs = [CanonicalLeg(EVT, "ml", None, "tie", "F5"), I1_OVER]
    assert svc.price_on_demand("novig", GAME, legs) is None
    assert calls == []
