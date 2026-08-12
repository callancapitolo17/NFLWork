"""Issue #84: SGPService.price_on_demand must fail CLOSED on non-FG legs.

Every book's ``resolve_legs`` dispatches on ``market_type`` alone against
FULL-GAME structures — pre-#85 an F5-period leg reaching a resolver would
silently resolve to the FG selection and return a confidently wrong price
(the worst failure mode in this pipeline). Until #85 plumbs period-aware
structures through all six books, any leg set containing a non-FG leg must
price as None (book declines) WITHOUT touching the book's hooks, so the
on-demand flight completes zero-book and the maker declines
``live_fetch_timeout`` / ``too_few_books`` — never quotes a wrong number.
"""
from kalshi_common.legset import CanonicalLeg
from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import GameRef, ResolvedLeg

EVT = "26AUG122210KCLAD"
GAME = GameRef(game_id="g1", home_team="Los Angeles Dodgers",
               away_team="Kansas City Royals", commence_time=None)


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


def test_f5_legs_decline_without_touching_book_hooks():
    calls = []
    svc = SGPService(books=("novig",),
                     on_demand_hooks={"novig": _tracking_hooks(calls)})
    legs = [CanonicalLeg(EVT, "spread", -1.5, "home", "F5"),
            CanonicalLeg(EVT, "total", 6.5, "over", "F5")]
    assert svc.price_on_demand("novig", GAME, legs) is None
    assert calls == []                 # no fetch, no resolve, no price


def test_mixed_fg_f5_leg_set_also_declines():
    calls = []
    svc = SGPService(books=("novig",),
                     on_demand_hooks={"novig": _tracking_hooks(calls)})
    legs = [CanonicalLeg(EVT, "spread", -1.5, "home"),          # FG
            CanonicalLeg(EVT, "total", 6.5, "over", "F5")]      # F5
    assert svc.price_on_demand("novig", GAME, legs) is None
    assert calls == []


def test_fg_legs_still_price_through_hooks():
    calls = []
    svc = SGPService(books=("novig",),
                     on_demand_hooks={"novig": _tracking_hooks(calls)})
    legs = [CanonicalLeg(EVT, "spread", -1.5, "home"),
            CanonicalLeg(EVT, "total", 8.5, "over")]
    res = svc.price_on_demand("novig", GAME, legs)
    assert res is not None
    assert len(calls) > 0
