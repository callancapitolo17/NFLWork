"""Issues #84/#85: SGPService.price_on_demand period guard.

#84 failed CLOSED on every non-FG leg because the six book resolvers were
FG-only. #85 plumbed period-aware structures through all six books, so F5
spread/total legs now flow to the hooks. What stays closed:

* F5 MONEYLINE legs (the KXMLBF5 winner family) — books' F5 ML is
  push-refund / 3-way while Kalshi's team markets resolve NO on a tie after
  5; mapping one onto the other misprices by ~P(tie). Closed until #86
  lands the conversion math. (classify_subcombo already routes these
  combos "unpriceable"; this is the service-level backstop.)
* Any UNKNOWN period — a future period value must never reach resolvers
  that don't know it.

Both closures must decline WITHOUT touching the book's hooks: zero wire
calls, flight completes zero-book, maker declines live_fetch_timeout /
too_few_books — never a confidently wrong price.
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


def _service(calls):
    return SGPService(books=("novig",),
                      on_demand_hooks={"novig": _tracking_hooks(calls)})


def test_f5_spread_total_legs_price_through_hooks():
    calls = []
    svc = _service(calls)
    legs = [CanonicalLeg(EVT, "spread", -0.5, "home", "F5"),
            CanonicalLeg(EVT, "total", 4.5, "over", "F5")]
    res = svc.price_on_demand("novig", GAME, legs)
    assert res is not None
    assert len(calls) > 0


def test_mixed_fg_f5_leg_set_prices_through_hooks():
    calls = []
    svc = _service(calls)
    legs = [CanonicalLeg(EVT, "spread", -1.5, "home"),          # FG
            CanonicalLeg(EVT, "total", 4.5, "over", "F5")]      # F5
    res = svc.price_on_demand("novig", GAME, legs)
    assert res is not None
    assert len(calls) > 0


def test_f5_moneyline_declines_without_touching_book_hooks():
    calls = []
    svc = _service(calls)
    legs = [CanonicalLeg(EVT, "ml", None, "home", "F5"),
            CanonicalLeg(EVT, "total", 4.5, "over", "F5")]
    assert svc.price_on_demand("novig", GAME, legs) is None
    assert calls == []                 # no fetch, no resolve, no price


def test_unknown_period_declines_without_touching_book_hooks():
    calls = []
    svc = _service(calls)
    legs = [CanonicalLeg(EVT, "total", 2.5, "over", "F3"),
            CanonicalLeg(EVT, "total", 8.5, "over")]
    assert svc.price_on_demand("novig", GAME, legs) is None
    assert calls == []


def test_fg_legs_still_price_through_hooks():
    calls = []
    svc = _service(calls)
    legs = [CanonicalLeg(EVT, "spread", -1.5, "home"),
            CanonicalLeg(EVT, "total", 8.5, "over")]
    res = svc.price_on_demand("novig", GAME, legs)
    assert res is not None
    assert len(calls) > 0
