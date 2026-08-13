"""Issues #84/#85/#86: SGPService.price_on_demand period guard.

#84 failed CLOSED on every non-FG leg because the six book resolvers were
FG-only. #85 plumbed period-aware structures through all six books. #86
opens F5-winner TEAM legs (books' push-refund F5 ML converts to Kalshi's
unconditional tie-loses space downstream) with the adapter labeling its
semantics on the result. What stays closed, declining WITHOUT touching the
book's hooks (zero wire calls — flight completes zero-book, maker declines
live_fetch_timeout / too_few_books — never a confidently wrong price):

* F5-winner TIE-side legs — no book prices a bare F5-tie leg in an SGP.
* F5 ML legs at a book whose hooks carry NO ``f5_ml_semantics`` label —
  unlabeled conditional fairs must never leak into consensus.
* Any UNKNOWN period — a future period value must never reach resolvers
  that don't know it.
"""
from kalshi_common.legset import CanonicalLeg
from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import GameRef, ResolvedLeg

EVT = "26AUG122210KCLAD"
GAME = GameRef(game_id="g1", home_team="Los Angeles Dodgers",
               away_team="Kansas City Royals", commence_time=None)


def _tracking_hooks(calls, *, f5_ml_semantics=None, tie_prob=None):
    def price(client, refs, event):
        calls.append(list(refs))
        return 2.0

    hooks = {
        "match_event": lambda client, game: calls.append("match") or "EV",
        "build_structure": lambda client, ev, game: {"fake": "structure"},
        "resolve": lambda structure, legs, home, away: [
            ResolvedLeg(ref=f"r{j}", opposite_ref=f"o{j}",
                        single_decimal=1.9, opposite_decimal=1.9)
            for j in range(len(legs))],
        "price": price,
    }
    if f5_ml_semantics is not None:
        hooks["f5_ml_semantics"] = f5_ml_semantics
    if tie_prob is not None:
        hooks["tie_prob"] = tie_prob
    return hooks


def _service(calls, **hook_kwargs):
    return SGPService(
        books=("novig",),
        on_demand_hooks={"novig": _tracking_hooks(calls, **hook_kwargs)})


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


def test_f5_winner_team_leg_prices_through_labeled_hooks():
    """#86: an F5-winner TEAM leg flows to the hooks when the adapter labels
    its F5 ML semantics; the result carries the label + the book's own tie
    anchor for the engine-side conversion."""
    calls = []
    svc = _service(calls, f5_ml_semantics="push_2way",
                   tie_prob=lambda structure: 0.14)
    legs = [CanonicalLeg(EVT, "ml", None, "home", "F5"),
            CanonicalLeg(EVT, "total", 4.5, "over", "F5")]
    res = svc.price_on_demand("novig", GAME, legs)
    assert res is not None
    assert len(calls) > 0
    assert res.f5_semantics == "push_2way"
    assert res.p_tie_book == 0.14


def test_f5_winner_without_tie_prob_hook_still_prices_unanchored():
    """Books with no 3-way market (DK/NV/CZR) label semantics but carry no
    anchor — p_tie_book None; the engine converts them with the flight's
    median anchor (or excludes them when no anchor landed at all)."""
    calls = []
    svc = _service(calls, f5_ml_semantics="push_2way")
    legs = [CanonicalLeg(EVT, "ml", None, "home", "F5"),
            CanonicalLeg(EVT, "total", 4.5, "over", "F5")]
    res = svc.price_on_demand("novig", GAME, legs)
    assert res is not None
    assert res.f5_semantics == "push_2way"
    assert res.p_tie_book is None


def test_f5_winner_without_semantics_label_declines_zero_wire():
    """A book whose hooks carry no f5_ml_semantics must never price an
    F5-winner set: an unlabeled conditional fair in consensus is exactly
    the ~P(tie) pick-off #86 exists to close."""
    calls = []
    svc = _service(calls)
    legs = [CanonicalLeg(EVT, "ml", None, "home", "F5"),
            CanonicalLeg(EVT, "total", 4.5, "over", "F5")]
    assert svc.price_on_demand("novig", GAME, legs) is None
    assert calls == []                 # no fetch, no resolve, no price


def test_f5_tie_side_leg_declines_zero_wire_even_when_labeled():
    calls = []
    svc = _service(calls, f5_ml_semantics="push_2way",
                   tie_prob=lambda structure: 0.14)
    for side in ("tie", "not_tie"):
        legs = [CanonicalLeg(EVT, "ml", None, side, "F5"),
                CanonicalLeg(EVT, "total", 4.5, "over", "F5")]
        assert svc.price_on_demand("novig", GAME, legs) is None
    assert calls == []


def test_fg_set_carries_no_f5_semantics_even_when_hooks_labeled():
    """The label rides ONLY F5-winner flights — an FG/F5 spread-total set at
    a labeled book must not be tagged (the engine would else convert a fair
    that was never conditional)."""
    calls = []
    svc = _service(calls, f5_ml_semantics="push_2way",
                   tie_prob=lambda structure: 0.14)
    legs = [CanonicalLeg(EVT, "spread", -0.5, "home", "F5"),
            CanonicalLeg(EVT, "total", 4.5, "over", "F5")]
    res = svc.price_on_demand("novig", GAME, legs)
    assert res is not None
    assert res.f5_semantics is None
    assert res.p_tie_book is None


def test_f5_winner_broken_tie_prob_hook_prices_unanchored():
    """A raising tie_prob hook must cost the anchor, not the book."""
    calls = []

    def boom(structure):
        raise ValueError("bad structure")

    svc = _service(calls, f5_ml_semantics="push_2way", tie_prob=boom)
    legs = [CanonicalLeg(EVT, "ml", None, "home", "F5"),
            CanonicalLeg(EVT, "total", 4.5, "over", "F5")]
    res = svc.price_on_demand("novig", GAME, legs)
    assert res is not None
    assert res.f5_semantics == "push_2way"
    assert res.p_tie_book is None


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
