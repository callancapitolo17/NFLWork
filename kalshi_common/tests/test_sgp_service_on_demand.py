"""Tests for SGPService.price_on_demand (Phase 2, Task 11).

All tests inject fake per-book hooks via the `on_demand_hooks` test seam
(mirroring the existing `runners` seam) — no HTTP, no clients. Hook shape
per book:

    {"match_event":     f(client, game) -> event | None,
     "build_structure": f(client, event, game) -> structure | None,
     "resolve":         f(structure, legs, home, away) -> [ResolvedLeg] | None,
     "price":           f(client, refs, event) -> float | None}

Cell ordering contract under test (shared with legset.enumerate_partition):
cell index i in range(2^N); bit j of i (LSB = leg 0) means leg j is FLIPPED
to its opposite ref; cell 0 = the target and is priced FIRST.
"""
import pytest

from kalshi_common import fair_value
from kalshi_common.legset import CanonicalLeg
from kalshi_common.sgp_service import SGPService, ON_DEMAND_VIG_FALLBACK
from mlb_sgp._shared import GameRef, OnDemandBookResult, ResolvedLeg

EVT = "KXMLBGAME-26JUL09NYYBOS"
GAME = GameRef(game_id="g1", home_team="Boston Red Sox",
               away_team="New York Yankees", commence_time=None)


def _legs(n=3):
    all_legs = [
        CanonicalLeg(EVT, "spread", -1.5, "home"),
        CanonicalLeg(EVT, "total", 8.5, "over"),
        CanonicalLeg(EVT, "ml", None, "home"),
        CanonicalLeg(EVT, "total", 9.5, "under"),
    ]
    return all_legs[:n]


def _resolved(n=3, two_sided=True, singles=True):
    """n ResolvedLegs with refs r0..rN / o0..oN and optional single odds."""
    out = []
    for j in range(n):
        out.append(ResolvedLeg(
            ref=f"r{j}",
            opposite_ref=(f"o{j}" if two_sided else None),
            single_decimal=(1.9 if singles else None),
            opposite_decimal=(1.9 if (singles and two_sided) else None),
        ))
    return out


def _hooks(resolved, price_fn, event="EV"):
    return {
        "match_event": lambda client, game: event,
        "build_structure": lambda client, ev, game: {"fake": "structure"},
        "resolve": lambda structure, legs, home, away: resolved,
        "price": price_fn,
    }


def _svc(book, hooks):
    return SGPService(books=(book,), on_demand_hooks={book: hooks})


# --------------------------------------------------------------------- #
# Route A (partition)                                                    #
# --------------------------------------------------------------------- #

# 8 joint fair probs summing to 1.0, uniform 1.30 overround (< 1 + 0.12*3).
CELL_FAIRS = [0.30, 0.10, 0.10, 0.05, 0.15, 0.10, 0.10, 0.10]
CELL_DECS = [1.0 / (1.30 * p) for p in CELL_FAIRS]


def test_route_a_three_legs_prices_eight_cells_in_bitmask_order():
    calls = []

    def price(client, refs, event):
        calls.append(list(refs))
        return CELL_DECS[len(calls) - 1]

    svc = _svc("novig", _hooks(_resolved(3), price))
    res = svc.price_on_demand("novig", GAME, _legs(3))

    assert isinstance(res, OnDemandBookResult)
    assert res.book == "novig"
    assert res.route == "partition"
    assert res.n_cells_priced == 8
    assert res.latency_sec >= 0.0
    assert res.fair == pytest.approx(fair_value.devig_partition(CELL_DECS, 3))
    # Bitmask ordering: cell i flips leg j to o{j} iff bit j of i is set;
    # cell 0 (all chosen refs) priced FIRST.
    assert len(calls) == 8
    for i, refs in enumerate(calls):
        assert refs == [f"o{j}" if (i >> j) & 1 else f"r{j}" for j in range(3)]


def test_route_a_later_cell_missing_falls_to_route_b_without_repricing_target():
    calls = []

    def price(client, refs, event):
        calls.append(list(refs))
        if len(calls) == 1:
            return 3.0          # cell 0 = target
        return None             # cell 1 missing -> abandon partition

    svc = _svc("novig", _hooks(_resolved(2), price))
    res = svc.price_on_demand("novig", GAME, _legs(2))

    assert res is not None
    assert res.route == "transfer"
    assert res.n_cells_priced == 1              # only the target cell priced
    # Abandoned after cell 1; target NOT re-priced (cell 0 reused).
    assert calls == [["r0", "r1"], ["o0", "r1"]]
    fair = fair_value.devig_two_way(1.9, 1.9)[0]
    expected = fair_value.fair_by_correlation_transfer(
        3.0, [(1.0 / 1.9, fair)] * 2)
    assert res.fair == pytest.approx(expected)


def test_route_a_cell0_missing_reprices_target_for_route_b():
    calls = []

    def price(client, refs, event):
        calls.append(list(refs))
        if len(calls) == 1:
            return None          # cell 0 itself failed
        return 3.0               # Route B target re-price

    svc = _svc("novig", _hooks(_resolved(2), price))
    res = svc.price_on_demand("novig", GAME, _legs(2))

    assert res is not None
    assert res.route == "transfer"
    assert res.n_cells_priced == 1
    # Same refs both times: cell 0 attempt, then the Route B target call.
    assert calls == [["r0", "r1"], ["r0", "r1"]]


def test_partition_gate_failure_falls_to_route_b_reusing_cell0():
    # All 4 cells price, but the overround is insane -> devig_partition None
    # -> Route B without any extra wire call (target = cell 0).
    insane = [3.0, 1.05, 1.05, 1.05]
    calls = []

    def price(client, refs, event):
        calls.append(list(refs))
        return insane[len(calls) - 1]

    svc = _svc("novig", _hooks(_resolved(2), price))
    res = svc.price_on_demand("novig", GAME, _legs(2))

    assert res is not None
    assert res.route == "transfer"
    assert res.n_cells_priced == 4
    assert len(calls) == 4                       # no 5th (re-price) call


# --------------------------------------------------------------------- #
# Route B (correlation transfer)                                         #
# --------------------------------------------------------------------- #

def test_one_sided_leg_skips_route_a_and_haircuts_missing_opposite():
    calls = []

    def price(client, refs, event):
        calls.append(list(refs))
        return 3.0

    resolved = [
        ResolvedLeg(ref="r0", opposite_ref="o0",
                    single_decimal=1.9, opposite_decimal=1.9),
        ResolvedLeg(ref="r1", opposite_ref=None,     # one-sided at the book
                    single_decimal=1.9, opposite_decimal=None),
    ]
    svc = _svc("novig", _hooks(resolved, price))
    res = svc.price_on_demand("novig", GAME, _legs(2))

    assert res is not None
    assert res.route == "transfer"
    assert calls == [["r0", "r1"]]               # target only, no cells
    fair0 = fair_value.devig_two_way(1.9, 1.9)[0]
    vig = ON_DEMAND_VIG_FALLBACK["novig"]
    fair1 = (1.0 / 1.9) / (1.0 + vig)            # haircut fallback
    expected = fair_value.fair_by_correlation_transfer(
        3.0, [(1.0 / 1.9, fair0), (1.0 / 1.9, fair1)])
    assert res.fair == pytest.approx(expected)


def test_dk_style_no_singles_prices_each_side_as_one_leg_call():
    # 4 legs -> Route A ineligible (N > 3); no structure odds (DK) -> each
    # leg's two sides priced via 1-leg calls, then devig_two_way.
    calls = []

    def price(client, refs, event):
        calls.append(list(refs))
        if len(refs) == 4:
            return 3.5                          # the target SGP
        (ref,) = refs
        return 1.35 if ref.startswith("r") else 3.4

    svc = _svc("draftkings", _hooks(_resolved(4, singles=False), price))
    res = svc.price_on_demand("draftkings", GAME, _legs(4))

    assert res is not None
    assert res.route == "transfer"
    assert res.n_cells_priced == 1               # target only
    # #50 item 4: the 1-leg calls fire concurrently, so only the target's
    # position is ordered; the 8 single calls are an unordered multiset.
    assert calls[0] == ["r0", "r1", "r2", "r3"]
    expected_singles = sorted(
        [[f"r{j}"] for j in range(4)] + [[f"o{j}"] for j in range(4)])
    assert sorted(calls[1:]) == expected_singles
    fair = fair_value.devig_two_way(1.35, 3.4)[0]
    expected = fair_value.fair_by_correlation_transfer(
        3.5, [(1.0 / 1.35, fair)] * 4)
    assert res.fair == pytest.approx(expected)


def test_dk_style_one_sided_leg_haircuts_after_one_leg_call():
    calls = []

    def price(client, refs, event):
        calls.append(list(refs))
        if len(refs) == 2:
            return 3.0
        (ref,) = refs
        return {"r0": 1.9, "o0": 2.0, "r1": 1.9}[ref]

    resolved = [
        ResolvedLeg(ref="r0", opposite_ref="o0",
                    single_decimal=None, opposite_decimal=None),
        ResolvedLeg(ref="r1", opposite_ref=None,     # one-sided AND no odds
                    single_decimal=None, opposite_decimal=None),
    ]
    svc = _svc("draftkings", _hooks(resolved, price))
    res = svc.price_on_demand("draftkings", GAME, _legs(2))

    assert res is not None
    assert res.route == "transfer"
    # target first, then r0/o0/r1 concurrently (#50) — never the missing o1.
    assert calls[0] == ["r0", "r1"]
    assert sorted(calls[1:]) == [["o0"], ["r0"], ["r1"]]
    fair0 = fair_value.devig_two_way(1.9, 2.0)[0]
    vig = ON_DEMAND_VIG_FALLBACK["draftkings"]
    fair1 = (1.0 / 1.9) / (1.0 + vig)
    expected = fair_value.fair_by_correlation_transfer(
        3.0, [(1.0 / 1.9, fair0), (1.0 / 1.9, fair1)])
    assert res.fair == pytest.approx(expected)


def test_leg_with_no_usable_single_returns_none():
    def price(client, refs, event):
        if len(refs) == 2:
            return 3.0
        return None                              # 1-leg calls all fail

    resolved = [
        ResolvedLeg(ref="r0", opposite_ref=None,
                    single_decimal=1.9, opposite_decimal=None),
        ResolvedLeg(ref="r1", opposite_ref=None,     # no odds anywhere
                    single_decimal=None, opposite_decimal=None),
    ]
    svc = _svc("draftkings", _hooks(resolved, price))
    assert svc.price_on_demand("draftkings", GAME, _legs(2)) is None


# --------------------------------------------------------------------- #
# Failure accounting: transport vs clean-None                            #
# --------------------------------------------------------------------- #

def test_transport_exception_increments_failures_and_returns_none():
    def price(client, refs, event):
        raise RuntimeError("connection reset")

    svc = _svc("novig", _hooks(_resolved(2), price))
    assert svc.price_on_demand("novig", GAME, _legs(2)) is None
    assert svc._state["novig"].failures == 1


def test_match_event_exception_is_a_transport_failure():
    hooks = _hooks(_resolved(2), lambda c, r, e: 3.0)
    def boom(client, game):
        raise RuntimeError("dns failure")
    hooks["match_event"] = boom
    svc = _svc("novig", hooks)
    assert svc.price_on_demand("novig", GAME, _legs(2)) is None
    assert svc._state["novig"].failures == 1


def test_clean_resolve_miss_is_not_a_failure():
    hooks = _hooks(None, lambda c, r, e: 3.0)    # resolve -> None
    svc = _svc("novig", hooks)
    assert svc.price_on_demand("novig", GAME, _legs(2)) is None
    assert svc._state["novig"].failures == 0
    assert svc._state["novig"].last_success is None


def test_event_match_miss_is_clean_none():
    hooks = _hooks(_resolved(2), lambda c, r, e: 3.0, event=None)
    svc = _svc("novig", hooks)
    assert svc.price_on_demand("novig", GAME, _legs(2)) is None
    assert svc._state["novig"].failures == 0


def test_unpriceable_target_is_clean_none():
    # Book returns None for every combo price: clean "won't price this".
    svc = _svc("novig", _hooks(_resolved(2), lambda c, r, e: None))
    assert svc.price_on_demand("novig", GAME, _legs(2)) is None
    assert svc._state["novig"].failures == 0


def test_success_touches_neither_last_success_nor_failures():
    svc = _svc("novig", _hooks(
        _resolved(2), lambda c, r, e: CELL_DECS[0]))
    svc._state["novig"].failures = 2             # pre-existing sweep strikes
    res = svc.price_on_demand("novig", GAME, _legs(2))
    assert res is not None
    assert svc._state["novig"].last_success is None   # belongs to the sweep
    assert svc._state["novig"].failures == 2          # untouched on success


def test_unknown_book_returns_none():
    svc = _svc("novig", _hooks(_resolved(2), lambda c, r, e: 3.0))
    assert svc.price_on_demand("pinnacle", GAME, _legs(2)) is None


def test_empty_legs_returns_none():
    svc = _svc("novig", _hooks(_resolved(2), lambda c, r, e: 3.0))
    assert svc.price_on_demand("novig", GAME, []) is None


# --------------------------------------------------------------------- #
# Plumbing: event threading + MGM fixture_id                             #
# --------------------------------------------------------------------- #

def test_matched_event_is_threaded_to_price_hook():
    sentinel = object()
    seen = []

    def price(client, refs, event):
        seen.append(event)
        return CELL_DECS[len(seen) - 1]

    svc = _svc("novig", _hooks(_resolved(2), price, event=sentinel))
    res = svc.price_on_demand("novig", GAME, _legs(2))
    assert res is not None
    assert seen and all(e is sentinel for e in seen)


def test_mgm_real_price_hook_threads_fixture_id():
    """The real betmgm hook must call price_selection_set(client, refs,
    fixture_id) with the matched Event's event_id (3-arg MGM contract)."""
    from types import SimpleNamespace

    class FakeMGMClient:
        def __init__(self):
            self.calls = []

        def price_picks(self, fixture_id, picks):
            self.calls.append((fixture_id, picks))
            return {"decimal": 3.2}

    svc = SGPService(books=("betmgm",))
    st = svc._state["betmgm"]
    st.client = FakeMGMClient()
    hooks = svc._book_on_demand_hooks("betmgm")
    ev = SimpleNamespace(event_id="F1", home_team="Boston Red Sox",
                         away_team="New York Yankees")
    dec = hooks["price"](st.client, [("m1", "o1"), ("m2", "o2")], ev)
    assert dec == 3.2
    assert st.client.calls == [("F1", [("m1", "o1"), ("m2", "o2")])]


# --------------------------------------------------------------------- #
# Issue #50 item 4: DK's Route-B 1-leg calls fire concurrently            #
# --------------------------------------------------------------------- #

def test_route_b_one_leg_calls_fire_concurrently_not_serially():
    """DK's structure carries no single-leg odds, so Route B pays 2 extra
    wire calls per leg. Serially that is 8 × RTT on a 4-leg combo — the
    straggler behavior #50 item 4 removes. At least two 1-leg calls must
    be in flight at the same moment."""
    import threading
    import time as _time

    state = {"active": 0, "max_active": 0}
    lock = threading.Lock()

    def price(client, refs, event):
        if len(refs) == 4:
            return 3.5                          # the target SGP
        with lock:
            state["active"] += 1
            state["max_active"] = max(state["max_active"], state["active"])
        _time.sleep(0.15)                       # a realistic 1-leg RTT
        with lock:
            state["active"] -= 1
        (ref,) = refs
        return 1.35 if ref.startswith("r") else 3.4

    svc = _svc("draftkings", _hooks(_resolved(4, singles=False), price))
    t0 = _time.monotonic()
    res = svc.price_on_demand("draftkings", GAME, _legs(4))
    wall = _time.monotonic() - t0

    assert res is not None
    assert res.route == "transfer"
    # 8 one-leg calls at 0.15s each: serial ≥ 1.2s, parallel ≪ that.
    assert state["max_active"] >= 2, "1-leg calls ran strictly serially"
    assert wall < 1.0, f"Route B singles took {wall:.2f}s — serial wire calls"
    # the math is unchanged by the parallelism
    fair = fair_value.devig_two_way(1.35, 3.4)[0]
    expected = fair_value.fair_by_correlation_transfer(
        3.5, [(1.0 / 1.35, fair)] * 4)
    assert res.fair == pytest.approx(expected)


def test_service_init_makes_canonical_match_resolvable():
    """betmgm/caesars/novig matchers import canonical_match (lives in
    "Answer Keys"); the sweep only resolves it via legacy-scraper import
    side effects, so the service must insert the path deterministically."""
    import sys
    from kalshi_common import sgp_service as mod
    SGPService = mod.SGPService
    SGPService()
    assert any(p.endswith("Answer Keys") for p in sys.path)
    import importlib
    assert importlib.util.find_spec("canonical_match") is not None
