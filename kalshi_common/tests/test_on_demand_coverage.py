"""Issue #81 — on-demand coverage counter: `scrape_done`'s observability heir.

With the maker's sweep gone there is no periodic per-book row count, so a
silently-dead book would be visible only through #37 alerts (which can be
disabled). `OnDemandCoverageCounter` is the replacement record: an
UNCONDITIONAL per-book outcome tally over the on_demand path, owned by
SGPService (never by the alerter — `build_alerter` returns None when
alerting is off, and research coverage must not die with it), drained
periodically by the maker's tick into an `on_demand_coverage` research
event.

Written FIRST, failing on pre-#81 code, per repo rule.
"""
from kalshi_common.legset import CanonicalLeg
from kalshi_common.sgp_health import OnDemandCoverageCounter
from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import GameRef, ResolvedLeg

GAME = GameRef(game_id="g1", home_team="Boston Red Sox",
               away_team="New York Yankees", commence_time=None)

# Four joint fairs summing to 1.0 with a uniform 1.2 overround.
CELL_FAIRS = [0.40, 0.20, 0.25, 0.15]
CELL_DECS = [1.0 / (1.2 * p) for p in CELL_FAIRS]


def _legs():
    evt = "KXMLBGAME-26JUL09NYYBOS"
    return [CanonicalLeg(evt, "spread", -1.5, "home"),
            CanonicalLeg(evt, "total", 8.5, "over")]


def _resolved():
    return [ResolvedLeg(ref=f"r{j}", opposite_ref=f"o{j}",
                        single_decimal=1.9, opposite_decimal=1.9)
            for j in range(2)]


def _hooks(price_fn):
    return {
        "match_event": lambda client, game: "EV",
        "build_structure": lambda client, ev, game: {"fake": "s"},
        "resolve": lambda structure, legs, home, away: _resolved(),
        "price": price_fn,
    }


# ------------------------------------------------------------------ #
# Counter semantics                                                   #
# ------------------------------------------------------------------ #

def test_counts_on_demand_path_only():
    c = OnDemandCoverageCounter()
    c.observe(book="draftkings", path="on_demand", outcome="ok")
    c.observe(book="draftkings", path="on_demand", outcome="ok")
    c.observe(book="novig", path="on_demand", outcome="transport_error")
    # sweep/warming outcomes are NOT coverage — the heir replaces the
    # sweep's book counts with the path that actually prices quotes.
    c.observe(book="draftkings", path="sweep", outcome="ok")
    c.observe(book="caesars", path="warming", outcome="ok")
    assert c.snapshot_and_reset() == {
        "draftkings": {"ok": 2},
        "novig": {"transport_error": 1},
    }


def test_snapshot_resets():
    c = OnDemandCoverageCounter()
    c.observe(book="fanduel", path="on_demand", outcome="empty")
    assert c.snapshot_and_reset() == {"fanduel": {"empty": 1}}
    assert c.snapshot_and_reset() == {}


def test_snapshot_is_a_copy_not_a_view():
    c = OnDemandCoverageCounter()
    c.observe(book="fanduel", path="on_demand", outcome="ok")
    snap = c.snapshot_and_reset()
    snap["fanduel"]["ok"] = 999
    c.observe(book="fanduel", path="on_demand", outcome="ok")
    assert c.snapshot_and_reset() == {"fanduel": {"ok": 1}}


# ------------------------------------------------------------------ #
# Service wiring                                                      #
# ------------------------------------------------------------------ #

def test_service_owns_counter_and_feeds_it_per_fetch():
    """Every price_on_demand outcome lands in svc.coverage — with NO
    book_health configured (coverage is unconditional, alerting is not)."""
    book = "draftkings"
    svc = SGPService(books=(book,),
                     on_demand_hooks={book: _hooks(
                         lambda client, refs, event: CELL_DECS[0])})
    assert svc.price_on_demand(book, GAME, _legs()) is not None
    snap = svc.coverage.snapshot_and_reset()
    assert snap == {book: {"ok": 1}}


def test_service_counts_failures_too():
    book = "novig"

    def price_none(client, refs, event):
        return None                      # book declines every cell

    svc = SGPService(books=(book,),
                     on_demand_hooks={book: _hooks(price_none)})
    assert svc.price_on_demand(book, GAME, _legs()) is None
    snap = svc.coverage.snapshot_and_reset()
    assert list(snap) == [book]
    assert sum(snap[book].values()) == 1
    assert "ok" not in snap[book]
