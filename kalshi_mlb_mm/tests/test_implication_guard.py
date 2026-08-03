"""Issue #66: same-game spread+moneyline leg sets are refused end-to-end.

The classifier returning "unpriceable" is only half the claim. What matters
operationally is that (a) the RFQ falls out of scope with a NAMED reason that
reaches quote_decisions, and (b) the refusal happens BEFORE any book is
queried — a wire call made and then discarded is latency and rate-limit budget
spent for nothing, and Novig's 403 threshold is measured in call count.
"""
import pandas as pd
import pytest

from kalshi_common import legset
from kalshi_mlb_mm import router

EVT = "KXMLBGAME-25JUN271905NYYBOS"      # away=NYY, home=BOS
EVT_B = "KXMLBGAME-25JUN272005SDLAD"     # away=SD, home=LAD


def _spread(team, n, side):
    return {"market_ticker": f"KXMLBSPREAD-25JUN271905NYYBOS-{team}{n}",
            "event_ticker": EVT, "side": side}


def _total(n, side):
    return {"market_ticker": f"KXMLBTOTAL-25JUN271905NYYBOS-{n}",
            "event_ticker": EVT, "side": side}


def _ml(team, side):
    return {"market_ticker": f"KXMLBGAME-25JUN271905NYYBOS-{team}",
            "event_ticker": EVT, "side": side}


SPREAD_ML = [_spread("BOS", 2, "yes"), _ml("BOS", "yes")]
THREE_LEG = [_spread("BOS", 2, "yes"), _total(9, "yes"), _ml("BOS", "yes")]


# ---------------------------------------------------------------- #
# (a) the RFQ leaves scope with a named reason                      #
# ---------------------------------------------------------------- #

@pytest.mark.parametrize("legs", [SPREAD_ML, THREE_LEG],
                         ids=["spread_ml", "spread_total_ml"])
def test_maker_declines_and_names_the_reason(legs):
    from kalshi_mlb_mm import main
    canon = legset.parse_legs(legs)
    assert canon is not None                      # legs are well-formed MLB
    assert main._priceable(canon) is False
    assert main._priceable_in_phase1(canon) is False
    assert main._out_of_scope_reason(legs, canon) == \
        "out_of_scope_spread_ml_implication"


def test_named_reason_has_a_monitor_gloss():
    """A reason with no gloss renders blank in kalshi_mlb_monitor, which is
    how a deliberate refusal ends up looking like a silent drop."""
    from kalshi_mlb_monitor import bots
    assert bots.gloss("out_of_scope_spread_ml_implication") != ""


def test_other_out_of_scope_reasons_are_unchanged():
    from kalshi_mlb_mm import main
    assert main._out_of_scope_reason([], None) == "out_of_scope_unparseable"
    bad = [{"market_ticker": "KXMLBPLAYER-foo", "event_ticker": EVT, "side": "yes"}]
    assert main._out_of_scope_reason(bad, None) == "out_of_scope_non_mlb"
    one = legset.parse_legs([_ml("BOS", "yes")])
    assert main._out_of_scope_reason([_ml("BOS", "yes")], one) == \
        "out_of_scope_lone_single"


def test_in_scope_shapes_are_still_in_scope():
    """Non-regression: the two grid routes and a cross-game pair keep quoting."""
    from kalshi_mlb_mm import main
    grid_st = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes")])
    grid_mt = legset.parse_legs([_ml("BOS", "yes"), _total(9, "yes")])
    cross = legset.parse_legs([
        _spread("BOS", 2, "yes"),
        {"market_ticker": "KXMLBGAME-25JUN272005SDLAD-LAD",
         "event_ticker": EVT_B, "side": "yes"}])
    assert main._priceable(grid_st) is True
    assert main._priceable(grid_mt) is True
    assert main._priceable(cross) is True, \
        "a spread and an ML in DIFFERENT games are independent — still quotable"


# ---------------------------------------------------------------- #
# (b) refused BEFORE any wire call                                  #
# ---------------------------------------------------------------- #

class _Tripwire:
    """Stands in for the on-demand lookup, which is the only thing that can
    cause a book fetch. Any call at all is a failure of the guard."""

    def __init__(self):
        self.calls = []

    def __call__(self, leg_set_hash):
        self.calls.append(leg_set_hash)
        raise AssertionError(
            "on_demand lookup reached for a refused leg set — the refusal "
            "must happen before any book query")


@pytest.mark.parametrize("legs", [SPREAD_ML, THREE_LEG],
                         ids=["spread_ml", "spread_total_ml"])
def test_router_never_consults_the_on_demand_feed(legs):
    tripwire = _Tripwire()
    canon = legset.parse_legs(legs)
    cons, reason = router.subcombo_consensus(
        "g1", canon, sgp_df=None, min_books=2, sigma_z_max=0.07,
        on_demand_fairs=tripwire)
    assert cons is None and reason == "unpriceable"
    assert tripwire.calls == []


def test_combo_fair_detail_never_consults_the_feed():
    tripwire = _Tripwire()
    detail, reason = router.combo_fair_detail(
        THREE_LEG, pd.DataFrame(), lambda gl: "g1", 2, 0.07,
        on_demand_fairs=tripwire)
    assert detail is None and reason == "unpriceable"
    assert tripwire.calls == []


class _TripwireEngine:
    """Every method that could touch a book raises. The guard must mean none
    of them is reachable for a refused leg set."""

    def lookup(self, h):
        raise AssertionError("lookup reached for a refused leg set")

    def lookup_results(self, h):
        raise AssertionError("lookup_results reached for a refused leg set")

    def landed_at(self, h):
        raise AssertionError("landed_at reached for a refused leg set")

    def ensure_fetch(self, h, game, legs):
        raise AssertionError("ensure_fetch reached — this issues book traffic")

    def refetch_now(self, jobs, deadline_sec):
        raise AssertionError("refetch_now reached for a refused leg set")


class _NoQuoteGW:
    def submit_quote(self, *a):
        raise AssertionError("must not quote a refused leg set")


def test_discovery_tick_makes_no_book_call_and_logs_the_named_reason(
        monkeypatch, tmp_path):
    """The live path, with scope computed rather than pre-seeded: a refused
    RFQ must never reach _ENGINE.ensure_fetch (which issues the six books'
    SGP calls), and the named reason must land in quote_decisions."""
    import importlib
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "guard.duckdb")
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()
    # Non-empty: _discovery_tick returns at the book-freshness gate before it
    # ever evaluates scope if _SGP_ODDS is empty. Content is irrelevant here —
    # a refused leg set never reaches a grid lookup.
    monkeypatch.setattr(main, "_SGP_ODDS", pd.DataFrame(
        {"game_id": ["g"], "combo": ["Home Spread + Over"], "period": ["FG"],
         "bookmaker": ["dk"], "sgp_decimal": [2.0], "fetch_time": [None],
         "spread_line": [-1.5], "total_line": [8.5]}))
    monkeypatch.setattr(main, "_today_fills", lambda: [])
    monkeypatch.setattr(main, "_today_fills_by_game", lambda: [])
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)
    monkeypatch.setattr(main, "_ENGINE", _TripwireEngine())
    monkeypatch.setattr(main, "_SCOPE_CACHE", {})     # force a real scope decision
    monkeypatch.setattr(main, "_OD_RESULT_EMITTED", {})

    class Src:
        def poll(self):
            return [{"id": "r-guard", "market_ticker": "COMBO-3LEG",
                     "contracts": 1}]

        def get_market(self, t):
            return {"mve_selected_legs": THREE_LEG}

    main._discovery_tick(Src(), _NoQuoteGW(), dry_run=False)

    with db.connect(read_only=True) as con:
        decision = con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "ORDER BY observed_at DESC LIMIT 1").fetchone()
        seen = con.execute(
            "SELECT in_scope, last_decision FROM seen_rfqs "
            "WHERE rfq_id='r-guard'").fetchone()
    assert decision == ("skipped", "out_of_scope_spread_ml_implication")
    assert seen == (False, "out_of_scope_spread_ml_implication")
