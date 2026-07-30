"""Minimal per-book drives of ``price_sgps`` for the counter/tripwire tests.

NOT a test module — imported by ``test_parse_tripwires.py`` and
``test_library_logging.py``.

Each drive runs one book's real sweep orchestration end to end with the
network faked at the client boundary, in two flavors:

  healthy=True   the book offers the target line  -> rows, no tripwire
  healthy=False  the book's structure parses to NOTHING (a mangled runner
                 name / renamed market) while events still match — the exact
                 shape of the FanDuel alt-total paren regression.

The parsers themselves are covered by their own fixture tests; these drives
exist to prove the COUNTERS and the TRIPWIRE, so per-book structure is stubbed
at the smallest seam that keeps ``price_sgps``'s own control flow real.

``counters=`` (issue #38) forwards a caller-owned ``FetchCounters`` into the
sweep so a test can inspect what the real orchestration counted — that is the
seam ``SGPService`` uses to put sweep counters on a health row. Default None
keeps every pre-#38 caller unchanged.
"""
from __future__ import annotations

import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(MLB_SGP_DIR))

from mlb_sgp._shared import TargetLine  # noqa: E402

CT = datetime(2026, 6, 8, 23, 0, tzinfo=timezone.utc)
SPREAD = -1.5
TOTAL = 8.5


def target(game_id: str = "g1") -> TargetLine:
    return TargetLine(game_id=game_id, home_team="New York Yankees",
                      away_team="Boston Red Sox", commence_time=CT,
                      period="FG", spread=SPREAD, total=TOTAL)


# ------------------------------------------------------------------ #
# Caesars                                                             #
# ------------------------------------------------------------------ #

def _czr_leg(decimal: float) -> dict:
    return {"price": {"decimal": decimal}}


def drive_caesars(monkeypatch, *, healthy=True, combo_decimal=2.5,
                  leg_decimal=1.9, price_exc=None, n_targets=1,
                  counters=None):
    from mlb_sgp import caesars

    event = SimpleNamespace(event_id="e1", home_team="New York Yankees",
                            away_team="Boston Red Sox", start_time=CT)
    monkeypatch.setattr(caesars, "_match_events",
                        lambda events, targets: {"g1": event})
    if healthy:
        parsed = {
            "FG": {"spreads": {SPREAD: {"home": _czr_leg(leg_decimal),
                                        "away": _czr_leg(leg_decimal)}},
                   "totals": {TOTAL: {"over": _czr_leg(leg_decimal),
                                      "under": _czr_leg(leg_decimal)}},
                   "moneyline": None},
            "F5": {"spreads": {}, "totals": {}, "moneyline": None},
        }
    else:
        # Structure fetched fine; our parser recognized no line in it.
        parsed = {"FG": {"spreads": {}, "totals": {}, "moneyline": None},
                  "F5": {"spreads": {}, "totals": {}, "moneyline": None}}
    monkeypatch.setattr(caesars, "parse_markets", lambda ev: parsed)

    client = MagicMock()
    client.list_events.return_value = [event]
    client.fetch_event.return_value = {"id": "e1"}
    if price_exc is not None:
        client.price_combo.side_effect = price_exc
    elif combo_decimal is None:
        # How CZR/MGM report BOTH "book declines this combo" and a swallowed
        # price-stage transport failure: the client returns None, so the
        # orchestrator never sees an exception.
        client.price_combo.return_value = None
    else:
        client.price_combo.return_value = {"decimal": combo_decimal,
                                           "american": 150}
    # n_targets>1 walks the same game at several totals, which is how a real
    # slate reaches PriceCallTally's 8-price-call verdict floor. TargetLine is
    # a frozen dataclass (NOT a namedtuple), so replace() is the way to vary
    # it — an earlier cut used _replace and silently produced N identical
    # targets instead.
    targets = [replace(target(), total=TOTAL + i) for i in range(n_targets)]
    if n_targets > 1 and healthy:
        for i in range(n_targets):
            parsed["FG"]["totals"][TOTAL + i] = {
                "over": _czr_leg(leg_decimal), "under": _czr_leg(leg_decimal)}
    return caesars.price_sgps(targets, periods=("FG",), client=client,
                              counters=counters)


# ------------------------------------------------------------------ #
# BetMGM                                                              #
# ------------------------------------------------------------------ #

def drive_betmgm(monkeypatch, *, healthy=True, combo_decimal=2.5,
                 leg_decimal=1.9, counters=None):
    from mlb_sgp import betmgm

    event = SimpleNamespace(event_id="e1", home_team="New York Yankees",
                            away_team="Boston Red Sox", start_time=CT)
    monkeypatch.setattr(betmgm, "_match_events",
                        lambda events, targets: {"g1": event})
    # MGM legs are (market_id, option_id, decimal) tuples.
    leg = ("m1", "o1", leg_decimal)
    if healthy:
        parsed = {"FG": {"spreads": {SPREAD: {"home": leg, "away": leg}},
                         "totals": {TOTAL: {"over": leg, "under": leg}},
                         "moneyline": None}}
    else:
        parsed = {"FG": {"spreads": {}, "totals": {}, "moneyline": None}}
    monkeypatch.setattr(betmgm, "parse_markets",
                        lambda markets, home, away: parsed)

    client = MagicMock()
    client.list_events.return_value = [event]
    client.fetch_markets.return_value = [{"raw": True}]
    client.price_picks.return_value = {"decimal": combo_decimal}
    return betmgm.price_sgps([target()], periods=("FG",), client=client,
                             counters=counters)


# ------------------------------------------------------------------ #
# FanDuel                                                             #
# ------------------------------------------------------------------ #

def drive_fanduel(monkeypatch, *, healthy=True, price_exc=None,
                  counters=None):
    import scraper_fanduel_sgp as legacy
    from mlb_sgp import fanduel

    if healthy:
        spreads = {("home", SPREAD): ("SM0", "SH0"),
                   ("away", -SPREAD): ("SM0", "SA0")}
        totals = {("O", TOTAL): ("TM0", "TO0"), ("U", TOTAL): ("TM0", "TU0")}
    else:
        # FD renamed its alt-total runners; every line failed the parse, so
        # the runners dict comes back EMPTY while the event still matched.
        spreads, totals = {}, {}
    runners = {"fg": {"spreads": spreads, "totals": totals},
               "f5": {"spreads": {}, "totals": {}}}

    monkeypatch.setattr(legacy, "match_events", lambda evs, tgt: [
        {"game_id": "g1", "fd_event_id": "e1",
         "fd_home": "New York Yankees", "fd_away": "Boston Red Sox"}])
    def _price(session, smid, ssid, tmid, tsid, verbose=False):
        if price_exc is not None:
            raise price_exc
        return {"decimal": 3.4, "american": 240}

    monkeypatch.setattr(legacy, "price_combo", _price)
    fetchers = {
        "fetch_fd_events": lambda session: [{"id": "e1"}],
        "fetch_event_runners": lambda session, eid, h, a: runners,
    }
    return fanduel.price_sgps([target()], periods=("FG",), client=MagicMock(),
                              parallelism=2, fetchers=fetchers,
                              counters=counters)


# ------------------------------------------------------------------ #
# DraftKings                                                          #
# ------------------------------------------------------------------ #

def drive_draftkings(monkeypatch, *, healthy=True, counters=None):
    import scraper_draftkings_sgp as legacy
    from mlb_sgp import draftkings

    # DK selection ids carry their market number in the prefix
    # (_MARKET_NUM_RE: "0HC<num>" / "0OU<num>"), and _price_combo drops any
    # leg whose market number is not in `canonical` — so the fake ids have to
    # look real or nothing prices.
    if healthy:
        # DK keys spreads (sign, abs_line, participant); "N" = negative side.
        spreads = {("N", abs(SPREAD), "1"): ["0HC1#home"],
                   ("P", abs(SPREAD), "3"): ["0HC1#away"]}
        totals = {("O", TOTAL): ["0OU2#over"], ("U", TOTAL): ["0OU2#under"]}
    else:
        spreads, totals = {}, {}
    sel_ids = {"fg": {"spreads": spreads, "totals": totals,
                      "canonical": {"1", "2"}, "moneyline": {}},
               "f5": {"spreads": {}, "totals": {}, "canonical": set(),
                      "moneyline": {}}}

    monkeypatch.setattr(legacy, "match_events", lambda evs, tgt: [
        {"game_id": "g1", "dk_event_id": "e1",
         "dk_home": "New York Yankees", "dk_away": "Boston Red Sox"}])
    monkeypatch.setattr(legacy, "calculate_sgp",
                        lambda *a, **kw: {"trueOdds": 3.4})
    fetchers = {
        "fetch_dk_events": lambda session: [{"id": "e1"}],
        "fetch_main_market_nums": lambda session, eid: {"fg": "1"},
        "fetch_selection_ids": lambda session, eid, nums, verbose=False: sel_ids,
    }
    return draftkings.price_sgps([target()], periods=("FG",),
                                 client=MagicMock(), parallelism=2,
                                 fetchers=fetchers, counters=counters)


# ------------------------------------------------------------------ #
# Novig                                                               #
# ------------------------------------------------------------------ #

def _nv_market(mtype, strike, offered=True):
    return {"type": mtype, "strike": strike, "is_consensus": True,
            "outcomes": [] if not offered else [{"id": "o1"}]}


def drive_novig(monkeypatch, *, healthy=True, counters=None):
    import scraper_novig_sgp as legacy
    from mlb_sgp import novig

    monkeypatch.setattr(legacy, "match_events", lambda evs, lines: [
        {"game_id": "g1", "nv_event_id": "e1", "nv_home": "New York Yankees",
         "nv_away": "Boston Red Sox", "nv_home_sym": "NYY",
         "nv_away_sym": "BOS", "fg_total_line": TOTAL}])
    leg = {"id": "u1", "available": 0.52}
    legs = {"fg": {"home_spread": leg, "away_spread": leg,
                   "over": leg, "under": leg,
                   "home_ml": leg, "away_ml": leg},
            "f5": {}}
    markets = [_nv_market("SPREAD", SPREAD), _nv_market("TOTAL", TOTAL)]
    monkeypatch.setattr(legacy, "fetch_event_legs",
                        lambda session, game, verbose=False, **kw: (legs, markets))
    monkeypatch.setattr(legacy, "submit_parlay",
                        lambda session, ids, verbose=False:
                        ({"decimal": 2.5, "american": 150}, False))
    offered = {"spreads": {SPREAD}, "totals": {TOTAL}} if healthy else \
              {"spreads": set(), "totals": set()}
    monkeypatch.setattr(novig, "_extract_offered_lines_nv",
                        lambda markets, period_key:
                        offered if period_key == "fg"
                        else {"spreads": set(), "totals": set()})

    client = MagicMock()
    client.list_events.return_value = [
        SimpleNamespace(event_id="e1", home_team="New York Yankees",
                        away_team="Boston Red Sox", home_sym="NYY",
                        away_sym="BOS", start_time=CT)]
    return novig.price_sgps([target()], periods=("FG",), client=client,
                            counters=counters)


# ------------------------------------------------------------------ #
# ProphetX                                                            #
# ------------------------------------------------------------------ #

def drive_prophetx(monkeypatch, *, healthy=True, counters=None):
    import scraper_prophetx_sgp as legacy
    from mlb_sgp import prophetx

    monkeypatch.setattr(legacy, "match_events", lambda evs, lines: [
        {"game_id": "g1", "px_event_id": "e1", "px_home": "New York Yankees",
         "px_away": "Boston Red Sox", "px_home_competitor_id": 1,
         "px_away_competitor_id": 3}])
    monkeypatch.setattr(legacy, "_verify_competitor_ids",
                        lambda markets, home_id, away_id: True)

    def _sel(sel_id, line, competitor_id=None, name=None):
        return {"id": sel_id, "line": line, "competitorId": competitor_id,
                "name": name, "odds": -110}

    px_markets = [
        {"market_id": "m_spread", "name": "Run Line",
         "selections": [[_sel("s_home", SPREAD, competitor_id=1),
                         _sel("s_away", -SPREAD, competitor_id=3)]],
         "outcomes": [], "ml_selections": []},
        {"market_id": "m_total", "name": "Total Runs",
         "selections": [[_sel("t_over", TOTAL, name="Over"),
                         _sel("t_under", TOTAL, name="Under")]],
         "outcomes": [], "ml_selections": []},
    ]
    market_objs = [SimpleNamespace(market_id=m["market_id"], name=m["name"],
                                   selections=[{"selections": m["selections"]}],
                                   outcomes=m["outcomes"],
                                   ml_selections=m["ml_selections"])
                   for m in px_markets]
    offered = {"spreads": {SPREAD}, "totals": {TOTAL}} if healthy else \
              {"spreads": set(), "totals": set()}
    monkeypatch.setattr(prophetx, "_extract_offered_lines_px",
                        lambda markets, period_key, home_id:
                        offered if period_key == "fg"
                        else {"spreads": set(), "totals": set()})

    client = MagicMock()
    client.list_events.return_value = [
        SimpleNamespace(event_id="e1", home_team="New York Yankees",
                        away_team="Boston Red Sox", home_id=1, away_id=3,
                        start_time=CT)]
    client.fetch_event_markets.return_value = market_objs
    client.submit_parlay_rfq.return_value = ({"odds": 150}, False)
    return prophetx.price_sgps([target()], periods=("FG",), client=client,
                               parallelism=2, counters=counters)


DRIVES = {
    "caesars": drive_caesars,
    "betmgm": drive_betmgm,
    "fanduel": drive_fanduel,
    "draftkings": drive_draftkings,
    "novig": drive_novig,
    "prophetx": drive_prophetx,
}
