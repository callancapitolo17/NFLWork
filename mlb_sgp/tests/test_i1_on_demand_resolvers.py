"""Issue #87: I1 (1st-inning 0.5 total) legs through every book's plumbing.

RFI legs re-encode at parse as ``CanonicalLeg(gid, "total", 0.5,
"over"/"under", "I1")``, so each book needs an "I1" bucket in its on-demand
structure. This file tests the STRUCTURE BUILDERS (the fail-first surface —
resolvers already pick buckets generically by ``leg.period``):

- DK ``fetch_selection_ids``: an "i1" bucket keyed off the parlays payload's
  ``"Runs - 1st Inning"`` market name (live 2026-08-13: mnum 85860052 with
  O/U selections at 0.5 and 1.5; the line-range heuristic CANNOT classify it
  — its <=2.0 totals were deliberately skipped as inning noise for FG/F5).
- FD ``fetch_event_runners``: ``"1st Inning 0.5 Runs"`` market — runners are
  bare "Over"/"Under" with handicap **0** (live 2026-08-13), so the line is
  FIXED at 0.5 from the market name, not parsed from the runner.
- MGM ``parse_markets``: ``"Will there be a run in the 1st inning?"``
  (Yes/No, live 2026-08-13) maps Yes->over / No->under at 0.5. MGM's
  "1st Inning Runs" market is 3-way (0/1/2+) and must NOT be read.
- Novig ``build_line_structure(period="i1")``: type ``FIRST_INNING_TOTAL``
  (live 2026-08-13, strike 0.5, "Over 0.5"/"Under 0.5" outcomes).
- PX ``_OD_MARKET_NAMES``: ``("total","I1") -> "1st Inning Total Runs"``
  (listed live 2026-08-13; board carried no priced lines, so resolution
  declines cleanly until PX carries liquidity — same posture as F5).
- CZR: deliberately NO "I1" bucket (name recon blocked by a WAF auth outage
  2026-08-13) — resolve declines the whole book cleanly.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from kalshi_common.legset import CanonicalLeg  # noqa: E402

HOME, AWAY = "Los Angeles Dodgers", "Milwaukee Brewers"

I1_OVER = CanonicalLeg("g1", "total", 0.5, "over", "I1")
I1_UNDER = CanonicalLeg("g1", "total", 0.5, "under", "I1")
FG_ML = CanonicalLeg("g1", "ml", None, "home")


# ---------------------------------------------------------------------------
# DraftKings — fetch_selection_ids builds the i1 bucket
# ---------------------------------------------------------------------------

class _FakeResp:
    def __init__(self, text, payload):
        self.status_code = 200
        self.text = text
        self._payload = payload

    def json(self):
        return self._payload


class _FakeSession:
    def __init__(self, resp):
        self._resp = resp

    def get(self, url, **kwargs):
        return self._resp


def _dk_parlays_stub():
    # Live shapes 2026-08-13 (PHI@MIN): I1 total mnum 85860052 with O/U at
    # 0.5 AND 1.5; an FG total mnum rides along to prove period separation.
    text = ("0OU85860052O50_1 0OU85860052U50_3 "
            "0OU85860052O150_1 0OU85860052U150_3 "
            "0OU85850300O850_1 0OU85850300U850_3")
    payload = {"data": {"markets": [
        {"id": "3_85860052", "name": "Runs - 1st Inning"},
        {"id": "3_85850300", "name": "Total"},
    ]}}
    return _FakeResp(text, payload)


def test_dk_fetch_selection_ids_builds_i1_totals_bucket():
    from scraper_draftkings_sgp import fetch_selection_ids
    out = fetch_selection_ids(_FakeSession(_dk_parlays_stub()), "evt1")
    i1 = out.get("i1")
    assert i1 is not None, "fetch_selection_ids must emit an 'i1' bucket"
    assert i1["totals"].get(("O", 0.5)) == ["0OU85860052O50_1"]
    assert i1["totals"].get(("U", 0.5)) == ["0OU85860052U50_3"]
    # The 1.5 alt rides along in the same bucket (harmless, Kalshi never
    # asks for it today).
    assert i1["totals"].get(("O", 1.5)) == ["0OU85860052O150_1"]


def test_dk_i1_market_stays_out_of_fg_and_f5_buckets():
    from scraper_draftkings_sgp import fetch_selection_ids
    out = fetch_selection_ids(_FakeSession(_dk_parlays_stub()), "evt1")
    for per in ("fg", "f5"):
        for (ou, line), sels in out[per]["totals"].items():
            assert line > 2.0 or not any("85860052" in s for s in sels), \
                f"I1 selections leaked into {per}: {sels}"


def test_dk_resolve_i1_leg_from_i1_bucket():
    from mlb_sgp.draftkings import resolve_legs
    structure = {
        "FG": {"spreads": {}, "totals": {("O", 8.5): ["0OU300O850_1"],
                                         ("U", 8.5): ["0OU300U850_1"]},
               "moneyline": {"1": ["0ML400_1"], "3": ["0ML400_3"]},
               "canonical": {"300", "400"}},
        "F5": None,
        "I1": {"spreads": {}, "totals": {("O", 0.5): ["0OU100O50_1"],
                                         ("U", 0.5): ["0OU100U50_3"]},
               "moneyline": {"1": [], "3": []}, "canonical": {"100"}},
    }
    out = resolve_legs(structure, [I1_OVER, FG_ML], HOME, AWAY)
    assert out is not None
    assert out[0].ref == "0OU100O50_1"
    assert out[0].opposite_ref == "0OU100U50_3"
    assert out[1].ref == "0ML400_1"


# ---------------------------------------------------------------------------
# FanDuel — fixed-line 0.5 market, runners are bare Over/Under (handicap 0)
# ---------------------------------------------------------------------------

def _fd_event_stub():
    return {"attachments": {"markets": {
        "m-i1": {"marketId": "m-i1", "marketName": "1st Inning 0.5 Runs",
                 "runners": [
                     {"selectionId": 11, "runnerName": "Over", "handicap": 0,
                      "winRunnerOdds": {"trueOdds": {"decimalOdds":
                                        {"decimalOdds": 1.75}}}},
                     {"selectionId": 12, "runnerName": "Under", "handicap": 0,
                      "winRunnerOdds": {"trueOdds": {"decimalOdds":
                                        {"decimalOdds": 2.1}}}},
                 ]},
        "m-tot": {"marketId": "m-tot", "marketName": "Total Runs",
                  "runners": [
                      {"selectionId": 21, "runnerName": "Over", "handicap": 8.5,
                       "winRunnerOdds": {"trueOdds": {"decimalOdds":
                                         {"decimalOdds": 1.9}}}},
                      {"selectionId": 22, "runnerName": "Under", "handicap": 8.5,
                       "winRunnerOdds": {"trueOdds": {"decimalOdds":
                                         {"decimalOdds": 1.9}}}},
                  ]},
    }}}


def test_fd_fetch_event_runners_builds_i1_bucket():
    from scraper_fanduel_sgp import fetch_event_runners
    out = fetch_event_runners(_FakeSession(_FakeResp("", _fd_event_stub())),
                              "fd-evt", HOME, AWAY)
    i1 = out.get("i1")
    assert i1 is not None, "fetch_event_runners must emit an 'i1' bucket"
    assert i1["totals"].get(("O", 0.5)) == ("m-i1", 11, 1.75)
    assert i1["totals"].get(("U", 0.5)) == ("m-i1", 12, 2.1)
    # FG untouched by the fixed-line path.
    assert out["fg"]["totals"].get(("O", 8.5)) == ("m-tot", 21, 1.9)
    # handicap 0 must NOT register a line-0.0 total anywhere.
    assert ("O", 0.0) not in out["fg"]["totals"]
    assert ("O", 0.0) not in i1["totals"]


def test_fd_resolve_i1_leg():
    from mlb_sgp.fanduel import resolve_legs
    structure = {
        "FG": {"spreads": {}, "totals": {}, "moneyline":
               {"home": ("m9", 91, 2.2), "away": ("m9", 92, 1.7)}},
        "F5": {"spreads": {}, "totals": {}, "moneyline": {}},
        "I1": {"spreads": {}, "totals": {("O", 0.5): ("m-i1", 11, 1.75),
                                         ("U", 0.5): ("m-i1", 12, 2.1)},
               "moneyline": {}},
    }
    out = resolve_legs(structure, [I1_OVER, FG_ML], HOME, AWAY)
    assert out is not None
    assert out[0].ref == ("m-i1", 11)
    assert out[0].opposite_ref == ("m-i1", 12)
    assert out[0].single_decimal == 1.75


# ---------------------------------------------------------------------------
# BetMGM — YRFI Yes/No market maps to over/under 0.5
# ---------------------------------------------------------------------------

def _mgm_markets_stub():
    def nm(v):
        return {"value": v}
    return [
        {"id": 1, "name": nm("Will there be a run in the 1st inning?"),
         "isBetBuilder": True,
         "options": [
             {"id": 10, "name": nm("Yes"), "price": {"odds": 2.15}},
             {"id": 11, "name": nm("No"), "price": {"odds": 1.67}},
         ]},
        # 3-way sibling (live 2026-08-13) — must NOT be read as a total.
        {"id": 2, "name": nm("1st Inning Runs"), "isBetBuilder": True,
         "options": [
             {"id": 20, "name": nm("0 Runs"), "price": {"odds": 1.67}},
             {"id": 21, "name": nm("1 Run"), "price": {"odds": 4.25}},
             {"id": 22, "name": nm("2 or more runs"), "price": {"odds": 4.0}},
         ]},
        {"id": 3, "name": nm("Totals"), "isBetBuilder": True,
         "options": [
             {"id": 30, "name": nm("Over 8,5"), "price": {"odds": 1.9}},
             {"id": 31, "name": nm("Under 8,5"), "price": {"odds": 1.9}},
         ]},
    ]


def test_mgm_parse_markets_builds_i1_totals_from_yrfi():
    from mlb_sgp.betmgm import parse_markets
    out = parse_markets(_mgm_markets_stub(), HOME, AWAY)
    i1 = out.get("I1")
    assert i1 is not None, "parse_markets must emit an 'I1' bucket"
    assert i1["totals"].get(0.5) == {"over": (1, 10, 2.15),
                                     "under": (1, 11, 1.67)}
    # FG regression + 3-way exclusion.
    assert out["FG"]["totals"].get(8.5) is not None
    assert 0.5 not in out["FG"]["totals"]
    for line, sides in i1["totals"].items():
        for _side, (mid, _oid, _dec) in sides.items():
            assert mid != 2, "3-way '1st Inning Runs' market must be ignored"


def test_mgm_resolve_i1_leg():
    from mlb_sgp.betmgm import resolve_legs
    structure = {
        "FG": {"spreads": {}, "totals": {},
               "moneyline": {"home": (9, 90, 2.25), "away": (9, 91, 1.65)}},
        "F5": {"spreads": {}, "totals": {}, "moneyline": None},
        "I1": {"spreads": {}, "totals": {0.5: {"over": (1, 10, 2.15),
                                               "under": (1, 11, 1.67)}},
               "moneyline": None},
    }
    out = resolve_legs(structure, [I1_UNDER, FG_ML], HOME, AWAY)
    assert out is not None
    assert out[0].ref == (1, 11)
    assert out[0].opposite_ref == (1, 10)
    assert out[1].ref == (9, 90)


# ---------------------------------------------------------------------------
# Novig — FIRST_INNING_TOTAL type at strike 0.5
# ---------------------------------------------------------------------------

def _nv_markets_stub():
    return [
        {"type": "FIRST_INNING_TOTAL", "strike": 0.5, "is_consensus": True,
         "outcomes": [
             {"id": "u-over", "description": "Over 0.5", "available": 0.525},
             {"id": "u-under", "description": "Under 0.5", "available": 0.48},
         ]},
        {"type": "TOTAL", "strike": 8.5, "is_consensus": True,
         "outcomes": [
             {"id": "u-fg-o", "description": "Over 8.5", "available": 0.5},
             {"id": "u-fg-u", "description": "Under 8.5", "available": 0.5},
         ]},
    ]


def test_novig_build_line_structure_i1_period():
    from mlb_sgp.novig import build_line_structure
    out = build_line_structure(_nv_markets_stub(), "LAD", "MIL", period="i1")
    assert out["over"].get(0.5, {}).get("id") == "u-over"
    assert out["under"].get(0.5, {}).get("id") == "u-under"
    # FG types must not leak into the i1 bucket.
    assert 8.5 not in out["over"]


def test_novig_i1_bucket_ignores_fg_and_fg_ignores_i1():
    from mlb_sgp.novig import build_line_structure
    fg = build_line_structure(_nv_markets_stub(), "LAD", "MIL", period="fg")
    assert 0.5 not in fg["over"]
    assert fg["over"].get(8.5, {}).get("id") == "u-fg-o"


def test_novig_resolve_i1_leg():
    from mlb_sgp.novig import resolve_legs
    structure = {
        "FG": {"home_spread": {}, "away_spread": {}, "over": {}, "under": {},
               "home_ml": {"id": "u-hml", "available": 0.6},
               "away_ml": {"id": "u-aml", "available": 0.42}},
        "F5": {"home_spread": {}, "away_spread": {}, "over": {}, "under": {},
               "home_ml": None, "away_ml": None},
        "I1": {"home_spread": {}, "away_spread": {},
               "over": {0.5: {"id": "u-over", "available": 0.525}},
               "under": {0.5: {"id": "u-under", "available": 0.48}},
               "home_ml": None, "away_ml": None},
    }
    out = resolve_legs(structure, [I1_OVER, FG_ML], HOME, AWAY)
    assert out is not None
    assert out[0].ref == "u-over"
    assert out[0].opposite_ref == "u-under"


# ---------------------------------------------------------------------------
# ProphetX — name alias registered; declines cleanly when board is bare
# ---------------------------------------------------------------------------

def _px_structure(markets):
    return {"event_id": "e1", "markets": markets, "home_id": 1, "away_id": 2}


def test_px_i1_total_name_alias_registered_and_resolves():
    from mlb_sgp.prophetx import resolve_legs
    markets = [
        {"id": 1500029961, "name": "1st Inning Total Runs",
         "marketLines": [{"selections": [[
             {"id": 2, "name": "Over 0.5", "line": 0.5,
              "lineID": "lid-over", "odds": -110},
             {"id": 1, "name": "Under 0.5", "line": 0.5,
              "lineID": "lid-under", "odds": -120},
         ]]}],
         "outcomes": [], "selections": []},
    ]
    out = resolve_legs(_px_structure(markets), [I1_OVER], HOME, AWAY)
    assert out is not None and len(out) == 1
    assert out[0].ref.line_id == "lid-over"
    assert out[0].opposite_ref.line_id == "lid-under"


def test_px_i1_leg_declines_cleanly_when_market_absent():
    """Live 2026-08-13: PX lists '1st Inning Total Runs' but with no priced
    lines — resolution must decline the whole book, never guess."""
    from mlb_sgp.prophetx import resolve_legs
    markets = [{"id": 9, "name": "Total Runs", "marketLines": [],
                "outcomes": [], "selections": []}]
    assert resolve_legs(_px_structure(markets), [I1_OVER], HOME, AWAY) is None


# ---------------------------------------------------------------------------
# Caesars — deliberately unmapped (WAF outage blocked name recon 2026-08-13)
# ---------------------------------------------------------------------------

def test_czr_i1_leg_declines_cleanly_without_i1_bucket():
    from mlb_sgp.caesars import resolve_legs
    structure = {
        "FG": {"spreads": {}, "totals": {}, "moneyline": None},
        "F5": {"spreads": {}, "totals": {}, "moneyline": None},
        # no "I1" key — parse_markets does not emit one
    }
    assert resolve_legs(structure, [I1_OVER], HOME, AWAY) is None
