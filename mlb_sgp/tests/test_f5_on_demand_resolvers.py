"""Issue #85: F5 legs resolve through every book's on-demand resolver.

Contract under test (new in #85): each book's ``resolve_legs`` receives the
NORMALIZED two-period structure ``{"FG": <bucket>, "F5": <bucket>}`` (keys are
``CanonicalLeg.period`` values) and picks the bucket per leg from
``leg.period``. A leg whose period bucket is missing/None fails the whole
book (clean decline, None) — never a silent FG price for an F5 leg.

ProphetX is name-driven (no period buckets): its ``_OD_MARKET_NAMES`` becomes
a ``(market_type, period) -> name aliases`` map. PX listed NO F5 markets on
the 2026-08-12 recon board, so the F5 names are the 2026-04-24 catalogue +
aliases; an F5 leg with no matching market name declines cleanly.

FG regression: the FG fixtures here are byte-copies of the per-book on-demand
test fixtures; expected refs are identical to the pre-#85 tests, proving FG
resolution is unchanged through the new structure shape.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from kalshi_common.legset import CanonicalLeg  # noqa: E402
from mlb_sgp._shared import american_to_decimal  # noqa: E402


HOME, AWAY = "New York Yankees", "Boston Red Sox"

F5_SPREAD = CanonicalLeg("g1", "spread", -0.5, "home", "F5")
F5_TOTAL = CanonicalLeg("g1", "total", 4.5, "over", "F5")
FG_SPREAD = CanonicalLeg("g1", "spread", -1.5, "home")
FG_TOTAL = CanonicalLeg("g1", "total", 8.5, "over")


# ---------------------------------------------------------------------------
# DraftKings
# ---------------------------------------------------------------------------
def _dk_structure():
    return {
        "FG": {
            "spreads": {
                ("N", 1.5, "1"): ["0HC100N150_1"],
                ("P", 1.5, "3"): ["0HC100P150_3"],
            },
            "totals": {
                ("O", 8.5): ["0OU200O850_1"],
                ("U", 8.5): ["0OU200U850_1"],
            },
            "moneyline": {"1": ["0ML300_1"], "3": ["0ML300_3"]},
            "canonical": {"100", "200", "300"},
        },
        "F5": {
            "spreads": {
                ("N", 0.5, "1"): ["0HC500N050_1"],
                ("P", 0.5, "3"): ["0HC500P050_3"],
            },
            "totals": {
                ("O", 4.5): ["0OU600O450_1"],
                ("U", 4.5): ["0OU600U450_1"],
            },
            "moneyline": {"1": ["0ML700_1"], "3": ["0ML700_3"]},
            "canonical": {"500", "600", "700"},
        },
    }


def test_dk_f5_spread_total_combo_resolves():
    from mlb_sgp.draftkings import resolve_legs
    out = resolve_legs(_dk_structure(), [F5_SPREAD, F5_TOTAL], HOME, AWAY)
    assert out is not None and len(out) == 2
    assert out[0].ref == "0HC500N050_1"
    assert out[0].opposite_ref == "0HC500P050_3"
    assert out[1].ref == "0OU600O450_1"
    assert out[1].opposite_ref == "0OU600U450_1"


def test_dk_mixed_fg_f5_set_resolves_each_from_its_bucket():
    from mlb_sgp.draftkings import resolve_legs
    out = resolve_legs(_dk_structure(), [FG_TOTAL, F5_SPREAD], HOME, AWAY)
    assert out is not None
    assert out[0].ref == "0OU200O850_1"       # FG bucket
    assert out[1].ref == "0HC500N050_1"       # F5 bucket


def test_dk_fg_resolution_unchanged_through_period_structure():
    from mlb_sgp.draftkings import resolve_legs
    out = resolve_legs(_dk_structure(), [FG_SPREAD, FG_TOTAL], HOME, AWAY)
    assert [rl.ref for rl in out] == ["0HC100N150_1", "0OU200O850_1"]


def test_dk_missing_f5_bucket_declines_whole_book():
    from mlb_sgp.draftkings import resolve_legs
    s = _dk_structure()
    s["F5"] = None                            # book posted no F5 markets
    assert resolve_legs(s, [F5_SPREAD], HOME, AWAY) is None
    # A mixed set fails as a whole too (all-or-nothing contract).
    assert resolve_legs(s, [FG_TOTAL, F5_SPREAD], HOME, AWAY) is None
    # ...but pure FG still resolves.
    assert resolve_legs(s, [FG_TOTAL], HOME, AWAY) is not None


def test_dk_f5_line_never_resolves_from_fg_bucket():
    """The wrong-price trap #84's guard existed for: an F5 leg at a line the
    FG bucket happens to carry must NOT resolve to the FG selection."""
    from mlb_sgp.draftkings import resolve_legs
    s = _dk_structure()
    s["F5"] = {"spreads": {}, "totals": {}, "moneyline": {"1": [], "3": []},
               "canonical": set()}
    leg = CanonicalLeg("g1", "spread", -1.5, "home", "F5")   # FG carries -1.5
    assert resolve_legs(s, [leg], HOME, AWAY) is None


# ---------------------------------------------------------------------------
# FanDuel
# ---------------------------------------------------------------------------
def _fd_structure():
    return {
        "FG": {
            "spreads": {("home", -1.5): ("m1", "s1", 2.2),
                        ("away", 1.5): ("m1", "s2", 1.7)},
            "totals": {("O", 8.5): ("m2", "s3", 1.95),
                       ("U", 8.5): ("m2", "s4", 1.87)},
            "moneyline": {"home": ("m3", "s5", 1.6),
                          "away": ("m3", "s6", 2.4)},
        },
        "F5": {
            "spreads": {("home", -0.5): ("m4", "s7", 1.9),
                        ("away", 0.5): ("m4", "s8", 1.9)},
            "totals": {("O", 4.5): ("m5", "s9", 2.05),
                       ("U", 4.5): ("m5", "s10", 1.79)},
            "moneyline": {"home": ("m6", "s11", 1.65),
                          "away": ("m6", "s12", 2.3)},
        },
    }


def test_fd_f5_spread_total_combo_resolves():
    from mlb_sgp.fanduel import resolve_legs
    out = resolve_legs(_fd_structure(), [F5_SPREAD, F5_TOTAL], HOME, AWAY)
    assert out is not None and len(out) == 2
    assert out[0].ref == ("m4", "s7")
    assert out[0].opposite_ref == ("m4", "s8")
    assert out[0].single_decimal == 1.9
    assert out[1].ref == ("m5", "s9")
    assert out[1].opposite_ref == ("m5", "s10")


def test_fd_mixed_fg_f5_and_fg_regression():
    from mlb_sgp.fanduel import resolve_legs
    out = resolve_legs(_fd_structure(), [FG_SPREAD, F5_TOTAL], HOME, AWAY)
    assert out is not None
    assert out[0].ref == ("m1", "s1")
    assert out[1].ref == ("m5", "s9")
    out = resolve_legs(_fd_structure(), [FG_SPREAD, FG_TOTAL], HOME, AWAY)
    assert [rl.ref for rl in out] == [("m1", "s1"), ("m2", "s3")]


def test_fd_missing_f5_bucket_declines():
    from mlb_sgp.fanduel import resolve_legs
    s = _fd_structure()
    s["F5"] = None
    assert resolve_legs(s, [F5_TOTAL], HOME, AWAY) is None
    assert resolve_legs(s, [FG_TOTAL], HOME, AWAY) is not None


# ---------------------------------------------------------------------------
# BetMGM
# ---------------------------------------------------------------------------
def _mgm_structure():
    return {
        "FG": {
            "spreads": {-1.5: {"home": ("om1", "o1", 2.25),
                               "away": ("om1", "o2", 1.67)}},
            "totals": {8.5: {"over": ("om2", "o3", 1.95),
                             "under": ("om2", "o4", 1.87)}},
            "moneyline": {"home": ("om3", "o5", 1.62),
                          "away": ("om3", "o6", 2.35)},
        },
        "F5": {
            "spreads": {-0.5: {"home": ("om4", "o7", 1.91),
                               "away": ("om4", "o8", 1.91)}},
            "totals": {4.5: {"over": ("om5", "o9", 2.0),
                             "under": ("om5", "o10", 1.83)}},
            "moneyline": None,
        },
    }


def test_mgm_f5_spread_total_combo_resolves():
    from mlb_sgp.betmgm import resolve_legs
    out = resolve_legs(_mgm_structure(), [F5_SPREAD, F5_TOTAL], HOME, AWAY)
    assert out is not None and len(out) == 2
    assert out[0].ref == ("om4", "o7")
    assert out[0].opposite_ref == ("om4", "o8")
    assert out[1].ref == ("om5", "o9")
    assert out[1].single_decimal == 2.0


def test_mgm_fg_regression_and_missing_f5_bucket():
    from mlb_sgp.betmgm import resolve_legs
    out = resolve_legs(_mgm_structure(), [FG_SPREAD, FG_TOTAL], HOME, AWAY)
    assert [rl.ref for rl in out] == [("om1", "o1"), ("om2", "o3")]
    s = _mgm_structure()
    s["F5"] = None
    assert resolve_legs(s, [F5_SPREAD], HOME, AWAY) is None
    assert resolve_legs(s, [FG_SPREAD], HOME, AWAY) is not None


# ---------------------------------------------------------------------------
# Caesars (refs ARE the /bets/details leg dicts)
# ---------------------------------------------------------------------------
def _czr_leg(name, dec):
    return {"name": name, "price": {"d": dec}}


def _czr_structure():
    return {
        "FG": {
            "spreads": {-1.5: {"home": _czr_leg("fg-h", 2.2),
                               "away": _czr_leg("fg-a", 1.7)}},
            "totals": {8.5: {"over": _czr_leg("fg-o", 1.95),
                             "under": _czr_leg("fg-u", 1.87)}},
            "moneyline": {"home": _czr_leg("fg-hml", 1.6),
                          "away": _czr_leg("fg-aml", 2.4)},
        },
        "F5": {
            "spreads": {-0.5: {"home": _czr_leg("f5-h", 1.9),
                               "away": _czr_leg("f5-a", 1.9)}},
            "totals": {4.5: {"over": _czr_leg("f5-o", 2.05),
                             "under": _czr_leg("f5-u", 1.79)}},
            "moneyline": None,
        },
    }


def test_czr_f5_spread_total_combo_resolves():
    from mlb_sgp.caesars import resolve_legs
    out = resolve_legs(_czr_structure(), [F5_SPREAD, F5_TOTAL], HOME, AWAY)
    assert out is not None and len(out) == 2
    assert out[0].ref["name"] == "f5-h"
    assert out[0].opposite_ref["name"] == "f5-a"
    assert out[1].ref["name"] == "f5-o"
    assert out[1].single_decimal == 2.05


def test_czr_fg_regression_and_missing_f5_bucket():
    from mlb_sgp.caesars import resolve_legs
    out = resolve_legs(_czr_structure(), [FG_SPREAD, FG_TOTAL], HOME, AWAY)
    assert [rl.ref["name"] for rl in out] == ["fg-h", "fg-o"]
    s = _czr_structure()
    s["F5"] = None
    assert resolve_legs(s, [F5_TOTAL], HOME, AWAY) is None
    assert resolve_legs(s, [FG_TOTAL], HOME, AWAY) is not None


# ---------------------------------------------------------------------------
# Novig (per-period buckets built by build_line_structure(period=...))
# ---------------------------------------------------------------------------
def _nv_leg(uuid, available):
    return {"id": uuid, "available": available}


def _nv_structure():
    return {
        "FG": {
            "home_spread": {-1.5: _nv_leg("hs-15", 0.60)},
            "away_spread": {1.5: _nv_leg("as+15", 0.44)},
            "over": {8.5: _nv_leg("ov85", 0.52)},
            "under": {8.5: _nv_leg("un85", 0.50)},
            "home_ml": _nv_leg("hml", 0.62),
            "away_ml": _nv_leg("aml", 0.42),
        },
        "F5": {
            "home_spread": {-0.5: _nv_leg("f5hs-05", 0.55)},
            "away_spread": {0.5: _nv_leg("f5as+05", 0.49)},
            "over": {4.5: _nv_leg("f5ov45", 0.51)},
            "under": {4.5: _nv_leg("f5un45", 0.51)},
            "home_ml": _nv_leg("f5hml", 0.58),
            "away_ml": _nv_leg("f5aml", 0.46),
        },
    }


def test_nv_f5_spread_total_combo_resolves():
    from mlb_sgp.novig import resolve_legs
    out = resolve_legs(_nv_structure(), [F5_SPREAD, F5_TOTAL], HOME, AWAY)
    assert out is not None and len(out) == 2
    assert out[0].ref == "f5hs-05"
    assert out[0].opposite_ref == "f5as+05"
    assert out[1].ref == "f5ov45"
    assert out[1].single_decimal == 1.0 / 0.51


def test_nv_fg_regression_and_missing_f5_bucket():
    from mlb_sgp.novig import resolve_legs
    out = resolve_legs(_nv_structure(), [FG_SPREAD, FG_TOTAL], HOME, AWAY)
    assert [rl.ref for rl in out] == ["hs-15", "ov85"]
    s = _nv_structure()
    s["F5"] = None
    assert resolve_legs(s, [F5_SPREAD], HOME, AWAY) is None
    assert resolve_legs(s, [FG_SPREAD], HOME, AWAY) is not None


def test_nv_build_line_structure_f5_reads_1h_market_types():
    """Period pass-through: SPREAD_1H/TOTAL_1H/MONEY_1H markets land in the
    f5 bucket (live-verified 2026-08-12: _1H strikes sit at F5 scale —
    TOTAL_1H 3.5/4.5 vs FG TOTAL 7.5-12.5, SPREAD_1H ±0.5)."""
    from mlb_sgp.novig import build_line_structure
    markets = [
        {"type": "SPREAD_1H", "strike": -0.5, "is_consensus": True,
         "outcomes": [
             {"id": "nvh", "description": "NYY -0.5", "available": 0.55,
              "competitor": {"symbol": "NYY"}},
             {"id": "nva", "description": "BOS +0.5", "available": 0.49,
              "competitor": {"symbol": "BOS"}},
         ]},
        {"type": "TOTAL_1H", "strike": 4.5, "is_consensus": True,
         "outcomes": [
             {"id": "nvo", "description": "Over 4.5", "available": 0.51},
             {"id": "nvu", "description": "Under 4.5", "available": 0.51},
         ]},
        # FG market at the same strike must NOT leak into the f5 bucket.
        {"type": "TOTAL", "strike": 4.5, "is_consensus": True,
         "outcomes": [
             {"id": "fgo", "description": "Over 4.5", "available": 0.9},
             {"id": "fgu", "description": "Under 4.5", "available": 0.12},
         ]},
    ]
    out = build_line_structure(markets, "NYY", "BOS", period="f5")
    assert out["home_spread"][-0.5]["id"] == "nvh"
    assert out["away_spread"][0.5]["id"] == "nva"
    assert out["over"][4.5]["id"] == "nvo"
    assert out["under"][4.5]["id"] == "nvu"


# ---------------------------------------------------------------------------
# ProphetX ((market_type, period) name map with aliases)
# ---------------------------------------------------------------------------
PX_HOME_ID, PX_AWAY_ID = 111, 222


def _px_f5_markets(spread_name="1st-5th Inning Spread",
                   total_name="1st-5th Inning Total Runs"):
    return [
        {
            "id": 40, "name": spread_name,
            "marketLines": [
                {"outcomes": [
                    {"id": 4001, "lineID": "L5-h-05", "line": -0.5,
                     "competitorId": PX_HOME_ID, "name": "Yankees -0.5",
                     "odds": -105},
                    {"id": 4002, "lineID": "L5-a+05", "line": 0.5,
                     "competitorId": PX_AWAY_ID, "name": "Red Sox +0.5",
                     "odds": -115},
                ]},
            ],
            "outcomes": [], "selections": [],
        },
        {
            "id": 50, "name": total_name,
            "marketLines": [
                {"outcomes": [
                    {"id": 5001, "lineID": "L5-o45", "line": 4.5,
                     "name": "Over 4.5", "odds": -110},
                    {"id": 5002, "lineID": "L5-u45", "line": 4.5,
                     "name": "Under 4.5", "odds": -110},
                ]},
            ],
            "outcomes": [], "selections": [],
        },
    ]


def _px_fg_markets():
    return [
        {
            "id": 10, "name": "Run Line",
            "marketLines": [
                {"outcomes": [
                    {"id": 1001, "lineID": "L-h-15", "line": -1.5,
                     "competitorId": PX_HOME_ID, "name": "Yankees -1.5",
                     "odds": -115},
                    {"id": 1002, "lineID": "L-a+15", "line": 1.5,
                     "competitorId": PX_AWAY_ID, "name": "Red Sox +1.5",
                     "odds": -105},
                ]},
            ],
            "outcomes": [], "selections": [],
        },
    ]


def _px_structure(markets):
    return {"event_id": "ev900", "markets": markets,
            "home_id": PX_HOME_ID, "away_id": PX_AWAY_ID}


def test_px_f5_spread_total_combo_resolves():
    from mlb_sgp.prophetx import resolve_legs
    s = _px_structure(_px_fg_markets() + _px_f5_markets())
    out = resolve_legs(s, [F5_SPREAD, F5_TOTAL], HOME, AWAY)
    assert out is not None and len(out) == 2
    assert out[0].ref.outcome_id == "4001"
    assert out[0].opposite_ref.outcome_id == "4002"
    assert out[1].ref.outcome_id == "5001"
    assert out[1].single_decimal == american_to_decimal(-110)


def test_px_f5_alias_names_dispatch():
    from mlb_sgp.prophetx import resolve_legs
    s = _px_structure(_px_f5_markets(spread_name="1st 5 Innings Spread",
                                     total_name="F5 Total Runs"))
    out = resolve_legs(s, [F5_SPREAD, F5_TOTAL], HOME, AWAY)
    assert out is not None
    assert out[0].ref.outcome_id == "4001"
    assert out[1].ref.outcome_id == "5001"


def test_px_no_f5_markets_declines_and_fg_regression():
    from mlb_sgp.prophetx import resolve_legs
    s = _px_structure(_px_fg_markets())
    # 2026-08-12 recon: PX board carries zero F5 markets — must decline.
    assert resolve_legs(s, [F5_SPREAD], HOME, AWAY) is None
    out = resolve_legs(s, [FG_SPREAD], HOME, AWAY)
    assert out[0].ref.outcome_id == "1001"


def test_px_f5_moneyline_has_no_known_name_declines():
    from mlb_sgp.prophetx import resolve_legs
    s = _px_structure(_px_fg_markets() + _px_f5_markets())
    leg = CanonicalLeg("g1", "ml", None, "home", "F5")
    assert resolve_legs(s, [leg], HOME, AWAY) is None


# ---------------------------------------------------------------------------
# Issue #86: F5-winner (ml, F5) legs — resolution + tie anchors per book
# ---------------------------------------------------------------------------
F5_ML = CanonicalLeg("g1", "ml", None, "home", "F5")


def test_fd_f5_moneyline_resolves_from_f5_bucket():
    from mlb_sgp.fanduel import resolve_legs
    out = resolve_legs(_fd_structure(), [F5_ML, F5_TOTAL], HOME, AWAY)
    assert out is not None and len(out) == 2
    assert out[0].ref == ("m6", "s11")
    assert out[0].opposite_ref == ("m6", "s12")


def test_nv_f5_moneyline_resolves_from_f5_bucket():
    from mlb_sgp.novig import resolve_legs
    out = resolve_legs(_nv_structure(), [F5_ML, F5_TOTAL], HOME, AWAY)
    assert out is not None and len(out) == 2
    assert out[0].ref == "f5hml"


def test_dk_f5_moneyline_resolves_from_f5_bucket():
    from mlb_sgp.draftkings import resolve_legs
    out = resolve_legs(_dk_structure(), [F5_ML, F5_TOTAL], HOME, AWAY)
    assert out is not None and len(out) == 2
    assert out[0].ref == "0ML700_1"
    assert out[0].opposite_ref == "0ML700_3"


def test_fd_tie_prob_devigs_the_result3_bucket():
    """FD's "First 5 Innings Result" runners (home/tie/away decimals) ride
    the F5 bucket as ``result3``; f5_tie_prob deviggs P(tie) off them."""
    from mlb_sgp.fanduel import f5_tie_prob
    s = _fd_structure()
    s["F5"]["result3"] = {"home": 1.72463768115942, "tie": 6.5, "away": 2.88}
    out = f5_tie_prob(s)
    assert out is not None
    assert 0.10 < out < 0.17


def test_fd_tie_prob_none_when_result3_missing_or_partial():
    from mlb_sgp.fanduel import f5_tie_prob
    assert f5_tie_prob(_fd_structure()) is None            # no result3 key
    s = _fd_structure()
    s["F5"]["result3"] = {"home": 1.72, "away": 2.88}      # tie runner absent
    assert f5_tie_prob(s) is None
    s["F5"] = None                                         # no F5 board
    assert f5_tie_prob(s) is None


def test_fd_market_map_carries_first_5_innings_result():
    from scraper_fanduel_sgp import _MARKET_MAP
    assert _MARKET_MAP.get("First 5 Innings Result") == ("f5", "result3", "main")


def test_mgm_parse_markets_captures_f5_1x2_and_tie_prob():
    """BetMGM "1st 5 Innings - 1X2" options land as F5 result3 decimals;
    f5_tie_prob deviggs them. Live recon 2026-08-12: Royals 2.65 / Tie 6.75
    / Dodgers 1.8 (sum 1.0811)."""
    from mlb_sgp.betmgm import parse_markets, f5_tie_prob

    def opt(name, oid, odds):
        return {"name": {"value": name}, "id": oid,
                "price": {"odds": odds}}

    markets = [{
        "name": {"value": "1st 5 Innings - 1X2"},
        "id": "om9",
        "isBetBuilder": True,
        "options": [opt("Royals", "x1", 2.65), opt("Tie", "x2", 6.75),
                    opt("Dodgers", "x3", 1.8)],
    }]
    s = parse_markets(markets, "Los Angeles Dodgers", "Kansas City Royals")
    assert s["F5"]["result3"] == {"home": 1.8, "tie": 6.75, "away": 2.65}
    out = f5_tie_prob(s)
    assert out is not None
    assert 0.10 < out < 0.18


def test_mgm_parse_markets_result3_absent_leaves_tie_prob_none():
    from mlb_sgp.betmgm import parse_markets, f5_tie_prob
    s = parse_markets([], "Los Angeles Dodgers", "Kansas City Royals")
    assert f5_tie_prob(s) is None


def test_czr_classifies_f5_moneyline_names():
    """CZR posts its F5 board near lineups; the ML name mirrors its posted
    "1st 5 Innings Run Line"/"Total Runs" family (and its live "1st 3
    Innings Money Line"). Both plain and dashed variants classify."""
    from mlb_sgp.caesars import parse_markets

    def mkt(name, sels):
        return {"name": name, "id": "m1", "selections": sels}

    sels = [{"type": "home", "name": "|LAD|", "price": {"d": 1.6}},
            {"type": "away", "name": "|KC|", "price": {"d": 2.3}}]
    for name in ("1st 5 Innings Money Line", "1st 5 innings - money line"):
        event = {"id": "e1", "competitionId": "c1",
                 "keyMarketGroups": [{"markets": [mkt(name, sels)]}]}
        s = parse_markets(event)
        assert s["F5"]["moneyline"] is not None, name
        assert s["F5"]["moneyline"]["home"]["price"]["d"] == 1.6


def test_czr_f5_moneyline_resolves_once_parsed():
    from mlb_sgp.caesars import resolve_legs
    s = _czr_structure()
    s["F5"]["moneyline"] = {"home": _czr_leg("f5-hml", 1.6),
                            "away": _czr_leg("f5-aml", 2.3)}
    out = resolve_legs(s, [F5_ML, F5_TOTAL], HOME, AWAY)
    assert out is not None and len(out) == 2
    assert out[0].ref["name"] == "f5-hml"


def test_dk_main_market_nums_accepts_bare_1st_5_innings_as_f5_ml(monkeypatch):
    """DK's live board (2026-08-12) names the F5 winner market bare
    "1st 5 Innings" — no "Moneyline - 1st 5 Innings" row existed. Both
    names must seed the F5 moneyline market_num (bare name only as
    fallback, never overriding an explicit Moneyline row)."""
    import scraper_draftkings_sgp as dk

    def fake_subcats(session, event_id, subcat_id, profile=None):
        if subcat_id == "4519":
            return [("m_1", "Run Line"), ("m_2", "Total"), ("m_3", "Moneyline")]
        return [("2_100", "Run Line - 1st 5 Innings"),
                ("3_200", "Total Runs - 1st 5 Innings"),
                ("359435579", "1st 5 Innings")]

    monkeypatch.setattr(dk, "_fetch_subcat_markets", fake_subcats)
    out = dk.fetch_main_market_nums(None, "evt")
    assert out["f5"]["moneyline"] == "359435579"


def test_dk_explicit_f5_moneyline_name_beats_bare_alias(monkeypatch):
    import scraper_draftkings_sgp as dk

    def fake_subcats(session, event_id, subcat_id, profile=None):
        if subcat_id == "4519":
            return []
        return [("359435579", "1st 5 Innings"),
                ("4_300", "Moneyline - 1st 5 Innings")]

    monkeypatch.setattr(dk, "_fetch_subcat_markets", fake_subcats)
    out = dk.fetch_main_market_nums(None, "evt")
    assert out["f5"]["moneyline"] == "300"
