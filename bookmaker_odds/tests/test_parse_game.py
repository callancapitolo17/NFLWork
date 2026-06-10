"""Tests for _parse_game — BKM Derivatives.line -> DuckDB records.

BKM packs the main and all alt run lines / alt totals into one
`Derivatives.line` array keyed by `index` (0 = main, non-zero = alts). The
parser now emits one record per usable line index: the main keeps the ML and
the original market label ("spreads"/"spreads_f5"/"spreads_h2"); non-zero
indices emit alt rows with market = "alternate_spreads"-prefixed label and
no ML (ML is per-game, not per-line).
"""
import sys
from pathlib import Path

# Allow `import scraper` when pytest runs from the repo root.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import scraper  # noqa: E402


def _game(*lines, vtm="Houston Astros", htm="Texas Rangers",
          idgm="123", gmdt="05/28", gmtm="20:05"):
    """Helper to build a minimal BKM game dict."""
    return {
        "Stat": "O",
        "vtm": vtm,
        "htm": htm,
        "idgm": idgm,
        "gmdt": gmdt,
        "gmtm": gmtm,
        "Derivatives": {"line": list(lines)},
    }


def _line(index, *, vsprdt=None, hsprdt=None, vsprdoddst=None, hsprdoddst=None,
          ovt=None, ovoddst=None, unoddst=None, voddst=None, hoddst=None):
    """Build a Derivatives.line entry. Status flags auto-set from fields."""
    line = {"index": str(index)}
    if vsprdt is not None:
        line.update({"s_sp": 1, "vsprdt": vsprdt, "hsprdt": hsprdt,
                     "vsprdoddst": vsprdoddst, "hsprdoddst": hsprdoddst})
    else:
        line["s_sp"] = 0
    if ovt is not None:
        line.update({"s_tot": 1, "ovt": ovt, "ovoddst": ovoddst, "unoddst": unoddst})
    else:
        line["s_tot"] = 0
    if voddst is not None:
        line.update({"s_ml": 1, "voddst": voddst, "hoddst": hoddst})
    else:
        line["s_ml"] = 0
    return line


def test_single_line_game_emits_one_main_record():
    """A game with only index 0 returns one record carrying spread + total + ML."""
    game = _game(_line(0,
                       vsprdt="-1.5", hsprdt="1.5", vsprdoddst="106", hsprdoddst="-125",
                       ovt="9.5", ovoddst="-102", unoddst="-115",
                       voddst="-150", hoddst="133"))
    out = scraper._parse_game(game, market="spreads", period="fg",
                              sport_key="baseball_mlb",
                              team_dict={}, canonical_games=[],
                              fetch_time="2026-05-29T03:00:00Z")
    assert len(out) == 1
    r = out[0]
    assert r["market"] == "spreads"
    assert r["period"] == "fg"
    assert r["away_spread"] == -1.5
    assert r["home_spread"] == 1.5
    assert r["total"] == 9.5
    assert r["away_ml"] == -150
    assert r["home_ml"] == 133


def test_multi_index_game_emits_main_plus_alts():
    """7-index Derivatives.line -> 1 main + 6 alt records.

    Mirrors the HOU @ TEX live recon: spreads -2.5..-5.5 paired with totals
    8.5..11.5 across indices -3..+3.
    """
    lines = []
    spreads = {-3: ("-2.5", "2.5", "-300", "223"),
               -2: ("-3",   "3",   "-236", "184"),
               -1: ("-3.5", "3.5", "-160", "130"),
                0: ("-4",   "4",   "-105", "-114"),
                1: ("-4.5", "4.5", "133",  "-163"),
                2: ("-5",   "5",   "201",  "-260"),
                3: ("-5.5", "5.5", "239",  "-326")}
    totals = {-3: ("8.5", "-105", "-115"), -2: ("9", "-105", "-115"),
              -1: ("9.5", "-105", "-115"),  0: ("10", "-105", "-115"),
               1: ("10.5", "-105", "-115"), 2: ("11", "-105", "-115"),
               3: ("11.5", "-105", "-115")}
    for idx in sorted(spreads):
        v, h, vp, hp = spreads[idx]
        ot, ov, un = totals[idx]
        # ML only on main per BKM convention
        ml = {"voddst": "-1076", "hoddst": "700"} if idx == 0 else {}
        lines.append(_line(idx, vsprdt=v, hsprdt=h, vsprdoddst=vp, hsprdoddst=hp,
                           ovt=ot, ovoddst=ov, unoddst=un, **ml))
    game = _game(*lines)

    out = scraper._parse_game(game, market="spreads", period="fg",
                              sport_key="baseball_mlb",
                              team_dict={}, canonical_games=[],
                              fetch_time="2026-05-29T03:00:00Z")

    assert len(out) == 7
    mains = [r for r in out if r["market"] == "spreads"]
    alts = [r for r in out if r["market"] == "alternate_spreads"]
    assert len(mains) == 1 and len(alts) == 6

    # Main has ML; alts do not.
    m = mains[0]
    assert m["away_spread"] == -4.0
    assert m["total"] == 10.0
    assert m["away_ml"] == -1076

    for a in alts:
        assert a["away_ml"] is None
        assert a["home_ml"] is None
        # Alt rows carry both spread and total at that index.
        assert a["away_spread"] is not None
        assert a["total"] is not None

    alt_spreads = sorted(r["away_spread"] for r in alts)
    assert alt_spreads == [-5.5, -5.0, -4.5, -3.5, -3.0, -2.5]
    alt_totals = sorted(r["total"] for r in alts)
    assert alt_totals == [8.5, 9.0, 9.5, 10.5, 11.0, 11.5]


def test_alt_market_label_includes_period_suffix():
    """For non-FG leagues, the alt market label preserves the period suffix."""
    game = _game(
        _line(0, vsprdt="-0.5", hsprdt="0.5", vsprdoddst="100", hsprdoddst="-120"),
        _line(1, vsprdt="-1",   hsprdt="1",   vsprdoddst="130", hsprdoddst="-150"),
    )
    out = scraper._parse_game(game, market="spreads_f5", period="F5",
                              sport_key="baseball_mlb",
                              team_dict={}, canonical_games=[],
                              fetch_time="2026-05-29T03:00:00Z")
    assert {r["market"] for r in out} == {"spreads_f5", "alternate_spreads_f5"}


def test_stat_not_open_returns_empty_list():
    """Games not in status 'O' (open) yield no records."""
    game = _game(_line(0, vsprdt="-1.5", hsprdt="1.5",
                       vsprdoddst="106", hsprdoddst="-125"))
    game["Stat"] = "C"
    out = scraper._parse_game(game, market="spreads", period="fg",
                              sport_key="baseball_mlb",
                              team_dict={}, canonical_games=[],
                              fetch_time="2026-05-29T03:00:00Z")
    assert out == []


def test_alt_index_with_no_usable_data_is_skipped():
    """If a non-zero index has neither spread nor total populated, it's dropped
    (we don't emit empty alt rows). Main is still emitted as before."""
    game = _game(
        _line(0, vsprdt="-1.5", hsprdt="1.5",
              vsprdoddst="106", hsprdoddst="-125",
              ovt="9.5", ovoddst="-102", unoddst="-115",
              voddst="-150", hoddst="133"),
        _line(1),  # s_sp=0, s_tot=0, s_ml=0 — nothing
    )
    out = scraper._parse_game(game, market="spreads", period="fg",
                              sport_key="baseball_mlb",
                              team_dict={}, canonical_games=[],
                              fetch_time="2026-05-29T03:00:00Z")
    assert len(out) == 1
    assert out[0]["market"] == "spreads"
