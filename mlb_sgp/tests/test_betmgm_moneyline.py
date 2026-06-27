"""BetMGM moneyline-in-SGP: classify, parse, and ML×total emission.

BetMGM names its moneyline market "Money Line" (FG) / "1st 5 innings - Money
Line" (F5). ML option names use SHORT team names ("Giants") while the canonical
home/away names are full ("San Francisco Giants"), so the parser matches
bidirectionally. price_sgps prices the 4 FG ML×total combos via the same
price_picks endpoint, tagging spread_line=None.
"""
import sys
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock

sys.path.insert(0, str(Path(__file__).parent.parent))

from mlb_sgp._shared import TargetLine  # noqa: E402


def test_classify_moneyline():
    from betmgm import _classify_market
    assert _classify_market("Money Line") == ("FG", "moneyline")
    assert _classify_market("1st 5 innings - Money Line") == ("F5", "moneyline")
    # Other-period / margin markets must NOT classify as the game moneyline.
    assert _classify_market("1st 7 Innings - Money Line") is None
    assert _classify_market("Winning Margin") is None


def _om(name, opts, bb=True):
    return {"name": {"value": name}, "id": abs(hash(name)) % 100000,
            "isBetBuilder": bb, "options": opts}


def _opt(name, oid, odds, attr=None):
    return {"name": {"value": name}, "id": oid, "attr": attr,
            "price": {"odds": odds}}


def _raw_markets():
    return [
        _om("Money Line", [_opt("Giants", 11, 2.5), _opt("Athletics", 12, 1.53)]),
        _om("Run Line Spread", [
            _opt("San Francisco Giants -1.5", 13, 2.4, "-1,5"),
            _opt("Athletics +1.5", 14, 1.6, "+1,5")]),
        _om("Totals", [_opt("Over 8.5", 15, 1.9), _opt("Under 8.5", 16, 1.9)]),
    ]


def test_parse_extracts_moneyline_bidirectional():
    from betmgm import parse_markets
    parsed = parse_markets(_raw_markets(), "San Francisco Giants", "Athletics")
    ml = parsed["FG"]["moneyline"]
    assert ml is not None
    # "Giants" matched the full home name; "Athletics" matched away exactly.
    assert ml["home"][1] == 11 and ml["home"][2] == 2.5
    assert ml["away"][1] == 12 and ml["away"][2] == 1.53


def test_price_sgps_emits_ml_rows(monkeypatch):
    import betmgm
    from betmgm_client import Event

    ev = Event(event_id="e1", home_team="San Francisco Giants",
               away_team="Athletics", start_time="")
    monkeypatch.setattr(betmgm, "_match_events", lambda events, targets: {"g1": ev})

    client = MagicMock()
    client.accessid.return_value = "acc"
    client.list_events.return_value = [ev]
    client.fetch_markets.return_value = _raw_markets()
    client.price_picks.return_value = {"decimal": 3.5, "american": 250}

    targets = [TargetLine("g1", "San Francisco Giants", "Athletics",
                          datetime(2026, 6, 25, 23, 0, tzinfo=timezone.utc),
                          "FG", -1.5, 8.5)]
    rows = betmgm.price_sgps(targets, client=client, verbose=False)

    ml_rows = [r for r in rows if "ML" in r.combo]
    assert {r.combo for r in ml_rows} == {
        "Home ML + Over", "Home ML + Under", "Away ML + Over", "Away ML + Under"}
    for r in ml_rows:
        assert r.spread_line is None
        assert r.total_line == 8.5
        assert r.bookmaker == "betmgm"
