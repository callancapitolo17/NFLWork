"""Caesars moneyline-in-SGP: classify, parse, and ML×total emission.

Caesars selections are typed home/away/over/under directly, so the moneyline
market reuses the exact mechanism proven for the run-line market — the only
differences are the market name and that ML markets carry no `line` (handled
before the line-None skip). price_sgps prices the 4 FG ML×total combos via the
same /bets/details endpoint, tagging spread_line=None.

NB: the exact Caesars ML market-name string was not confirmable live at build
time (WAF IP-throttle); the classifier accepts the common variants
{"moneyline", "money line"}. A wrong name fails closed (no rows, like a
throttled book), never bad data.
"""
import sys
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock

sys.path.insert(0, str(Path(__file__).parent.parent))

from mlb_sgp._shared import TargetLine  # noqa: E402


def test_classify_moneyline():
    from caesars import _classify
    assert _classify("Moneyline") == ("FG", "moneyline")
    assert _classify("Money Line") == ("FG", "moneyline")
    # Run line / total still classify; unrelated names do not.
    assert _classify("Run Line") == ("FG", "spread")
    assert _classify("1st 5 Innings Moneyline") is None


def _sel(t, sid, d):
    return {"id": sid, "type": t, "name": t, "price": {"d": d}}


def _mkt(name, line, sels, mid):
    return {"id": mid, "name": name, "line": line, "selections": sels}


def _event():
    return {"id": "e1", "competitionId": "c1", "keyMarketGroups": [{"markets": [
        _mkt("Moneyline", None,
             [_sel("home", "mh", 1.8), _sel("away", "ma", 2.1)], "m_ml"),
        _mkt("Run Line", -1.5,
             [_sel("home", "sh", 2.4), _sel("away", "sa", 1.6)], "m_rl"),
        _mkt("Total Runs", 8.5,
             [_sel("over", "ov", 1.9), _sel("under", "un", 1.9)], "m_tot"),
    ]}]}


def test_parse_extracts_moneyline():
    from caesars import parse_markets
    parsed = parse_markets(_event())
    ml = parsed["FG"]["moneyline"]
    assert ml is not None
    assert ml["home"]["selectionId"] == "mh"
    assert ml["away"]["selectionId"] == "ma"
    # Spread + total still parse alongside.
    assert -1.5 in parsed["FG"]["spreads"]
    assert 8.5 in parsed["FG"]["totals"]


def test_price_sgps_emits_ml_rows(monkeypatch):
    import caesars
    from caesars_client import Event

    ev = Event(event_id="e1", competition_id="c1",
               home_team="Tigers", away_team="Astros", start_time="")
    monkeypatch.setattr(caesars, "_match_events", lambda events, targets: {"g1": ev})

    client = MagicMock()
    client.list_events.return_value = [ev]
    client.fetch_event.return_value = _event()
    client.price_combo.return_value = {"decimal": 3.5, "american": 250}

    targets = [TargetLine("g1", "Tigers", "Astros",
                          datetime(2026, 6, 25, 23, 0, tzinfo=timezone.utc),
                          "FG", -1.5, 8.5)]
    rows = caesars.price_sgps(targets, client=client, verbose=False)

    ml_rows = [r for r in rows if "ML" in r.combo]
    assert {r.combo for r in ml_rows} == {
        "Home ML + Over", "Home ML + Under", "Away ML + Over", "Away ML + Under"}
    for r in ml_rows:
        assert r.spread_line is None
        assert r.total_line == 8.5
        assert r.bookmaker == "caesars"
