"""Phase 1: parse Odds API h2h/spreads/totals into mlb_singles_odds rows."""
from datetime import datetime, timezone

from kalshi_common.sgp_runner import _parse_singles_payload

_FT = datetime(2026, 6, 22, 20, 0, tzinfo=timezone.utc)

_PAYLOAD = [{
    "id": "game123",
    "home_team": "Los Angeles Angels",
    "away_team": "Texas Rangers",
    "bookmakers": [{
        "key": "draftkings",
        "markets": [
            {"key": "h2h", "outcomes": [
                {"name": "Los Angeles Angels", "price": 1.80},
                {"name": "Texas Rangers", "price": 2.10},
            ]},
            {"key": "spreads", "outcomes": [
                {"name": "Los Angeles Angels", "price": 2.05, "point": -1.5},
                {"name": "Texas Rangers", "price": 1.80, "point": 1.5},
            ]},
            {"key": "totals", "outcomes": [
                {"name": "Over", "price": 1.95, "point": 8.5},
                {"name": "Under", "price": 1.88, "point": 8.5},
            ]},
        ],
    }],
}]


def _rows():
    return _parse_singles_payload(_PAYLOAD, _FT)


def test_moneyline_rows_home_away():
    rows = _rows()
    ml = [r for r in rows if r["market"] == "moneyline"]
    assert {r["outcome"] for r in ml} == {"home", "away"}
    home = next(r for r in ml if r["outcome"] == "home")
    assert home["game_id"] == "game123"
    assert home["bookmaker"] == "draftkings"
    assert home["decimal"] == 1.80
    assert home["line"] is None


def test_spread_rows_home_perspective_line():
    rows = _rows()
    sp = [r for r in rows if r["market"] == "spread"]
    home = next(r for r in sp if r["outcome"] == "home")
    away = next(r for r in sp if r["outcome"] == "away")
    assert home["line"] == -1.5    # home favored by 1.5 (home-perspective)
    assert away["line"] == -1.5    # away row stores the SAME home-perspective line
    assert home["decimal"] == 2.05


def test_total_rows_over_under():
    rows = _rows()
    tot = [r for r in rows if r["market"] == "total"]
    over = next(r for r in tot if r["outcome"] == "over")
    under = next(r for r in tot if r["outcome"] == "under")
    assert over["line"] == 8.5 and under["line"] == 8.5
    assert over["decimal"] == 1.95


def test_unknown_team_name_dropped():
    bad = [{"id": "g", "home_team": "A", "away_team": "B",
            "bookmakers": [{"key": "dk", "markets": [
                {"key": "h2h", "outcomes": [
                    {"name": "Someone Else", "price": 2.0},
                    {"name": "B", "price": 1.8}]}]}]}]
    rows = _parse_singles_payload(bad, _FT)
    # The h2h pair is incomplete (one name unmatched) -> no moneyline rows emitted
    assert [r for r in rows if r["market"] == "moneyline"] == []
