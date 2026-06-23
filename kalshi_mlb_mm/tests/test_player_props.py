"""Phase 4: player-prop odds parsing + name normalization."""
from datetime import datetime, timezone

from kalshi_common.sgp_runner import _norm_player, _parse_props_event

_FT = datetime(2026, 6, 23, 20, 0, tzinfo=timezone.utc)


def test_norm_player_strips_accents_and_punct():
    # Kalshi titles and Odds API descriptions must normalize to the same key.
    assert _norm_player("Julio Rodríguez") == _norm_player("Julio Rodriguez")
    assert _norm_player("Vladimir Guerrero Jr.") == "vladimir guerrero"
    assert _norm_player("Shohei Ohtani") == "shohei ohtani"
    assert _norm_player("J.D. Martinez") == "jd martinez"


_EVENT = {
    "id": "gP",
    "home_team": "Seattle Mariners",
    "away_team": "Pittsburgh Pirates",
    "bookmakers": [
        {"key": "draftkings", "markets": [
            {"key": "batter_home_runs", "outcomes": [
                {"name": "Over", "description": "Cal Raleigh", "price": 3.5, "point": 0.5},
                {"name": "Under", "description": "Cal Raleigh", "price": 1.30, "point": 0.5},
            ]},
        ]},
        {"key": "fanduel", "markets": [
            {"key": "batter_home_runs", "outcomes": [
                {"name": "Over", "description": "Cal Raleigh", "price": 3.6, "point": 0.5},
                {"name": "Under", "description": "Cal Raleigh", "price": 1.28, "point": 0.5},
            ]},
            {"key": "pitcher_strikeouts", "outcomes": [
                {"name": "Over", "description": "Some Pitcher", "price": 1.9, "point": 5.5},
                {"name": "Under", "description": "Some Pitcher", "price": 1.9, "point": 5.5},
            ]},
        ]},
    ],
}


def _rows():
    return _parse_props_event(_EVENT, _FT)


def test_props_market_type_mapping():
    rows = _rows()
    types = {r["market_type"] for r in rows}
    assert "home_runs" in types
    assert "strikeouts" in types


def test_props_player_normalized_and_line():
    rows = _rows()
    hr = [r for r in rows if r["market_type"] == "home_runs"]
    assert all(r["player"] == "cal raleigh" for r in hr)
    assert all(r["line"] == 0.5 for r in hr)
    assert {r["outcome"] for r in hr} == {"over", "under"}
    assert {r["bookmaker"] for r in hr} == {"draftkings", "fanduel"}


def test_props_incomplete_pair_dropped():
    ev = {"id": "g", "home_team": "A", "away_team": "B", "bookmakers": [
        {"key": "dk", "markets": [{"key": "batter_hits", "outcomes": [
            {"name": "Over", "description": "X Y", "price": 2.0, "point": 1.5}]}]}]}
    rows = _parse_props_event(ev, _FT)
    assert rows == []   # only over present -> can't devig -> dropped


def test_props_game_id_carried():
    rows = _rows()
    assert all(r["game_id"] == "gP" for r in rows)
