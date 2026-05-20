"""Unit tests for resolve_leagues() — the dynamic league discovery logic.

These tests pin down the *disambiguation contract* introduced when we moved
BKM's scraper off hardcoded leagueIds:

  1. (sportId, region, leagueDescEn) is the primary key — sportId alone is
     ambiguous (BASKETBALL/sportId=NBA has both region=NBA and region=WNBA
     under it, so matching on (sportId, leagueDescEn) would pick one
     arbitrarily based on dict insertion order).
  2. Missing markets are reported with a clear diagnostic, not silently
     skipped — and the diagnostic shows what IS available at the requested
     (sportId, region) scope so the operator can fix the config.
  3. Empty / missing routing_info returns an empty list without crashing.

All tests are pure Python — no network, no DuckDB.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import scraper  # noqa: E402


def _routing(*leagues):
    """Build a minimal routing_info dict from a flat list of league entries."""
    return {"valid": True, "routedSports": [{"routedLeagues": list(leagues)}]}


def _league(sport_id, region, ldEn, league_id):
    return {"sportId": sport_id, "region": region,
            "leagueDescEn": ldEn, "leagueId": league_id}


# ---- The headline disambiguation case --------------------------------------

def test_nba_wnba_share_sportId_resolves_to_correct_region():
    """NBA and WNBA both have sportId='NBA' — region is the disambiguator.

    This is the bug the dynamic-discovery PoC caught: keying on just
    (sportId, leagueDescEn) lets WNBA's entry overwrite NBA's in the
    lookup dict (or vice versa) depending on insertion order.
    """
    routing = _routing(
        _league("NBA", "NBA",  "GAME LINES", "3"),
        _league("NBA", "WNBA", "GAME LINES", "12139"),
    )
    markets = [{"league_pattern": "GAME LINES", "market": "spreads", "period": "fg"}]

    nba = scraper.resolve_leagues(routing, "NBA", "NBA", markets)
    wnba = scraper.resolve_leagues(routing, "NBA", "WNBA", markets)

    assert nba == [{"id": "3", "name": "GAME LINES",
                    "market": "spreads", "period": "fg"}]
    assert wnba == [{"id": "12139", "name": "GAME LINES",
                     "market": "spreads", "period": "fg"}]


# ---- Happy path ------------------------------------------------------------

def test_resolves_multiple_markets_in_order():
    routing = _routing(
        _league("MLB", "MLB", "GAME LINES",    "5"),
        _league("MLB", "MLB", "1ST 5 INNINGS", "6"),
        _league("MLB", "MLB", "2ND HALVES",    "503"),
    )
    markets = [
        {"league_pattern": "GAME LINES",    "market": "spreads",    "period": "fg"},
        {"league_pattern": "1ST 5 INNINGS", "market": "spreads_f5", "period": "F5"},
        {"league_pattern": "2ND HALVES",    "market": "spreads_h2", "period": "H2"},
    ]
    out = scraper.resolve_leagues(routing, "MLB", "MLB", markets)

    assert [lg["id"] for lg in out] == ["5", "6", "503"]
    assert [lg["period"] for lg in out] == ["fg", "F5", "H2"]


def test_leagueId_coerced_to_string():
    """BKM returns leagueId as a string already, but the API contract is
    'str' regardless. If a future response uses an int, we should still
    pass a string to fetch_schedule.
    """
    routing = _routing(_league("MLB", "MLB", "GAME LINES", 5))  # int, not str
    markets = [{"league_pattern": "GAME LINES", "market": "spreads", "period": "fg"}]
    out = scraper.resolve_leagues(routing, "MLB", "MLB", markets)
    assert out[0]["id"] == "5"


# ---- Missing-market path (WARN-and-skip) -----------------------------------

def test_missing_market_returns_empty_and_warns(capsys):
    routing = _routing(_league("MLB", "MLB", "GAME LINES", "5"))
    markets = [{"league_pattern": "DOES NOT EXIST",
                "market": "spreads", "period": "fg"}]
    out = scraper.resolve_leagues(routing, "MLB", "MLB", markets)
    assert out == []
    captured = capsys.readouterr().out
    assert "WARNING" in captured
    assert "DOES NOT EXIST" in captured
    # And the diagnostic should list what WAS available at that scope.
    assert "GAME LINES" in captured


def test_unknown_sport_id_warns_with_no_leagues_diagnostic(capsys):
    """When the (sportId, region) scope has no entries at all, the warning
    should say so explicitly — not just 'no match for pattern'."""
    routing = _routing(_league("MLB", "MLB", "GAME LINES", "5"))
    markets = [{"league_pattern": "GAME LINES",
                "market": "spreads", "period": "fg"}]
    out = scraper.resolve_leagues(routing, "NCB", "NCAAB", markets)
    assert out == []
    captured = capsys.readouterr().out
    assert "no BKM leagues under" in captured
    assert "NCB" in captured and "NCAAB" in captured


# ---- Defensive paths -------------------------------------------------------

def test_none_routing_info_returns_empty():
    assert scraper.resolve_leagues(None, "MLB", "MLB", [
        {"league_pattern": "GAME LINES", "market": "spreads", "period": "fg"}
    ]) == []


def test_empty_routing_info_returns_empty():
    routing = {"valid": True, "routedSports": []}
    assert scraper.resolve_leagues(routing, "MLB", "MLB", [
        {"league_pattern": "GAME LINES", "market": "spreads", "period": "fg"}
    ]) == []


def test_partial_match_resolves_what_it_can(capsys):
    """If only some markets match, return the matches and warn for the rest."""
    routing = _routing(
        _league("MLB", "MLB", "GAME LINES", "5"),
        _league("MLB", "MLB", "1ST 5 INNINGS", "6"),
    )
    markets = [
        {"league_pattern": "GAME LINES",    "market": "spreads",    "period": "fg"},
        {"league_pattern": "2ND HALVES",    "market": "spreads_h2", "period": "H2"},
        {"league_pattern": "1ST 5 INNINGS", "market": "spreads_f5", "period": "F5"},
    ]
    out = scraper.resolve_leagues(routing, "MLB", "MLB", markets)
    assert [lg["id"] for lg in out] == ["5", "6"]
    assert "2ND HALVES" in capsys.readouterr().out


def test_malformed_league_entry_without_leagueId_is_ignored():
    """A catalog entry missing leagueId would crash later in fetch_schedule.
    The resolver should drop it instead of indexing a phantom match."""
    routing = _routing(
        # Malformed first (no leagueId)
        {"sportId": "MLB", "region": "MLB", "leagueDescEn": "GAME LINES"},
        # Good entry with the same name — should be the one we resolve.
        _league("MLB", "MLB", "GAME LINES", "5"),
    )
    markets = [{"league_pattern": "GAME LINES",
                "market": "spreads", "period": "fg"}]
    out = scraper.resolve_leagues(routing, "MLB", "MLB", markets)
    assert out == [{"id": "5", "name": "GAME LINES",
                    "market": "spreads", "period": "fg"}]
