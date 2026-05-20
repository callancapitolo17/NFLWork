"""Unit tests for resolve_leagues() — WZ's dynamic league discovery.

These tests pin down the resolver contract so future drift surfaces here
rather than as bad data in wagerzon.duckdb. Pure Python — no network.

The resolver matches declared markets by (IndexName, Description) against
WZ's ActiveLeaguesHelper catalog. The triple (IndexName, Description)
isn't strictly necessary today (WZ Description strings start with the
sport prefix like "MLB - ") but adding IndexName as a guard is cheap
and protects against future cross-sport collisions.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import scraper_v2  # noqa: E402


def _league(idx_name: str, desc: str, league_id, sport_id: str = "MLB",
            active: bool = True) -> dict:
    """Build a minimal catalog row matching WZ's ActiveLeaguesHelper shape."""
    return {
        "IdLeague": league_id,
        "Description": desc,
        "IdSport": sport_id,
        "IndexName": idx_name,
        "Active": active,
    }


# ---- Happy path ------------------------------------------------------------

def test_resolves_multiple_markets_in_order():
    catalog = [
        _league("MLB", "MLB - GAME LINES",         416),
        _league("MLB", "MLB - 1ST 5 FULL INNINGS", 417),
        _league("MLB", "MLB - 3 INNINGS LINE",     732),
    ]
    markets = [
        {"description": "MLB - GAME LINES",         "period": "fg", "kind": "lines"},
        {"description": "MLB - 1ST 5 FULL INNINGS", "period": "F5", "kind": "lines"},
        {"description": "MLB - 3 INNINGS LINE",     "period": "F3", "kind": "lines"},
    ]
    out = scraper_v2.resolve_leagues(catalog, "MLB", markets)
    assert [lg["id"] for lg in out] == ["416", "417", "732"]
    assert [lg["period"] for lg in out] == ["fg", "F5", "F3"]
    assert all(lg["kind"] == "lines" for lg in out)


def test_leagueId_coerced_to_string():
    """WZ returns IdLeague as int; downstream code uses string for the URL.
    The resolver should coerce so callers don't have to worry about it.
    """
    catalog = [_league("MLB", "MLB - GAME LINES", 416)]  # int
    markets = [{"description": "MLB - GAME LINES", "period": "fg", "kind": "lines"}]
    out = scraper_v2.resolve_leagues(catalog, "MLB", markets)
    assert out[0]["id"] == "416"


def test_index_name_disambiguates_cross_sport():
    """If a Description appeared under two IndexNames, resolving with one
    should not match the other. This isn't currently a real WZ scenario
    (each Description is prefixed with its sport) but the guard is cheap
    and protects against future WZ data shape changes.
    """
    catalog = [
        _league("MLB", "TEAM TOTALS", 418),
        _league("NBA", "TEAM TOTALS", 626),
    ]
    markets = [{"description": "TEAM TOTALS", "period": "fg", "kind": "team_total"}]
    nba = scraper_v2.resolve_leagues(catalog, "NBA", markets)
    mlb = scraper_v2.resolve_leagues(catalog, "MLB", markets)
    assert nba[0]["id"] == "626"
    assert mlb[0]["id"] == "418"


# ---- Missing-market WARN path ---------------------------------------------

def test_missing_market_returns_empty_and_warns(capsys):
    catalog = [_league("MLB", "MLB - GAME LINES", 416)]
    markets = [{"description": "MLB - DOES NOT EXIST",
                "period": "fg", "kind": "lines"}]
    out = scraper_v2.resolve_leagues(catalog, "MLB", markets)
    assert out == []
    log = capsys.readouterr().out
    assert "WARNING" in log
    assert "DOES NOT EXIST" in log
    # Should also list what WAS available under MLB.
    assert "MLB - GAME LINES" in log


def test_unknown_index_name_warns_with_off_season_hint(capsys):
    """When the catalog has no entries under the requested IndexName at all
    (e.g., off-season NFL), the warning should say so — not just 'pattern
    not found'.
    """
    catalog = [_league("MLB", "MLB - GAME LINES", 416)]
    markets = [{"description": "NFL - GAME LINES",
                "period": "fg", "kind": "lines"}]
    out = scraper_v2.resolve_leagues(catalog, "NFL", markets)
    assert out == []
    log = capsys.readouterr().out
    assert "no WZ leagues under IndexName='NFL'" in log
    assert "off-season" in log


def test_partial_match_resolves_what_it_can(capsys):
    """If some markets match and others don't, return the matches and warn
    individually for the rest. Off-season for a few markets shouldn't tank
    the whole sport.
    """
    catalog = [
        _league("MLB", "MLB - GAME LINES", 416),
        _league("MLB", "MLB - TEAM TOTALS", 418),
    ]
    markets = [
        {"description": "MLB - GAME LINES",       "period": "fg", "kind": "lines"},
        {"description": "MLB - 1ST 5 FULL INNINGS", "period": "F5", "kind": "lines"},
        {"description": "MLB - TEAM TOTALS",      "period": "fg", "kind": "team_total"},
    ]
    out = scraper_v2.resolve_leagues(catalog, "MLB", markets)
    assert [lg["id"] for lg in out] == ["416", "418"]
    log = capsys.readouterr().out
    assert "1ST 5 FULL INNINGS" in log
    assert "WARNING" in log


def test_long_pattern_list_truncated_in_diagnostic(capsys):
    """When the IndexName has more than 20 leagues, the hint should cap
    the listed patterns to avoid log spam (MLB has 25+ leagues today).
    """
    # 30 leagues under MLB so the cap (20) is exercised.
    catalog = [_league("MLB", f"MLB - FAKE {i}", 9000 + i) for i in range(30)]
    markets = [{"description": "MLB - DOES NOT EXIST",
                "period": "fg", "kind": "lines"}]
    scraper_v2.resolve_leagues(catalog, "MLB", markets)
    log = capsys.readouterr().out
    assert "first 20" in log


# ---- Defensive paths -------------------------------------------------------

def test_none_catalog_returns_empty(capsys):
    out = scraper_v2.resolve_leagues(None, "MLB", [
        {"description": "MLB - GAME LINES", "period": "fg", "kind": "lines"}
    ])
    assert out == []
    assert "ActiveLeaguesHelper unavailable" in capsys.readouterr().out


def test_empty_catalog_returns_empty():
    out = scraper_v2.resolve_leagues([], "MLB", [
        {"description": "MLB - GAME LINES", "period": "fg", "kind": "lines"}
    ])
    assert out == []


def test_malformed_catalog_row_without_idleague_ignored():
    """A row missing IdLeague would crash later in the URL build. The
    resolver should drop it instead of indexing a phantom match.
    """
    catalog = [
        # Malformed first (no IdLeague)
        {"IndexName": "MLB", "Description": "MLB - GAME LINES", "IdSport": "MLB"},
        # Good entry with the same desc — should be the one we resolve.
        _league("MLB", "MLB - GAME LINES", 416),
    ]
    markets = [{"description": "MLB - GAME LINES",
                "period": "fg", "kind": "lines"}]
    out = scraper_v2.resolve_leagues(catalog, "MLB", markets)
    assert out == [{"id": "416", "description": "MLB - GAME LINES",
                    "period": "fg", "kind": "lines"}]


def test_missing_period_defaults_to_fg():
    """Period is a per-market field; if a config entry forgets it, fall
    back to 'fg' rather than crashing. Same for kind."""
    catalog = [_league("MLB", "MLB - GAME LINES", 416)]
    markets = [{"description": "MLB - GAME LINES"}]  # no period, no kind
    out = scraper_v2.resolve_leagues(catalog, "MLB", markets)
    assert out[0]["period"] == "fg"
    assert out[0]["kind"] == "lines"
