"""Tests for Wagerzon scraper_v2 — 3-way market handling."""
import json
import sys
import tempfile
from pathlib import Path

import duckdb
import pytest

# pytest discovers this via conftest path; for now import directly
sys.path.insert(0, str(Path(__file__).parent.parent))


def test_init_database_adds_draw_ml_to_fresh_db(monkeypatch):
    """Fresh DB created with init_database() must have draw_ml column."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "wagerzon.duckdb"
        # Patch the module-level DB_PATH that init_database uses
        import scraper_v2
        monkeypatch.setattr(scraper_v2, "DB_PATH", db_path)

        scraper_v2.init_database("mlb")

        con = duckdb.connect(str(db_path))
        cols = [r[0] for r in con.execute("DESCRIBE mlb_odds").fetchall()]
        con.close()

        assert "draw_ml" in cols, f"draw_ml missing from fresh schema: {cols}"


def test_init_database_idempotent_ALTER_on_existing_db(monkeypatch):
    """Existing DB without draw_ml gets the column added; second call is no-op."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "wagerzon.duckdb"

        # Pre-create an old-schema mlb_odds (no draw_ml)
        con = duckdb.connect(str(db_path))
        con.execute("""
            CREATE TABLE mlb_odds (
                fetch_time TIMESTAMP, sport_key VARCHAR, game_id VARCHAR,
                game_date VARCHAR, game_time VARCHAR, away_team VARCHAR,
                home_team VARCHAR, market VARCHAR, period VARCHAR,
                away_spread FLOAT, away_spread_price INTEGER,
                home_spread FLOAT, home_spread_price INTEGER,
                total FLOAT, over_price INTEGER, under_price INTEGER,
                away_ml INTEGER, home_ml INTEGER, idgm INTEGER
            )
        """)
        con.close()

        import scraper_v2
        monkeypatch.setattr(scraper_v2, "DB_PATH", db_path)

        # First call adds draw_ml
        scraper_v2.init_database("mlb")
        con = duckdb.connect(str(db_path))
        cols_after_1 = [r[0] for r in con.execute("DESCRIBE mlb_odds").fetchall()]
        assert "draw_ml" in cols_after_1

        # Second call: idempotent, no error, no duplicate column
        scraper_v2.init_database("mlb")
        cols_after_2 = [r[0] for r in con.execute("DESCRIBE mlb_odds").fetchall()]
        con.close()
        assert cols_after_1 == cols_after_2, "Second init_database call changed schema"


def _make_3way_base():
    """Synthetic base dict matching what parse_odds builds before calling helpers."""
    return {
        "fetch_time": "2026-05-03 19:30:00",
        "sport_key": "baseball_mlb",
        "away_team": "Cleveland Guardians",
        "home_team": "Athletics",
        "game_date": "05/03",
        "game_time": "19:30",
        "idgm": 5635900,
    }


def test_parse_3way_line_emits_record_with_draw_ml():
    """A 3-way GameLine with all three prices populates the row correctly."""
    from scraper_v2 import parse_3way_line
    line = {
        "voddst": "105",   # away ML
        "hoddst": "130",   # home ML
        "vspoddst": "475", # draw price
        # ...other fields like vsprdt/vsprdoddst exist but are empty/irrelevant
    }
    rec = parse_3way_line(line, "test-game-1", "f5", "h2h_3way_1st_5_innings", _make_3way_base())

    assert rec is not None
    assert rec["market"] == "h2h_3way_1st_5_innings"
    assert rec["period"] == "f5"
    assert rec["away_ml"] == 105
    assert rec["home_ml"] == 130
    assert rec["draw_ml"] == 475
    # Spread/total fields must be NULL for 3-way rows
    assert rec["away_spread"] is None
    assert rec["home_spread"] is None
    assert rec["total"] is None
    assert rec["over_price"] is None
    assert rec["under_price"] is None
    # Base fields preserved
    assert rec["away_team"] == "Cleveland Guardians"
    assert rec["home_team"] == "Athletics"
    assert rec["game_id"] == "test-game-1"


def test_parse_3way_line_returns_none_when_all_three_prices_missing():
    """If voddst/hoddst/vspoddst are all empty, no row should be emitted."""
    from scraper_v2 import parse_3way_line
    line = {"voddst": "", "hoddst": "", "vspoddst": ""}
    rec = parse_3way_line(line, "test-game-2", "f5", "h2h_3way_1st_5_innings", _make_3way_base())
    assert rec is None


def test_parse_3way_line_handles_negative_prices():
    """Negative American odds (e.g. heavy favorite) parse correctly."""
    from scraper_v2 import parse_3way_line
    line = {"voddst": "-150", "hoddst": "+200", "vspoddst": "+500"}
    rec = parse_3way_line(line, "test-game-3", "f5", "h2h_3way_1st_5_innings", _make_3way_base())
    assert rec["away_ml"] == -150
    assert rec["home_ml"] == 200
    assert rec["draw_ml"] == 500


def test_parse_3way_line_emits_record_when_only_two_prices_present():
    """Partial price coverage (no draw posted yet) must still emit a row.

    Locks the `if all three are None` guard against a future refactor that
    accidentally flips AND to OR — which would silently drop rows where
    Wagerzon hasn't yet posted the third leg.
    """
    from scraper_v2 import parse_3way_line
    line = {"voddst": "105", "hoddst": "130", "vspoddst": ""}
    rec = parse_3way_line(line, "g", "f5", "h2h_3way_1st_5_innings", _make_3way_base())
    assert rec is not None
    assert rec["away_ml"] == 105
    assert rec["home_ml"] == 130
    assert rec["draw_ml"] is None
