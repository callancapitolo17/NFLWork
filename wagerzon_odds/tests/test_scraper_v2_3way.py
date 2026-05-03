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
