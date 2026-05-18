"""Unit tests for _resolve_wagerzon_play_idgm.

The Wagerzon scraper writes spread + total + ML for a given (game, period)
combined on ONE row keyed by market='spreads' (FG) or market='spreads_<period>'.
The resolver has to route totals/h2h lookups to those rows and read the
total/over_price/under_price/home_ml/away_ml columns from there.

These tests build a tiny fixture wagerzon.duckdb that mirrors the real
schema and verify the resolver handles every market kind correctly,
including the negative case (a spreads-only row must NOT satisfy a
totals lookup).
"""
import sys
from pathlib import Path
from unittest.mock import patch

import duckdb
import pytest

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))            # MLB Dashboard/
sys.path.insert(0, str(HERE.parents[3]))            # NFLWork/

import mlb_dashboard_server as srv


HOME = "Washington Nationals"
AWAY = "New York Mets"


def _build_fixture_db(db_path: Path):
    """Create a wagerzon.duckdb mirroring the production mlb_odds schema."""
    db_path.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(db_path))
    con.execute("""
        CREATE TABLE mlb_odds (
            fetch_time TIMESTAMP,
            sport_key VARCHAR,
            game_id VARCHAR,
            game_date VARCHAR,
            game_time VARCHAR,
            away_team VARCHAR,
            home_team VARCHAR,
            market VARCHAR,
            period VARCHAR,
            away_spread FLOAT,
            away_spread_price INTEGER,
            home_spread FLOAT,
            home_spread_price INTEGER,
            total FLOAT,
            over_price INTEGER,
            under_price INTEGER,
            away_ml INTEGER,
            home_ml INTEGER,
            idgm INTEGER,
            draw_ml INTEGER
        )
    """)
    rows = [
        # FG: spread + total + ML combined on a spreads row.
        ("2026-05-18 17:14:01", "mlb", "g1", "2026-05-18", "15:46",
         AWAY, HOME, "spreads", "fg",
         -1.5, 101, 1.5, -121,                # away spread / home spread
         10.0, -110, -110,                    # total / over / under
         140, -160,                           # away ml / home ml
         500001, None),
        # F7: combined row carries total + over/under prices. NO spread.
        ("2026-05-18 17:14:01", "mlb", "g1", "2026-05-18", "15:46",
         AWAY, HOME, "spreads_f7", "f7",
         0.0, None, 0.0, None,
         7.5, -125, -115,
         None, None,
         500002, None),
        # F3: combined row carries total. NO spread, NO ML.
        ("2026-05-18 17:14:01", "mlb", "g1", "2026-05-18", "15:46",
         AWAY, HOME, "spreads_f3", "f3",
         0.0, None, 0.0, None,
         3.5, 115, -155,
         None, None,
         500003, None),
        # F5 / H1: spread + total + ML combined.
        ("2026-05-18 17:14:01", "mlb", "g1", "2026-05-18", "15:46",
         AWAY, HOME, "spreads_h1", "h1",
         -0.5, -135, 0.5, 115,
         5.5, -120, 100,
         130, -150,
         500004, None),
        # Alt totals — its own row.
        ("2026-05-18 17:14:01", "mlb", "g1", "2026-05-18", "15:46",
         AWAY, HOME, "alternate_totals_fg", "fg",
         None, None, None, None,
         11.5, 140, -180,
         None, None,
         500005, None),
        # Alt spreads — its own row, no total.
        ("2026-05-18 17:14:01", "mlb", "g1", "2026-05-18", "15:46",
         AWAY, HOME, "alternate_spreads_fg", "fg",
         2.5, -180, -2.5, 140,
         None, None, None,
         None, None,
         500006, None),
        # NEGATIVE FIXTURE: a spreads_f7 row that ONLY carries spread (no total).
        # Totals lookup MUST NOT return this row.
        ("2026-05-18 17:14:01", "mlb", "g2", "2026-05-18", "15:46",
         "Atlanta Braves", "Miami Marlins", "spreads_f7", "f7",
         -1.5, -110, 1.5, -110,
         None, None, None,                   # total NULL, prices NULL
         None, None,
         500099, None),
    ]
    con.executemany(
        "INSERT INTO mlb_odds VALUES "
        "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    con.close()


@pytest.fixture
def repo_root(tmp_path, monkeypatch):
    """Patch _REPO_ROOT so the resolver reads our fixture db."""
    _build_fixture_db(tmp_path / "wagerzon_odds" / "wagerzon.duckdb")
    monkeypatch.setattr(srv, "_REPO_ROOT", tmp_path)
    return tmp_path


def test_f7_total_over_resolves(repo_root):
    got = srv._resolve_wagerzon_play_idgm({
        "market": "totals_1st_7_innings", "bet_on": "Over",
        "home_team": HOME, "away_team": AWAY, "line": 7.5,
    })
    assert got == {"idgm": 500002, "play": 2}


def test_f7_total_under_resolves(repo_root):
    got = srv._resolve_wagerzon_play_idgm({
        "market": "totals_1st_7_innings", "bet_on": "Under",
        "home_team": HOME, "away_team": AWAY, "line": 7.5,
    })
    assert got == {"idgm": 500002, "play": 3}


def test_f3_total_resolves(repo_root):
    got = srv._resolve_wagerzon_play_idgm({
        "market": "totals_1st_3_innings", "bet_on": "Over",
        "home_team": HOME, "away_team": AWAY, "line": 3.5,
    })
    assert got == {"idgm": 500003, "play": 2}


def test_f5_total_resolves(repo_root):
    got = srv._resolve_wagerzon_play_idgm({
        "market": "totals_1st_5_innings", "bet_on": "Over",
        "home_team": HOME, "away_team": AWAY, "line": 5.5,
    })
    assert got == {"idgm": 500004, "play": 2}


def test_f5_ml_home_resolves(repo_root):
    got = srv._resolve_wagerzon_play_idgm({
        "market": "h2h_1st_5_innings", "bet_on": HOME,
        "home_team": HOME, "away_team": AWAY, "line": None,
    })
    assert got == {"idgm": 500004, "play": 5}


def test_f5_ml_away_resolves(repo_root):
    got = srv._resolve_wagerzon_play_idgm({
        "market": "h2h_1st_5_innings", "bet_on": AWAY,
        "home_team": HOME, "away_team": AWAY, "line": None,
    })
    assert got == {"idgm": 500004, "play": 4}


def test_fg_total_resolves(repo_root):
    got = srv._resolve_wagerzon_play_idgm({
        "market": "totals", "bet_on": "Over",
        "home_team": HOME, "away_team": AWAY, "line": 10.0,
    })
    assert got == {"idgm": 500001, "play": 2}


def test_fg_h2h_home_resolves(repo_root):
    got = srv._resolve_wagerzon_play_idgm({
        "market": "h2h", "bet_on": HOME,
        "home_team": HOME, "away_team": AWAY, "line": None,
    })
    assert got == {"idgm": 500001, "play": 5}


def test_fg_spread_home_unchanged(repo_root):
    got = srv._resolve_wagerzon_play_idgm({
        "market": "spreads", "bet_on": HOME,
        "home_team": HOME, "away_team": AWAY, "line": 1.5,
    })
    assert got == {"idgm": 500001, "play": 1}


def test_alt_total_over_unchanged(repo_root):
    got = srv._resolve_wagerzon_play_idgm({
        "market": "alternate_totals_fg", "bet_on": "Over",
        "home_team": HOME, "away_team": AWAY, "line": 11.5,
    })
    assert got == {"idgm": 500005, "play": 2}


def test_alt_spread_unchanged(repo_root):
    got = srv._resolve_wagerzon_play_idgm({
        "market": "alternate_spreads_fg", "bet_on": HOME,
        "home_team": HOME, "away_team": AWAY, "line": -2.5,
    })
    assert got == {"idgm": 500006, "play": 1}


def test_spreads_only_row_does_not_satisfy_total(repo_root):
    """The Atlanta@Miami spreads_f7 row has total IS NULL — a totals lookup
    must NOT return it just because the market/team match."""
    got = srv._resolve_wagerzon_play_idgm({
        "market": "totals_1st_7_innings", "bet_on": "Over",
        "home_team": "Miami Marlins", "away_team": "Atlanta Braves",
        "line": 6.5,
    })
    assert got is None


def test_unknown_market_returns_none(repo_root):
    got = srv._resolve_wagerzon_play_idgm({
        "market": "team_totals_home_fg", "bet_on": "Over",
        "home_team": HOME, "away_team": AWAY, "line": 4.5,
    })
    assert got is None


def test_missing_wagerzon_db_returns_none(tmp_path, monkeypatch):
    monkeypatch.setattr(srv, "_REPO_ROOT", tmp_path)  # no DB on disk
    got = srv._resolve_wagerzon_play_idgm({
        "market": "totals_1st_7_innings", "bet_on": "Over",
        "home_team": HOME, "away_team": AWAY, "line": 7.5,
    })
    assert got is None
