import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "wagerzon_odds"))

import pytest
from unittest.mock import MagicMock
import parlay_pricer


def test_get_combined_parlay_price_builds_sel_with_per_leg_idgm():
    """Cross-game parlay must encode each leg's own idgm in the sel string."""
    legs = [
        {"idgm": 100001, "play": 1, "points": "-1.5", "odds": 110},   # game A home spread
        {"idgm": 100001, "play": 2, "points": "9.5", "odds": -105},   # game A over
        {"idgm": 100002, "play": 0, "points": "+1.5", "odds": -120},  # game B away spread
        {"idgm": 100002, "play": 3, "points": "7.5", "odds": -110},   # game B under
    ]
    expected_sel_parts = [
        "1_100001_-1.5_110",
        "2_100001_9.5_-105",
        "0_100002_+1.5_-120",
        "3_100002_7.5_-110",
    ]

    fake_session = MagicMock()
    fake_response = MagicMock()
    fake_response.json.return_value = {
        "result": {
            "details": [{"Win": 21000, "Risk": 1000}],
        }
    }
    fake_response.raise_for_status.return_value = None
    fake_session.post.return_value = fake_response

    result = parlay_pricer.get_combined_parlay_price(fake_session, legs, amount=1000)

    call_args = fake_session.post.call_args
    sel_value = call_args.kwargs["data"]["sel"]
    for part in expected_sel_parts:
        assert part in sel_value, f"Missing leg in sel: {part}"
    assert result is not None
    assert result["win"] == 21000
    assert result["decimal"] == pytest.approx(22.0, rel=0.01)


def test_get_combined_parlay_price_returns_none_on_error():
    legs = [
        {"idgm": 100001, "play": 1, "points": "-1.5", "odds": 110},
        {"idgm": 100002, "play": 0, "points": "+1.5", "odds": -120},
    ]
    fake_session = MagicMock()
    fake_response = MagicMock()
    fake_response.json.return_value = {
        "result": {
            "details": [{"Win": 0, "Risk": 0}],
            "ErrorMsgKey": "MAXPARLAYRISKEXCEED",
        }
    }
    fake_response.raise_for_status.return_value = None
    fake_session.post.return_value = fake_response

    result = parlay_pricer.get_combined_parlay_price(fake_session, legs, amount=10000)
    assert result is None


import duckdb
import importlib.util


def _load_migration(path):
    spec = importlib.util.spec_from_file_location("migration", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_migration_adds_combo_columns(tmp_path):
    db = tmp_path / "test.duckdb"
    # Seed a minimal placed_parlays table matching the production schema
    con = duckdb.connect(str(db))
    con.execute("""
        CREATE TABLE placed_parlays (
            parlay_hash VARCHAR PRIMARY KEY,
            game_id VARCHAR NOT NULL,
            home_team VARCHAR NOT NULL,
            away_team VARCHAR NOT NULL,
            game_time TIMESTAMP,
            combo VARCHAR NOT NULL,
            spread_line FLOAT,
            total_line FLOAT,
            fair_odds INTEGER,
            wz_odds INTEGER NOT NULL,
            edge_pct FLOAT,
            kelly_bet FLOAT NOT NULL,
            actual_size FLOAT,
            placed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            status VARCHAR DEFAULT 'pending'
        )
    """)
    con.close()

    migration_path = (
        Path(__file__).resolve().parents[1]
        / "migrations" / "001_combined_parlay_columns.py"
    )
    mod = _load_migration(migration_path)
    mod.run(str(db))

    con = duckdb.connect(str(db))
    cols = {row[0]: row[1] for row in con.execute("DESCRIBE placed_parlays").fetchall()}
    con.close()

    assert "is_combo" in cols
    assert "combo_leg_ids" in cols
    assert "parent_combo_id" in cols


def test_migration_is_idempotent(tmp_path):
    db = tmp_path / "test.duckdb"
    con = duckdb.connect(str(db))
    con.execute("""
        CREATE TABLE placed_parlays (
            parlay_hash VARCHAR PRIMARY KEY,
            game_id VARCHAR NOT NULL,
            home_team VARCHAR NOT NULL,
            away_team VARCHAR NOT NULL,
            combo VARCHAR NOT NULL,
            wz_odds INTEGER NOT NULL,
            kelly_bet FLOAT NOT NULL
        )
    """)
    con.close()

    migration_path = (
        Path(__file__).resolve().parents[1]
        / "migrations" / "001_combined_parlay_columns.py"
    )
    mod = _load_migration(migration_path)
    mod.run(str(db))
    mod.run(str(db))  # second run must not error

    con = duckdb.connect(str(db))
    cols = {row[0]: row[1] for row in con.execute("DESCRIBE placed_parlays").fetchall()}
    con.close()

    assert "is_combo" in cols
    assert "combo_leg_ids" in cols
    assert "parent_combo_id" in cols
