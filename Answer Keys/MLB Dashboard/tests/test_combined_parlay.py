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


def _seed_endpoint_dbs(tmp_path):
    """Set up isolated DBs and return paths."""
    test_dashboard_db = tmp_path / "mlb_dashboard.duckdb"
    test_mlb_db = tmp_path / "mlb.duckdb"

    con = duckdb.connect(str(test_dashboard_db))
    con.execute("""
        CREATE TABLE placed_parlays (
            parlay_hash VARCHAR PRIMARY KEY, game_id VARCHAR, home_team VARCHAR,
            away_team VARCHAR, combo VARCHAR, wz_odds INTEGER, kelly_bet FLOAT,
            is_combo BOOLEAN DEFAULT FALSE, combo_leg_ids VARCHAR,
            parent_combo_id INTEGER, status VARCHAR DEFAULT 'pending'
        )
    """)
    con.execute("""
        CREATE TABLE sizing_settings (param VARCHAR PRIMARY KEY, value FLOAT)
    """)
    con.execute("INSERT INTO sizing_settings VALUES ('parlay_bankroll', 1000), ('parlay_kelly_mult', 0.5)")
    con.close()

    con = duckdb.connect(str(test_mlb_db))
    con.execute("""
        CREATE TABLE mlb_parlay_opportunities (
            parlay_hash VARCHAR, fair_dec FLOAT, wz_dec FLOAT, idgm INTEGER,
            spread_line FLOAT, total_line FLOAT, spread_price INTEGER,
            total_price INTEGER, combo VARCHAR, game_id VARCHAR
        )
    """)
    con.execute("""
        INSERT INTO mlb_parlay_opportunities VALUES
        ('hash_a', 4.32, 4.55, 100001, -1.5, 9.5, 110, -105, 'Home -1.5 + Over 9.5', 'NYY@BOS_2026'),
        ('hash_b', 4.85, 5.10, 100002, -1.5, 7.5, 120, -110, 'Home -1.5 + Under 7.5', 'LAD@SD_2026')
    """)
    con.close()
    return test_dashboard_db, test_mlb_db


def test_price_combined_parlay_endpoint(monkeypatch, tmp_path):
    test_dashboard_db, test_mlb_db = _seed_endpoint_dbs(tmp_path)

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import mlb_dashboard_server as svr

    monkeypatch.setattr(svr, "DB_PATH", test_dashboard_db)
    monkeypatch.setattr(svr, "MLB_MM_DB", test_mlb_db)

    def fake_get_combined_parlay_price(session, legs, amount=10000):
        return {"win": 18000, "decimal": 19.0, "american": 1800, "amount": 1000}
    monkeypatch.setattr(svr, "wz_get_combined_parlay_price", fake_get_combined_parlay_price)
    monkeypatch.setattr(svr, "_get_wz_session", lambda: object())
    svr._COMBO_PRICE_CACHE.clear()

    client = svr.app.test_client()
    resp = client.post("/api/price-combined-parlay", json={
        "parlay_hash_a": "hash_a", "parlay_hash_b": "hash_b"
    })
    assert resp.status_code == 200, resp.get_data(as_text=True)
    body = resp.get_json()
    assert body["joint_fair_dec"] == pytest.approx(4.32 * 4.85, rel=0.001)
    assert body["wz_dec"] == 19.0
    assert "kelly_stake" in body
    assert "joint_edge" in body


def test_price_combined_parlay_caches_within_ttl(monkeypatch, tmp_path):
    """Re-pricing the same pair within TTL must NOT call WZ a second time."""
    test_dashboard_db, test_mlb_db = _seed_endpoint_dbs(tmp_path)

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import mlb_dashboard_server as svr
    monkeypatch.setattr(svr, "DB_PATH", test_dashboard_db)
    monkeypatch.setattr(svr, "MLB_MM_DB", test_mlb_db)

    call_count = {"n": 0}
    def fake_wz(session, legs, amount=10000):
        call_count["n"] += 1
        return {"win": 18000, "decimal": 19.0, "american": 1800, "amount": 1000}
    monkeypatch.setattr(svr, "wz_get_combined_parlay_price", fake_wz)
    monkeypatch.setattr(svr, "_get_wz_session", lambda: object())
    svr._COMBO_PRICE_CACHE.clear()

    client = svr.app.test_client()
    payload = {"parlay_hash_a": "hash_a", "parlay_hash_b": "hash_b"}
    r1 = client.post("/api/price-combined-parlay", json=payload)
    r2 = client.post("/api/price-combined-parlay", json=payload)
    assert r1.status_code == 200 and r2.status_code == 200
    assert call_count["n"] == 1, "WZ pricer should be called exactly once across two requests"


def test_price_combined_parlay_rejects_same_game(monkeypatch, tmp_path):
    test_dashboard_db, test_mlb_db = _seed_endpoint_dbs(tmp_path)

    # Make both rows have the same game_id
    con = duckdb.connect(str(test_mlb_db))
    con.execute("UPDATE mlb_parlay_opportunities SET game_id = 'SAME_GAME' WHERE parlay_hash IN ('hash_a', 'hash_b')")
    con.close()

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import mlb_dashboard_server as svr
    monkeypatch.setattr(svr, "DB_PATH", test_dashboard_db)
    monkeypatch.setattr(svr, "MLB_MM_DB", test_mlb_db)
    svr._COMBO_PRICE_CACHE.clear()

    client = svr.app.test_client()
    resp = client.post("/api/price-combined-parlay", json={
        "parlay_hash_a": "hash_a", "parlay_hash_b": "hash_b"
    })
    assert resp.status_code == 400
    assert "Same-game" in resp.get_json()["error"]


def test_price_combined_parlay_sends_correct_play_codes(monkeypatch, tmp_path):
    """Captures legs sent to WZ pricer; asserts play codes + points sign per the canonical encoding."""
    test_dashboard_db, test_mlb_db = _seed_endpoint_dbs(tmp_path)

    # Add an Away/Under combo and a Home/Over combo to ensure both branches are tested
    con = duckdb.connect(str(test_mlb_db))
    con.execute("DELETE FROM mlb_parlay_opportunities")
    con.execute("""
        INSERT INTO mlb_parlay_opportunities VALUES
        ('hash_home_over',  4.32, 4.55, 100001, -1.5, 9.5, 110, -105, 'Home -1.5 + Over 9.5',  'NYY@BOS_2026'),
        ('hash_away_under', 4.85, 5.10, 100002,  1.5, 7.5, 120, -110, 'Away +1.5 + Under 7.5', 'LAD@SD_2026')
    """)
    con.close()

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import mlb_dashboard_server as svr
    monkeypatch.setattr(svr, "DB_PATH", test_dashboard_db)
    monkeypatch.setattr(svr, "MLB_MM_DB", test_mlb_db)

    captured_legs = {}
    def capturing_fake(session, legs, amount=10000):
        captured_legs["legs"] = legs
        return {"win": 18000, "decimal": 19.0, "american": 1800, "amount": 1000}
    monkeypatch.setattr(svr, "wz_get_combined_parlay_price", capturing_fake)
    monkeypatch.setattr(svr, "_get_wz_session", lambda: object())
    svr._COMBO_PRICE_CACHE.clear()

    client = svr.app.test_client()
    resp = client.post("/api/price-combined-parlay", json={
        "parlay_hash_a": "hash_home_over", "parlay_hash_b": "hash_away_under"
    })
    assert resp.status_code == 200, resp.get_data(as_text=True)

    legs = captured_legs["legs"]
    assert len(legs) == 4

    # Leg 0: Home spread, leg A
    assert legs[0]["play"] == 1, "Home spread should be play=1"
    assert legs[0]["idgm"] == 100001
    assert legs[0]["points"] == "-1.5"

    # Leg 1: Over total, leg A — points must be NEGATIVE
    assert legs[1]["play"] == 2, "Over total should be play=2"
    assert legs[1]["idgm"] == 100001
    assert legs[1]["points"] == "-9.5", f"Over points must be negative, got {legs[1]['points']!r}"

    # Leg 2: Away spread, leg B
    assert legs[2]["play"] == 0, "Away spread should be play=0"
    assert legs[2]["idgm"] == 100002
    assert legs[2]["points"] == "1.5"

    # Leg 3: Under total, leg B — points must be POSITIVE
    assert legs[3]["play"] == 3, "Under total should be play=3"
    assert legs[3]["idgm"] == 100002
    assert legs[3]["points"] == "7.5", f"Under points must be positive, got {legs[3]['points']!r}"


def test_place_combined_parlay_creates_combo_row(monkeypatch, tmp_path):
    test_dashboard_db = tmp_path / "mlb_dashboard.duckdb"
    con = duckdb.connect(str(test_dashboard_db))
    con.execute("""
        CREATE TABLE placed_parlays (
            parlay_hash VARCHAR PRIMARY KEY, game_id VARCHAR, home_team VARCHAR,
            away_team VARCHAR, combo VARCHAR, wz_odds INTEGER, kelly_bet FLOAT,
            actual_size FLOAT, status VARCHAR DEFAULT 'pending',
            placed_at TIMESTAMP, is_combo BOOLEAN DEFAULT FALSE,
            combo_leg_ids VARCHAR, parent_combo_id INTEGER
        )
    """)
    con.close()

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import mlb_dashboard_server as svr
    monkeypatch.setattr(svr, "DB_PATH", test_dashboard_db)

    client = svr.app.test_client()
    resp = client.post("/api/place-combined-parlay", json={
        "combo_hash": "combo_hash_xyz",
        "parlay_hash_a": "hash_a",
        "parlay_hash_b": "hash_b",
        "wz_odds": 1810,
        "kelly_bet": 9.40,
        "actual_size": 9.40,
        "combo_label": "NYY @ BOS + LAD @ SD (4-leg)",
    })
    assert resp.status_code == 200, resp.get_data(as_text=True)

    con = duckdb.connect(str(test_dashboard_db))
    rows = con.execute("SELECT parlay_hash, is_combo, combo_leg_ids FROM placed_parlays").fetchall()
    con.close()
    assert len(rows) == 1
    assert rows[0][0] == "combo_hash_xyz"
    assert rows[0][1] is True
    import json as _json
    assert sorted(_json.loads(rows[0][2])) == ["hash_a", "hash_b"]


def test_place_combined_parlay_rejects_duplicate(monkeypatch, tmp_path):
    test_dashboard_db = tmp_path / "mlb_dashboard.duckdb"
    con = duckdb.connect(str(test_dashboard_db))
    con.execute("""
        CREATE TABLE placed_parlays (
            parlay_hash VARCHAR PRIMARY KEY, game_id VARCHAR, home_team VARCHAR,
            away_team VARCHAR, combo VARCHAR, wz_odds INTEGER, kelly_bet FLOAT,
            actual_size FLOAT, status VARCHAR DEFAULT 'pending',
            placed_at TIMESTAMP, is_combo BOOLEAN DEFAULT FALSE,
            combo_leg_ids VARCHAR, parent_combo_id INTEGER
        )
    """)
    con.close()

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import mlb_dashboard_server as svr
    monkeypatch.setattr(svr, "DB_PATH", test_dashboard_db)

    client = svr.app.test_client()
    payload = {
        "combo_hash": "combo_dup",
        "parlay_hash_a": "h_a", "parlay_hash_b": "h_b",
        "wz_odds": 1500, "kelly_bet": 10.0, "actual_size": 10.0,
        "combo_label": "test",
    }
    r1 = client.post("/api/place-combined-parlay", json=payload)
    r2 = client.post("/api/place-combined-parlay", json=payload)
    assert r1.status_code == 200
    assert r2.status_code == 409
    assert "already placed" in r2.get_json()["error"].lower()


def test_place_combined_parlay_missing_fields(monkeypatch, tmp_path):
    test_dashboard_db = tmp_path / "mlb_dashboard.duckdb"
    con = duckdb.connect(str(test_dashboard_db))
    con.execute("""
        CREATE TABLE placed_parlays (
            parlay_hash VARCHAR PRIMARY KEY, game_id VARCHAR, home_team VARCHAR,
            away_team VARCHAR, combo VARCHAR, wz_odds INTEGER, kelly_bet FLOAT,
            actual_size FLOAT, status VARCHAR DEFAULT 'pending',
            placed_at TIMESTAMP, is_combo BOOLEAN DEFAULT FALSE,
            combo_leg_ids VARCHAR, parent_combo_id INTEGER
        )
    """)
    con.close()

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import mlb_dashboard_server as svr
    monkeypatch.setattr(svr, "DB_PATH", test_dashboard_db)

    client = svr.app.test_client()
    resp = client.post("/api/place-combined-parlay", json={"combo_hash": "x"})
    assert resp.status_code == 400
    assert "Missing fields" in resp.get_json()["error"]
