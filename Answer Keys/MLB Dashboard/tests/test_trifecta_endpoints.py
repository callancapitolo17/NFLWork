"""Tests for /api/place-trifecta and /api/remove-trifecta endpoints."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import duckdb
import pytest

# Make the server module importable
DASHBOARD_DIR = Path(__file__).parent.parent
MLB_DIR       = DASHBOARD_DIR.parent
NFLWORK_DIR   = MLB_DIR.parent
sys.path.insert(0, str(DASHBOARD_DIR))


@pytest.fixture
def server_with_temp_dbs(tmp_path, monkeypatch):
    """Spin up the Flask server with isolated temp DBs."""
    dash_db = tmp_path / "mlb_dashboard.duckdb"
    mlb_db  = tmp_path / "mlb.duckdb"

    # Seed the trifecta opportunities table the server reads from
    con = duckdb.connect(str(mlb_db))
    con.execute("""
        CREATE TABLE mlb_trifecta_opportunities (
            trifecta_hash VARCHAR, game_id VARCHAR, game VARCHAR,
            game_time TIMESTAMP, target_team VARCHAR, prop_type VARCHAR,
            side VARCHAR, description VARCHAR, n_samples INTEGER,
            model_odds INTEGER, dk_odds INTEGER, fair_odds INTEGER,
            book_odds INTEGER, edge_pct DOUBLE, kelly_bet DOUBLE
        )
    """)
    con.execute("""
        INSERT INTO mlb_trifecta_opportunities VALUES
        ('hashA', 'g1', 'NYY @ BOS', NULL, 'Yankees', 'TRIPLE-PLAY', 'away',
         'YANKEES — SCR 1ST, 1H & GM', 500, 350, 380, 365, 500, 8.5, 5.0)
    """)
    con.close()

    # Patch DB paths BEFORE importing the server (which reads them at import time)
    monkeypatch.setenv("MLB_DASHBOARD_DB", str(dash_db))
    monkeypatch.setenv("MLB_DB_PATH", str(mlb_db))

    # Force re-import so the patched env is picked up
    if "mlb_dashboard_server" in sys.modules:
        del sys.modules["mlb_dashboard_server"]
    import mlb_dashboard_server as server
    monkeypatch.setattr(server, "DB_PATH", dash_db)
    monkeypatch.setattr(server, "MLB_DB", mlb_db)
    server.init_db()

    server.app.config["TESTING"] = True
    return server.app.test_client()


def test_place_trifecta_inserts_row(server_with_temp_dbs):
    client = server_with_temp_dbs
    r = client.post("/api/place-trifecta",
                    json={"trifecta_hash": "hashA", "actual_wager": 7.0})
    assert r.status_code == 200, r.get_data(as_text=True)
    body = json.loads(r.get_data(as_text=True))
    assert body.get("success") is True or body.get("ok") is True


def test_place_trifecta_idempotent_on_duplicate(server_with_temp_dbs):
    client = server_with_temp_dbs
    client.post("/api/place-trifecta",
                json={"trifecta_hash": "hashA", "actual_wager": 7.0})
    r = client.post("/api/place-trifecta",
                    json={"trifecta_hash": "hashA", "actual_wager": 7.0})
    # Second call must NOT crash; either returns ok:true (no-op) or 409.
    assert r.status_code in (200, 409)


def test_place_trifecta_unknown_hash_returns_404(server_with_temp_dbs):
    client = server_with_temp_dbs
    r = client.post("/api/place-trifecta",
                    json={"trifecta_hash": "no_such_hash", "actual_wager": 7.0})
    assert r.status_code == 404


def test_remove_trifecta_deletes_row(server_with_temp_dbs):
    client = server_with_temp_dbs
    client.post("/api/place-trifecta",
                json={"trifecta_hash": "hashA", "actual_wager": 7.0})
    r = client.post("/api/remove-trifecta", json={"trifecta_hash": "hashA"})
    assert r.status_code == 200
    body = json.loads(r.get_data(as_text=True))
    assert body.get("success") is True or body.get("ok") is True


def test_place_trifecta_missing_hash_returns_400(server_with_temp_dbs):
    client = server_with_temp_dbs
    r = client.post("/api/place-trifecta", json={"actual_wager": 7.0})
    assert r.status_code == 400
