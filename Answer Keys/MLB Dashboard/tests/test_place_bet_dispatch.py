"""Tests for /api/place-bet dispatcher.

Three branches:
  - wagerzon -> call single_placer.place_single, return its response
  - hoop88 / bfa / betonlineag -> spawn Playwright via internal helper
  - everything else -> 400 with manual-log message
"""
import json
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))            # MLB Dashboard/
sys.path.insert(0, str(HERE.parents[3]))            # NFLWork/

import mlb_dashboard_server as srv


@pytest.fixture
def client():
    srv.app.config["TESTING"] = True
    with srv.app.test_client() as c:
        yield c


def _wz_bet(**overrides):
    payload = dict(
        bet_hash="abc123",
        bookmaker_key="wagerzon",
        account="primary",
        kelly_bet=42, actual_size=42,
        idgm=100001, play=1, line=-1.5,
        american_odds=110, wz_odds_at_place=110,
        game_id="g1", market="spreads_1st_5_innings",
        bet_on="Home", home_team="NYY", away_team="BOS",
    )
    payload.update(overrides)
    return payload


def _mock_no_existing_bet():
    """Return a context manager that patches duckdb.connect so the 409 duplicate
    guard sees no existing row (fetchone() -> None).  Use in WZ-path tests that
    patch _insert_placement_breadcrumb but don't want to hit a real DB file."""
    mock_con = MagicMock()
    mock_con.execute.return_value.fetchone.return_value = None
    mock_con.__enter__ = lambda s: s
    mock_con.__exit__ = MagicMock(return_value=False)
    return patch("mlb_dashboard_server.duckdb.connect", return_value=mock_con)


def test_dispatch_wagerzon_calls_single_placer(client):
    fake_result = {"status": "placed", "ticket_number": "T-9876",
                   "balance_after": 153.50,
                   "error_msg": None, "error_msg_key": None}
    with _mock_no_existing_bet(), \
         patch("mlb_dashboard_server.single_placer.place_single",
               return_value=fake_result) as mock_place, \
         patch("mlb_dashboard_server._insert_placement_breadcrumb"), \
         patch("mlb_dashboard_server._finalize_placement"):
        resp = client.post("/api/place-bet", json=_wz_bet())
    assert resp.status_code == 200
    body = resp.get_json()
    assert body["status"] == "placed"
    assert body["ticket_number"] == "T-9876"
    mock_place.assert_called_once()


def test_dispatch_hoop88_routes_to_playwright(client):
    with patch("mlb_dashboard_server._spawn_playwright_placer") as mock_pw, \
         patch("mlb_dashboard_server._insert_placement_breadcrumb"):
        mock_pw.return_value = {"success": True,
                                 "status": "playwright_launched",
                                 "message": "Browser launching..."}
        resp = client.post("/api/place-bet",
                           json=_wz_bet(bookmaker_key="hoop88"))
    assert resp.status_code == 200
    mock_pw.assert_called_once()


def test_dispatch_bfa_routes_to_playwright(client):
    with patch("mlb_dashboard_server._spawn_playwright_placer") as mock_pw, \
         patch("mlb_dashboard_server._insert_placement_breadcrumb"):
        mock_pw.return_value = {"success": True,
                                 "status": "playwright_launched",
                                 "message": "Browser launching..."}
        resp = client.post("/api/place-bet",
                           json=_wz_bet(bookmaker_key="bfa"))
    assert resp.status_code == 200
    mock_pw.assert_called_once()


def test_dispatch_betonlineag_routes_to_playwright(client):
    with patch("mlb_dashboard_server._spawn_playwright_placer") as mock_pw, \
         patch("mlb_dashboard_server._insert_placement_breadcrumb"):
        mock_pw.return_value = {"success": True,
                                 "status": "playwright_launched",
                                 "message": "Browser launching..."}
        resp = client.post("/api/place-bet",
                           json=_wz_bet(bookmaker_key="betonlineag"))
    assert resp.status_code == 200
    mock_pw.assert_called_once()


def test_dispatch_unsupported_book_returns_400(client):
    resp = client.post("/api/place-bet",
                       json=_wz_bet(bookmaker_key="draftkings"))
    assert resp.status_code == 400
    body = resp.get_json()
    assert "manual log" in (body.get("error") or "").lower() or \
           "draftkings" in (body.get("error") or "").lower()


def test_wagerzon_dispatch_requires_account(client):
    resp = client.post("/api/place-bet",
                       json=_wz_bet(account=None))
    assert resp.status_code == 400
    body = resp.get_json()
    assert "account" in (body.get("error") or "").lower()


def test_wagerzon_dispatch_writes_breadcrumb_before_call(client):
    fake_result = {"status": "placed", "ticket_number": "T-1",
                   "balance_after": 100.0,
                   "error_msg": None, "error_msg_key": None}
    with _mock_no_existing_bet(), \
         patch("mlb_dashboard_server.single_placer.place_single",
               return_value=fake_result), \
         patch("mlb_dashboard_server._insert_placement_breadcrumb") as mock_bc, \
         patch("mlb_dashboard_server._finalize_placement"):
        client.post("/api/place-bet", json=_wz_bet())
    mock_bc.assert_called_once()
    # breadcrumb's status arg must be "placing"
    call_args = mock_bc.call_args
    # Accept either kwarg style or positional
    status_arg = call_args.kwargs.get("status")
    if status_arg is None and len(call_args.args) >= 4:
        # positional: (bet_hash, account, bet_meta, status="placing")
        status_arg = call_args.args[3] if len(call_args.args) > 3 else None
    assert status_arg == "placing" or "placing" in str(call_args)


def test_wagerzon_failure_finalizes_with_error_status(client):
    """If single_placer returns price_moved, finalize must record that status."""
    fake_result = {"status": "price_moved", "ticket_number": None,
                   "balance_after": None,
                   "error_msg": "Odds drifted", "error_msg_key": "drift"}
    with _mock_no_existing_bet(), \
         patch("mlb_dashboard_server.single_placer.place_single",
               return_value=fake_result), \
         patch("mlb_dashboard_server._insert_placement_breadcrumb"), \
         patch("mlb_dashboard_server._finalize_placement") as mock_fin:
        resp = client.post("/api/place-bet", json=_wz_bet())
    assert resp.status_code == 200
    assert resp.get_json()["status"] == "price_moved"
    mock_fin.assert_called_once()


def test_log_bet_records_without_book_call(client):
    """POST /api/log-bet inserts a placed row without contacting WZ or Playwright."""
    with patch("mlb_dashboard_server._insert_placement_breadcrumb") as mock_bc, \
         patch("mlb_dashboard_server.single_placer.place_single") as mock_wz, \
         patch("mlb_dashboard_server._spawn_playwright_placer") as mock_pw:
        resp = client.post("/api/log-bet",
                           json=_wz_bet(bookmaker_key="draftkings"))
    assert resp.status_code == 200
    body = resp.get_json()
    assert body.get("status") == "placed" or body.get("success") is True
    mock_wz.assert_not_called()
    mock_pw.assert_not_called()
    mock_bc.assert_called_once()


def test_log_bet_works_for_unsupported_book(client):
    """log-bet doesn't care about bookmaker_key — it just records."""
    with patch("mlb_dashboard_server._insert_placement_breadcrumb"):
        resp = client.post("/api/log-bet",
                           json=_wz_bet(bookmaker_key="pinnacle"))
    assert resp.status_code == 200


def test_missing_bet_hash_returns_400(client):
    resp = client.post("/api/place-bet",
                       json={"bookmaker_key": "wagerzon", "account": "primary"})
    assert resp.status_code == 400
    body = resp.get_json()
    assert "bet_hash" in (body.get("error") or "").lower()


def test_log_bet_missing_bet_hash_returns_400(client):
    resp = client.post("/api/log-bet",
                       json={"bookmaker_key": "wagerzon"})
    assert resp.status_code == 400


def test_dispatch_wagerzon_returns_409_on_in_flight_duplicate(client):
    """If placed_bets already has a 'placing' row for this bet_hash,
    /api/place-bet returns 409 instead of re-calling single_placer."""
    with patch("mlb_dashboard_server.duckdb.connect") as mock_connect, \
         patch("mlb_dashboard_server.single_placer.place_single") as mock_place:
        mock_con = MagicMock()
        mock_con.execute.return_value.fetchone.return_value = ("placing",)
        mock_connect.return_value = mock_con
        resp = client.post("/api/place-bet", json=_wz_bet())
    assert resp.status_code == 409
    mock_place.assert_not_called()


def test_integration_wagerzon_breadcrumb_writes_to_real_db(tmp_path, monkeypatch):
    """Real DuckDB insertion test — catches schema mismatches the unit tests
    miss. Builds a fresh placed_bets table with the production NOT NULL
    columns (model_prob, model_ev) and verifies the dispatcher actually
    INSERTs without violating constraints."""
    import duckdb
    import mlb_dashboard_server as srv

    # Create a temporary DB mirroring the production placed_bets schema.
    db_path = tmp_path / "test_dashboard.duckdb"
    con = duckdb.connect(str(db_path))
    con.execute("""
        CREATE TABLE placed_bets (
            bet_hash VARCHAR PRIMARY KEY,
            game_id VARCHAR,
            home_team VARCHAR,
            away_team VARCHAR,
            game_time TIMESTAMP,
            market VARCHAR,
            bet_on VARCHAR,
            line FLOAT,
            model_prob FLOAT NOT NULL,
            model_ev FLOAT NOT NULL,
            recommended_size FLOAT,
            actual_size FLOAT,
            odds INTEGER,
            bookmaker VARCHAR,
            placed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            status VARCHAR DEFAULT 'pending',
            ticket_number VARCHAR,
            error_msg VARCHAR,
            error_msg_key VARCHAR,
            wz_odds_at_place INTEGER,
            account VARCHAR
        )
    """)
    con.close()

    monkeypatch.setattr(srv, "DB_PATH", db_path)
    fake_result = {"status": "placed", "ticket_number": "T-1",
                   "balance_after": 100.0,
                   "error_msg": None, "error_msg_key": None}
    monkeypatch.setattr(srv.single_placer, "place_single",
                        lambda **kwargs: fake_result)
    # Bypass the WZ scraper lookup — supply idgm/play directly.
    monkeypatch.setattr(srv, "_resolve_wagerzon_play_idgm",
                        lambda bet: {"idgm": 100001, "play": 1})

    srv.app.config["TESTING"] = True
    with srv.app.test_client() as client:
        bet = dict(
            bet_hash="real_insert_test",
            bookmaker_key="wagerzon",
            account="primary",
            kelly_bet=42, actual_size=42,
            line=-1.5, american_odds=110, wz_odds_at_place=110,
            game_id="g1", market="spreads_1st_5_innings",
            bet_on="Home", home_team="NYY", away_team="BOS",
            model_prob=0.57, model_ev=0.06,
        )
        resp = client.post("/api/place-bet", json=bet)

    assert resp.status_code == 200, resp.get_json()

    con = duckdb.connect(str(db_path), read_only=True)
    row = con.execute(
        "SELECT model_prob, model_ev, status, ticket_number "
        "FROM placed_bets WHERE bet_hash = ?",
        ["real_insert_test"]
    ).fetchone()
    con.close()
    assert row is not None
    assert abs(row[0] - 0.57) < 1e-6
    assert abs(row[1] - 0.06) < 1e-6
    assert row[2] == "placed"
    assert row[3] == "T-1"
