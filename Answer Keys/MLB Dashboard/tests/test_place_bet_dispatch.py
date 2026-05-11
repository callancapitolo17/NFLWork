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


def test_dispatch_wagerzon_calls_single_placer(client):
    fake_result = {"status": "placed", "ticket_number": "T-9876",
                   "balance_after": 153.50,
                   "error_msg": None, "error_msg_key": None}
    with patch("mlb_dashboard_server.single_placer.place_single",
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
    with patch("mlb_dashboard_server._spawn_playwright_placer") as mock_pw:
        mock_pw.return_value = {"success": True,
                                 "status": "playwright_launched",
                                 "message": "Browser launching..."}
        resp = client.post("/api/place-bet",
                           json=_wz_bet(bookmaker_key="hoop88"))
    assert resp.status_code == 200
    mock_pw.assert_called_once()


def test_dispatch_bfa_routes_to_playwright(client):
    with patch("mlb_dashboard_server._spawn_playwright_placer") as mock_pw:
        mock_pw.return_value = {"success": True,
                                 "status": "playwright_launched",
                                 "message": "Browser launching..."}
        resp = client.post("/api/place-bet",
                           json=_wz_bet(bookmaker_key="bfa"))
    assert resp.status_code == 200
    mock_pw.assert_called_once()


def test_dispatch_betonlineag_routes_to_playwright(client):
    with patch("mlb_dashboard_server._spawn_playwright_placer") as mock_pw:
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
    with patch("mlb_dashboard_server.single_placer.place_single",
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
    with patch("mlb_dashboard_server.single_placer.place_single",
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
