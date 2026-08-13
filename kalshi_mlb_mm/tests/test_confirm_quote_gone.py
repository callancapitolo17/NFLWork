"""Phantom-open-quotes fix: expired-RFQ quotes 404 on GET and must close.

When an RFQ expires, Kalshi DELETES the maker's quote object —
GET /communications/quotes/{qid} returns HTTP 404
{'error': {'code': 'not_found'}} (live-verified 2026-08-12). Before the
fix the confirm tick only matched status in ("cancelled", "expired",
"executed"); with no quote object st=None matched nothing and the row
wedged as status='open' forever, eating MAX_OPEN_QUOTES and both
exposure caps until the bot silently stopped quoting.

Fail-safe contract under test: ONLY an explicit 404 with a parsed JSON
body closes the row. 5xx, non-dict bodies, and raised connection errors
must leave a live quote open. Scaffolding follows
test_confirm_singles_veto.py.
"""
import importlib
from datetime import datetime, timezone

import pytest


class MustNotTouchGW:
    """The gone-quote path must neither confirm nor cancel anything."""

    def confirm(self, qid):
        raise AssertionError("gone-quote handling must not confirm")

    def cancel(self, qid):
        raise AssertionError("gone-quote handling must not cancel")


def _setup(monkeypatch, tmp_path, db_name, api_response):
    """One open quote in live_quotes + a stubbed quote-GET response.
    `api_response` is either the (status, body, headers) tuple to return
    or an Exception instance to raise (models timeouts/connection drops,
    which auth_client.api does NOT catch)."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_common import auth_client
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
            "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
            "status, submitted_at, closed_at, leg_prices_json) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["q-gone", "r-gone", "COMBO-GONE", "game1", 0.28, 0.68,
             None, 0.30, 0.30, "open", now, None, None])

    def fake_api(method, path, *a, **k):
        if isinstance(api_response, Exception):
            raise api_response
        return api_response
    monkeypatch.setattr(auth_client, "api", fake_api)
    return main, db


def _row_state(db):
    with db.connect(read_only=True) as con:
        return con.execute(
            "SELECT status, closed_at FROM live_quotes "
            "WHERE quote_id='q-gone'").fetchone()


def test_404_closes_quote_as_expired(monkeypatch, tmp_path):
    """THE regression: Kalshi deleted the quote with its expired RFQ ->
    the row must close as 'expired' with closed_at set, not wedge open."""
    main, db = _setup(monkeypatch, tmp_path, "gone_404.duckdb",
                      (404, {"error": {"code": "not_found"}}, {}))
    main._confirm_tick(MustNotTouchGW(), dry_run=False)
    status, closed_at = _row_state(db)
    assert status == "expired"
    assert closed_at is not None


def test_500_leaves_quote_open(monkeypatch, tmp_path):
    """Transient server error must NOT close a live quote."""
    main, db = _setup(monkeypatch, tmp_path, "gone_500.duckdb",
                      (500, {"error": {"code": "internal"}}, {}))
    main._confirm_tick(MustNotTouchGW(), dry_run=False)
    status, closed_at = _row_state(db)
    assert status == "open"
    assert closed_at is None


def test_404_with_non_dict_body_leaves_quote_open(monkeypatch, tmp_path):
    """A 404 whose body didn't parse as JSON (proxy/gateway HTML page) is
    not the verified Kalshi gone-quote shape — fail safe, stay open."""
    main, db = _setup(monkeypatch, tmp_path, "gone_404_text.duckdb",
                      (404, "<html>not found</html>", {}))
    main._confirm_tick(MustNotTouchGW(), dry_run=False)
    status, _ = _row_state(db)
    assert status == "open"


def test_connection_error_leaves_quote_open(monkeypatch, tmp_path):
    """auth_client.api raises on timeouts/connection drops; the main loop
    catches it upstream. The row must still be open afterward."""
    main, db = _setup(monkeypatch, tmp_path, "gone_conn.duckdb",
                      ConnectionError("boom"))
    with pytest.raises(ConnectionError):
        main._confirm_tick(MustNotTouchGW(), dry_run=False)
    status, _ = _row_state(db)
    assert status == "open"


def test_200_with_unknown_status_leaves_quote_open(monkeypatch, tmp_path):
    """A readable quote in a non-terminal, non-accepted state (e.g. still
    'open' on Kalshi's side) is untouched."""
    main, db = _setup(monkeypatch, tmp_path, "gone_open.duckdb",
                      (200, {"quote": {"status": "open"}}, {}))
    main._confirm_tick(MustNotTouchGW(), dry_run=False)
    status, _ = _row_state(db)
    assert status == "open"
