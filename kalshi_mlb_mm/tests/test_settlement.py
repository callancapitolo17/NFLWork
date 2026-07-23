"""Regression tests for the settlement sweep (issue #12, audit M-1).

These fail on pre-#12 code by construction: nothing populated
fills.realized_pnl and there was no `settlements` table. All Kalshi API
traffic is stubbed via auth_client.api.
"""
import importlib
import json
import logging
from datetime import datetime, timezone

import pytest


# ---------------------------------------------------------------------------
# Harness
# ---------------------------------------------------------------------------

def _setup_db(monkeypatch, tmp_path, name):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / f"{name}.duckdb")
    importlib.reload(db)
    db.init_database()
    return db


def _insert_fill(db, fill_id, ticker, side_held, contracts, price, fee,
                 reconciled=True, realized_pnl=None):
    with db.connect() as con:
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
            "game_id, side_held, contracts, price, fee, model_fair_at_quote, "
            "book_fair_at_quote, blended_fair_at_quote, fair_at_confirm, "
            "realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            [fill_id, f"q-{fill_id}", f"r-{fill_id}", ticker, "g1", side_held,
             contracts, price, fee, None, 0.55, 0.55, 0.56, realized_pnl,
             datetime.now(timezone.utc), reconciled])


def _get_pnl(db, fill_id):
    with db.connect(read_only=True) as con:
        return con.execute("SELECT realized_pnl FROM fills WHERE fill_id=?",
                           [fill_id]).fetchone()[0]


def _settlement_rows(db):
    with db.connect(read_only=True) as con:
        return con.execute(
            "SELECT combo_market_ticker, result, raw_payload FROM settlements "
            "ORDER BY combo_market_ticker").fetchall()


class ApiStub:
    """Programmable auth_client.api stand-in. Records every call."""

    def __init__(self, market_result="yes", market_status="settled",
                 market_http=200, portfolio_fills_body=None, raise_exc=None):
        self.market_result = market_result
        self.market_status = market_status
        self.market_http = market_http
        self.portfolio_fills_body = portfolio_fills_body or {"fills": []}
        self.raise_exc = raise_exc
        self.calls = []

    def __call__(self, method, path, body=None, timeout=30):
        self.calls.append((method, path))
        if self.raise_exc is not None:
            raise self.raise_exc
        if path.startswith("/markets/"):
            if self.market_http != 200:
                return self.market_http, {"error": "boom"}, {}
            ticker = path.rsplit("/", 1)[-1]
            return 200, {"market": {"ticker": ticker,
                                    "status": self.market_status,
                                    "result": self.market_result,
                                    "settled_time": "2026-07-22T03:00:00Z"}}, {}
        if path.startswith("/portfolio/fills"):
            return 200, self.portfolio_fills_body, {}
        if path.startswith("/portfolio/settlements"):
            return 200, {"settlements": []}, {}
        raise AssertionError(f"unexpected path {path}")


@pytest.fixture()
def emitted(monkeypatch):
    from kalshi_mlb_mm import settlement
    events = []
    monkeypatch.setattr(settlement.research, "emit",
                        lambda et, **kw: events.append((et, kw)))
    return events


def _patch_api(monkeypatch, stub):
    from kalshi_mlb_mm import settlement
    monkeypatch.setattr(settlement.auth_client, "api", stub)


# ---------------------------------------------------------------------------
# Signed P&L math — hand-computed expected values (acceptance criteria)
# ---------------------------------------------------------------------------

def test_yes_win_pnl(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "yes_win")
    # Hold YES 10 @ $0.55, fee $0.01/ct, market resolves YES:
    # pnl = 10 * ((1 - 0.55) - 0.01) = 10 * 0.44 = 4.40
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01)
    _patch_api(monkeypatch, ApiStub(market_result="yes"))
    settlement.settlement_sweep_tick()
    assert _get_pnl(db, "f1") == pytest.approx(4.40)
    # First-ever settlement also triggers the one-time fee verification.
    assert [e for e, _ in emitted] == ["settlement_recorded", "fee_verification"]


def test_yes_loss_pnl(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "yes_loss")
    # Hold YES 10 @ $0.55, fee $0.01/ct, market resolves NO:
    # pnl = 10 * (-0.55 - 0.01) = -5.60
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01)
    _patch_api(monkeypatch, ApiStub(market_result="no"))
    settlement.settlement_sweep_tick()
    assert _get_pnl(db, "f1") == pytest.approx(-5.60)


def test_no_win_pnl(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "no_win")
    # Hold NO 5 @ $0.48, fee $0.01/ct, market resolves NO:
    # pnl = 5 * ((1 - 0.48) - 0.01) = 5 * 0.51 = 2.55
    _insert_fill(db, "f1", "TKR-A", "no", 5, 0.48, 0.01)
    _patch_api(monkeypatch, ApiStub(market_result="no"))
    settlement.settlement_sweep_tick()
    assert _get_pnl(db, "f1") == pytest.approx(2.55)


def test_no_loss_pnl(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "no_loss")
    # Hold NO 5 @ $0.48, fee $0.01/ct, market resolves YES:
    # pnl = 5 * (-0.48 - 0.01) = -2.45
    _insert_fill(db, "f1", "TKR-A", "no", 5, 0.48, 0.01)
    _patch_api(monkeypatch, ApiStub(market_result="yes"))
    settlement.settlement_sweep_tick()
    assert _get_pnl(db, "f1") == pytest.approx(-2.45)


def test_multiple_fills_same_ticker_all_settle(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "multi")
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01)   # +4.40
    _insert_fill(db, "f2", "TKR-A", "no", 5, 0.48, 0.01)     # loses: -2.45
    _patch_api(monkeypatch, ApiStub(market_result="yes"))
    settlement.settlement_sweep_tick()
    assert _get_pnl(db, "f1") == pytest.approx(4.40)
    assert _get_pnl(db, "f2") == pytest.approx(-2.45)
    assert len([e for e, _ in emitted if e == "settlement_recorded"]) == 2
    assert len(_settlement_rows(db)) == 1


# ---------------------------------------------------------------------------
# Unsettled / failure tolerance / retry
# ---------------------------------------------------------------------------

def test_unsettled_market_untouched(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "unsettled")
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01)
    _patch_api(monkeypatch, ApiStub(market_result="", market_status="active"))
    settlement.settlement_sweep_tick()
    assert _get_pnl(db, "f1") is None
    assert _settlement_rows(db) == []
    assert emitted == []


def test_terminal_nonbinary_result_untouched_and_logged(monkeypatch, tmp_path,
                                                        emitted, caplog):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "voided")
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01)
    _patch_api(monkeypatch, ApiStub(market_result="void", market_status="settled"))
    with caplog.at_level(logging.WARNING, logger="kalshi_mlb_mm.settlement"):
        settlement.settlement_sweep_tick()
    assert _get_pnl(db, "f1") is None
    assert "settlement_unrecognized_result" in caplog.text


def test_api_http_error_retried_next_sweep(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "http_err")
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01)
    _patch_api(monkeypatch, ApiStub(market_http=500))
    settlement.settlement_sweep_tick()          # API down — no crash
    assert _get_pnl(db, "f1") is None
    _patch_api(monkeypatch, ApiStub(market_result="yes"))
    settlement.settlement_sweep_tick()          # next sweep succeeds
    assert _get_pnl(db, "f1") == pytest.approx(4.40)


def test_api_exception_never_raises(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "api_exc")
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01)
    _patch_api(monkeypatch, ApiStub(raise_exc=RuntimeError("connection reset")))
    settlement.settlement_sweep_tick()          # must not propagate
    assert _get_pnl(db, "f1") is None


def test_no_unsettled_fills_no_api_calls(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import settlement
    _setup_db(monkeypatch, tmp_path, "empty")
    stub = ApiStub()
    _patch_api(monkeypatch, stub)
    settlement.settlement_sweep_tick()
    assert stub.calls == []


# ---------------------------------------------------------------------------
# Idempotency + reconciled gate
# ---------------------------------------------------------------------------

def test_idempotent_rerun(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "idem")
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01)
    _patch_api(monkeypatch, ApiStub(market_result="yes"))
    settlement.settlement_sweep_tick()
    settlement.settlement_sweep_tick()
    assert _get_pnl(db, "f1") == pytest.approx(4.40)
    assert len(_settlement_rows(db)) == 1
    # settlement_recorded emitted once — second sweep found nothing to settle.
    assert len([e for e, _ in emitted if e == "settlement_recorded"]) == 1


def test_unreconciled_fill_excluded_until_reconciled(monkeypatch, tmp_path,
                                                     emitted):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "unrec")
    # Reconcile may still rewrite side/size — settling early would bake in a
    # potentially wrong P&L forever.
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01, reconciled=False)
    stub = ApiStub(market_result="yes")
    _patch_api(monkeypatch, stub)
    settlement.settlement_sweep_tick()
    assert _get_pnl(db, "f1") is None
    assert stub.calls == []          # not even polled
    with db.connect() as con:
        con.execute("UPDATE fills SET reconciled=TRUE WHERE fill_id='f1'")
    settlement.settlement_sweep_tick()
    assert _get_pnl(db, "f1") == pytest.approx(4.40)


def test_mixed_reconciled_state_settles_only_reconciled(monkeypatch, tmp_path,
                                                        emitted):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "mixed")
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01, reconciled=True)
    _insert_fill(db, "f2", "TKR-A", "yes", 3, 0.60, 0.01, reconciled=False)
    _patch_api(monkeypatch, ApiStub(market_result="yes"))
    settlement.settlement_sweep_tick()
    assert _get_pnl(db, "f1") == pytest.approx(4.40)
    assert _get_pnl(db, "f2") is None    # waits for reconcile, next sweep
    settlement.settlement_sweep_tick()   # still unreconciled — still waits
    assert _get_pnl(db, "f2") is None
    with db.connect() as con:
        con.execute("UPDATE fills SET reconciled=TRUE WHERE fill_id='f2'")
    settlement.settlement_sweep_tick()
    # 3 * ((1 - 0.60) - 0.01) = 1.17
    assert _get_pnl(db, "f2") == pytest.approx(1.17)


# ---------------------------------------------------------------------------
# settlements audit table
# ---------------------------------------------------------------------------

def test_settlements_row_written_with_raw_payload(monkeypatch, tmp_path,
                                                  emitted):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "audit")
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01)
    _patch_api(monkeypatch, ApiStub(market_result="yes"))
    settlement.settlement_sweep_tick()
    rows = _settlement_rows(db)
    assert len(rows) == 1
    ticker, result, raw = rows[0]
    assert ticker == "TKR-A"
    assert result == "yes"
    payload = json.loads(raw)
    assert payload["status"] == "settled"     # raw Kalshi payload round-trips
    assert payload["result"] == "yes"


# ---------------------------------------------------------------------------
# Fee verification (issue item 3)
# ---------------------------------------------------------------------------

def test_fee_mismatch_logs_loudly(monkeypatch, tmp_path, emitted, caplog):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "feemis")
    # Our records: fee $0.01/ct * 10 = $0.10 expected. Kalshi says $0.50.
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01)
    stub = ApiStub(market_result="yes",
                   portfolio_fills_body={"fills": [{"ticker": "TKR-A",
                                                    "fee_cost": "0.50"}]})
    _patch_api(monkeypatch, stub)
    with caplog.at_level(logging.WARNING, logger="kalshi_mlb_mm.settlement"):
        settlement.settlement_sweep_tick()
    assert "FEE_MISMATCH" in caplog.text
    fee_events = [kw for e, kw in emitted if e == "fee_verification"]
    assert len(fee_events) == 1
    assert fee_events[0]["payload"]["verdict"] == "MISMATCH"


def test_fee_match_no_warning(monkeypatch, tmp_path, emitted, caplog):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "feeok")
    # Expected $0.10; Kalshi reports exactly $0.10.
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01)
    _patch_api(monkeypatch, ApiStub(
        market_result="yes",
        portfolio_fills_body={"fills": [{"ticker": "TKR-A", "fee_cost": "0.10"}]}))
    with caplog.at_level(logging.WARNING, logger="kalshi_mlb_mm.settlement"):
        settlement.settlement_sweep_tick()
    assert "FEE_MISMATCH" not in caplog.text
    fee_events = [kw for e, kw in emitted if e == "fee_verification"]
    assert fee_events[0]["payload"]["verdict"] == "matched"


def test_fee_verification_runs_once_ever(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "feeonce")
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01)
    _patch_api(monkeypatch, ApiStub(market_result="yes"))
    settlement.settlement_sweep_tick()
    # Second settled ticker, later sweep: no second verification.
    _insert_fill(db, "f2", "TKR-B", "yes", 4, 0.50, 0.01)
    _patch_api(monkeypatch, ApiStub(market_result="yes"))
    settlement.settlement_sweep_tick()
    fee_events = [kw for e, kw in emitted if e == "fee_verification"]
    assert len(fee_events) == 1


def test_fee_field_absent_logs_unverified(monkeypatch, tmp_path, emitted,
                                          caplog):
    from kalshi_mlb_mm import settlement
    db = _setup_db(monkeypatch, tmp_path, "feenone")
    _insert_fill(db, "f1", "TKR-A", "yes", 10, 0.55, 0.01)
    _patch_api(monkeypatch, ApiStub(
        market_result="yes",
        portfolio_fills_body={"fills": [{"ticker": "TKR-A", "count": 10}]}))
    with caplog.at_level(logging.WARNING, logger="kalshi_mlb_mm.settlement"):
        settlement.settlement_sweep_tick()
    assert "fee_unverified" in caplog.text
    # P&L still lands — fee verification is advisory, never blocking.
    assert _get_pnl(db, "f1") == pytest.approx(4.40)


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def test_realized_pnl_usd_formula():
    from kalshi_mlb_mm.settlement import realized_pnl_usd
    assert realized_pnl_usd("yes", 10, 0.55, 0.01, "yes") == pytest.approx(4.40)
    assert realized_pnl_usd("yes", 10, 0.55, 0.01, "no") == pytest.approx(-5.60)
    assert realized_pnl_usd("no", 5, 0.48, 0.01, "no") == pytest.approx(2.55)
    assert realized_pnl_usd("no", 5, 0.48, 0.01, "yes") == pytest.approx(-2.45)
    # Phantom fill (contracts=0 after reconcile) settles to exactly 0.
    assert realized_pnl_usd("yes", 0, 0.55, 0.01, "no") == 0.0


def test_collect_fee_fields_nested():
    from kalshi_mlb_mm.settlement import _collect_fee_fields
    found = _collect_fee_fields(
        {"a": {"fills": [{"fee_cost": "0.10", "count": 3},
                         {"maker_fee": 0.02}]},
         "fee_total": "not_a_number"})
    assert found == {"a.fills[0].fee_cost": 0.10, "a.fills[1].maker_fee": 0.02}
