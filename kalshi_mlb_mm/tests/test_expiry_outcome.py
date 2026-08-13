"""Regression tests for the expired-quote outcome labeler.

These fail on pre-labeler code by construction: there was no
`quote_expiry_outcomes` table and nothing read the combo trade tape. All
Kalshi API traffic is stubbed via auth_client.api.

Timezone note: live_quotes.submitted_at/closed_at come back from DuckDB as
America/Los_Angeles tz-aware datetimes; the tape's created_time is UTC. The
harness inserts LA-aware timestamps on purpose so a naive or wrong-zone
comparison shows up as a mislabeled quote (LA and UTC differ by 7-8h — far
larger than any test window).
"""
import importlib
from datetime import datetime, timedelta, timezone
from zoneinfo import ZoneInfo

import pytest

LA = ZoneInfo("America/Los_Angeles")

# Fixed reference instant, defined in UTC, far enough in the past that every
# quote's grace window has elapsed. Windows are built relative to it.
CLOSED_UTC = datetime(2026, 8, 13, 1, 0, 0, tzinfo=timezone.utc)
SUBMITTED_UTC = CLOSED_UTC - timedelta(seconds=30)


def _iso_z(dt: datetime) -> str:
    return dt.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------------------
# Harness (mirrors test_settlement.py)
# ---------------------------------------------------------------------------

def _setup_db(monkeypatch, tmp_path, name):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / f"{name}.duckdb")
    importlib.reload(db)
    db.init_database()
    return db


def _insert_quote(db, quote_id, ticker, status="expired",
                  submitted_at=None, closed_at=None):
    """Insert a live_quotes row with LA-aware timestamps (what DuckDB round-
    trips on this machine)."""
    submitted_at = (submitted_at or SUBMITTED_UTC).astimezone(LA)
    closed_at = closed_at.astimezone(LA) if closed_at is not None \
        else CLOSED_UTC.astimezone(LA)
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
            "game_id, yes_bid, no_bid, status, submitted_at, closed_at) "
            "VALUES (?,?,?,?,?,?,?,?,?)",
            [quote_id, f"r-{quote_id}", ticker, "g1", 0.45, 0.49, status,
             submitted_at, closed_at])


def _outcome_rows(db):
    with db.connect(read_only=True) as con:
        return con.execute(
            "SELECT quote_id, combo_market_ticker, quote_status, label, "
            "n_trades_in_window, first_trade_after_close_sec "
            "FROM quote_expiry_outcomes ORDER BY quote_id").fetchall()


def _label_of(db, quote_id):
    rows = [r for r in _outcome_rows(db) if r[0] == quote_id]
    return rows[0][3] if rows else None


class TapeStub:
    """Programmable auth_client.api stand-in serving /markets/trades.

    Reproduces the live-verified MVE-combo quirk: yes_price/no_price/count
    are None on every print — only created_time and taker_side carry data.
    """

    def __init__(self, trade_times=None, http=200, raise_exc=None,
                 trades_by_ticker=None):
        self.trade_times = trade_times or []
        self.trades_by_ticker = trades_by_ticker
        self.http = http
        self.raise_exc = raise_exc
        self.calls = []

    def __call__(self, method, path, body=None, timeout=30):
        self.calls.append((method, path, timeout))
        if self.raise_exc is not None:
            raise self.raise_exc
        assert path.startswith("/markets/trades?ticker="), path
        if self.http != 200:
            return self.http, {"error": "boom"}, {}
        ticker = path.split("ticker=")[1].split("&")[0]
        times = (self.trades_by_ticker.get(ticker, [])
                 if self.trades_by_ticker is not None else self.trade_times)
        trades = [{"trade_id": f"t{i}", "ticker": ticker,
                   "created_time": _iso_z(t), "taker_side": "yes",
                   "yes_price": None, "no_price": None, "count": None}
                  for i, t in enumerate(times)]
        return 200, {"trades": trades, "cursor": ""}, {}


@pytest.fixture()
def emitted(monkeypatch):
    from kalshi_mlb_mm import expiry_outcome
    events = []
    monkeypatch.setattr(expiry_outcome.research, "emit",
                        lambda et, **kw: events.append((et, kw)))
    return events


def _patch_api(monkeypatch, stub):
    from kalshi_mlb_mm import expiry_outcome
    monkeypatch.setattr(expiry_outcome.auth_client, "api", stub)


# ---------------------------------------------------------------------------
# Pure labeling rule
# ---------------------------------------------------------------------------

def test_label_in_window_trade_is_competitor():
    from kalshi_mlb_mm import expiry_outcome as eo
    label, n, after = eo.label_quote(
        SUBMITTED_UTC, CLOSED_UTC,
        [SUBMITTED_UTC + timedelta(seconds=10)], grace_sec=300)
    assert label == eo.LABEL_COMPETITOR_TRADED
    assert n == 1
    assert after is None


def test_label_window_bounds_inclusive():
    from kalshi_mlb_mm import expiry_outcome as eo
    for boundary in (SUBMITTED_UTC, CLOSED_UTC):
        label, n, _ = eo.label_quote(SUBMITTED_UTC, CLOSED_UTC, [boundary],
                                     grace_sec=300)
        assert label == eo.LABEL_COMPETITOR_TRADED
        assert n == 1


def test_label_grace_window_trade_is_shortly_after():
    from kalshi_mlb_mm import expiry_outcome as eo
    label, n, after = eo.label_quote(
        SUBMITTED_UTC, CLOSED_UTC,
        [CLOSED_UTC + timedelta(seconds=120)], grace_sec=300)
    assert label == eo.LABEL_TRADED_SHORTLY_AFTER
    assert n == 0
    assert after == pytest.approx(120.0)


def test_label_trade_past_grace_is_no_trade_with_diagnostic():
    from kalshi_mlb_mm import expiry_outcome as eo
    label, n, after = eo.label_quote(
        SUBMITTED_UTC, CLOSED_UTC,
        [CLOSED_UTC + timedelta(seconds=900)], grace_sec=300)
    assert label == eo.LABEL_NO_TRADE
    assert n == 0
    assert after == pytest.approx(900.0)    # diagnostic still recorded


def test_label_empty_tape_is_no_trade():
    from kalshi_mlb_mm import expiry_outcome as eo
    label, n, after = eo.label_quote(SUBMITTED_UTC, CLOSED_UTC, [],
                                     grace_sec=300)
    assert label == eo.LABEL_NO_TRADE
    assert n == 0
    assert after is None


def test_label_trade_before_window_ignored():
    from kalshi_mlb_mm import expiry_outcome as eo
    label, n, after = eo.label_quote(
        SUBMITTED_UTC, CLOSED_UTC,
        [SUBMITTED_UTC - timedelta(seconds=60)], grace_sec=300)
    assert label == eo.LABEL_NO_TRADE
    assert n == 0
    assert after is None


def test_label_in_window_beats_grace_and_counts_all():
    from kalshi_mlb_mm import expiry_outcome as eo
    label, n, after = eo.label_quote(
        SUBMITTED_UTC, CLOSED_UTC,
        [SUBMITTED_UTC + timedelta(seconds=5),
         SUBMITTED_UTC + timedelta(seconds=20),
         CLOSED_UTC + timedelta(seconds=60)], grace_sec=300)
    assert label == eo.LABEL_COMPETITOR_TRADED
    assert n == 2
    assert after == pytest.approx(60.0)


# ---------------------------------------------------------------------------
# Trade-time parsing (None-price quirk rides through the sweep tests below)
# ---------------------------------------------------------------------------

def test_parse_trade_time_z_suffix_and_fractional():
    from kalshi_mlb_mm import expiry_outcome as eo
    parsed = eo.parse_trade_time("2026-08-13T01:00:00Z")
    assert parsed == CLOSED_UTC
    parsed = eo.parse_trade_time("2026-08-13T01:00:00.123456Z")
    assert parsed is not None and parsed.tzinfo is not None


def test_parse_trade_time_garbage_is_none():
    from kalshi_mlb_mm import expiry_outcome as eo
    assert eo.parse_trade_time(None) is None
    assert eo.parse_trade_time("not-a-time") is None


# ---------------------------------------------------------------------------
# Sweep end-to-end: tz conversion, quirk tolerance, dedupe, statuses
# ---------------------------------------------------------------------------

def test_sweep_labels_competitor_across_timezones(monkeypatch, tmp_path,
                                                  emitted):
    """LA-aware quote window vs UTC tape print inside it — a naive or
    wrong-zone comparison misses by hours and mislabels as no_trade."""
    from kalshi_mlb_mm import expiry_outcome as eo
    db = _setup_db(monkeypatch, tmp_path, "tz")
    _insert_quote(db, "q1", "TKR-A")
    _patch_api(monkeypatch,
               TapeStub(trade_times=[SUBMITTED_UTC + timedelta(seconds=10)]))
    eo.expiry_outcome_sweep_tick()
    assert _label_of(db, "q1") == eo.LABEL_COMPETITOR_TRADED
    (event_type, kw), = emitted
    assert event_type == "quote_expiry_outcome"
    assert kw["quote_id"] == "q1"
    assert kw["payload"]["label"] == eo.LABEL_COMPETITOR_TRADED
    assert kw["payload"]["n_trades_in_window"] == 1


def test_sweep_no_trade_and_none_price_quirk(monkeypatch, tmp_path, emitted):
    """Tape prints with yes_price/no_price/count == None must still parse —
    the labeler only consumes created_time (+ taker_side in the payload)."""
    from kalshi_mlb_mm import expiry_outcome as eo
    db = _setup_db(monkeypatch, tmp_path, "quirk")
    _insert_quote(db, "q1", "TKR-A")
    # Only print is well past grace → no_trade, despite None-valued prices.
    _patch_api(monkeypatch,
               TapeStub(trade_times=[CLOSED_UTC + timedelta(seconds=3600)]))
    eo.expiry_outcome_sweep_tick()
    assert _label_of(db, "q1") == eo.LABEL_NO_TRADE
    row = _outcome_rows(db)[0]
    assert row[4] == 0                                # n_trades_in_window
    assert row[5] == pytest.approx(3600.0)            # first_after_close_sec


def test_sweep_grace_label(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import expiry_outcome as eo
    db = _setup_db(monkeypatch, tmp_path, "grace")
    _insert_quote(db, "q1", "TKR-A")
    _patch_api(monkeypatch,
               TapeStub(trade_times=[CLOSED_UTC + timedelta(seconds=90)]))
    eo.expiry_outcome_sweep_tick()
    assert _label_of(db, "q1") == eo.LABEL_TRADED_SHORTLY_AFTER


def test_sweep_dedupe_labels_once(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import expiry_outcome as eo
    db = _setup_db(monkeypatch, tmp_path, "dedupe")
    _insert_quote(db, "q1", "TKR-A")
    stub = TapeStub()
    _patch_api(monkeypatch, stub)
    eo.expiry_outcome_sweep_tick()
    eo.expiry_outcome_sweep_tick()
    assert len(_outcome_rows(db)) == 1
    assert len(emitted) == 1
    assert len(stub.calls) == 1        # second sweep found nothing pending


def test_sweep_one_fetch_per_distinct_ticker(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import expiry_outcome as eo
    db = _setup_db(monkeypatch, tmp_path, "onefetch")
    _insert_quote(db, "q1", "TKR-A")
    _insert_quote(db, "q2", "TKR-A",
                  submitted_at=SUBMITTED_UTC - timedelta(minutes=10),
                  closed_at=CLOSED_UTC - timedelta(minutes=10))
    _insert_quote(db, "q3", "TKR-B")
    stub = TapeStub()
    _patch_api(monkeypatch, stub)
    eo.expiry_outcome_sweep_tick()
    assert len(_outcome_rows(db)) == 3
    assert len(stub.calls) == 2        # TKR-A once, TKR-B once


def test_sweep_waits_for_grace_window(monkeypatch, tmp_path, emitted):
    """A quote closed 60s ago (grace 300) is not yet decidable — labeled on
    a later sweep once the grace window has elapsed."""
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import expiry_outcome as eo
    db = _setup_db(monkeypatch, tmp_path, "wait")
    now_utc = datetime.now(timezone.utc)
    _insert_quote(db, "q1", "TKR-A",
                  submitted_at=now_utc - timedelta(seconds=90),
                  closed_at=now_utc - timedelta(seconds=60))
    stub = TapeStub()
    _patch_api(monkeypatch, stub)
    eo.expiry_outcome_sweep_tick()
    assert _outcome_rows(db) == []
    assert stub.calls == []            # not even fetched
    monkeypatch.setattr(cfg, "EXPIRY_OUTCOME_GRACE_SEC", 30)
    eo.expiry_outcome_sweep_tick()     # grace now elapsed relative to close
    assert _label_of(db, "q1") == eo.LABEL_NO_TRADE


def test_sweep_ignores_open_and_filled_quotes(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import expiry_outcome as eo
    db = _setup_db(monkeypatch, tmp_path, "statuses")
    _insert_quote(db, "q1", "TKR-A", status="open", closed_at=CLOSED_UTC)
    _insert_quote(db, "q2", "TKR-A", status="executed")
    _insert_quote(db, "q3", "TKR-A", status="voided")
    stub = TapeStub()
    _patch_api(monkeypatch, stub)
    eo.expiry_outcome_sweep_tick()
    assert _outcome_rows(db) == []
    assert stub.calls == []


def test_sweep_cancelled_gated_by_knob(monkeypatch, tmp_path, emitted):
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import expiry_outcome as eo
    db = _setup_db(monkeypatch, tmp_path, "cancelled")
    _insert_quote(db, "q1", "TKR-A", status="cancelled")
    _patch_api(monkeypatch, TapeStub())
    eo.expiry_outcome_sweep_tick()
    assert _outcome_rows(db) == []                     # off by default
    monkeypatch.setattr(cfg, "EXPIRY_OUTCOME_INCLUDE_CANCELLED", True)
    eo.expiry_outcome_sweep_tick()
    assert _label_of(db, "q1") == eo.LABEL_NO_TRADE
    assert _outcome_rows(db)[0][2] == "cancelled"      # quote_status recorded


# ---------------------------------------------------------------------------
# Failure tolerance / bounds
# ---------------------------------------------------------------------------

def test_sweep_api_error_retried_next_sweep(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import expiry_outcome as eo
    db = _setup_db(monkeypatch, tmp_path, "apierr")
    _insert_quote(db, "q1", "TKR-A")
    _patch_api(monkeypatch, TapeStub(http=500))
    eo.expiry_outcome_sweep_tick()             # tape down — no label, no crash
    assert _outcome_rows(db) == []
    _patch_api(monkeypatch, TapeStub())
    eo.expiry_outcome_sweep_tick()             # next sweep succeeds
    assert _label_of(db, "q1") == eo.LABEL_NO_TRADE


def test_sweep_api_exception_never_raises(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import expiry_outcome as eo
    db = _setup_db(monkeypatch, tmp_path, "apiexc")
    _insert_quote(db, "q1", "TKR-A")
    _patch_api(monkeypatch, TapeStub(raise_exc=RuntimeError("conn reset")))
    eo.expiry_outcome_sweep_tick()             # must not propagate
    assert _outcome_rows(db) == []


def test_sweep_ticker_failure_isolated(monkeypatch, tmp_path, emitted):
    """A failing ticker must not block labeling of the others."""
    from kalshi_mlb_mm import expiry_outcome as eo
    db = _setup_db(monkeypatch, tmp_path, "isolate")
    _insert_quote(db, "q1", "TKR-BAD")
    _insert_quote(db, "q2", "TKR-OK")

    good = TapeStub()

    def flaky(method, path, body=None, timeout=30):
        if "TKR-BAD" in path:
            raise RuntimeError("boom")
        return good(method, path, body=body, timeout=timeout)

    _patch_api(monkeypatch, flaky)
    eo.expiry_outcome_sweep_tick()
    assert _label_of(db, "q2") == eo.LABEL_NO_TRADE
    assert _label_of(db, "q1") is None         # retried next sweep


def test_sweep_ticker_cap_bounds_fetches(monkeypatch, tmp_path, emitted):
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import expiry_outcome as eo
    db = _setup_db(monkeypatch, tmp_path, "cap")
    for i in range(4):
        _insert_quote(db, f"q{i}", f"TKR-{i}")
    monkeypatch.setattr(cfg, "EXPIRY_OUTCOME_MAX_TICKERS_PER_SWEEP", 2)
    stub = TapeStub()
    _patch_api(monkeypatch, stub)
    eo.expiry_outcome_sweep_tick()
    assert len(stub.calls) == 2
    assert len(_outcome_rows(db)) == 2
    eo.expiry_outcome_sweep_tick()             # remainder next sweep
    assert len(_outcome_rows(db)) == 4


def test_api_calls_use_bounded_timeout(monkeypatch, tmp_path, emitted):
    """The sweep shares the trading loop's thread — every tape call must use
    a short explicit timeout so a hung call can't starve the confirm loop."""
    from kalshi_mlb_mm import expiry_outcome as eo
    db = _setup_db(monkeypatch, tmp_path, "timeouts")
    _insert_quote(db, "q1", "TKR-A")
    stub = TapeStub()
    _patch_api(monkeypatch, stub)
    eo.expiry_outcome_sweep_tick()
    assert stub.calls
    assert all(timeout <= 10 for _, _, timeout in stub.calls)


def test_no_pending_quotes_no_api_calls(monkeypatch, tmp_path, emitted):
    from kalshi_mlb_mm import expiry_outcome as eo
    _setup_db(monkeypatch, tmp_path, "empty")
    stub = TapeStub()
    _patch_api(monkeypatch, stub)
    eo.expiry_outcome_sweep_tick()
    assert stub.calls == []
