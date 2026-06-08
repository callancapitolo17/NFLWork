"""Tests for the pre-launch hardening bundle: N7/N8/N9/N10/N11/N12."""

# ---------------------------------------------------------------------------
# N7 — per-combo cap counts in-flight live_quotes' worst-case exposure.
# Pre-seed 4 open live_quotes on the ticker (each worth MAX_RFQ_CONTRACTS=$5
# worst-case). With MAX_COMBO_EXPOSURE_USD=10, the cap check must fire even
# though fills is empty (inflight + this-quote-max already exceeds cap).
# ---------------------------------------------------------------------------
def test_n7_inflight_quotes_trigger_per_combo_cap(monkeypatch, tmp_path):
    from datetime import datetime, timezone
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "n7.duckdb")
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    # MAX_RFQ_CONTRACTS=2, MAX_COMBO_EXPOSURE_USD=10 → 4 open quotes × 2 = 8,
    # plus this-quote-max=2 → 10 which is NOT > 10.  Use 5+1=6 > 10 via cap=5.
    # Simpler: 4 open quotes × MAX_RFQ_CONTRACTS (5 default) = 20 >> cap=10.
    monkeypatch.setattr(cfg, "MAX_COMBO_EXPOSURE_USD", 10.0)
    monkeypatch.setattr(cfg, "MAX_RFQ_CONTRACTS", 2)
    # With MAX_RFQ_CONTRACTS=2: 4 in-flight × $2 worst = $8, + this-quote $2 = $10.
    # 10 > 10 is False, so change to MAX_RFQ_CONTRACTS=3 → 4×3+3=15>10=True.
    monkeypatch.setattr(cfg, "MAX_RFQ_CONTRACTS", 3)

    import importlib
    importlib.reload(db)
    db.init_database()

    import pandas as pd
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["g7"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, m: True)
    monkeypatch.setattr(main, "_today_fills", lambda: [])

    legs = [{"market_ticker": "KXMLBSPREAD-N7", "event_ticker": "EVT-N7",
             "side": "yes", "count": 1},
            {"market_ticker": "KXMLBTOTAL-N7", "event_ticker": "EVT-N7",
             "side": "yes", "count": 1}]
    monkeypatch.setattr(main, "_SCOPE_CACHE", {"COMBO-N7": (True, "g7", legs)})
    monkeypatch.setattr(main, "_resolve_game_and_lines",
                        lambda ticker, legs: ("g7", -1.5, 8.5))

    # Pre-seed 4 open live_quotes on COMBO-N7 — fills table is EMPTY.
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        for i in range(4):
            con.execute(
                "INSERT INTO live_quotes VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                [f"qid-n7-{i}", f"r-n7-{i}", "COMBO-N7", "g7",
                 0.50, 0.43, None, 0.55, 0.55, "open", now, None])

    submit_calls = []

    class GW:
        def submit_quote(self, *a):
            submit_calls.append(a)
            return "qid-new-n7"

    class Src:
        def poll(self):
            # 5th RFQ on the same combo
            return [{"id": "r-n7-new", "market_ticker": "COMBO-N7", "contracts": 3}]

        def get_market(self, t):
            return {}

    main._discovery_tick(Src(), GW(), dry_run=False)

    assert submit_calls == [], f"N7: must not submit when in-flight quotes exceed combo cap"
    with db.connect(read_only=True) as con:
        d = con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "WHERE rfq_id='r-n7-new' ORDER BY observed_at DESC LIMIT 1").fetchone()
    assert d is not None, "N7: expected a skip decision row"
    assert d[0] == "skipped" and d[1] == "per_combo_cap", (
        f"N7: expected per_combo_cap, got {d}")


# ---------------------------------------------------------------------------
# N8 — unreconciled fills are counted conservatively (GREATEST(contracts,
# MAX_RFQ_CONTRACTS)) in cap queries.  Pre-seed one unreconciled fill with
# contracts=1 on a game; the daily/per-game cap is just above 1*price but
# below MAX_RFQ_CONTRACTS*price, so only the conservative counting triggers.
# ---------------------------------------------------------------------------
def test_n8_unreconciled_fill_counted_conservatively_in_today_fills(monkeypatch, tmp_path):
    from datetime import datetime, timezone
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "n8.duckdb")
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")

    # BANKROLL=100, DAILY_EXPOSURE_CAP_PCT=1.0 → daily cap = 100.
    # MAX_RFQ_CONTRACTS=20. One unreconciled fill: contracts=1, price=0.50
    # → real exposure = 0.50 but conservative = 20 * 0.50 = 10.0.
    # Make BANKROLL=10 and DAILY_CAP_PCT=1.0 → cap=10.
    # Conservative exposure = 10, new quote MAX = 20*1 = 20. Total > 10 → skip.
    # But we test via _today_fills() indirectly — we need to confirm the
    # per_game_cap or daily_cap fires due to conservative counting.
    monkeypatch.setattr(cfg, "BANKROLL", 5.0)
    monkeypatch.setattr(cfg, "DAILY_EXPOSURE_CAP_PCT", 1.0)  # cap = $5
    monkeypatch.setattr(cfg, "MAX_RFQ_CONTRACTS", 10)        # conservative = 10 contracts

    import importlib
    importlib.reload(db)
    db.init_database()

    import pandas as pd
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["g8"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, m: True)

    legs = [{"market_ticker": "KXMLBSPREAD-N8", "event_ticker": "EVT-N8",
             "side": "yes", "count": 1},
            {"market_ticker": "KXMLBTOTAL-N8", "event_ticker": "EVT-N8",
             "side": "yes", "count": 1}]
    monkeypatch.setattr(main, "_SCOPE_CACHE", {"COMBO-N8": (True, "g8", legs)})
    monkeypatch.setattr(main, "_resolve_game_and_lines",
                        lambda ticker, legs: ("g8", -1.5, 8.5))

    # Pre-seed one UNRECONCILED fill with contracts=1, price=0.50.
    # Real exposure = $0.50 (well under $5 cap).
    # Conservative exposure = 10 × $0.50 = $5.00 (at cap → daily_cap fires).
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
            "game_id, side_held, contracts, price, fee, model_fair_at_quote, "
            "book_fair_at_quote, blended_fair_at_quote, fair_at_confirm, "
            "realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["fill-n8", "qid-n8", "r-n8-prior", "COMBO-N8", "g8",
             "yes", 1, 0.50, 0.0, None, 0.5, 0.5, 0.5, None, now, False])

    # Verify _today_fills returns conservative amounts when reconciled=False.
    fills = main._today_fills()
    assert len(fills) == 1, f"Expected 1 fill, got {len(fills)}"
    effective_price = fills[0]["price"]
    # Conservative: GREATEST(1, 10) × 0.50 = 10 × 0.50 = 5.0
    assert effective_price == 5.0, (
        f"N8: expected conservative price=5.0, got {effective_price}")

    submit_calls = []

    class GW:
        def submit_quote(self, *a):
            submit_calls.append(a)
            return "qid-n8-new"

    class Src:
        def poll(self):
            return [{"id": "r-n8-new", "market_ticker": "COMBO-N8", "contracts": 1}]

        def get_market(self, t):
            return {}

    main._discovery_tick(Src(), GW(), dry_run=False)

    # Cap should fire: conservative fill already at $5 = daily cap.
    assert submit_calls == [], f"N8: must skip when unreconciled fill conservatively fills cap"
    with db.connect(read_only=True) as con:
        d = con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "WHERE rfq_id='r-n8-new' ORDER BY observed_at DESC LIMIT 1").fetchone()
    assert d is not None, "N8: expected a skip decision row"
    assert d[0] == "skipped" and d[1] in ("daily_cap", "per_game_cap", "per_combo_cap"), (
        f"N8: expected a cap reason, got {d}")


# ---------------------------------------------------------------------------
# N9 — POST must NOT be retried on 429; GET that returns 429 then 200 should
# retry and return 200.
# ---------------------------------------------------------------------------
import io
import json
import urllib.error


def _configure_auth_stub():
    from kalshi_common import auth_client
    auth_client._API_KEY_ID = "stub-key"
    auth_client._PRIVATE_KEY_PATH = "/tmp/none"
    auth_client._BASE_URL = "https://example.invalid/trade-api/v2"
    auth_client._sign_request = lambda _pk, _ts, _m, _p: "fake-sig"


class _FakeResponse:
    def __init__(self, status: int, body: dict):
        self.status = status
        self._payload = json.dumps(body).encode()
        self.headers = {}

    def read(self):
        return self._payload

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _http_error(status: int, body: dict) -> urllib.error.HTTPError:
    payload = json.dumps(body).encode()
    return urllib.error.HTTPError(
        url="https://example.invalid", code=status, msg="x",
        hdrs=None, fp=io.BytesIO(payload))


def test_n9_post_not_retried_on_429(monkeypatch):
    """N9: POST must return immediately on 429, no retry, no sleep."""
    from kalshi_common import auth_client
    _configure_auth_stub()

    sleeps = []
    monkeypatch.setattr(auth_client.time, "sleep", lambda s: sleeps.append(s))

    calls = {"n": 0}

    def fake_urlopen(req, timeout=30):
        calls["n"] += 1
        raise _http_error(429, {"error": "rate_limited"})

    monkeypatch.setattr(auth_client.urllib.request, "urlopen", fake_urlopen)

    status, body, _ = auth_client.api("POST", "/communications/rfqs/r1/quotes")
    assert status == 429, f"N9: expected 429 back, got {status}"
    assert calls["n"] == 1, f"N9: POST must not retry, expected 1 call, got {calls['n']}"
    assert sleeps == [], f"N9: POST must not sleep on retry, got sleeps={sleeps}"


def test_n9_get_retries_429_then_200(monkeypatch):
    """N9: GET that returns 429 then 200 should retry and return 200 (unchanged)."""
    from kalshi_common import auth_client
    _configure_auth_stub()

    monkeypatch.setattr(auth_client.time, "sleep", lambda s: None)

    calls = {"n": 0}

    def fake_urlopen(req, timeout=30):
        calls["n"] += 1
        if calls["n"] == 1:
            raise _http_error(429, {"error": "rate_limited"})
        return _FakeResponse(200, {"ok": True})

    monkeypatch.setattr(auth_client.urllib.request, "urlopen", fake_urlopen)

    status, body, _ = auth_client.api("GET", "/exchange/status")
    assert status == 200, f"N9: expected 200 after retry, got {status}"
    assert calls["n"] == 2, f"N9: expected 2 calls (1 fail + 1 success), got {calls['n']}"


def test_n9_post_503_no_retry(monkeypatch):
    """N9: POST 503 must also not be retried."""
    from kalshi_common import auth_client
    _configure_auth_stub()

    calls = {"n": 0}

    def fake_urlopen(req, timeout=30):
        calls["n"] += 1
        raise _http_error(503, {"error": "unavailable"})

    monkeypatch.setattr(auth_client.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(auth_client.time, "sleep", lambda s: None)

    status, body, _ = auth_client.api("POST", "/some/endpoint")
    assert status == 503
    assert calls["n"] == 1, f"N9: POST 503 must not retry, got {calls['n']} calls"


# ---------------------------------------------------------------------------
# N10 — phantom fill: delta=0 when expected was non-zero → mark contracts=0,
# reconciled=TRUE, log [phantom_fill].
# ---------------------------------------------------------------------------
def test_n10_phantom_fill_marked_zero_contracts(monkeypatch, capsys, tmp_path):
    from datetime import datetime, timezone
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "n10.duckdb")
    import importlib
    importlib.reload(db)
    db.init_database()

    ticker = "COMBO-N10"
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, game_id, "
            "side_held, contracts, price, fee, model_fair_at_quote, book_fair_at_quote, "
            "blended_fair_at_quote, fair_at_confirm, realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["fill-n10", "qid-n10", "r-n10", ticker, "gN10",
             "yes", 5, 0.55, 0.01, None, 0.55, 0.55, 0.55, None, now, False])

    # Kalshi says position=0 (delta=0) but we recorded side='yes', contracts=5
    # (expected_signed=5 != 0). This is a phantom fill.
    # No other fills on this ticker, so prior_total=0, delta = 0 - 0 = 0.
    monkeypatch.setattr(main, "_get_position_contracts", lambda ticker, **kw: 0)

    main._reconcile_sweep_tick()

    captured = capsys.readouterr().out
    assert "[phantom_fill]" in captured, f"N10: expected [phantom_fill] in output: {captured!r}"

    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT contracts, reconciled FROM fills WHERE fill_id='fill-n10'").fetchone()
    assert row is not None
    contracts, reconciled = row
    assert int(contracts) == 0, f"N10: expected contracts=0 (phantom), got {contracts}"
    assert reconciled is True, "N10: phantom fill must be marked reconciled=TRUE"


# ---------------------------------------------------------------------------
# N11 — max-age fallback: fill older than MAX_RECONCILE_AGE_SEC with positions
# API down → mark reconciled=TRUE with recorded values, log
# [reconcile_max_age_fallback].
# ---------------------------------------------------------------------------
def test_n11_max_age_fallback_marks_reconciled(monkeypatch, capsys, tmp_path):
    from datetime import datetime, timedelta, timezone
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "n11.duckdb")
    monkeypatch.setattr(cfg, "MAX_RECONCILE_AGE_SEC", 300)
    import importlib
    importlib.reload(db)
    db.init_database()

    ticker = "COMBO-N11"
    # Fill inserted 400 seconds ago (older than MAX_RECONCILE_AGE_SEC=300).
    old_ts = datetime.now(timezone.utc) - timedelta(seconds=400)
    with db.connect() as con:
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, game_id, "
            "side_held, contracts, price, fee, model_fair_at_quote, book_fair_at_quote, "
            "blended_fair_at_quote, fair_at_confirm, realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["fill-n11", "qid-n11", "r-n11", ticker, "gN11",
             "yes", 3, 0.55, 0.01, None, 0.55, 0.55, 0.55, None, old_ts, False])

    # Positions API is persistently down.
    monkeypatch.setattr(main, "_get_position_contracts", lambda ticker, **kw: None)

    main._reconcile_sweep_tick()

    captured = capsys.readouterr().out
    assert "[reconcile_max_age_fallback]" in captured, (
        f"N11: expected [reconcile_max_age_fallback] in output: {captured!r}")

    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT contracts, side_held, reconciled FROM fills WHERE fill_id='fill-n11'").fetchone()
    assert row is not None
    contracts, side_held, reconciled = row
    # Recorded values must be preserved (best-effort fallback).
    assert int(contracts) == 3, f"N11: recorded contracts should be preserved, got {contracts}"
    assert side_held == "yes", f"N11: recorded side should be preserved, got {side_held}"
    assert reconciled is True, "N11: fill must be marked reconciled=TRUE after max-age fallback"


def test_n11_young_fill_stays_unreconciled_when_api_down(monkeypatch, capsys, tmp_path):
    """N11: A fill younger than MAX_RECONCILE_AGE_SEC stays reconciled=FALSE when API is down."""
    from datetime import datetime, timezone
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "n11_young.duckdb")
    monkeypatch.setattr(cfg, "MAX_RECONCILE_AGE_SEC", 300)
    import importlib
    importlib.reload(db)
    db.init_database()

    ticker = "COMBO-N11Y"
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, game_id, "
            "side_held, contracts, price, fee, model_fair_at_quote, book_fair_at_quote, "
            "blended_fair_at_quote, fair_at_confirm, realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["fill-n11y", "qid-n11y", "r-n11y", ticker, "gN11Y",
             "yes", 3, 0.55, 0.01, None, 0.55, 0.55, 0.55, None, now, False])

    monkeypatch.setattr(main, "_get_position_contracts", lambda ticker, **kw: None)

    main._reconcile_sweep_tick()

    captured = capsys.readouterr().out
    # Should see position_reconcile_unavailable (young fill, retry later).
    assert "position_reconcile_unavailable" in captured, (
        f"N11: young fill should log unavailable, got: {captured!r}")
    # Must NOT have been marked reconciled.
    with db.connect(read_only=True) as con:
        rec = con.execute(
            "SELECT reconciled FROM fills WHERE fill_id='fill-n11y'").fetchone()[0]
    assert rec is False, "N11: young fill must stay reconciled=FALSE when API is down"


# ---------------------------------------------------------------------------
# N12 — void-rate halt edge-triggered notification.
# ---------------------------------------------------------------------------
def test_n12_void_rate_halt_triggers_notify_on_transition(monkeypatch, tmp_path):
    """N12: notify.halt('void_rate') called exactly once when halt becomes True."""
    from datetime import datetime, timezone
    import uuid as _uuid
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.notify as notify_mod
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "n12.duckdb")
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    import importlib
    importlib.reload(db)
    db.init_database()

    import pandas as pd
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["g"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))

    # Ensure halt starts from False state.
    monkeypatch.setattr(main, "_VOID_HALT_ACTIVE", False)

    # Record halt/resume calls.
    halt_calls = []
    resume_calls = []
    monkeypatch.setattr(notify_mod, "halt", lambda reason, detail="": halt_calls.append((reason, detail)))
    monkeypatch.setattr(notify_mod, "resume", lambda reason, detail="": resume_calls.append((reason, detail)))

    # Pre-seed 4 voids + 1 confirmed → 80% void rate > 25% threshold.
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        for d in ("voided_last_look", "voided_no_fresh_books",
                  "voided_blend_failed", "voided_no_legs", "confirmed"):
            con.execute(
                "INSERT INTO quote_decisions (decision_id, rfq_id, quote_id, "
                "combo_market_ticker, game_id, decision, reason, model_fair, "
                "book_fair, blended_fair, yes_bid, no_bid, observed_at) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                [str(_uuid.uuid4()), "r-n12", "qid-n12", "COMBO-N12", "gN12",
                 d, None, None, None, None, None, None, now])

    class Src:
        def poll(self):
            return []

        def get_market(self, t):
            return {}

    class GW:
        def submit_quote(self, *a):
            raise AssertionError("must not submit when halted")

    # First tick: False → True transition.
    main._discovery_tick(Src(), GW(), dry_run=False)

    assert len(halt_calls) == 1, f"N12: expected 1 halt call on False→True, got {halt_calls}"
    assert halt_calls[0][0] == "void_rate", (
        f"N12: expected halt reason='void_rate', got '{halt_calls[0][0]}'")
    assert "void_rate=" in halt_calls[0][1], (
        f"N12: expected void_rate detail in halt call, got '{halt_calls[0][1]}'")
    assert resume_calls == [], f"N12: no resume on first halt, got {resume_calls}"

    # Second tick: still in halt (True → True) → no new notification.
    halt_calls.clear()
    resume_calls.clear()
    main._discovery_tick(Src(), GW(), dry_run=False)
    assert halt_calls == [], f"N12: no re-notification when already halted, got {halt_calls}"
    assert resume_calls == [], f"N12: no resume when still halted, got {resume_calls}"


def test_n12_void_rate_resume_triggers_notify_on_recovery(monkeypatch, tmp_path):
    """N12: notify.resume called exactly once when halt goes from True to False."""
    from datetime import datetime, timezone
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.notify as notify_mod
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "n12r.duckdb")
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    import importlib
    importlib.reload(db)
    db.init_database()

    import pandas as pd
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["g"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))

    # Start in halt-active state (True).
    monkeypatch.setattr(main, "_VOID_HALT_ACTIVE", True)

    halt_calls = []
    resume_calls = []
    monkeypatch.setattr(notify_mod, "halt", lambda reason, detail="": halt_calls.append((reason, detail)))
    monkeypatch.setattr(notify_mod, "resume", lambda reason, detail="": resume_calls.append((reason, detail)))

    # No voids in quote_decisions → void rate = 0 → halt is False.
    # quote_decisions table is empty (fresh db).

    class Src:
        def poll(self):
            return []

        def get_market(self, t):
            return {}

    class GW:
        def submit_quote(self, *a):
            raise AssertionError("unexpected submit")

    # Tick: True → False transition.
    main._discovery_tick(Src(), GW(), dry_run=False)

    assert len(resume_calls) == 1, (
        f"N12: expected 1 resume call on True→False, got {resume_calls}")
    assert resume_calls[0][0] == "void_rate", (
        f"N12: expected resume reason='void_rate', got '{resume_calls[0][0]}'")
    assert halt_calls == [], f"N12: no halt call on recovery, got {halt_calls}"
    # _VOID_HALT_ACTIVE should now be False.
    assert main._VOID_HALT_ACTIVE is False, "N12: _VOID_HALT_ACTIVE must be False after recovery"
