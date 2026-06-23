def test_main_imports():
    from kalshi_mlb_mm import main          # must import without error
    assert hasattr(main, "main_loop")


def test_discovery_tick_skips_out_of_scope(monkeypatch, tmp_path):
    import kalshi_mlb_mm.config as cfg
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "t.duckdb")
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    import importlib, kalshi_mlb_mm.db as db
    importlib.reload(db)
    db.init_database()
    from kalshi_mlb_mm import main
    # v1 hardening replaced the samples-staleness gate with a book-freshness
    # gate. Provide a non-empty _SGP_ODDS so the gate passes and we reach the
    # out-of-scope logging path the test is exercising.
    import pandas as pd
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["g"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))

    class Src:
        def poll(self):
            return [{"id": "r1", "market_ticker": "OTHER-1", "contracts": 1}]

        def get_market(self, t):
            return {"mve_selected_legs": [{"market_ticker": "FOO"}]}

    class GW:
        def submit_quote(self, *a):
            raise AssertionError("should not submit out-of-scope")

    main._discovery_tick(Src(), GW(), dry_run=True)
    with db.connect(read_only=True) as con:
        d = con.execute("SELECT decision FROM quote_decisions").fetchone()
    assert d[0] == "skipped"


# ---------------------------------------------------------------------------
# FIX I1/M1 — dedup: an in-scope RFQ already open at the same price must NOT
# cause submit_quote to be called again on the next tick.
# ---------------------------------------------------------------------------
def test_discovery_dedup_no_resubmit_when_price_unchanged(monkeypatch, tmp_path):
    """An OPEN live_quote for rfq_id 'r1' with matching yes/no bid means the
    discovery tick should skip re-submission (hysteresis branch → continue).
    We assert no new 'quoted' decision row is written and submit_quote is not
    called.
    """
    from datetime import datetime, timezone
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "dedup.duckdb")

    import importlib
    importlib.reload(db)
    db.init_database()

    # v1 hardening: book-freshness gate replaced the samples gate. Set _SGP_ODDS
    # non-empty so the discovery tick proceeds.
    import pandas as pd
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["game1"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))
    # Bypass kill-switch: make sure KILL_FILE doesn't exist (tmp_path is clean).
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")

    # Quote that pricing would produce for blended=0.55, TARGET_ROI=0.05.
    # Rather than replicating the math exactly we monkeypatch pricing.quote to
    # return a fixed Quote so we know exactly what yes/no bids to pre-seed.
    from kalshi_mlb_mm.pricing import Quote
    import kalshi_mlb_mm.pricing as pricing_mod
    fixed_quote = Quote(yes_bid=0.500, no_bid=0.430)
    monkeypatch.setattr(pricing_mod, "quote", lambda fair, roi: fixed_quote)

    # Monkeypatch the helpers that require real DBs / network.
    monkeypatch.setattr(main, "_today_fills", lambda: [])
    monkeypatch.setattr(main, "_combo_games", lambda legs: ["game1"])
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    # Make tipoff_ok pass (commence_time is None → normally fails; override).
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)
    # General pricer returns (fair, agreeing_books).
    monkeypatch.setattr(main, "_price_combo", lambda legs: (0.55, 3))

    # Scope cache: make the ticker appear in-scope with a minimal legs list.
    legs = [{"market_ticker": "KXMLBSPREAD-ABC", "event_ticker": "EVT-ABC",
             "side": "yes", "count": 1},
            {"market_ticker": "KXMLBTOTAL-ABC", "event_ticker": "EVT-ABC",
             "side": "yes", "count": 1}]
    monkeypatch.setattr(main, "_SCOPE_CACHE",
                        {"COMBO-1": (True, legs)})

    # Pre-insert the open live_quote at exactly the fixed bid prices.
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            ["qid-existing", "r1", "COMBO-1", "game1",
             fixed_quote.yes_bid, fixed_quote.no_bid,
             0.55, 0.55, 0.55, "open",
             datetime.now(timezone.utc), None])

    submit_calls = []

    class GW:
        def submit_quote(self, rid, yb, nb):
            submit_calls.append((rid, yb, nb))
            return "qid-new"

    class Src:
        def poll(self):
            return [{"id": "r1", "market_ticker": "COMBO-1", "contracts": 1}]

        def get_market(self, t):
            return {}

    main._discovery_tick(Src(), GW(), dry_run=False)

    # submit_quote must NOT have been called (hysteresis skip)
    assert submit_calls == [], f"Expected no resubmit, got: {submit_calls}"
    # No new 'quoted' decision row
    with db.connect(read_only=True) as con:
        quoted = con.execute(
            "SELECT COUNT(*) FROM quote_decisions WHERE decision='quoted'").fetchone()[0]
    assert quoted == 0, f"Expected 0 quoted decisions, got {quoted}"


# ---------------------------------------------------------------------------
# FIX I3 — daily cap: an RFQ that would otherwise be quoted is skipped when
# today's fills have already consumed the full daily exposure cap.
# ---------------------------------------------------------------------------
def test_discovery_skips_when_daily_cap_exhausted(monkeypatch, tmp_path):
    """With _today_fills returning an amount >= daily_exposure_cap_usd(), the
    discovery tick must log 'skipped' with reason='daily_cap' and must NOT call
    submit_quote.
    """
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "cap.duckdb")

    import importlib
    importlib.reload(db)
    db.init_database()

    # Bypass guards that would short-circuit before reaching the cap check.
    # v1 hardening: book-freshness gate replaced samples-staleness — set
    # _SGP_ODDS non-empty so we reach the cap check.
    import pandas as pd
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["game2"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)

    # Scope cache: in-scope combo.
    legs = [{"market_ticker": "KXMLBSPREAD-XYZ", "event_ticker": "EVT-XYZ",
             "side": "yes", "count": 1},
            {"market_ticker": "KXMLBTOTAL-XYZ", "event_ticker": "EVT-XYZ",
             "side": "yes", "count": 1}]
    monkeypatch.setattr(main, "_SCOPE_CACHE",
                        {"COMBO-2": (True, legs)})

    # _combo_games must return a valid game_id so we reach the cap check.
    monkeypatch.setattr(main, "_combo_games", lambda legs: ["game2"])
    monkeypatch.setattr(main, "_price_combo", lambda legs: (0.55, 3))

    # Make daily cap exhausted: return fills totalling far above cap.
    big_fill = [{"game_id": "game2", "price": 99999.0}]
    monkeypatch.setattr(main, "_today_fills", lambda: big_fill)

    submit_calls = []

    class GW:
        def submit_quote(self, rid, yb, nb):
            submit_calls.append((rid, yb, nb))
            return "qid-cap"

    class Src:
        def poll(self):
            return [{"id": "r2", "market_ticker": "COMBO-2", "contracts": 1}]

        def get_market(self, t):
            return {}

    main._discovery_tick(Src(), GW(), dry_run=False)

    # Must not have submitted
    assert submit_calls == [], f"Expected no submit (daily cap), got: {submit_calls}"
    # Decision logged as 'skipped' with reason='daily_cap'
    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "WHERE rfq_id='r2' ORDER BY observed_at DESC LIMIT 1").fetchone()
    assert row is not None, "Expected a skipped decision to be logged"
    assert row[0] == "skipped", f"Expected 'skipped', got '{row[0]}'"
    assert row[1] == "daily_cap", f"Expected reason='daily_cap', got '{row[1]}'"


# ---------------------------------------------------------------------------
# v1 hardening (Change A) — book-freshness gate. With _SGP_ODDS None or empty,
# _discovery_tick must return BEFORE polling RFQs (no submit_quote possible,
# no quote_decisions written, no source.poll() call).
# ---------------------------------------------------------------------------
def test_discovery_returns_early_when_books_empty(monkeypatch, tmp_path):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "books_empty.duckdb")
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    import importlib
    importlib.reload(db)
    db.init_database()

    # Set _SGP_ODDS to None: book-freshness gate must fire.
    monkeypatch.setattr(main, "_SGP_ODDS", None)

    polled = []

    class Src:
        def poll(self):
            polled.append(1)
            return []

        def get_market(self, t):
            raise AssertionError("must not call get_market when books empty")

    class GW:
        def submit_quote(self, *a):
            raise AssertionError("must not submit when books empty")

    main._discovery_tick(Src(), GW(), dry_run=False)
    assert polled == [], "source.poll must not be called when books are stale/missing"
    with db.connect(read_only=True) as con:
        n = con.execute("SELECT COUNT(*) FROM quote_decisions").fetchone()[0]
    assert n == 0, "no decision rows should be written before the gate"

    # Also with an empty DataFrame:
    import pandas as pd
    monkeypatch.setattr(main, "_SGP_ODDS", pd.DataFrame())
    main._discovery_tick(Src(), GW(), dry_run=False)
    assert polled == []


# ---------------------------------------------------------------------------
# v1 hardening (Change C) — risk sweep cancels open quotes when current book
# consensus has drifted more than BOOK_MOVE_CB_THRESHOLD from book_fair-at-quote.
# Catches gradual drift the per-tick circuit breaker misses.
# ---------------------------------------------------------------------------
def test_risk_sweep_cancels_on_drift_since_quote(monkeypatch, tmp_path):
    from datetime import datetime, timezone
    import json
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "drift.duckdb")
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    import importlib
    importlib.reload(db)
    db.init_database()

    # books are fresh (so books_stale path doesn't fire) and tipoff is far away.
    import pandas as pd
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["g1"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)

    # Book_fair_at_q was 0.40; now combo consensus is 0.45 (drift 0.05 > 0.03).
    monkeypatch.setattr(main, "_combo_games", lambda legs: ["g1"])
    monkeypatch.setattr(main, "_price_combo", lambda legs: (0.45, 3))

    # Pre-seed an open live_quote referencing rfq r-drift, and a seen_rfqs row
    # with legs_json so the sweep can compute spread/total lines.
    legs = [{"market_ticker": "KXMLBSPREAD-G1", "event_ticker": "EVT-G1",
             "side": "yes", "count": 1, "no_sub_title": "-1.5"},
            {"market_ticker": "KXMLBTOTAL-G1", "event_ticker": "EVT-G1",
             "side": "yes", "count": 1, "no_sub_title": "Over 8.5"}]
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            ["qid-drift", "r-drift", "COMBO-G1", "g1",
             0.40, 0.45, None, 0.40, 0.40, "open",
             datetime.now(timezone.utc), None])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs "
            "(rfq_id, market_ticker, in_scope, game_id, legs_json, "
            "first_seen_at, last_decision) VALUES (?,?,?,?,?,?,?)",
            ["r-drift", "COMBO-G1", True, "g1", json.dumps(legs),
             datetime.now(timezone.utc), "quoted"])

    cancelled = []

    class GW:
        def cancel(self, qid):
            cancelled.append(qid)
            return True

    main._risk_sweep_tick(GW())

    assert "qid-drift" in cancelled, "drift > threshold should trigger cancel"
    with db.connect(read_only=True) as con:
        st = con.execute("SELECT status FROM live_quotes WHERE quote_id='qid-drift'").fetchone()[0]
    assert st == "cancelled", f"expected status='cancelled', got '{st}'"


# ---------------------------------------------------------------------------
# v1 hardening (Change D) — last-look hard-fail when fresh book fairs are
# unavailable at confirm time. Must VOID, not silently fall back to prev_fair.
# ---------------------------------------------------------------------------
def test_confirm_voids_when_no_fresh_books(monkeypatch, tmp_path):
    from datetime import datetime, timezone
    import json
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main
    from kalshi_common import auth_client

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "voidbooks.duckdb")
    import importlib
    importlib.reload(db)
    db.init_database()

    legs = [{"market_ticker": "KXMLBSPREAD-G2", "event_ticker": "EVT-G2",
             "side": "yes", "count": 1, "no_sub_title": "-1.5"},
            {"market_ticker": "KXMLBTOTAL-G2", "event_ticker": "EVT-G2",
             "side": "yes", "count": 1, "no_sub_title": "Over 8.5"}]
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            ["qid-nofresh", "r-nofresh", "COMBO-G2", "g2",
             0.50, 0.43, None, 0.55, 0.55, "open",
             datetime.now(timezone.utc), None])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs "
            "(rfq_id, market_ticker, in_scope, game_id, legs_json, "
            "first_seen_at, last_decision) VALUES (?,?,?,?,?,?,?)",
            ["r-nofresh", "COMBO-G2", True, "g2", json.dumps(legs),
             datetime.now(timezone.utc), "quoted"])

    # The quote came back 'accepted' from the API.
    def fake_api(method, path, *a, **kw):
        if path.startswith("/communications/quotes/"):
            return 200, {"quote": {"status": "accepted", "accepted_side": "yes",
                                   "contracts": 1}}, None
        raise AssertionError(f"unexpected api call: {method} {path}")
    monkeypatch.setattr(auth_client, "api", fake_api)

    # Can't re-price (no fresh books / failed consensus) — _price_combo None.
    monkeypatch.setattr(main, "_price_combo", lambda legs: None)

    confirm_calls = []

    class GW:
        def confirm(self, qid):
            confirm_calls.append(qid)
            return True

        def cancel(self, qid):
            return True

    main._confirm_tick(GW(), dry_run=False)

    assert confirm_calls == [], "confirm must NOT be called when it can't re-price"
    with db.connect(read_only=True) as con:
        st = con.execute(
            "SELECT status FROM live_quotes WHERE quote_id='qid-nofresh'").fetchone()[0]
        n_fills = con.execute("SELECT COUNT(*) FROM fills").fetchone()[0]
        decision = con.execute(
            "SELECT decision FROM quote_decisions WHERE quote_id='qid-nofresh' "
            "ORDER BY observed_at DESC LIMIT 1").fetchone()
    assert st == "voided", f"expected 'voided', got '{st}'"
    assert n_fills == 0, "no fills row should be written when voided"
    assert decision is not None and decision[0] == "voided_blend_failed"


# ---------------------------------------------------------------------------
# N5 — confirm path no longer reconciles inline. The fill row is INSERTed
# immediately with EXPECTED side/size and reconciled=FALSE; no positions API
# call is made in the hot path. (Previously this could take up to 15s/fill via
# retries — a burst of accepts blew the 30s confirm window.)
# ---------------------------------------------------------------------------
def test_confirm_records_fill_fast_with_reconciled_false(monkeypatch, tmp_path):
    from datetime import datetime, timezone
    import json
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main
    from kalshi_common import auth_client

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "confirm_fast.duckdb")
    import importlib
    importlib.reload(db)
    db.init_database()

    legs = [{"market_ticker": "KXMLBSPREAD-G3", "event_ticker": "EVT-G3",
             "side": "yes", "count": 1, "no_sub_title": "-1.5"},
            {"market_ticker": "KXMLBTOTAL-G3", "event_ticker": "EVT-G3",
             "side": "yes", "count": 1, "no_sub_title": "Over 8.5"}]
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            ["qid-fast", "r-fast", "COMBO-G3", "g3",
             0.45, 0.45, None, 0.55, 0.55, "open",
             datetime.now(timezone.utc), None])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs "
            "(rfq_id, market_ticker, in_scope, game_id, legs_json, "
            "first_seen_at, last_decision) VALUES (?,?,?,?,?,?,?)",
            ["r-fast", "COMBO-G3", True, "g3", json.dumps(legs),
             datetime.now(timezone.utc), "quoted"])

    monkeypatch.setattr(main, "_price_combo", lambda legs: (0.56, 3))
    monkeypatch.setattr(main, "_combo_games", lambda legs: ["g3"])

    # Track API calls; if confirm path ever hits /portfolio/positions, we fail.
    api_paths = []

    def fake_api(method, path, *a, **kw):
        api_paths.append(path)
        if path.startswith("/communications/quotes/"):
            return 200, {"quote": {"status": "accepted", "accepted_side": "no",
                                   "contracts": 5}}, None
        if path.startswith("/portfolio/positions"):
            raise AssertionError(
                "confirm path must NOT call /portfolio/positions (N5: moved to sweep)")
        raise AssertionError(f"unexpected api call: {method} {path}")
    monkeypatch.setattr(auth_client, "api", fake_api)

    # If the hot path ever called _get_position_contracts directly, this would
    # raise — confirm path must not call it.
    def _no_positions(*a, **kw):
        raise AssertionError("confirm hot path must not call _get_position_contracts")
    monkeypatch.setattr(main, "_get_position_contracts", _no_positions)

    class GW:
        def confirm(self, qid):
            return True

        def cancel(self, qid):
            return True

    main._confirm_tick(GW(), dry_run=False)

    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT side_held, contracts, reconciled FROM fills "
            "WHERE quote_id='qid-fast'").fetchone()
    assert row is not None, "fill row should be written immediately on confirm"
    side, contracts, reconciled = row
    # Expected side: accepted_side='no' means taker took our no_bid, so we hold 'yes'.
    assert side == "yes", f"expected expected-side='yes', got '{side}'"
    assert int(contracts) == 5, f"expected expected-contracts=5, got {contracts}"
    assert reconciled is False, f"new fill must be reconciled=FALSE, got {reconciled}"


# ---------------------------------------------------------------------------
# N5 — _reconcile_sweep_tick: picks up unreconciled fills and UPDATEs them to
# match Kalshi's reported position (delta logic). Sets reconciled=TRUE.
# ---------------------------------------------------------------------------
def test_reconcile_sweep_corrects_mismatched_fill(monkeypatch, caplog, tmp_path):
    from datetime import datetime, timezone
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "sweep.duckdb")
    import importlib
    importlib.reload(db)
    db.init_database()

    ticker = "COMBO-SWEEP"
    # Pre-seed a prior reconciled fill: 5 YES on this ticker.
    # And an unreconciled fill that the bot recorded as 5 YES expected, but
    # Kalshi actually shows aggregate=2 → delta = 2-5 = -3 → side='no', ct=3.
    with db.connect() as con:
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, game_id, "
            "side_held, contracts, price, fee, model_fair_at_quote, book_fair_at_quote, "
            "blended_fair_at_quote, fair_at_confirm, realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["fill-prior", "qid-prior", "r-prior", ticker, "gSW",
             "yes", 5, 0.55, 0.01, None, 0.55, 0.55, 0.55, None,
             datetime.now(timezone.utc), True])
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, game_id, "
            "side_held, contracts, price, fee, model_fair_at_quote, book_fair_at_quote, "
            "blended_fair_at_quote, fair_at_confirm, realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["fill-new", "qid-new", "r-new", ticker, "gSW",
             "yes", 5, 0.55, 0.01, None, 0.55, 0.55, 0.55, None,
             datetime.now(timezone.utc), False])

    # Kalshi reports aggregate position = 2 (so the new fill must be 'no', 3 ct).
    monkeypatch.setattr(main, "_get_position_contracts", lambda ticker, **kw: 2)

    main._reconcile_sweep_tick()

    assert "position_mismatch" in caplog.text, f"expected mismatch log, got: {caplog.text!r}"

    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT side_held, contracts, reconciled FROM fills "
            "WHERE fill_id='fill-new'").fetchone()
    assert row is not None
    side, contracts, reconciled = row
    assert side == "no", f"expected corrected side='no', got '{side}'"
    assert int(contracts) == 3, f"expected corrected contracts=3, got {contracts}"
    assert reconciled is True, "fill must be marked reconciled=TRUE after sweep"

    # The prior fill must NOT have been re-touched (still reconciled=True, side/ct intact).
    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT side_held, contracts, reconciled FROM fills "
            "WHERE fill_id='fill-prior'").fetchone()
    assert row == ("yes", 5, True)


# ---------------------------------------------------------------------------
# N5 — sweep tolerates API outage: if positions API returns None, the fill
# stays reconciled=FALSE for the next sweep to retry.
# ---------------------------------------------------------------------------
def test_reconcile_sweep_skips_when_api_down(monkeypatch, caplog, tmp_path):
    from datetime import datetime, timezone
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "sweep_down.duckdb")
    import importlib
    importlib.reload(db)
    db.init_database()

    ticker = "COMBO-DOWN"
    with db.connect() as con:
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, game_id, "
            "side_held, contracts, price, fee, model_fair_at_quote, book_fair_at_quote, "
            "blended_fair_at_quote, fair_at_confirm, realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["fill-down", "qid-down", "r-down", ticker, "gDN",
             "yes", 5, 0.55, 0.01, None, 0.55, 0.55, 0.55, None,
             datetime.now(timezone.utc), False])

    monkeypatch.setattr(main, "_get_position_contracts", lambda ticker, **kw: None)
    main._reconcile_sweep_tick()

    assert "position_reconcile_unavailable" in caplog.text

    with db.connect(read_only=True) as con:
        rec = con.execute(
            "SELECT reconciled FROM fills WHERE fill_id='fill-down'").fetchone()[0]
    assert rec is False, "API-down: fill must stay reconciled=FALSE for retry"


# ---------------------------------------------------------------------------
# H4 — global void-rate halt: too many recent voids triggers halt; no
# downstream RFQ processing happens.
# ---------------------------------------------------------------------------
def test_discovery_halts_on_high_void_rate(monkeypatch, tmp_path):
    from datetime import datetime, timezone
    import uuid as _uuid
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "void_halt.duckdb")
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    import importlib
    importlib.reload(db)
    db.init_database()

    # books fresh — so the gate before void-rate doesn't short-circuit.
    import pandas as pd
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["g"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))

    # Pre-seed 4 voids + 1 confirmed in the last hour → 80% void rate > 25%.
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        for d in ("voided_last_look", "voided_no_fresh_books",
                  "voided_blend_failed", "voided_no_legs", "confirmed"):
            con.execute(
                "INSERT INTO quote_decisions (decision_id, rfq_id, quote_id, "
                "combo_market_ticker, game_id, decision, reason, model_fair, "
                "book_fair, blended_fair, yes_bid, no_bid, observed_at) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                [str(_uuid.uuid4()), "r-vd", "qid-vd", "COMBO-VD", "gVD",
                 d, None, None, None, None, None, None, now])

    polled = []

    class Src:
        def poll(self):
            polled.append(1)
            return [{"id": "r-new", "market_ticker": "COMBO-XX", "contracts": 1}]

        def get_market(self, t):
            raise AssertionError("must not reach get_market when halted")

    class GW:
        def submit_quote(self, *a):
            raise AssertionError("must not submit when halted")

    main._discovery_tick(Src(), GW(), dry_run=False)

    # poll should NOT have been called — halt fires before poll.
    assert polled == [], "discovery must short-circuit before polling on void-rate halt"
    with db.connect(read_only=True) as con:
        d = con.execute(
            "SELECT decision FROM quote_decisions WHERE decision='halted_high_void_rate'"
        ).fetchone()
    assert d is not None, "halted_high_void_rate decision must be logged"


# ---------------------------------------------------------------------------
# H4 — per-creator fill halt: a counterparty with N+ recent fills gets skipped.
# ---------------------------------------------------------------------------
def test_discovery_skips_creator_with_too_many_fills(monkeypatch, tmp_path):
    from datetime import datetime, timezone
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "creator_halt.duckdb")
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
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, m_: True)

    # Pre-seed: creator 'abc' has 10 fills (>= PER_CREATOR_FILL_HALT default).
    # Each fill is tied to a seen_rfqs row with creator_id='abc'.
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        for i in range(10):
            rid = f"r-prior-{i}"
            con.execute(
                "INSERT OR REPLACE INTO seen_rfqs "
                "(rfq_id, market_ticker, in_scope, game_id, legs_json, "
                "first_seen_at, last_decision, creator_id) "
                "VALUES (?,?,?,?,?,?,?,?)",
                [rid, f"COMBO-PR-{i}", True, "gPR", None, now, "confirmed", "abc"])
            con.execute(
                "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
                "game_id, side_held, contracts, price, fee, model_fair_at_quote, "
                "book_fair_at_quote, blended_fair_at_quote, fair_at_confirm, "
                "realized_pnl, filled_at, reconciled) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                [f"fill-{i}", f"qid-{i}", rid, f"COMBO-PR-{i}", "gPR",
                 "yes", 1, 0.5, 0.0, None, 0.5, 0.5, 0.5, None, now, True])

    # In-scope ticker so we reach the halt check.
    legs = [{"market_ticker": "KXMLBSPREAD-XX", "event_ticker": "EVT-XX",
             "side": "yes", "count": 1},
            {"market_ticker": "KXMLBTOTAL-XX", "event_ticker": "EVT-XX",
             "side": "yes", "count": 1}]
    monkeypatch.setattr(main, "_SCOPE_CACHE", {"COMBO-NEW": (True, legs)})
    monkeypatch.setattr(main, "_combo_games", lambda legs: ["gNEW"])

    submit_calls = []

    class GW:
        def submit_quote(self, *a):
            submit_calls.append(a)
            return "qid-x"

    class Src:
        def poll(self):
            # RFQ from the same creator 'abc' that's already over-farmed us.
            return [{"id": "r-new", "market_ticker": "COMBO-NEW",
                     "contracts": 1, "creator_user_id": "abc"}]

        def get_market(self, t):
            return {}

    main._discovery_tick(Src(), GW(), dry_run=False)

    assert submit_calls == [], "must NOT submit to over-farming creator"
    with db.connect(read_only=True) as con:
        d = con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "WHERE rfq_id='r-new' ORDER BY observed_at DESC LIMIT 1").fetchone()
    assert d is not None and d[0] == "skipped" and d[1] == "skipped_creator_halt", (
        f"expected skipped_creator_halt, got {d}")


# ---------------------------------------------------------------------------
# H8 — per-combo exposure cap. Pre-seed fills totalling >= cap on one ticker;
# discovery on that ticker is skipped with reason='per_combo_cap'.
# The per-combo cap now runs AFTER pricing (needs the quote), so stubs for
# _book_fairs / _commence_time are required to reach it.
# Cap is $50 (MAX_COMBO_EXPOSURE_USD default). Pre-seed $51 in reconciled fills.
# ---------------------------------------------------------------------------
def test_discovery_skips_when_combo_exposure_capped(monkeypatch, tmp_path):
    from datetime import datetime, timezone
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "combo_cap.duckdb")
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    monkeypatch.setattr(cfg, "BANKROLL", 500.0)
    monkeypatch.setattr(cfg, "MAX_FILL_EXPOSURE_PCT", 0.10)   # cap = $50
    monkeypatch.setattr(cfg, "MAX_COMBO_EXPOSURE_USD", 50.0)
    import importlib
    importlib.reload(db)
    db.init_database()

    import pandas as pd
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["g"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, m_: True)
    monkeypatch.setattr(main, "_today_fills", lambda: [])
    # Per-combo cap now runs after pricing — need a price stub so pricing runs.
    monkeypatch.setattr(main, "_price_combo", lambda legs: (0.55, 3))
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)

    legs = [{"market_ticker": "KXMLBSPREAD-CC", "event_ticker": "EVT-CC",
             "side": "yes", "count": 1},
            {"market_ticker": "KXMLBTOTAL-CC", "event_ticker": "EVT-CC",
             "side": "yes", "count": 1}]
    monkeypatch.setattr(main, "_SCOPE_CACHE", {"COMBO-CAP": (True, legs)})
    monkeypatch.setattr(main, "_combo_games", lambda legs: ["gCAP"])
    monkeypatch.setattr(main, "_PREV_BOOK_FAIR", {})

    # Pre-seed fills totaling $51 on this ticker (102 contracts × $0.50 = $51.00).
    # reconciled=True → counted as price*contracts = $51, which exceeds cap $50.
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
            "game_id, side_held, contracts, price, fee, model_fair_at_quote, "
            "book_fair_at_quote, blended_fair_at_quote, fair_at_confirm, "
            "realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["fill-cap", "qid-cap", "r-cap", "COMBO-CAP", "gCAP",
             "yes", 102, 0.50, 0.0, None, 0.5, 0.5, 0.5, None, now, True])

    class GW:
        def submit_quote(self, *a):
            raise AssertionError("must not submit when per-combo cap hit")

    class Src:
        def poll(self):
            return [{"id": "r-cap", "market_ticker": "COMBO-CAP", "contracts": 1}]

        def get_market(self, t):
            return {}

    main._discovery_tick(Src(), GW(), dry_run=False)

    with db.connect(read_only=True) as con:
        d = con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "WHERE rfq_id='r-cap' ORDER BY observed_at DESC LIMIT 1").fetchone()
    assert d is not None and d[0] == "skipped" and d[1] == "per_combo_cap", (
        f"expected per_combo_cap, got {d}")


# ---------------------------------------------------------------------------
# H9 — discovery side: a combo with an active cooldown row is skipped with
# reason='in_cooldown'.
# ---------------------------------------------------------------------------
def test_discovery_skips_when_combo_in_cooldown(monkeypatch, tmp_path):
    from datetime import datetime, timedelta, timezone
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "cooldown_disc.duckdb")
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
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, m_: True)
    monkeypatch.setattr(main, "_today_fills", lambda: [])

    legs = [{"market_ticker": "KXMLBSPREAD-CD", "event_ticker": "EVT-CD",
             "side": "yes", "count": 1},
            {"market_ticker": "KXMLBTOTAL-CD", "event_ticker": "EVT-CD",
             "side": "yes", "count": 1}]
    monkeypatch.setattr(main, "_SCOPE_CACHE", {"COMBO-CD": (True, legs)})
    monkeypatch.setattr(main, "_combo_games", lambda legs: ["gCD"])
    monkeypatch.setattr(main, "_price_combo", lambda legs: (0.55, 3))

    now = datetime.now(timezone.utc)
    with db.connect() as con:
        con.execute(
            "INSERT OR REPLACE INTO combo_cooldown (combo_market_ticker, cooled_until) "
            "VALUES (?, ?)", ["COMBO-CD", now + timedelta(seconds=30)])

    class GW:
        def submit_quote(self, *a):
            raise AssertionError("must not submit while in cooldown")

    class Src:
        def poll(self):
            return [{"id": "r-cd", "market_ticker": "COMBO-CD", "contracts": 1}]

        def get_market(self, t):
            return {}

    main._discovery_tick(Src(), GW(), dry_run=False)

    with db.connect(read_only=True) as con:
        d = con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "WHERE rfq_id='r-cd' ORDER BY observed_at DESC LIMIT 1").fetchone()
    assert d is not None and d[0] == "skipped" and d[1] == "in_cooldown", (
        f"expected in_cooldown, got {d}")


# ---------------------------------------------------------------------------
# H9 — confirm side: a successful confirm INSERTs a combo_cooldown row with
# cooled_until ~= now + COMBO_COOLDOWN_SEC.
# ---------------------------------------------------------------------------
def test_confirm_arms_combo_cooldown(monkeypatch, tmp_path):
    from datetime import datetime, timedelta, timezone
    import json
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main
    from kalshi_common import auth_client

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "cooldown_arm.duckdb")
    import importlib
    importlib.reload(db)
    db.init_database()

    legs = [{"market_ticker": "KXMLBSPREAD-AR", "event_ticker": "EVT-AR",
             "side": "yes", "count": 1, "no_sub_title": "-1.5"},
            {"market_ticker": "KXMLBTOTAL-AR", "event_ticker": "EVT-AR",
             "side": "yes", "count": 1, "no_sub_title": "Over 8.5"}]
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            ["qid-arm", "r-arm", "COMBO-ARM", "gAR",
             0.45, 0.45, None, 0.55, 0.55, "open",
             datetime.now(timezone.utc), None])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs "
            "(rfq_id, market_ticker, in_scope, game_id, legs_json, "
            "first_seen_at, last_decision) VALUES (?,?,?,?,?,?,?)",
            ["r-arm", "COMBO-ARM", True, "gAR", json.dumps(legs),
             datetime.now(timezone.utc), "quoted"])

    monkeypatch.setattr(main, "_price_combo", lambda legs: (0.56, 3))
    monkeypatch.setattr(main, "_combo_games", lambda legs: ["g3"])

    def fake_api(method, path, *a, **kw):
        if path.startswith("/communications/quotes/"):
            return 200, {"quote": {"status": "accepted", "accepted_side": "no",
                                   "contracts": 1}}, None
        raise AssertionError(f"unexpected api call: {method} {path}")
    monkeypatch.setattr(auth_client, "api", fake_api)

    class GW:
        def confirm(self, qid):
            return True

        def cancel(self, qid):
            return True

    # DuckDB TIMESTAMP is naive — it stores aware datetimes as naive LOCAL
    # (production stays consistent because reads also come back naive-local).
    # The test asserts the delta between insert-now and cooled_until matches
    # COMBO_COOLDOWN_SEC within a small tolerance, which is tz-invariant.
    before_naive = datetime.now().replace(tzinfo=None)
    main._confirm_tick(GW(), dry_run=False)
    after_naive = datetime.now().replace(tzinfo=None)

    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT cooled_until FROM combo_cooldown WHERE combo_market_ticker='COMBO-ARM'"
        ).fetchone()
    assert row is not None, "combo_cooldown row must be created on successful confirm"
    cooled_until = row[0]  # naive local from DuckDB
    expected_min = before_naive + timedelta(seconds=cfg.COMBO_COOLDOWN_SEC - 1)
    expected_max = after_naive + timedelta(seconds=cfg.COMBO_COOLDOWN_SEC + 1)
    assert expected_min <= cooled_until <= expected_max, (
        f"cooled_until {cooled_until} not within [{expected_min}, {expected_max}]")

