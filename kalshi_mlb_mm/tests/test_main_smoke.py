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
    monkeypatch.setattr(main, "_resolve_game_and_lines",
                        lambda ticker, legs: ("game1", -1.5, 8.5))
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    # Make tipoff_ok pass (commence_time is None → normally fails; override).
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)
    monkeypatch.setattr(main, "_book_fairs", lambda gid, sl, tl: {"DK": 0.55, "FD": 0.55, "PX": 0.55})

    # Patch blended_fair to return a stable value in range (signature now
    # (legs, gid, bf) → (book_med, blended), model removed in v1 hardening).
    import kalshi_mlb_mm.fairs as fairs_mod
    monkeypatch.setattr(fairs_mod, "blended_fair",
                        lambda legs, gid, bf: (0.55, 0.55))

    # Scope cache: make the ticker appear in-scope with a minimal legs list.
    legs = [{"market_ticker": "KXMLBSPREAD-ABC", "event_ticker": "EVT-ABC",
             "side": "yes", "count": 1},
            {"market_ticker": "KXMLBTOTAL-ABC", "event_ticker": "EVT-ABC",
             "side": "yes", "count": 1}]
    monkeypatch.setattr(main, "_SCOPE_CACHE",
                        {"COMBO-1": (True, "game1", legs)})

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
                        {"COMBO-2": (True, "game2", legs)})

    # _resolve_game_and_lines must return a valid game_id so we reach the cap check.
    monkeypatch.setattr(main, "_resolve_game_and_lines",
                        lambda ticker, legs: ("game2", -1.5, 8.5))

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

    # Book_fair_at_q was 0.40; now consensus median is 0.45 (drift 0.05 > 0.03).
    monkeypatch.setattr(main, "_book_fairs",
                        lambda gid, sl, tl: {"DK": 0.44, "FD": 0.45, "PX": 0.46})

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
            "INSERT OR REPLACE INTO seen_rfqs VALUES (?,?,?,?,?,?,?)",
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
            "INSERT OR REPLACE INTO seen_rfqs VALUES (?,?,?,?,?,?,?)",
            ["r-nofresh", "COMBO-G2", True, "g2", json.dumps(legs),
             datetime.now(timezone.utc), "quoted"])

    # The quote came back 'accepted' from the API.
    def fake_api(method, path, *a, **kw):
        if path.startswith("/communications/quotes/"):
            return 200, {"quote": {"status": "accepted", "accepted_side": "yes",
                                   "contracts": 1}}, None
        raise AssertionError(f"unexpected api call: {method} {path}")
    monkeypatch.setattr(auth_client, "api", fake_api)

    # No fresh books available — _book_fairs returns {}.
    monkeypatch.setattr(main, "_book_fairs", lambda gid, sl, tl: {})

    confirm_calls = []

    class GW:
        def confirm(self, qid):
            confirm_calls.append(qid)
            return True

        def cancel(self, qid):
            return True

    main._confirm_tick(GW(), dry_run=False)

    assert confirm_calls == [], "confirm must NOT be called when no fresh books"
    with db.connect(read_only=True) as con:
        st = con.execute(
            "SELECT status FROM live_quotes WHERE quote_id='qid-nofresh'").fetchone()[0]
        n_fills = con.execute("SELECT COUNT(*) FROM fills").fetchone()[0]
        decision = con.execute(
            "SELECT decision FROM quote_decisions WHERE quote_id='qid-nofresh' "
            "ORDER BY observed_at DESC LIMIT 1").fetchone()
    assert st == "voided", f"expected 'voided', got '{st}'"
    assert n_fills == 0, "no fills row should be written when voided"
    assert decision is not None and decision[0] == "voided_no_fresh_books"


# ---------------------------------------------------------------------------
# v1 hardening (Change E) — H6 position reconciliation: when /portfolio/positions
# reports a different side or size than the bot expected from the quote response,
# trust Kalshi and record the actual side/size in the fills row.
# ---------------------------------------------------------------------------
def test_confirm_reconciles_position_mismatch(monkeypatch, capsys, tmp_path):
    from datetime import datetime, timezone
    import json
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_mlb_mm import main
    from kalshi_common import auth_client

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "reconcile.duckdb")
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
            # Quote priced so 'yes' side is accepted and last-look will pass.
            # yes_bid=0.45, no_bid=0.45 → if accepted_side='no' (taker took our no_bid),
            # side_held becomes 'yes' and price = (1 - yb) = 0.55. To pass last-look
            # with current_fair=0.60 and prev_fair=0.55 (drift 0.05 fails tol 0.02);
            # so make current_fair=0.56 and prev_fair=0.55 — drift 0.01 < 0.02, and
            # p_yes=0.56 > price 0.55 + fee → +EV.
            ["qid-recon", "r-recon", "COMBO-G3", "g3",
             0.45, 0.45, None, 0.55, 0.55, "open",  # last 0.55 == blended_fair (prev_fair)
             datetime.now(timezone.utc), None])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs VALUES (?,?,?,?,?,?,?)",
            ["r-recon", "COMBO-G3", True, "g3", json.dumps(legs),
             datetime.now(timezone.utc), "quoted"])

    # Stub _book_fairs to return a current consensus near prev_fair so last-look passes.
    monkeypatch.setattr(main, "_book_fairs",
                        lambda gid, sl, tl: {"dk": 0.56, "fd": 0.56, "px": 0.56})

    def fake_api(method, path, *a, **kw):
        if path.startswith("/communications/quotes/"):
            # Quote accepted on the 'no' side → side_held='yes', contracts=5 expected.
            return 200, {"quote": {"status": "accepted", "accepted_side": "no",
                                   "contracts": 5}}, None
        if path.startswith("/portfolio/positions"):
            # But Kalshi actually shows position_fp = -3.00 (short YES = long NO, 3 ct).
            return 200, {"market_positions": [
                {"ticker": "COMBO-G3", "position_fp": "-3.00"}
            ]}, None
        raise AssertionError(f"unexpected api call: {method} {path}")
    monkeypatch.setattr(auth_client, "api", fake_api)

    class GW:
        def confirm(self, qid):
            return True

        def cancel(self, qid):
            return True

    main._confirm_tick(GW(), dry_run=False)

    # Warning was printed on stdout
    captured = capsys.readouterr().out
    assert "position_mismatch" in captured, f"expected warning, got: {captured!r}"

    # Fills row reflects the corrected side ('no') and size (3) — Kalshi-trusted.
    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT side_held, contracts FROM fills WHERE quote_id='qid-recon'"
        ).fetchone()
    assert row is not None, "expected a fills row to be written after confirm"
    side, contracts = row
    assert side == "no", f"expected reconciled side='no', got '{side}'"
    assert int(contracts) == 3, f"expected reconciled contracts=3, got {contracts}"

