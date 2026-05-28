def test_main_imports():
    from kalshi_mlb_mm import main          # must import without error
    assert hasattr(main, "main_loop")


def test_discovery_tick_skips_out_of_scope(monkeypatch, tmp_path):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.risk as risk
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "t.duckdb")
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    # FIX I4 added a staleness gate at the top of _discovery_tick; bypass it here
    # so this test continues to exercise the out-of-scope logging path.
    monkeypatch.setattr(risk, "staleness_ok", lambda *a, **kw: True)
    import importlib, kalshi_mlb_mm.db as db
    importlib.reload(db)
    db.init_database()
    from kalshi_mlb_mm import main

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

    # We need the staleness gate and kill-switch to pass.
    # Bypass staleness gate: patch staleness_ok to always return True.
    monkeypatch.setattr(risk, "staleness_ok", lambda *a, **kw: True)
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
    monkeypatch.setattr(main, "_book_fairs", lambda gid, sl, tl: {"DK": 0.55, "FD": 0.55})

    import pandas as pd
    # Provide a minimal samples df so blended_fair won't short-circuit.
    fake_samples = pd.DataFrame({
        "sim_idx": range(100),
        "home_margin": [1.0] * 100,
        "total_final_score": [9.0] * 100,
        "home_margin_f5": [0.5] * 100,
        "total_f5": [4.5] * 100,
    })
    monkeypatch.setattr(main, "_SAMPLES", {"game1": fake_samples})

    # Patch blended_fair to return a stable value in range.
    import kalshi_mlb_mm.fairs as fairs_mod
    monkeypatch.setattr(fairs_mod, "blended_fair",
                        lambda legs, gid, samples, bf: (0.55, 0.55, 0.55))

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
    monkeypatch.setattr(risk, "staleness_ok", lambda *a, **kw: True)
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
