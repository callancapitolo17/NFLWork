def test_main_imports():
    from kalshi_mlb_mm import main          # must import without error
    assert hasattr(main, "main_loop")


def test_discovery_tick_skips_out_of_scope(monkeypatch, tmp_path):
    import kalshi_mlb_mm.config as cfg
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "t.duckdb")
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
