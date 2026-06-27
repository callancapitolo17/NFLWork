import importlib

def test_init_and_session(tmp_path, monkeypatch):
    import kalshi_mlb_mm.config as cfg
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "t.duckdb")
    import kalshi_mlb_mm.db as db
    importlib.reload(db)
    db.init_database()
    sid = db.start_session(pid=1, dry_run=True)
    db.end_session(sid)
    with db.connect(read_only=True) as con:
        row = con.execute("SELECT dry_run, ended_at FROM sessions WHERE session_id=?", [sid]).fetchone()
    assert row[0] is True and row[1] is not None
