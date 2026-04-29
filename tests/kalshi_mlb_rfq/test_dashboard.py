from kalshi_mlb_rfq import dashboard, db


def test_dashboard_runs_without_error(tmp_path, monkeypatch, capsys):
    p = tmp_path / "test.duckdb"
    monkeypatch.setattr(db, "DB_PATH", p)
    db.init_database()
    dashboard.main()
    out = capsys.readouterr().out
    assert "Kalshi MLB RFQ" in out
