import datetime
from unabated_edge import storage, config


def test_snapshot_flag_and_research(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MARKET_DB_PATH", tmp_path / "m.duckdb")
    monkeypatch.setattr(config, "RESEARCH_DB_PATH", tmp_path / "r.duckdb")
    storage.init()
    now = datetime.datetime.now(datetime.timezone.utc)
    storage.snapshot_lines("soccer", [{"ts": now, "event_id": 9, "market_source_id": 7, "bet_type": "bt1", "side": "1", "price": -150.0, "points": None}])
    storage.log_flagged({"ts": now, "sport": "soccer", "event_id": 9, "market_ticker": "X", "outcome": "home",
                         "fair_prob": 0.6, "yes_ask": 0.5, "ev_pct": 0.1, "kelly_contracts": 3, "dry_run": True})
    storage.emit("candidate_priced", "soccer", event_id=9, ev_pct=0.1)
    storage.flush()
    with storage.connect(config.MARKET_DB_PATH, read_only=True) as c:
        assert c.execute("SELECT count(*) FROM line_snapshots").fetchone()[0] == 1
        assert c.execute("SELECT count(*) FROM flagged_edges").fetchone()[0] == 1
    with storage.connect(config.RESEARCH_DB_PATH, read_only=True) as c:
        assert c.execute("SELECT count(*) FROM research_events").fetchone()[0] == 1
