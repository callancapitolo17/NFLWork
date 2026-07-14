import datetime
from unabated_edge import config
from unabated_edge.maker import store
from unabated_edge.storage import connect


def test_maker_config_defaults():
    assert config.MAKER_MODE == "off"
    assert config.MATCH_CAP_PCT == 0.40
    assert config.GLOBAL_CAP_PCT == 0.75
    assert config.MAX_QUOTE_PCT == 0.30
    assert config.QUOTE_PULL_MIN == 3.0


def test_store_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")
    store.init()
    now = datetime.datetime.now(datetime.timezone.utc)
    store.log_quote(now, "soccer", 1, "T-O25", "yes", "rest", 0.40, 30, 0.42, 0.02, False, "quote", "oid1")
    store.log_fill(now, "soccer", 1, "oid1", "T-O25", "yes", 0.40, 5.0, 0.01, -2.0, "tr1")
    store.log_fill(now, "soccer", 1, "oid1", "T-O25", "yes", 0.40, 5.0, 0.01, -2.0, "tr1")  # dup trade_id
    store.log_ledger(now, "soccer", 1, -2.0, [-2.0] * 11, 4)
    store.log_settlement(now, "soccer", "T-O25", 3.5)
    with connect(config.MAKER_DB_PATH, read_only=True) as c:
        assert c.execute("SELECT count(*) FROM maker_quotes").fetchone()[0] == 1
        assert c.execute("SELECT count(*) FROM maker_fills").fetchone()[0] == 1   # deduped
        assert c.execute("SELECT count(*) FROM ledger_snapshots").fetchone()[0] == 1
        assert c.execute("SELECT settled_pnl FROM maker_pnl").fetchone()[0] == 3.5
