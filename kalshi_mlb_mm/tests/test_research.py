"""Tests for kalshi_mlb_mm/research.py — the research event firehose."""
import duckdb
import pytest

import kalshi_mlb_mm.research as research


@pytest.fixture(autouse=True)
def isolated_research_db(tmp_path, monkeypatch):
    """Point research at a temp DB and reset module state per test."""
    db = tmp_path / "research.duckdb"
    monkeypatch.setattr(research.config, "RESEARCH_DB_PATH", db)
    monkeypatch.setattr(research.config, "RESEARCH_BUFFER_MAX", 5000)
    monkeypatch.setattr(research.config, "RESEARCH_FLUSH_WARN_RATE_LIMIT_SEC", 60)
    research._BUFFER.clear()
    research._SESSION_ID = None
    research._LAST_FLUSH_WARN = 0.0
    # Keep FLUSH_WARN_INTERVAL_SEC in sync with the monkeypatched config value
    research.FLUSH_WARN_INTERVAL_SEC = 60.0
    yield db


def test_init_creates_events_table(isolated_research_db):
    research.init_research_db()
    con = duckdb.connect(str(isolated_research_db), read_only=True)
    cols = {r[0] for r in con.execute(
        "SELECT column_name FROM information_schema.columns "
        "WHERE table_name='events'").fetchall()}
    con.close()
    assert {"event_id", "session_id", "event_type", "ts",
            "ticker", "rfq_id", "quote_id", "payload"}.issubset(cols)


def test_init_is_idempotent(isolated_research_db):
    research.init_research_db()
    research.init_research_db()  # must not raise


def test_emit_appends_to_buffer(isolated_research_db):
    research.set_session_id("sess-1")
    research.emit("rfq_received", ticker="KXMLB-X", rfq_id="r1",
                  payload=dict(rfq_keys=["id", "market_ticker"], rfq_raw={}))
    assert len(research._BUFFER) == 1
    ev = research._BUFFER[0]
    assert ev["event_type"] == "rfq_received"
    assert ev["session_id"] == "sess-1"
    assert ev["ticker"] == "KXMLB-X"
    assert ev["rfq_id"] == "r1"
    assert ev["payload"]["rfq_keys"] == ["id", "market_ticker"]


def test_emit_never_raises_on_bad_payload(isolated_research_db):
    class Weird:
        pass
    research.emit("rfq_received", payload=dict(obj=Weird()))  # must not raise
    assert len(research._BUFFER) == 1


def test_emit_respects_buffer_cap(isolated_research_db, monkeypatch):
    monkeypatch.setattr(research.config, "RESEARCH_BUFFER_MAX", 3)
    for i in range(5):
        research.emit("quote_priced", payload=dict(i=i))
    assert len(research._BUFFER) == 3
    # Only the 3 newest remain
    assert [e["payload"]["i"] for e in research._BUFFER] == [2, 3, 4]


def test_flush_writes_batch_and_clears(isolated_research_db):
    research.init_research_db()
    research.set_session_id("sess-1")
    research.emit("rfq_received", ticker="T1", rfq_id="r1",
                  payload=dict(rfq_keys=[], rfq_raw={}))
    research.emit("quote_priced", ticker="T1", rfq_id="r1",
                  payload=dict(blended_fair=0.55))
    n_flushed = research.flush()
    assert n_flushed == 2
    assert research._BUFFER == []
    con = duckdb.connect(str(isolated_research_db), read_only=True)
    n = con.execute("SELECT count(*) FROM events").fetchone()[0]
    fair = con.execute(
        "SELECT payload->>'blended_fair' FROM events "
        "WHERE event_type='quote_priced'").fetchone()[0]
    con.close()
    assert n == 2
    assert float(fair) == 0.55


def test_flush_returns_zero_on_empty_buffer(isolated_research_db):
    n = research.flush()
    assert n == 0


def test_flush_swallows_db_error_and_keeps_running(isolated_research_db,
                                                   monkeypatch, caplog):
    research.set_session_id("sess-1")
    research.emit("rfq_received", payload=dict(rfq_keys=[], rfq_raw={}))

    def boom(*a, **k):
        raise duckdb.IOException("locked")

    monkeypatch.setattr(research, "_connect", boom)
    result = research.flush()   # MUST NOT raise
    assert result == 0
    assert "research flush failed" in caplog.text.lower()
    assert len(research._BUFFER) == 1   # retained for retry on next tick


def test_flush_never_raises_on_unserializable_payload(isolated_research_db):
    """A payload object whose __repr__ raises must not escape flush() into
    the trading loop (row-building incl. json.dumps is inside the try)."""
    research.init_research_db()
    research.set_session_id("sess-1")

    class BrokenRepr:
        def __repr__(self):
            raise RuntimeError("boom")

    research.emit("rfq_received", payload=dict(obj=BrokenRepr()))
    research.flush()   # MUST NOT raise despite the broken repr


def test_set_session_id_stamps_subsequent_emits(isolated_research_db):
    research.init_research_db()
    research.set_session_id("sess-abc")
    research.emit("decision", payload=dict(decision="quoted"))
    research.flush()
    con = duckdb.connect(str(isolated_research_db), read_only=True)
    sid = con.execute("SELECT session_id FROM events LIMIT 1").fetchone()[0]
    con.close()
    assert sid == "sess-abc"


def test_prune_deletes_only_old_rows(isolated_research_db):
    research.init_research_db()
    con = duckdb.connect(str(isolated_research_db))
    con.execute(
        "INSERT INTO events (event_id, event_type, ts, payload) VALUES "
        "('old', 't', now() - INTERVAL 100 DAY, '{}'), "
        "('new', 't', now() - INTERVAL 1 DAY, '{}')")
    con.close()
    deleted = research.prune_research(days=90)
    assert deleted == 1
    con = duckdb.connect(str(isolated_research_db), read_only=True)
    ids = {r[0] for r in con.execute("SELECT event_id FROM events").fetchall()}
    con.close()
    assert ids == {"new"}
