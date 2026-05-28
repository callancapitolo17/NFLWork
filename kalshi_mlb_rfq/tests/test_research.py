"""Tests for kalshi_mlb_rfq/research.py — the research event firehose."""
import duckdb
import pytest

import kalshi_mlb_rfq.research as research


@pytest.fixture(autouse=True)
def isolated_research_db(tmp_path, monkeypatch):
    """Point research at a temp DB and reset module state per test."""
    db = tmp_path / "research.duckdb"
    monkeypatch.setattr(research.config, "RESEARCH_DB_PATH", db)
    research._BUFFER.clear()
    research._SESSION_ID = None
    research._LAST_FLUSH_WARN = 0.0
    yield db


def test_init_creates_events_table(isolated_research_db):
    research.init_research_db()
    con = duckdb.connect(str(isolated_research_db), read_only=True)
    cols = {r[0] for r in con.execute(
        "SELECT column_name FROM information_schema.columns "
        "WHERE table_name='events'").fetchall()}
    con.close()
    assert {"event_id", "session_id", "event_type", "ts",
            "game_id", "combo_ticker", "rfq_id", "quote_id",
            "payload"}.issubset(cols)


def test_init_is_idempotent(isolated_research_db):
    research.init_research_db()
    research.init_research_db()  # must not raise


def test_emit_appends_to_buffer(isolated_research_db):
    research.set_session("sess-1")
    research.emit("candidate_evaluated", game_id="g1",
                  combo_ticker="KX-X", model_fair=0.4, edge=0.03)
    assert len(research._BUFFER) == 1
    ev = research._BUFFER[0]
    assert ev["event_type"] == "candidate_evaluated"
    assert ev["session_id"] == "sess-1"
    assert ev["game_id"] == "g1"
    payload = ev["payload"]
    assert payload["model_fair"] == 0.4 and payload["edge"] == 0.03


def test_emit_never_raises_on_bad_payload(isolated_research_db):
    class Weird:
        pass
    research.emit("candidate_evaluated", obj=Weird())  # must not raise
    assert len(research._BUFFER) == 1


def test_emit_respects_buffer_cap(isolated_research_db, monkeypatch):
    monkeypatch.setattr(research.config, "RESEARCH_BUFFER_MAX", 3)
    for i in range(5):
        research.emit("candidate_evaluated", i=i)
    assert len(research._BUFFER) == 3
    assert [e["payload"]["i"] for e in research._BUFFER] == [2, 3, 4]


def test_flush_writes_batch_and_clears(isolated_research_db):
    research.init_research_db()
    research.set_session("sess-1")
    research.emit("candidate_evaluated", game_id="g1", model_fair=0.4)
    research.emit("gate_evaluated", game_id="g1", gate="tipoff")
    research.flush()
    assert research._BUFFER == []
    con = duckdb.connect(str(isolated_research_db), read_only=True)
    n = con.execute("SELECT count(*) FROM events").fetchone()[0]
    fair = con.execute(
        "SELECT payload->>'model_fair' FROM events "
        "WHERE event_type='candidate_evaluated'").fetchone()[0]
    con.close()
    assert n == 2
    assert float(fair) == 0.4


def test_flush_swallows_db_error_and_keeps_running(isolated_research_db,
                                                   monkeypatch, caplog):
    research.set_session("sess-1")
    research.emit("candidate_evaluated", game_id="g1")
    def boom(*a, **k):
        raise duckdb.IOException("locked")
    monkeypatch.setattr(research, "_connect", boom)
    research.flush()   # MUST NOT raise
    assert "research flush failed" in caplog.text.lower()


def test_flush_on_empty_buffer_is_noop(isolated_research_db):
    research.flush()   # must not raise, must not create rows


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
