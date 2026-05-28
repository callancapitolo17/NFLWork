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
