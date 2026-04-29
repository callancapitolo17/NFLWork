import duckdb
import pytest

from kalshi_mlb_rfq import db


@pytest.fixture
def tmpdb(tmp_path, monkeypatch):
    p = tmp_path / "test.duckdb"
    monkeypatch.setattr(db, "DB_PATH", p)
    db.init_database()
    return p


def test_init_creates_all_tables(tmpdb):
    con = duckdb.connect(str(tmpdb), read_only=True)
    try:
        names = {row[0] for row in con.execute(
            "SELECT table_name FROM information_schema.tables "
            "WHERE table_schema='main'"
        ).fetchall()}
    finally:
        con.close()
    assert names == {
        "combo_cache", "live_rfqs", "quote_log", "fills",
        "positions", "sessions", "combo_cooldown", "reference_lines",
    }


def test_session_round_trip(tmpdb):
    sid = db.start_session(pid=99, dry_run=True, version="0.1.0")
    assert sid
    db.end_session(sid)
    con = duckdb.connect(str(tmpdb), read_only=True)
    try:
        ended_at = con.execute(
            "SELECT ended_at FROM sessions WHERE session_id=?", [sid]
        ).fetchone()[0]
    finally:
        con.close()
    assert ended_at is not None
