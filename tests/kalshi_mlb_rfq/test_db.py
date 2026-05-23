from datetime import datetime, timezone

import duckdb
import pytest

from kalshi_mlb_rfq import db


@pytest.fixture
def tmpdb(tmp_path, monkeypatch):
    p = tmp_path / "test.duckdb"
    monkeypatch.setattr(db, "DB_PATH", p)
    db.init_database()
    return p


# ---- v2 side-column migration ---------------------------------------------

PRE_V2_SCHEMA = """
CREATE TABLE positions (
    combo_market_ticker  VARCHAR PRIMARY KEY,
    game_id              VARCHAR NOT NULL,
    net_contracts        DOUBLE NOT NULL,
    weighted_price       DOUBLE NOT NULL,
    legs_json            VARCHAR NOT NULL,
    updated_at           TIMESTAMP NOT NULL
);
CREATE INDEX idx_positions_game ON positions(game_id);

CREATE TABLE combo_cooldown (
    leg_set_hash  VARCHAR PRIMARY KEY,
    game_id       VARCHAR NOT NULL,
    cooled_until  TIMESTAMP NOT NULL,
    reason        VARCHAR
);
"""


@pytest.fixture
def pre_v2_db(tmp_path, monkeypatch):
    """Builds a pre-v2 DB with sample rows, then points db.DB_PATH at it.
    Caller invokes db.init_database() to trigger the migration."""
    p = tmp_path / "test.duckdb"
    con = duckdb.connect(str(p))
    try:
        con.execute(PRE_V2_SCHEMA)
        now = datetime.now(timezone.utc)
        con.execute(
            "INSERT INTO positions VALUES (?, ?, ?, ?, ?, ?)",
            ["COMBO-X", "game-1", 5.0, 0.45, "[]", now],
        )
        con.execute(
            "INSERT INTO combo_cooldown VALUES (?, ?, ?, ?)",
            ["hash-abc", "game-1", now, "post_accept"],
        )
    finally:
        con.close()
    monkeypatch.setattr(db, "DB_PATH", p)
    return p


def test_v2_migration_adds_side_with_yes_backfill(pre_v2_db):
    db.init_database()
    con = duckdb.connect(str(pre_v2_db), read_only=True)
    try:
        # positions: side column exists, existing row backfilled to 'yes'
        rows = con.execute(
            "SELECT combo_market_ticker, side, net_contracts "
            "FROM positions ORDER BY combo_market_ticker"
        ).fetchall()
        assert rows == [("COMBO-X", "yes", 5.0)]

        # PK is now (combo_market_ticker, side) — we can insert a NO row
        # for the same ticker without violating uniqueness
        con.close()
        con = duckdb.connect(str(pre_v2_db))
        con.execute(
            "INSERT INTO positions VALUES (?, ?, ?, ?, ?, ?, ?)",
            ["COMBO-X", "no", "game-1", 3.0, 0.55, "[]",
             datetime.now(timezone.utc)],
        )
        sides = {r[0] for r in con.execute(
            "SELECT side FROM positions WHERE combo_market_ticker='COMBO-X'"
        ).fetchall()}
        assert sides == {"yes", "no"}

        # combo_cooldown: same shape
        cooldown_sides = {r[0] for r in con.execute(
            "SELECT side FROM combo_cooldown WHERE leg_set_hash='hash-abc'"
        ).fetchall()}
        assert cooldown_sides == {"yes"}
    finally:
        con.close()


def test_v2_migration_is_idempotent(pre_v2_db):
    db.init_database()
    db.init_database()  # second call must not fail or duplicate rows
    con = duckdb.connect(str(pre_v2_db), read_only=True)
    try:
        n = con.execute("SELECT COUNT(*) FROM positions").fetchone()[0]
        assert n == 1  # original row, not duplicated
    finally:
        con.close()


def test_v2_migration_backs_up_db(pre_v2_db):
    db.init_database()
    # Backup file should sit next to the DB with .bak.v2_side_columns.<ts> suffix
    backups = list(pre_v2_db.parent.glob(f"{pre_v2_db.name}.bak.v2_side_columns.*"))
    assert len(backups) == 1
    # Backup is the *pre-v2* schema — no side column
    con = duckdb.connect(str(backups[0]), read_only=True)
    try:
        cols = {r[0] for r in con.execute(
            "SELECT column_name FROM information_schema.columns "
            "WHERE table_name='positions'"
        ).fetchall()}
    finally:
        con.close()
    assert "side" not in cols


def test_v2_migration_adds_quote_log_hedge_columns(pre_v2_db):
    db.init_database()
    con = duckdb.connect(str(pre_v2_db), read_only=True)
    try:
        cols = {r[0] for r in con.execute(
            "SELECT column_name FROM information_schema.columns "
            "WHERE table_name='quote_log'"
        ).fetchall()}
    finally:
        con.close()
    expected = {
        "chosen_side", "ev_yes_pct", "ev_no_pct",
        "hedge_added", "hedge_original_side", "hedge_original_price",
        "hedge_new_price", "hedge_current_fair", "hedge_projected_net",
    }
    assert expected.issubset(cols)


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
