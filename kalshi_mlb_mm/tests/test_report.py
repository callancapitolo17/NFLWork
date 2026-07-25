"""Regression tests for the daily state-of-the-maker report (issue #14).

These fail on pre-#14 code by construction: kalshi_mlb_mm/report.py did not
exist. Fixture DBs are seeded in tmp_path with the production schema helpers
(db.init_database / research.init_research_db / mlb_sgp.db.ensure_table) so
the tests exercise the real DDL, never a hand-rolled copy.
"""
import importlib
from datetime import datetime, timedelta, timezone

import duckdb
import pytest


NOW = datetime(2026, 7, 24, 18, 0, 0, tzinfo=timezone.utc)


# ---------------------------------------------------------------------------
# Harness
# ---------------------------------------------------------------------------

def _setup_dbs(monkeypatch, tmp_path):
    """Create the three maker DBs (state/research/market) fresh in tmp_path.

    Returns (state_path, research_path, market_path) as strings.
    """
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.research as research
    import mlb_sgp.db as sgp_db

    state = tmp_path / "state.duckdb"
    research_path = tmp_path / "research.duckdb"
    market = tmp_path / "market.duckdb"

    monkeypatch.setattr(cfg, "DB_PATH", state)
    importlib.reload(db)
    db.init_database()

    monkeypatch.setattr(cfg, "RESEARCH_DB_PATH", research_path)
    research.init_research_db()

    sgp_db.ensure_table(db_path=str(market))
    return str(state), str(research_path), str(market)


# ---------------------------------------------------------------------------
# Task 1: plumbing — fresh-DB safety, missing files, lock sentinel
# ---------------------------------------------------------------------------

def test_empty_dbs_report_runs(monkeypatch, tmp_path):
    """Fresh schemas, zero rows: the report renders instead of crashing."""
    from kalshi_mlb_mm import report
    state, research, market = _setup_dbs(monkeypatch, tmp_path)
    out = report.build_report(state, research, market, now=NOW)
    assert "State of the Maker" in out
    assert "no fills yet" in out


def test_missing_db_files_report_runs(tmp_path):
    """DB files that don't exist degrade to a note, never an exception."""
    from kalshi_mlb_mm import report
    out = report.build_report(
        str(tmp_path / "nope_state.duckdb"),
        str(tmp_path / "nope_research.duckdb"),
        str(tmp_path / "nope_market.duckdb"),
        now=NOW,
    )
    assert "State of the Maker" in out
    assert "not found" in out


def test_read_missing_table_returns_empty(monkeypatch, tmp_path):
    """A DB that exists but lacks the table (fresh install) reads as empty."""
    from kalshi_mlb_mm import report
    path = tmp_path / "bare.duckdb"
    duckdb.connect(str(path)).close()
    rows = report._read(str(path), "SELECT decision FROM quote_decisions")
    assert rows == []


def test_read_locked_returns_sentinel(monkeypatch, tmp_path):
    """A held write lock surfaces as LOCKED, not as empty data or a crash."""
    from kalshi_mlb_mm import report
    path = tmp_path / "locked.duckdb"
    writer = duckdb.connect(str(path))  # exclusive write lock
    try:
        monkeypatch.setattr(report, "_RETRY_SLEEP_SEC", 0.0)
        rows = report._read(str(path), "SELECT 1")
        assert rows is report.LOCKED
    finally:
        writer.close()


def test_read_missing_file_returns_missing(tmp_path):
    from kalshi_mlb_mm import report
    rows = report._read(str(tmp_path / "ghost.duckdb"), "SELECT 1")
    assert rows is report.MISSING


# ---------------------------------------------------------------------------
# Task 2: RFQ funnel + path split (grid / on-demand / out-of-scope)
# ---------------------------------------------------------------------------

def _insert_rfq(db, rfq_id, in_scope, last_decision, first_seen_at):
    with db.connect() as con:
        con.execute(
            "INSERT INTO seen_rfqs (rfq_id, market_ticker, in_scope, game_id, "
            "legs_json, first_seen_at, last_decision, creator_id) "
            "VALUES (?,?,?,?,?,?,?,?)",
            [rfq_id, f"T-{rfq_id}", in_scope, "g1", "[]",
             first_seen_at, last_decision, None])


def _insert_decision(db, rfq_id, decision, reason, observed_at,
                     blended_fair=None, yes_bid=None, no_bid=None):
    with db.connect() as con:
        con.execute(
            "INSERT INTO quote_decisions (decision_id, rfq_id, quote_id, "
            "combo_market_ticker, game_id, decision, reason, model_fair, "
            "book_fair, blended_fair, yes_bid, no_bid, observed_at) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            [f"d-{rfq_id}-{decision}-{reason}", rfq_id, None, f"T-{rfq_id}",
             "g1", decision, reason, None, None, blended_fair, yes_bid,
             no_bid, observed_at])


def test_funnel_math(monkeypatch, tmp_path):
    """Seen → in-scope → quoted → accepted → confirmed + path shares."""
    from kalshi_mlb_mm import report
    import kalshi_mlb_mm.db as db
    state, research, market = _setup_dbs(monkeypatch, tmp_path)
    recent = NOW - timedelta(hours=1)
    old = NOW - timedelta(days=8)  # outside every window

    _insert_rfq(db, "r0", False, "out_of_scope", old)
    _insert_rfq(db, "r1", False, "out_of_scope", recent)
    _insert_rfq(db, "r2", False, "out_of_scope_lone_single", recent)
    _insert_rfq(db, "r3", True, "quoted", recent)      # on-demand path
    _insert_rfq(db, "r4", True, "confirmed", recent)   # grid path
    _insert_rfq(db, "r5", True, "no_fair", recent)     # grid path, skipped

    _insert_decision(db, "r1", "skipped", "out_of_scope", recent)
    _insert_decision(db, "r2", "skipped", "out_of_scope_lone_single", recent)
    _insert_decision(db, "r3", "skipped", "on_demand_pending", recent)
    _insert_decision(db, "r3", "quoted", None, recent)
    _insert_decision(db, "r4", "quoted", None, recent)
    _insert_decision(db, "r4", "confirmed", None, recent)
    _insert_decision(db, "r5", "skipped", "no_fair", recent)

    counts = report.funnel_counts(state, NOW - timedelta(hours=24))
    assert counts["seen"] == 5
    assert counts["in_scope"] == 3
    assert counts["quoted"] == 2
    assert counts["accepted"] == 1
    assert counts["confirmed"] == 1
    assert counts["filled"] == 0
    assert counts["paths"]["out_of_scope"] == {"rfqs": 2, "quoted": 0}
    assert counts["paths"]["on_demand"] == {"rfqs": 1, "quoted": 1}
    assert counts["paths"]["grid"] == {"rfqs": 2, "quoted": 1}
    reasons = {(d, r): n for d, r, n in counts["decisions"]}
    assert reasons[("skipped", "no_fair")] == 1
    assert reasons[("skipped", "out_of_scope_lone_single")] == 1


def test_funnel_empty(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    state, research, market = _setup_dbs(monkeypatch, tmp_path)
    counts = report.funnel_counts(state, NOW - timedelta(hours=24))
    assert counts["seen"] == 0
    assert counts["paths"]["grid"] == {"rfqs": 0, "quoted": 0}
