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


# ---------------------------------------------------------------------------
# Task 3: on-demand coverage — fetch success rate + latency
# ---------------------------------------------------------------------------

_EVENT_SEQ = iter(range(10_000))


def _insert_event(research_path, event_type, ts, payload):
    import json
    con = duckdb.connect(research_path)
    try:
        con.execute(
            "INSERT INTO events (event_id, session_id, event_type, ts, "
            "ticker, rfq_id, quote_id, payload) VALUES (?,?,?,?,?,?,?,?)",
            [f"e{next(_EVENT_SEQ)}", "s-test", event_type, ts,
             None, None, None, json.dumps(payload)])
    finally:
        con.close()


def test_on_demand_stats(monkeypatch, tmp_path):
    """Request→result pairing on leg_set_hash: success rate + latencies."""
    from kalshi_mlb_mm import report
    state, research, market = _setup_dbs(monkeypatch, tmp_path)
    t0 = NOW - timedelta(hours=2)

    for h in ("h1", "h2", "h3"):
        _insert_event(research, "on_demand_requested", t0,
                      {"leg_set_hash": h, "game_id": "g1", "n_legs": 3})
    _insert_event(research, "on_demand_result", t0 + timedelta(seconds=5),
                  {"leg_set_hash": "h1",
                   "books": {"draftkings": {"fair": 0.2, "route": "exact",
                                            "latency_sec": 1.2},
                             "fanduel": {"fair": 0.21, "route": "transfer",
                                         "latency_sec": 2.0}}})
    _insert_event(research, "on_demand_result", t0 + timedelta(seconds=3),
                  {"leg_set_hash": "h2",
                   "books": {"draftkings": {"fair": 0.4, "route": "exact",
                                            "latency_sec": 0.8}}})

    stats = report.on_demand_stats(research, NOW - timedelta(days=7))
    assert stats["flights_requested"] == 3
    assert stats["hashes_requested"] == 3
    assert stats["hashes_landed"] == 2
    assert stats["success_rate"] == pytest.approx(2 / 3)
    assert stats["wall_p50"] == pytest.approx(4.0)
    per_book = {b: (n, p50) for b, n, p50, _p95 in stats["per_book"]}
    assert per_book["draftkings"] == (2, pytest.approx(1.0))
    assert per_book["fanduel"] == (1, pytest.approx(2.0))


def test_on_demand_stats_empty(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    state, research, market = _setup_dbs(monkeypatch, tmp_path)
    stats = report.on_demand_stats(research, NOW - timedelta(days=7))
    assert stats["flights_requested"] == 0
    assert stats["success_rate"] is None
    assert stats["per_book"] == []


# ---------------------------------------------------------------------------
# Task 4: quotable universe (market DB, naive-UTC timestamps)
# ---------------------------------------------------------------------------

def _insert_sgp_odds(market_path, game_id, combo, bookmaker, fetch_time,
                     spread_line, total_line):
    con = duckdb.connect(market_path)
    try:
        con.execute(
            "INSERT INTO mlb_sgp_odds (game_id, combo, period, bookmaker, "
            "sgp_decimal, sgp_american, fetch_time, source, spread_line, "
            "total_line) VALUES (?,?,?,?,?,?,?,?,?,?)",
            [game_id, combo, "full_game", bookmaker, 2.5, 150, fetch_time,
             "scraper", spread_line, total_line])
    finally:
        con.close()


def test_universe_stats(monkeypatch, tmp_path):
    """Per-day combos passing MIN_AGREEING_BOOKS + per-book participation."""
    from kalshi_mlb_mm import report
    state, research, market = _setup_dbs(monkeypatch, tmp_path)
    now_naive = NOW.replace(tzinfo=None)
    day1 = now_naive - timedelta(days=1)
    day2 = now_naive

    # Combo X (spread x total) both days; combo Y (ml x total, NULL spread)
    # day1 only with a single book, so it never passes consensus.
    for book in ("draftkings", "fanduel", "novig"):
        _insert_sgp_odds(market, "g1", "sxt", book, day1, -1.5, 8.5)
    _insert_sgp_odds(market, "g1", "mxt", "draftkings", day1, None, 9.5)
    for book in ("draftkings", "fanduel"):
        _insert_sgp_odds(market, "g1", "sxt", book, day2, -1.5, 8.5)

    stats = report.universe_stats(market, now_naive - timedelta(days=7),
                                  now_naive)
    per_day = {d.strftime("%Y-%m-%d") if hasattr(d, "strftime") else str(d):
               (passing, total) for d, passing, total in stats["per_day"]}
    assert per_day[day1.strftime("%Y-%m-%d")] == (1, 2)
    assert per_day[day2.strftime("%Y-%m-%d")] == (1, 1)
    per_book = {b: (days, combos) for b, days, combos, _age in
                stats["per_book"]}
    assert per_book["draftkings"] == (2, 2)
    assert per_book["fanduel"] == (2, 1)
    assert per_book["novig"] == (1, 1)


def test_universe_empty(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    state, research, market = _setup_dbs(monkeypatch, tmp_path)
    now_naive = NOW.replace(tzinfo=None)
    stats = report.universe_stats(market, now_naive - timedelta(days=7),
                                  now_naive)
    assert stats["per_day"] == []
    assert stats["per_book"] == []


# ---------------------------------------------------------------------------
# Task 5: staleness — book-data age at quote time / quote age at accept
# ---------------------------------------------------------------------------

def test_staleness_stats(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    state, research, market = _setup_dbs(monkeypatch, tmp_path)
    t0 = NOW - timedelta(hours=3)

    _insert_event(research, "scrape_done", t0, {"book_counts": {}})
    _insert_event(research, "scrape_done", t0 + timedelta(seconds=120),
                  {"book_counts": {}})
    # ages vs latest prior scrape: 30s and 30s
    _insert_event(research, "quote_priced", t0 + timedelta(seconds=30),
                  {"blended_fair": 0.5})
    _insert_event(research, "quote_priced", t0 + timedelta(seconds=150),
                  {"blended_fair": 0.6})
    # a quote before any scrape_done: unknown age, not a crash
    _insert_event(research, "quote_priced", t0 - timedelta(seconds=10),
                  {"blended_fair": 0.7})
    _insert_event(research, "confirm_singles_check",
                  t0 + timedelta(seconds=200), {"quote_age_sec": 5.0})
    _insert_event(research, "confirm_singles_check",
                  t0 + timedelta(seconds=300), {"quote_age_sec": 45.0})

    stats = report.staleness_stats(research, NOW - timedelta(days=7))
    assert stats["quote_age"]["n"] == 2
    assert stats["quote_age"]["p50"] == pytest.approx(30.0)
    assert stats["quote_age_unknown"] == 1
    assert stats["accept_age"]["n"] == 2
    assert stats["accept_age"]["p50"] == pytest.approx(25.0)
    assert stats["accept_age"]["max"] == pytest.approx(45.0)


def test_staleness_empty(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    state, research, market = _setup_dbs(monkeypatch, tmp_path)
    stats = report.staleness_stats(research, NOW - timedelta(days=7))
    assert stats["quote_age"]["n"] == 0
    assert stats["accept_age"]["n"] == 0
