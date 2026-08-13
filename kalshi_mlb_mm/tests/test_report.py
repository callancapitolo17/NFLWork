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
    """Create the two maker DBs the report reads (state/research) fresh in
    tmp_path (#81: the report no longer touches the market DB — the
    quotable-universe section reads on_demand_result events).

    Returns (state_path, research_path) as strings.
    """
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.research as research

    state = tmp_path / "state.duckdb"
    research_path = tmp_path / "research.duckdb"

    monkeypatch.setattr(cfg, "DB_PATH", state)
    importlib.reload(db)
    db.init_database()

    monkeypatch.setattr(cfg, "RESEARCH_DB_PATH", research_path)
    research.init_research_db()
    return str(state), str(research_path)


# ---------------------------------------------------------------------------
# Task 1: plumbing — fresh-DB safety, missing files, lock sentinel
# ---------------------------------------------------------------------------

def test_empty_dbs_report_runs(monkeypatch, tmp_path):
    """Fresh schemas, zero rows: the report renders instead of crashing."""
    from kalshi_mlb_mm import report
    state, research = _setup_dbs(monkeypatch, tmp_path)
    out = report.build_report(state, research, now=NOW)
    assert "State of the Maker" in out
    assert "no fills yet" in out


def test_missing_db_files_report_runs(tmp_path):
    """DB files that don't exist degrade to a note, never an exception."""
    from kalshi_mlb_mm import report
    out = report.build_report(
        str(tmp_path / "nope_state.duckdb"),
        str(tmp_path / "nope_research.duckdb"),
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
#
# The funnel is derived ONLY from quote_decisions: seen_rfqs is NOT a log of
# all flow — in-scope RFQs only get a row inside the live quote transaction
# (main.py:1124), so in-scope-but-skipped and dry-run RFQs never appear there.
# ---------------------------------------------------------------------------

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
    state, research = _setup_dbs(monkeypatch, tmp_path)
    recent = NOW - timedelta(hours=1)
    old = NOW - timedelta(days=8)  # outside every window

    _insert_decision(db, "r0", "skipped", "out_of_scope", old)
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
    state, research = _setup_dbs(monkeypatch, tmp_path)
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
    state, research = _setup_dbs(monkeypatch, tmp_path)
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
    state, research = _setup_dbs(monkeypatch, tmp_path)
    stats = report.on_demand_stats(research, NOW - timedelta(days=7))
    assert stats["flights_requested"] == 0
    assert stats["success_rate"] is None
    assert stats["per_book"] == []


# ---------------------------------------------------------------------------
# Task 4: quotable universe (#81: on_demand_result flights, research DB)
# ---------------------------------------------------------------------------

def _od_result(hash_, books):
    return {"leg_set_hash": hash_, "dropped_books": [],
            "books": {b: {"fair": 0.3, "route": "partition",
                          "latency_sec": 1.0} for b in books}}


def test_universe_stats(monkeypatch, tmp_path):
    """Per-day combos passing MIN_AGREEING_BOOKS + per-book participation,
    computed from the firehose's on_demand_result flights."""
    from kalshi_mlb_mm import report
    state, research = _setup_dbs(monkeypatch, tmp_path)
    day1 = NOW - timedelta(days=1)
    day2 = NOW

    # Combo h_sxt lands both days; h_mxt day1 only with a single book, so
    # it never passes consensus.
    _insert_event(research, "on_demand_result", day1,
                  _od_result("h_sxt", ["draftkings", "fanduel", "novig"]))
    _insert_event(research, "on_demand_result", day1 + timedelta(minutes=1),
                  _od_result("h_mxt", ["draftkings"]))
    _insert_event(research, "on_demand_result", day2,
                  _od_result("h_sxt", ["draftkings", "fanduel"]))

    stats = report.universe_stats(research, NOW - timedelta(days=7), NOW)
    per_day = {str(d): (passing, total)
               for d, passing, total in stats["per_day"]}
    assert per_day[str(day1.date())] == (1, 2)
    assert per_day[str(day2.date())] == (1, 1)
    per_book = {b: (days, combos) for b, days, combos, _age in
                stats["per_book"]}
    assert per_book["draftkings"] == (2, 2)
    assert per_book["fanduel"] == (2, 1)
    assert per_book["novig"] == (1, 1)


def test_universe_empty(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    state, research = _setup_dbs(monkeypatch, tmp_path)
    stats = report.universe_stats(research, NOW - timedelta(days=7), NOW)
    assert stats["per_day"] == []
    assert stats["per_book"] == []


# ---------------------------------------------------------------------------
# Task 5: staleness — book-data age at quote time / quote age at accept
# ---------------------------------------------------------------------------

def test_staleness_stats(monkeypatch, tmp_path):
    """#81: quote-time data age comes from quote_priced's live_games trace
    (max sub-fetch age); a quote without the trace counts as unknown."""
    from kalshi_mlb_mm import report
    state, research = _setup_dbs(monkeypatch, tmp_path)
    t0 = NOW - timedelta(hours=3)

    _insert_event(research, "quote_priced", t0 + timedelta(seconds=30),
                  {"blended_fair": 0.5,
                   "live_games": {"h1": {"age_sec": 2.0},
                                  "h2": {"age_sec": 30.0}}})
    _insert_event(research, "quote_priced", t0 + timedelta(seconds=150),
                  {"blended_fair": 0.6,
                   "live_games": {"h1": {"age_sec": 30.0}}})
    # engine hiccup — no live_games trace: unknown age, not a crash
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
    state, research = _setup_dbs(monkeypatch, tmp_path)
    stats = report.staleness_stats(research, NOW - timedelta(days=7))
    assert stats["quote_age"]["n"] == 0
    assert stats["accept_age"]["n"] == 0


# ---------------------------------------------------------------------------
# Task 6: demand curve — quotes vs accepts by margin x fair band
# ---------------------------------------------------------------------------

def test_demand_stats(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    import kalshi_mlb_mm.db as db
    state, research = _setup_dbs(monkeypatch, tmp_path)
    recent = NOW - timedelta(hours=1)

    # One quoted decision: yes side margin .60-.57=.03 at fair .60;
    # no side margin .40-.375=.025 at fair .40.
    _insert_decision(db, "r1", "quoted", None, recent,
                     blended_fair=0.60, yes_bid=0.57, no_bid=0.375)
    # One accept of the yes side at the quoted price.
    with db.connect() as con:
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, "
            "combo_market_ticker, game_id, side_held, contracts, price, fee, "
            "model_fair_at_quote, book_fair_at_quote, blended_fair_at_quote, "
            "fair_at_confirm, realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["f1", "q1", "r1", "T-r1", "g1", "yes", 10.0, 0.57, 0.01,
             None, 0.60, 0.60, 0.61, None, recent, True])

    stats = report.demand_stats(state, NOW - timedelta(days=7))
    assert stats["cells"][("3-4c", "[.50,.75)")] == {"quotes": 1, "accepts": 1}
    assert stats["cells"][("2-3c", "[.25,.50)")] == {"quotes": 1, "accepts": 0}
    assert stats["quote_sides"] == 2
    assert stats["accepts"] == 1


def test_demand_empty(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    state, research = _setup_dbs(monkeypatch, tmp_path)
    stats = report.demand_stats(state, NOW - timedelta(days=7))
    assert stats["cells"] == {}
    assert stats["quote_sides"] == 0


# ---------------------------------------------------------------------------
# Task 7: settlement P&L (markouts descoped with #13)
# ---------------------------------------------------------------------------

def _insert_fill(db, fill_id, contracts, price, fee, fair_at_quote,
                 realized_pnl):
    with db.connect() as con:
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, "
            "combo_market_ticker, game_id, side_held, contracts, price, fee, "
            "model_fair_at_quote, book_fair_at_quote, blended_fair_at_quote, "
            "fair_at_confirm, realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            [fill_id, f"q-{fill_id}", f"r-{fill_id}", f"T-{fill_id}", "g1",
             "yes", contracts, price, fee, None, fair_at_quote,
             fair_at_quote, fair_at_quote + 0.01, realized_pnl,
             NOW - timedelta(hours=2), True])


def test_pnl_stats_settled(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    import kalshi_mlb_mm.db as db
    state, research = _setup_dbs(monkeypatch, tmp_path)

    # Won yes-fill: 10 contracts at 55c, 1c/ct fee, fair at quote 58c.
    # realized = 10 * ((1 - 0.55) - 0.01) = 4.40 (settlement.py math)
    _insert_fill(db, "f1", 10.0, 0.55, 0.01, 0.58, 4.40)
    _insert_fill(db, "f2", 5.0, 0.30, 0.01, 0.33, None)  # not settled yet
    with db.connect() as con:
        con.execute(
            "INSERT INTO settlements (combo_market_ticker, result, "
            "settled_at, raw_payload, recorded_at) VALUES (?,?,?,?,?)",
            ["T-f1", "yes", NOW - timedelta(hours=1), "{}",
             NOW - timedelta(hours=1)])

    stats = report.pnl_stats(state)
    assert stats["fills"] == 2
    assert stats["settled_fills"] == 1
    assert stats["unsettled_fills"] == 1
    assert stats["realized_pnl_total"] == pytest.approx(4.40)
    assert stats["avg_quoted_margin_per_ct"] == pytest.approx(0.03, abs=1e-9)
    assert stats["realized_per_ct"] == pytest.approx(0.44)
    assert stats["by_result"] == [("yes", 1, pytest.approx(4.40))]


def test_pnl_stats_fills_but_no_settlements(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    import kalshi_mlb_mm.db as db
    state, research = _setup_dbs(monkeypatch, tmp_path)
    _insert_fill(db, "f1", 10.0, 0.55, 0.01, 0.58, None)
    stats = report.pnl_stats(state)
    assert stats["fills"] == 1
    assert stats["settled_fills"] == 0
    assert stats["realized_pnl_total"] is None
    out = report.build_report(state, research, now=NOW)
    assert "none settled yet" in out


def test_pnl_stats_empty(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    state, research = _setup_dbs(monkeypatch, tmp_path)
    stats = report.pnl_stats(state)
    assert stats["fills"] == 0


# ---------------------------------------------------------------------------
# Task 8: health — voids, halts, phantom fills, reconcile outcomes
# ---------------------------------------------------------------------------

def test_health_stats(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    import kalshi_mlb_mm.db as db
    state, research = _setup_dbs(monkeypatch, tmp_path)
    recent = NOW - timedelta(hours=1)

    _insert_decision(db, "r1", "confirmed", None, recent)
    _insert_decision(db, "r2", "confirmed", None, recent)
    _insert_decision(db, "r3", "voided_singles_moved", None, recent)
    _insert_decision(db, "r4", "sweep_cancel", "tipoff", recent)
    _insert_decision(db, "r5", "circuit_breaker",
                     "book_move_pulled_3_quotes", recent)
    _insert_decision(db, None, "halted_high_void_rate", None, recent)

    # Phantom fill: reconcile found no exchange record → contracts zeroed.
    with db.connect() as con:
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, "
            "combo_market_ticker, game_id, side_held, contracts, price, fee, "
            "model_fair_at_quote, book_fair_at_quote, blended_fair_at_quote, "
            "fair_at_confirm, realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["ph1", "q-ph1", "r-ph1", "T-ph1", "g1", "yes", 0.0, 0.5, 0.0,
             None, None, None, None, None, recent, True])

    _insert_event(research, "halt_event", recent,
                  {"halt_type": "void_rate", "transition": "fire"})
    _insert_event(research, "reconcile_done", recent, {"outcome": "phantom"})
    _insert_event(research, "reconcile_done", recent, {"outcome": "matched"})

    stats = report.health_stats(state, research, NOW - timedelta(days=7))
    assert stats["confirmed"] == 2
    assert stats["voids_by_decision"] == {"voided_singles_moved": 1}
    assert stats["void_rate"] == pytest.approx(1 / 3)
    assert stats["sweep_cancels_by_reason"] == {"tipoff": 1}
    assert stats["circuit_breakers"] == 1
    assert stats["void_rate_halts"] == 1
    assert stats["halt_transitions"] == {"fire": 1}
    assert stats["phantom_fills"] == 1
    assert stats["reconcile_outcomes"] == {"matched": 1, "phantom": 1}


def test_health_empty(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    state, research = _setup_dbs(monkeypatch, tmp_path)
    stats = report.health_stats(state, research, NOW - timedelta(days=7))
    assert stats["void_rate"] is None
    assert stats["phantom_fills"] == 0
