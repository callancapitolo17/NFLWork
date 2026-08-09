"""Issue #81 — remove the maker's sweep entirely.

What each block pins (written FIRST, failing on pre-#81 code):

1. Zombie-gate: the sweep's symbols are DELETED from the source, not
   flagged off — `_SGP_ODDS`, `_refresh_sgp`, `_book_fairs`, `scrape_done`,
   `sgp_cycle` gone from main.py; `SGP_REFRESH_SEC` and
   `SGP_SCRAPER_TIMEOUT_SEC` gone from the maker's config.py (the taker
   keeps its own copies). Mirrors #56/#57's source-scan style.
2. Target-line refresh is its own loop arm: `_target_line_tick` runs the
   extracted `target_line_cycle` (zero book requests) and re-homes the O-1
   game-metadata cache invalidation — including on enumeration failure.
3. Coverage heir: `_coverage_summary_tick` drains the service's
   unconditional on-demand coverage counter into an `on_demand_coverage`
   research event; an idle window (empty snapshot) emits nothing.
4. Config: new keyword-checked knobs `TARGET_LINE_REFRESH_SEC` and
   `COVERAGE_SUMMARY_SEC`; warming stays strictly faster than the
   target-line cadence.
5. report.py's "quotable universe" reads the research firehose's
   `on_demand_result` events — the live-flight record — instead of the
   now-unwritten `mlb_sgp_odds`.

Clock control is injected, never slept (repo correctness bar).
"""
import inspect
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest


def _src(module) -> str:
    return Path(inspect.getsourcefile(module)).read_text()


# --------------------------------------------------------------------------- #
# 1. Zombie-gate: sweep symbols are deleted, not flagged                       #
# --------------------------------------------------------------------------- #

def test_sweep_symbols_deleted_from_main():
    from kalshi_mlb_mm import main
    src = _src(main)
    for symbol in ("_SGP_ODDS", "_refresh_sgp", "_book_fairs",
                   "scrape_done", "sgp_cycle"):
        assert symbol not in src, (
            f"{symbol} must be deleted from main.py, not left behind "
            "(repo rule: delete dead code, git is the archive)")


def test_sweep_knobs_deleted_from_config():
    import kalshi_mlb_mm.config as cfg
    src = _src(cfg)
    for knob in ("SGP_REFRESH_SEC", "SGP_SCRAPER_TIMEOUT_SEC"):
        assert knob not in src, (
            f"{knob} must be deleted from the maker's config.py "
            "(the taker keeps its own copy)")


def test_target_line_arm_wired_in_main():
    from kalshi_mlb_mm import main
    src = _src(main)
    assert "target_line_cycle" in src
    assert "TARGET_LINE_REFRESH_SEC" in src
    assert "COVERAGE_SUMMARY_SEC" in src


# --------------------------------------------------------------------------- #
# 2. Target-line tick: refresh + O-1 cache-invalidation heir                   #
# --------------------------------------------------------------------------- #

def test_target_line_tick_refreshes_and_invalidates_caches(monkeypatch):
    from kalshi_common import sgp_runner
    from kalshi_mlb_mm import main

    calls = []
    monkeypatch.setattr(sgp_runner, "target_line_cycle",
                        lambda **kw: calls.append(kw) or [])
    monkeypatch.setitem(main._COMMENCE_CACHE, "g-old", None)
    monkeypatch.setitem(main._RESOLVE_CACHE, "h-old", "g-old")

    main._target_line_tick()

    assert len(calls) == 1
    assert calls[0]["both_teams"] is True     # away-margin tickers need both
    assert main._COMMENCE_CACHE == {}
    assert main._RESOLVE_CACHE == {}


def test_target_line_tick_invalidates_even_on_enumeration_failure(monkeypatch):
    """Kalshi API down mid-tick: the stale caches must still drop (the old
    `_refresh_sgp` invalidated before its early returns for the same
    reason), and the failure propagates for the loop arm to log."""
    from kalshi_common import sgp_runner
    from kalshi_mlb_mm import main

    def boom(**kw):
        raise RuntimeError("kalshi api down")

    monkeypatch.setattr(sgp_runner, "target_line_cycle", boom)
    monkeypatch.setitem(main._COMMENCE_CACHE, "g-old", None)

    with pytest.raises(RuntimeError):
        main._target_line_tick()
    assert main._COMMENCE_CACHE == {}


# --------------------------------------------------------------------------- #
# 3. Coverage summary: scrape_done's observability heir                        #
# --------------------------------------------------------------------------- #

class _FakeCoverage:
    def __init__(self, snap):
        self.snap = snap
        self.calls = 0

    def snapshot_and_reset(self):
        self.calls += 1
        out, self.snap = self.snap, {}
        return out


class _FakeSvc:
    def __init__(self, snap):
        self.coverage = _FakeCoverage(snap)


def test_coverage_summary_tick_emits_and_drains(monkeypatch):
    from kalshi_mlb_mm import main

    emitted = []
    monkeypatch.setattr(main.research, "emit",
                        lambda ev, **kw: emitted.append((ev, kw)))
    svc = _FakeSvc({"draftkings": {"ok": 12},
                    "novig": {"transport_error": 3}})

    main._coverage_summary_tick(sgp_service=svc)

    assert svc.coverage.calls == 1
    assert len(emitted) == 1
    event, kw = emitted[0]
    assert event == "on_demand_coverage"
    assert kw["payload"]["books"] == {"draftkings": {"ok": 12},
                                      "novig": {"transport_error": 3}}


def test_coverage_summary_tick_skips_idle_window(monkeypatch):
    from kalshi_mlb_mm import main

    emitted = []
    monkeypatch.setattr(main.research, "emit",
                        lambda ev, **kw: emitted.append((ev, kw)))
    main._coverage_summary_tick(sgp_service=_FakeSvc({}))
    assert emitted == []


def test_coverage_summary_tick_param_keyword_only():
    from kalshi_mlb_mm import main
    sig = inspect.signature(main._coverage_summary_tick)
    assert sig.parameters["sgp_service"].kind is inspect.Parameter.KEYWORD_ONLY


# --------------------------------------------------------------------------- #
# 4. Config knobs                                                              #
# --------------------------------------------------------------------------- #

def test_new_cadence_knobs_defaults():
    import kalshi_mlb_mm.config as cfg
    assert cfg.TARGET_LINE_REFRESH_SEC == 300
    assert cfg.COVERAGE_SUMMARY_SEC == 300
    # Warming must outrun target-line refresh: new games get structures
    # warmed within a warming period of appearing in mlb_target_lines.
    assert cfg.STRUCTURE_WARM_SEC < cfg.TARGET_LINE_REFRESH_SEC


# --------------------------------------------------------------------------- #
# 5. report.py universe reads the firehose flight record                       #
# --------------------------------------------------------------------------- #

NOW = datetime(2026, 8, 9, 18, 0, 0, tzinfo=timezone.utc)


def _research_db(monkeypatch, tmp_path):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.research as research
    path = tmp_path / "research.duckdb"
    monkeypatch.setattr(cfg, "RESEARCH_DB_PATH", path)
    research.init_research_db()
    return str(path)


def _insert_result(db_path, event_id, ts, hash_, books):
    import duckdb
    payload = json.dumps({"leg_set_hash": hash_, "dropped_books": [],
                          "books": {b: {"fair": 0.3, "route": "partition",
                                        "n_cells": 4, "latency_sec": 1.0}
                                    for b in books}})
    con = duckdb.connect(db_path)
    try:
        con.execute(
            "INSERT INTO events (event_id, session_id, event_type, ts, "
            "ticker, rfq_id, quote_id, payload) VALUES (?,?,?,?,?,?,?,?)",
            [event_id, "s1", "on_demand_result", ts, "T", "r1", None, payload])
    finally:
        con.close()


def test_universe_stats_from_on_demand_results(monkeypatch, tmp_path):
    import kalshi_mlb_mm.config as cfg
    from kalshi_mlb_mm import report

    monkeypatch.setattr(cfg, "MIN_AGREEING_BOOKS", 2)
    db = _research_db(monkeypatch, tmp_path)

    day1 = NOW - timedelta(days=2)
    # day1: h1 lands with 2 books (passing), h2 with 1 (not passing)
    _insert_result(db, "e1", day1, "h1", ["draftkings", "fanduel"])
    _insert_result(db, "e2", day1 + timedelta(minutes=5), "h2", ["draftkings"])
    # day2 (= NOW's day): h1 lands thin first, then with 3 books — a hash
    # passes if ANY landing that day reached the gate.
    _insert_result(db, "e3", NOW - timedelta(hours=2), "h1", ["draftkings"])
    _insert_result(db, "e4", NOW - timedelta(hours=1), "h1",
                   ["draftkings", "fanduel", "novig"])

    stats = report.universe_stats(db, NOW - timedelta(days=7), NOW)

    assert stats["unavailable"] is None
    assert stats["min_books"] == 2
    per_day = {str(day): (passing, total)
               for day, passing, total in stats["per_day"]}
    assert per_day[str(day1.date())] == (1, 2)
    assert per_day[str(NOW.date())] == (1, 1)

    per_book = {row[0]: row[1:] for row in stats["per_book"]}
    assert per_book["draftkings"][:2] == (2, 2)   # 2 days, 2 distinct hashes
    assert per_book["fanduel"][:2] == (2, 1)
    assert per_book["novig"][:2] == (1, 1)
    assert per_book["novig"][2] == 60             # minutes since last landing


def test_universe_stats_empty_research_db(monkeypatch, tmp_path):
    from kalshi_mlb_mm import report
    db = _research_db(monkeypatch, tmp_path)
    stats = report.universe_stats(db, NOW - timedelta(days=7), NOW)
    assert stats["per_day"] == []
    assert stats["unavailable"] is None
