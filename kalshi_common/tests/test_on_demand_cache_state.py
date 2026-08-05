"""Issue #50 item 1: explicit cold-vs-warm tagging of on-demand health rows.

"Cold" means the call paid a wire structure fetch (events listing and/or the
per-event structure); "warm" means everything structural came from the TTL
caches and only price calls hit the wire. The tag is a dedicated
``cache_state`` column on ``sgp_fetch_health`` so cold/warm percentiles are
one GROUP BY away — never inferred from durations.

Detection plumbing under test:
  - ``TTLCache.get_or_fetch`` grows a keyword-only ``miss_cb`` fired only
    when the fetch actually runs;
  - the on-demand hooks wire ``miss_cb`` to a new ``structure_fetches``
    counter;
  - ``price_on_demand`` derives cache_state from that counter and passes it
    to the health recorder.

The real-hooks test uses BetMGM with a fake client — the same monkeypatch
pattern as test_sgp_service_counters.py::test_real_hooks_thread_counters...
"""
import inspect
from types import SimpleNamespace

import duckdb
import pytest

from kalshi_common.legset import CanonicalLeg
from kalshi_common.sgp_health import FetchHealthRecorder
from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import GameRef, ResolvedLeg, TTLCache

EVT = "KXMLBGAME-26AUG05NYYBOS"
GAME = GameRef(game_id="g1", home_team="Boston Red Sox",
               away_team="New York Yankees", commence_time=None)

# Four joint fairs summing to 1.0 with a uniform 1.2 overround (2-leg Route A).
CELL_FAIRS = [0.40, 0.20, 0.25, 0.15]
CELL_DECS = [1.0 / (1.2 * p) for p in CELL_FAIRS]


def _legs(n=2):
    return [CanonicalLeg(EVT, "spread", -1.5, "home"),
            CanonicalLeg(EVT, "total", 8.5, "over")][:n]


def _resolved(n=2):
    return [ResolvedLeg(ref=f"r{j}", opposite_ref=f"o{j}",
                        single_decimal=1.9, opposite_decimal=1.9)
            for j in range(n)]


# --------------------------------------------------------------------- #
# TTLCache.miss_cb                                                       #
# --------------------------------------------------------------------- #

def test_ttlcache_miss_cb_fires_on_miss_not_on_hit():
    cache = TTLCache(ttl_sec=300)
    misses = []
    assert cache.get_or_fetch("k", lambda: "v1",
                              miss_cb=lambda: misses.append(1)) == "v1"
    assert misses == [1]                      # cold: fetch ran
    assert cache.get_or_fetch("k", lambda: "v2",
                              miss_cb=lambda: misses.append(1)) == "v1"
    assert misses == [1]                      # warm: served from cache


def test_ttlcache_miss_cb_fires_again_after_expiry():
    clock = SimpleNamespace(t=1000.0)
    cache = TTLCache(ttl_sec=10, now_fn=lambda: clock.t)
    misses = []
    cache.get_or_fetch("k", lambda: "v", miss_cb=lambda: misses.append(1))
    clock.t += 11
    cache.get_or_fetch("k", lambda: "v", miss_cb=lambda: misses.append(1))
    assert misses == [1, 1]


def test_ttlcache_miss_cb_is_keyword_only():
    sig = inspect.signature(TTLCache.get_or_fetch)
    assert sig.parameters["miss_cb"].kind is inspect.Parameter.KEYWORD_ONLY
    # bind() proof: a third positional argument must not reach miss_cb
    with pytest.raises(TypeError):
        sig.bind(object(), "key", lambda: None, lambda: None)


# --------------------------------------------------------------------- #
# End to end: real hooks, fake client -> cold row then warm row          #
# --------------------------------------------------------------------- #

class FakeMGMClient:
    """Counts wire-level structure calls so the test can prove which
    price_on_demand call paid them."""

    def __init__(self):
        self.structure_calls = 0

    def list_events(self, profile=None):
        self.structure_calls += 1
        return ["raw-event"]

    def fetch_markets(self, event_id, profile=None):
        self.structure_calls += 1
        return [{"raw": "market"}]


def _patched_mgm_service(monkeypatch, db_path):
    from mlb_sgp import betmgm as mod
    svc = SGPService(books=("betmgm",), health_db_path=db_path)
    st = svc._state["betmgm"]
    st.client = FakeMGMClient()
    monkeypatch.setattr(svc, "_ensure_client", lambda b: st)
    ev = SimpleNamespace(event_id="F1", home_team="Boston Red Sox",
                         away_team="New York Yankees")
    monkeypatch.setattr(mod, "_match_events",
                        lambda events, targets: {GAME.game_id: ev})
    monkeypatch.setattr(mod, "parse_markets",
                        lambda markets, home, away: {"FG": {"fake": "s"}})
    monkeypatch.setattr(
        mod, "resolve_legs",
        lambda structure, legs, home, away, counters=None: _resolved(2))
    calls = {"n": 0}

    def fake_price(client, refs, fixture_id, counters=None):
        dec = CELL_DECS[calls["n"] % 4]
        calls["n"] += 1
        return dec
    monkeypatch.setattr(mod, "price_selection_set", fake_price)
    return svc, st


def test_first_fetch_records_cold_then_cached_fetch_records_warm(
        tmp_path, monkeypatch):
    db = str(tmp_path / "health.duckdb")
    svc, st = _patched_mgm_service(monkeypatch, db)

    res1 = svc.price_on_demand("betmgm", GAME, _legs())
    assert res1 is not None
    paid_cold = st.client.structure_calls
    assert paid_cold > 0                      # first call hit the wire

    res2 = svc.price_on_demand("betmgm", GAME, _legs())
    assert res2 is not None
    assert st.client.structure_calls == paid_cold   # served from cache

    assert svc.health.flush() == 2
    con = duckdb.connect(db, read_only=True)
    try:
        rows = con.execute(
            "SELECT cache_state FROM sgp_fetch_health "
            "WHERE path = 'on_demand' ORDER BY fetched_at").fetchall()
    finally:
        con.close()
    assert [r[0] for r in rows] == ["cold", "warm"]


def test_cold_warm_counter_rides_on_the_result_counters(
        tmp_path, monkeypatch):
    """structure_fetches is a first-class #35 counter: >0 on the cold call,
    0 on the warm repeat — the health column is derived from it."""
    svc, _ = _patched_mgm_service(monkeypatch, str(tmp_path / "h.duckdb"))
    res1 = svc.price_on_demand("betmgm", GAME, _legs())
    res2 = svc.price_on_demand("betmgm", GAME, _legs())
    assert res1.counters.structure_fetches > 0
    assert res2.counters.structure_fetches == 0


# --------------------------------------------------------------------- #
# Schema: migration + sweep rows                                         #
# --------------------------------------------------------------------- #

OLD_SCHEMA_SQL = """
CREATE TABLE sgp_fetch_health (
    fetched_at     TIMESTAMPTZ,
    book           VARCHAR,
    path           VARCHAR,
    outcome        VARCHAR,
    rows_or_prices INTEGER,
    duration_sec   DOUBLE,
    error_class    VARCHAR,
    counters_json  VARCHAR
);
"""


def test_flush_migrates_a_pre_50_table_and_keeps_old_rows(tmp_path):
    db = str(tmp_path / "old.duckdb")
    con = duckdb.connect(db)
    con.execute(OLD_SCHEMA_SQL)
    con.execute(
        "INSERT INTO sgp_fetch_health VALUES "
        "(now(), 'fanduel', 'on_demand', 'ok', 4, 0.8, NULL, NULL)")
    con.close()

    rec = FetchHealthRecorder(db)
    rec.record(book="betmgm", path="on_demand", outcome="ok",
               rows_or_prices=4, duration_sec=1.2, cache_state="warm")
    assert rec.flush() == 1

    con = duckdb.connect(db, read_only=True)
    try:
        rows = con.execute(
            "SELECT book, cache_state FROM sgp_fetch_health "
            "ORDER BY book").fetchall()
    finally:
        con.close()
    # legacy row survives with NULL cache_state; new row carries the tag
    assert rows == [("betmgm", "warm"), ("fanduel", None)]


def test_sweep_rows_carry_null_cache_state(tmp_path):
    db = str(tmp_path / "sweep.duckdb")
    svc = SGPService(books=("fanduel",),
                     runners={"fanduel": lambda targets: []},
                     health_db_path=db)
    svc.refresh([])                            # writes + flushes one sweep row
    con = duckdb.connect(db, read_only=True)
    try:
        rows = con.execute(
            "SELECT path, cache_state FROM sgp_fetch_health").fetchall()
    finally:
        con.close()
    assert rows == [("sweep", None)]


def test_recorder_warns_but_never_raises_on_unknown_cache_state(
        tmp_path, caplog):
    """Same contract as outcome/path: a vocabulary typo is a WARNING that
    poisons queries, never an exception that breaks a fetch."""
    import logging
    rec = FetchHealthRecorder(str(tmp_path / "h.duckdb"))
    with caplog.at_level(logging.WARNING, logger="kalshi_common.sgp_health"):
        rec.record(book="fanduel", path="on_demand", outcome="ok",
                   cache_state="lukewarm")
    assert rec.pending() == 1                  # row still buffered
    assert any("lukewarm" in r.message for r in caplog.records)
