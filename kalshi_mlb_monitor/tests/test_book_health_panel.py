"""Monitor's per-book feed-health readout (issue #37, fix-spec item 5).

The panel exists so "DK has been dead for a week" is visible at a glance.
The ways it can lie:

  * painting a book red because it answered 'empty' (off-day / thin slate),
  * hiding a book that died days ago because the 24h window has no rows,
  * reporting LOCKED when the bot simply never wrote the table,
  * disagreeing with the threshold the bot actually alerts on.

Note this DB is a real DuckDB file in tmp_path — never a symlink into a
worktree, because DuckDB writes its WAL next to the PATH.
"""
from __future__ import annotations

from datetime import datetime, timedelta, timezone

import duckdb
import pytest

from kalshi_common import health_queries as HQ
from kalshi_common.book_health import DEFAULT_STREAK_THRESHOLD
from kalshi_common.sgp_health import SCHEMA_SQL
from kalshi_mlb_monitor import queries as Q

NOW = datetime.now(timezone.utc)


@pytest.fixture
def market_db(tmp_path, monkeypatch):
    path = tmp_path / "market.duckdb"
    con = duckdb.connect(str(path))
    con.execute(SCHEMA_SQL)
    con.close()
    monkeypatch.setitem(Q.BOTS["maker"], "market_db", str(path))
    return path


def _insert(path, rows):
    """rows: (book, path, minutes_ago, outcome, error_class)"""
    con = duckdb.connect(str(path))
    try:
        con.executemany(
            "INSERT INTO sgp_fetch_health (fetched_at, book, path, outcome, "
            "rows_or_prices, duration_sec, error_class, counters_json) "
            "VALUES (?, ?, ?, ?, 1, 1.0, ?, NULL)",
            [(NOW - timedelta(minutes=m), b, p, o, e)
             for (b, p, m, o, e) in rows])
    finally:
        con.close()


# ------------------------------------------------------------------ #
# The shared query loader                                             #
# ------------------------------------------------------------------ #

def test_loader_returns_the_shipped_streak_query():
    sql = HQ.named_query("current_failure_streak")
    assert "outcome IN ('ok', 'empty')" in sql
    assert "failing_streak" in sql


def test_unknown_query_name_raises():
    """An empty string would read downstream as 'no unhealthy books' — the
    silent-failure shape this epic removes."""
    with pytest.raises(KeyError):
        HQ.named_query("no_such_query")


# ------------------------------------------------------------------ #
# book_health()                                                       #
# ------------------------------------------------------------------ #

def test_missing_table_is_none_not_locked(tmp_path, monkeypatch):
    """A bot that hasn't run since #38 shipped has no table. That is 'no
    history', not 'the live data is locked'."""
    path = tmp_path / "empty.duckdb"
    duckdb.connect(str(path)).close()
    monkeypatch.setitem(Q.BOTS["maker"], "market_db", str(path))
    assert Q.book_health("maker") is None


def test_healthy_book_is_green_and_empty_counts_as_healthy(market_db):
    _insert(market_db, [("fanduel", "sweep", m, "empty", None)
                        for m in (30, 20, 10, 5, 1)])
    df = Q.book_health("maker")
    row = df[df.book == "fanduel"].iloc[0]
    assert row["failing_streak"] == 0
    assert row["status"] == "ok"
    # 100% ANSWERED even though 0% were 'ok' — the panel must not show a
    # green dot next to "0%", which would teach the reader to distrust it.
    assert row["pct_answered"] == 100.0


def test_dead_book_is_red_with_its_error_class(market_db):
    _insert(market_db, [("draftkings", "sweep", m, "transport_error",
                         "BookTransportError:events:403")
                        for m in (30, 20, 10)])
    row = Q.book_health("maker").set_index("book").loc["draftkings"]
    assert row["failing_streak"] == 3
    assert row["status"] == "dead"
    assert row["newest_error_class"] == "BookTransportError:events:403"


def test_partial_streak_is_amber(market_db):
    _insert(market_db, [("novig", "sweep", 30, "ok", None),
                        ("novig", "sweep", 10, "timeout", "FutureTimeout")])
    row = Q.book_health("maker").set_index("book").loc["novig"]
    assert row["failing_streak"] == 1
    assert row["status"] == "degraded"


def test_book_that_died_days_ago_still_appears(market_db):
    """It has no rows in the 24h mix. Dropping it would make the longest
    outage the most invisible one."""
    old = 60 * 24 * 3        # 3 days ago, inside the streak query's 7d window
    _insert(market_db, [("caesars", "sweep", old + i, "transport_error",
                         "BookTransportError:auth:403") for i in range(4)])
    df = Q.book_health("maker")
    assert "caesars" in set(df["book"])
    row = df.set_index("book").loc["caesars"]
    assert row["status"] == "dead"
    assert row["n_fetches"] == 0          # nothing in the 24h window


def test_paths_are_reported_separately(market_db):
    _insert(market_db,
            [("betmgm", "sweep", 5, "ok", None),
             ("betmgm", "on_demand", 5, "transport_error", "X"),
             ("betmgm", "on_demand", 4, "transport_error", "X"),
             ("betmgm", "on_demand", 3, "transport_error", "X")])
    df = Q.book_health("maker").set_index(["book", "path"])
    assert df.loc[("betmgm", "sweep"), "status"] == "ok"
    assert df.loc[("betmgm", "on_demand"), "status"] == "dead"


def test_panel_threshold_matches_the_bot_alert_threshold():
    assert Q.BOOK_ALERT_STREAK == DEFAULT_STREAK_THRESHOLD
    assert Q._health_status(DEFAULT_STREAK_THRESHOLD) == "dead"
    assert Q._health_status(DEFAULT_STREAK_THRESHOLD - 1) == "degraded"
    assert Q._health_status(0) == "ok"


def test_both_bots_have_a_market_db_path():
    for key in ("maker", "taker"):
        assert Q.BOTS[key]["market_db"].endswith("_market.duckdb"), key
