"""DuckDB mirror of the leg surface (issue #96).

Side effects: creates and writes ``kalshi_mlb_mm/kalshi_mlb_mm_surface.duckdb``
— its OWN file and write lock. Deliberately NOT the market DB
(``kalshi_mlb_mm_market.duckdb``), which the maker's pricing path reads: a
writer on a 20s cadence has no business contending with it.

This mirror is for research, the monitor, and #96's acceptance queries. The
quote path reads ``LegSurface`` in memory and never opens this file.

Tables:
  mlb_leg_surface     current state, one row per (book, leg). Replaced per
                      book, so duplicate keys are structurally impossible.
  surface_refresh_log one row per (book, pass) — achieved cadence and the
                      countable exclusion reasons.
  mlb_surface_slate   current slate + leg counts, for debugging a book that
                      is priced nothing.
"""
from __future__ import annotations

import logging
import time
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
from random import random

import duckdb

from kalshi_mlb_mm import config

log = logging.getLogger(__name__)

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS mlb_leg_surface (
    book             VARCHAR     NOT NULL,
    -- Kalshi event-ticker suffix (26AUG252138CLELAA): unique per
    -- doubleheader game, and exactly what CanonicalLeg.game_id carries.
    game_id          VARCHAR     NOT NULL,
    game_start_time  TIMESTAMPTZ NOT NULL,
    period           VARCHAR     NOT NULL,   -- 'FG' | 'F5' | 'I1'
    market_type      VARCHAR     NOT NULL,   -- 'ml' | 'spread' | 'total'
    line             DOUBLE,                 -- signed home-perspective; NULL for ml
    side             VARCHAR     NOT NULL,
    fair_prob        DOUBLE      NOT NULL,
    -- raw prices retained for debugging
    raw_decimal      DOUBLE      NOT NULL,
    raw_decimal_opp  DOUBLE      NOT NULL,
    raw_overround    DOUBLE      NOT NULL,   -- pre-devig implied sum
    route            VARCHAR     NOT NULL,   -- 'structure' | 'singles'
    -- when the BOOK DATA was fetched, not when the row was devigged or
    -- written: one structure fetch stamps every leg it served, so #99's age
    -- gate reads the price's true age.
    built_at         TIMESTAMPTZ NOT NULL
);
-- No PRIMARY KEY: DuckDB PKs reject NULL and `line` is legitimately NULL for
-- moneyline (the established ml x total convention — NULL, not a sentinel).
-- Uniqueness is structural instead: each flush DELETEs one book's rows and
-- re-inserts that book's in-memory slice, which is a dict on the key.
CREATE TABLE IF NOT EXISTS surface_refresh_log (
    book              VARCHAR     NOT NULL,
    route             VARCHAR     NOT NULL,
    started_at        TIMESTAMPTZ NOT NULL,
    duration_sec      DOUBLE      NOT NULL,
    games_attempted   INTEGER     NOT NULL,
    games_priced      INTEGER     NOT NULL,
    rungs_priced      INTEGER     NOT NULL,
    legs_written      INTEGER     NOT NULL,
    -- Fixed, small exclusion vocabulary: a reason is a column, so counting
    -- one is plain SQL. Per-rung exclusion ROWS were rejected — FanDuel
    -- alone one-sides ~51 deep-alt rungs per slate, which at this cadence is
    -- thousands of noise rows an hour.
    n_crossed         INTEGER     NOT NULL,  -- implied sum < 1.0
    n_overround       INTEGER     NOT NULL,  -- outside the vig envelope
    n_one_sided       INTEGER     NOT NULL,  -- book posts only one side
    n_unresolved      INTEGER     NOT NULL,  -- book does not offer the rung
    n_game_unmatched  INTEGER     NOT NULL,  -- book does not list the game
    n_game_ambiguous  INTEGER     NOT NULL,  -- 2+ book games matched; fail closed
    error_class       VARCHAR                -- transport failure, else NULL
);
CREATE TABLE IF NOT EXISTS mlb_surface_slate (
    game_id          VARCHAR     NOT NULL,
    home_team        VARCHAR     NOT NULL,
    away_team        VARCHAR     NOT NULL,
    game_start_time  TIMESTAMPTZ NOT NULL,
    n_legs           INTEGER     NOT NULL,
    discovered_at    TIMESTAMPTZ NOT NULL
);
"""

_SURFACE_COLS = ("book, game_id, game_start_time, period, market_type, line, "
                 "side, fair_prob, raw_decimal, raw_decimal_opp, "
                 "raw_overround, route, built_at")


@contextmanager
def connect(read_only: bool = False, retries: int = 10):
    """Retrying connect (the maker's db.connect pattern). A lock contention
    is transient; a raise here would kill an ingest worker."""
    last_err = None
    for attempt in range(retries):
        try:
            con = duckdb.connect(str(config.SURFACE_DB), read_only=read_only)
            try:
                yield con
            finally:
                con.close()
            return
        except duckdb.IOException as e:
            last_err = e
            time.sleep(0.05 * (2 ** attempt) + random() * 0.05)
    raise last_err


def init_database() -> None:
    with connect() as con:
        con.execute(SCHEMA_SQL)


def _aware(dt: datetime | None) -> datetime | None:
    """Naive values in this codebase are wall-clock UTC (the mlb_target_lines
    convention); TIMESTAMPTZ needs the offset made explicit or DuckDB reads
    them as local time."""
    if dt is None:
        return None
    return dt.replace(tzinfo=timezone.utc) if dt.tzinfo is None else dt


def replace_slice_rows(book: str, route: str, rows: list) -> int:
    """Atomically swap one (book, route) slice of mlb_leg_surface.

    Keyed by route as well as book because FanDuel publishes through both,
    each owning a disjoint set of (market_type, period) pairs.

    One connection for the whole slice (~1,000 rows), never one per row — a
    connect costs ~17ms.
    """
    values = [[r.book, r.game_id, _aware(r.game_start_time), r.period,
               r.market_type, r.line, r.side, r.fair_prob, r.raw_decimal,
               r.raw_decimal_opp, r.raw_overround, r.route, _aware(r.built_at)]
              for r in rows]
    with connect() as con:
        con.execute("BEGIN TRANSACTION")
        try:
            con.execute("DELETE FROM mlb_leg_surface "
                        "WHERE book = ? AND route = ?", [book, route])
            if values:
                con.executemany(
                    f"INSERT INTO mlb_leg_surface ({_SURFACE_COLS}) VALUES "
                    "(?,?,?,?,?,?,?,?,?,?,?,?,?)", values)
            con.execute("COMMIT")
        except Exception:
            con.execute("ROLLBACK")
            raise
    return len(values)


def record_pass(pass_result) -> None:
    """One surface_refresh_log row. ``pass_result`` is a runner.PassResult."""
    p = pass_result
    with connect() as con:
        con.execute(
            "INSERT INTO surface_refresh_log (book, route, started_at, "
            "duration_sec, games_attempted, games_priced, rungs_priced, "
            "legs_written, n_crossed, n_overround, n_one_sided, n_unresolved, "
            "n_game_unmatched, n_game_ambiguous, error_class) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            [p.book, p.route, _aware(p.started_at), p.duration_sec,
             p.games_attempted, p.games_priced, p.rungs_priced,
             p.legs_written, p.counts.crossed, p.counts.overround,
             p.counts.one_sided, p.counts.unresolved,
             p.counts.game_unmatched, p.counts.game_ambiguous,
             p.error_class])


def record_slate(games: list) -> None:
    """Replace mlb_surface_slate with the current slate."""
    now = datetime.now(timezone.utc)
    values = [[g.game_id, g.home_team, g.away_team, _aware(g.start_utc),
               len(g.legs), now] for g in games]
    with connect() as con:
        con.execute("BEGIN TRANSACTION")
        try:
            con.execute("DELETE FROM mlb_surface_slate")
            if values:
                con.executemany(
                    "INSERT INTO mlb_surface_slate (game_id, home_team, "
                    "away_team, game_start_time, n_legs, discovered_at) "
                    "VALUES (?,?,?,?,?,?)", values)
            con.execute("COMMIT")
        except Exception:
            con.execute("ROLLBACK")
            raise


def prune_refresh_log(retention_hours: float) -> int:
    """Drop refresh-log rows older than the retention window. The log is the
    only unbounded table here (one row per book per pass, ~15k/day at a 20s
    cadence)."""
    cutoff = datetime.now(timezone.utc) - timedelta(hours=retention_hours)
    with connect() as con:
        # DuckDB returns the deleted row count from DELETE, so this stays one
        # statement — a pair of COUNT(*)s to report a delta would scan the
        # whole table twice on every housekeeping tick.
        row = con.execute(
            "DELETE FROM surface_refresh_log WHERE started_at < ?",
            [cutoff]).fetchone()
    return int(row[0]) if row else 0
