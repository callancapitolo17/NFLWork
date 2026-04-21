"""Dashboard query functions — separated from Dash callbacks for testability.

These are pure(ish) functions that read from DuckDB and return plain Python
dicts/lists. No Dash imports here, so they can be unit-tested with a fresh
temp DuckDB file in each test.

Lock-contention safety
----------------------
The cron writer (``run.py --mode scrape``) can briefly hold the DuckDB file
lock while it ingests a batch of odds. If a Dash callback tries to read at
the exact same moment, DuckDB raises an ``IOException`` like::

    IO Error: Could not set lock on file ... Conflicting lock is held

That's a transient condition — the next tick of the Dash interval (15 sec on
the pre-draft cadence, 3 sec on draft day) will succeed. Rather than letting
the exception bubble up and 500 the whole callback, the query functions here
catch lock-family errors and return a ``QueryLocked`` sentinel. Callbacks
interpret the sentinel as ``PreventUpdate`` so the UI simply keeps showing
the last successful render.

Non-lock exceptions (schema errors, typos, etc.) still propagate — we only
short-circuit on lock/IO/conflict errors, which have a predictable substring
signature.
"""

import statistics
from typing import List, Dict, Any, Optional
from datetime import datetime

import duckdb

from nfl_draft.lib.db import read_connection


# Any row older than MAX_AGE_HOURS is excluded from Cross-Book Grid /
# ev_candidates. This prevents venues that have silently stopped scraping
# (e.g. FD regression 2026-04-19) from polluting the grid with 24h-old
# prices that get flagged as "edges" when in reality the feed just died.
# Tune via this constant — 2h is a conservative default for pre-draft
# cadence (scrapers run every few minutes).
MAX_AGE_HOURS = 2


class QueryLocked:
    """Sentinel returned by query functions when DuckDB lock contention
    prevents a read. Callers should treat this as 'no new data available —
    don't re-render'. A singleton instance (``_LOCKED``) is used; callers
    should check with ``isinstance(result, QueryLocked)``.
    """


_LOCKED = QueryLocked()


def _is_lock_error(err: BaseException) -> bool:
    """True if this exception matches the DuckDB lock-contention signature.

    DuckDB surfaces lock errors as IOException with messages like
    'Could not set lock on file' / 'Conflicting lock is held'. We also
    include 'attached' / 'file handle conflict' for the in-process collision
    mode that happens when two readers race to open the same file.
    """
    msg = str(err).lower()
    return any(marker in msg for marker in (
        "lock on file",
        "conflicting lock",
        "could not set lock",
        "i/o error",
        "io error",
        "file handle conflict",
    ))


def _safe_read(fn):
    """Invoke ``fn``; return its result, or ``_LOCKED`` on a lock error.

    Any other exception propagates unchanged so genuine bugs (missing table,
    bad SQL, etc.) still surface loudly.
    """
    try:
        return fn()
    except (duckdb.IOException, duckdb.Error) as e:
        if _is_lock_error(e):
            return _LOCKED
        raise


def cross_book_grid(threshold_pp: float = 10.0):
    """Return one row per market with all-venue prices, median, and outlier flags.

    For each market, take the most-recent devig_prob per book, compute the
    median across all posting books, and flag any book whose probability
    differs from that median by >= threshold_pp percentage points.

    A market with only one posting venue has no median to compare against,
    so it gets no flags (outlier_count = 0) — still listed for completeness.

    Returns ``QueryLocked`` sentinel instead of raising on DuckDB lock
    contention; see module docstring.
    """
    def _query() -> List[Dict[str, Any]]:
        with read_connection() as con:
            rows = con.execute(
                f"""
                WITH latest AS (
                  SELECT market_id, book, devig_prob,
                         ROW_NUMBER() OVER (PARTITION BY market_id, book ORDER BY fetched_at DESC) AS rn
                  FROM draft_odds
                  WHERE fetched_at > NOW() - INTERVAL '{MAX_AGE_HOURS} hours'
                )
                SELECT market_id, book, devig_prob FROM latest WHERE rn = 1
                """
            ).fetchall()

        by_market: Dict[str, Dict[str, float]] = {}
        for market_id, book, prob in rows:
            by_market.setdefault(market_id, {})[book] = prob

        threshold = threshold_pp / 100.0
        output: List[Dict[str, Any]] = []
        for market_id, books in by_market.items():
            if len(books) < 2:
                output.append({
                    "market_id": market_id,
                    "books": books,
                    "median": list(books.values())[0] if books else None,
                    "flags": {},
                    "outlier_count": 0,
                })
                continue
            median = statistics.median(books.values())
            flags = {b: abs(p - median) >= threshold for b, p in books.items()}
            output.append({
                "market_id": market_id,
                "books": books,
                "median": median,
                "flags": flags,
                "outlier_count": sum(flags.values()),
            })
        return output

    return _safe_read(_query)


def ev_candidates(threshold_pp: float = 10.0):
    """Flat list of flagged (market, venue) outliers, sorted by |delta| desc.

    Each row represents a single book/market combo that sits more than
    threshold_pp points away from the cross-book median. Signed delta
    preserves direction: positive = book is higher than consensus (bet NO),
    negative = book is lower than consensus (bet YES).

    Delegates to ``cross_book_grid``; propagates its ``QueryLocked``
    sentinel when the read fails.
    """
    grid = cross_book_grid(threshold_pp)
    if isinstance(grid, QueryLocked):
        return grid
    out: List[Dict[str, Any]] = []
    for m in grid:
        if m["median"] is None or not m["flags"]:
            continue
        for book, flagged in m["flags"].items():
            if not flagged:
                continue
            delta = m["books"][book] - m["median"]
            out.append({
                "market_id": m["market_id"],
                "book": book,
                "book_prob": m["books"][book],
                "median": m["median"],
                "delta": delta,
            })
    out.sort(key=lambda r: abs(r["delta"]), reverse=True)
    return out


def trade_tape(limit: int = 200, large_threshold_usd: float = 500.0, min_size_usd: float = 0.0):
    """Recent Kalshi trades with human-readable market context.

    LEFT JOIN against ``market_info`` (populated by the legacy Kalshi
    fetcher) to surface the market's title + subtitle. When the JOIN
    misses (ticker predates the market_info backfill), title/subtitle
    come back ``None`` and the dashboard falls back to the raw ticker.

    Args:
        limit: max rows to return.
        large_threshold_usd: fills >= this flagged ``is_large`` for
            highlighting. Not a filter — all rows still returned.
        min_size_usd: HARD filter floor — only trades with
            ``notional_usd >= min_size_usd`` are returned. Default 0
            (no filter) preserves existing behavior.

    Returns ``QueryLocked`` sentinel on lock contention.
    """
    def _query() -> List[Dict[str, Any]]:
        with read_connection() as con:
            rows = con.execute(
                """
                SELECT t.trade_id, t.ticker, t.side, t.price_cents, t.count,
                       t.notional_usd, t.traded_at,
                       t.notional_usd >= ? AS is_large,
                       mi.title AS market_title,
                       mi.subtitle AS market_subtitle
                FROM kalshi_trades t
                LEFT JOIN market_info mi ON mi.ticker = t.ticker
                WHERE t.notional_usd >= ?
                ORDER BY t.traded_at DESC
                LIMIT ?
                """,
                [large_threshold_usd, min_size_usd, limit],
            ).fetchall()
        cols = ["trade_id", "ticker", "side", "price_cents", "count",
                "notional_usd", "traded_at", "is_large",
                "market_title", "market_subtitle"]
        return [dict(zip(cols, r)) for r in rows]

    return _safe_read(_query)


def bet_log_rows():
    """All logged bets, most-recent first.

    Returns ``QueryLocked`` sentinel on lock contention.
    """
    def _query() -> List[Dict[str, Any]]:
        with read_connection() as con:
            rows = con.execute(
                """
                SELECT bet_id, market_id, book, side, american_odds, stake_usd, taken_at, note
                FROM draft_bets
                ORDER BY taken_at DESC
                """
            ).fetchall()
        cols = ["bet_id", "market_id", "book", "side", "american_odds", "stake_usd", "taken_at", "note"]
        return [dict(zip(cols, r)) for r in rows]

    return _safe_read(_query)


def all_market_ids():
    """Every known market_id from either the market catalog or the odds feed.

    Returns ``QueryLocked`` sentinel on lock contention.
    """
    def _query() -> List[str]:
        with read_connection() as con:
            rows = con.execute(
                """
                SELECT DISTINCT market_id FROM draft_markets
                UNION
                SELECT DISTINCT market_id FROM draft_odds
                """
            ).fetchall()
        return sorted({r[0] for r in rows if r[0]})

    return _safe_read(_query)


def latest_max_fetched_at(table: str):
    """Return MAX(fetched_at) for a table. Used by the cheap-poll guard to
    decide whether the tab needs to re-render on the interval tick.

    Returns ``QueryLocked`` sentinel on lock contention (callers should
    treat this like None — no new data info available).
    """
    def _query() -> Optional[datetime]:
        with read_connection() as con:
            result = con.execute(f"SELECT MAX(fetched_at) FROM {table}").fetchone()
        return result[0] if result and result[0] else None

    return _safe_read(_query)
