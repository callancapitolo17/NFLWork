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


# Any row older than MAX_AGE_MINUTES is excluded from Cross-Book Grid,
# +EV Candidates, and the Kalshi tooltip. This prevents venues that have
# silently stopped scraping (e.g. FD regression 2026-04-19) from polluting
# the grid with hours-old prices that get flagged as "edges" when the feed
# has actually died. 20 minutes = pre-draft scrape cadence (15 min from
# crontab.pre) + 5-minute cushion for a late cron run; it is also ~10x the
# draft-day cadence (2 min from crontab.draft), so dead venues drop out
# within ~20 min on draft day instead of the prior ~2h.
MAX_AGE_MINUTES = 20


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

    Median is ``statistics.median(devig_prob across posting venues)`` — the
    true cross-book fair. Flags compare each venue's **raw take price**
    (``implied_prob``) to that median. A flagged cell is a direct +EV signal:
    the venue is offering a price that differs from the cross-venue fair by
    more than ``threshold_pp`` percentage points.

    Signed delta (``implied_prob - median_fair``) preserves direction:
      * delta > 0: price above fair -> YES is overpriced -> bet NO
      * delta < 0: price below fair -> YES is underpriced -> bet YES

    Kalshi is treated uniformly: its ``implied_prob`` is ``yes_ask / 100``
    (the take price), so the same formula applies. Rows with
    ``implied_prob`` NULL (one-sided Kalshi with no ask) participate in
    the median via ``devig_prob`` but can't be flagged.

    A market with only one posting venue has no median to compare against,
    so it gets no flags — still listed for completeness.

    Returns ``QueryLocked`` sentinel instead of raising on DuckDB lock
    contention; see module docstring.
    """
    def _query() -> List[Dict[str, Any]]:
        with read_connection() as con:
            rows = con.execute(
                f"""
                WITH latest AS (
                  SELECT market_id, book, american_odds, implied_prob, devig_prob,
                         ROW_NUMBER() OVER (PARTITION BY market_id, book ORDER BY fetched_at DESC) AS rn
                  FROM draft_odds
                  WHERE fetched_at > NOW() - INTERVAL '{MAX_AGE_MINUTES} minutes'
                )
                SELECT market_id, book, american_odds, implied_prob, devig_prob
                FROM latest WHERE rn = 1
                """
            ).fetchall()

        by_market: Dict[str, Dict[str, Dict[str, Any]]] = {}
        for market_id, book, american_odds, implied_prob, devig_prob in rows:
            by_market.setdefault(market_id, {})[book] = {
                "american_odds": american_odds,
                "implied_prob": implied_prob,
                "devig_prob": devig_prob,
            }

        threshold = threshold_pp / 100.0
        output: List[Dict[str, Any]] = []
        for market_id, books in by_market.items():
            # Median uses devig_prob (fair). Pull every non-null devig for the median.
            fair_probs = [r["devig_prob"] for r in books.values() if r["devig_prob"] is not None]

            if len(books) < 2 or not fair_probs:
                # Only one posting venue -> no median -> no flags.
                output.append({
                    "market_id": market_id,
                    "books": books,
                    "median": fair_probs[0] if fair_probs else None,
                    "flags": {},
                    "outlier_count": 0,
                })
                continue

            median = statistics.median(fair_probs)
            flags: Dict[str, bool] = {}
            for book, r in books.items():
                # Uniform rule: flag if |raw take price - median fair| >= threshold.
                take = r["implied_prob"]
                flags[book] = (take is not None and abs(take - median) >= threshold)
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

    Each row represents a book/market combo where the book's **raw take
    price** sits more than threshold_pp from the cross-book median fair.
    Signed delta (``implied_prob - median_fair``) preserves direction and
    equals the EV in percentage points:
      * positive -> price above fair -> bet NO
      * negative -> price below fair -> bet YES

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
            book_prob = m["books"][book]["implied_prob"]
            delta = book_prob - m["median"]
            out.append({
                "market_id": m["market_id"],
                "book": book,
                "book_prob": book_prob,
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
        # Aggregation rationale
        # ---------------------
        # Kalshi reports each order-match as its own ``trade_id``. A single
        # logical trade (e.g. a 500-contract taker order that matched against
        # 11 resting makers) shows up as 11 rows with the same
        # ``(ticker, traded_at, price_cents, side)``. Without aggregation, the
        # Trade Tape renders those 11 fills as 11 separate rows — visually
        # implying 11 trades when there was really just one.
        #
        # Taker-side pricing
        # ------------------
        # We store ``price_cents`` as the YES-side price regardless of the
        # taker's ``side`` (see ``parse_trades_response``). For a NO taker
        # that means the cost per contract is ``100 - price_cents`` — if we
        # display the raw value, a 1-cent NO fill reads as "99¢" which is
        # wrong. Flip it at read time so the dashboard always shows what the
        # taker actually paid. The raw yes_price stays in the DB (no schema
        # migration).
        with read_connection() as con:
            rows = con.execute(
                """
                WITH aggregated AS (
                    SELECT
                        t.ticker,
                        t.traded_at,
                        t.side,
                        t.price_cents AS yes_price_cents,
                        CASE WHEN t.side = 'yes' THEN t.price_cents
                             ELSE 100 - t.price_cents
                        END AS taker_price_cents,
                        SUM(t.count) AS total_count,
                        -- Cast to DOUBLE so downstream Python gets floats
                        -- (integer * 0.01 yields DECIMAL in DuckDB, which
                        -- crosses the FFI as decimal.Decimal and breaks
                        -- arithmetic in tests / callbacks).
                        CAST(SUM(
                            t.count
                            * (CASE WHEN t.side = 'yes' THEN t.price_cents
                                    ELSE 100 - t.price_cents END)
                            * 0.01
                        ) AS DOUBLE) AS taker_notional_usd,
                        MIN(t.trade_id) AS trade_id,
                        MIN(t.fetched_at) AS fetched_at
                    FROM kalshi_trades t
                    GROUP BY t.ticker, t.traded_at, t.side, t.price_cents
                )
                -- market_info can contain duplicate rows per ticker
                -- (legacy data issue in the production DB). A naive LEFT JOIN
                -- against it multiplies each aggregated trade by the dup
                -- factor. De-duplicate to one row per ticker before joining.
                , market_info_dedup AS (
                    SELECT ticker, ANY_VALUE(title) AS title,
                           ANY_VALUE(subtitle) AS subtitle
                    FROM market_info
                    GROUP BY ticker
                )
                SELECT
                    a.trade_id,
                    a.ticker,
                    a.side,
                    a.taker_price_cents AS price_cents,
                    a.total_count AS count,
                    a.taker_notional_usd AS notional_usd,
                    a.traded_at,
                    a.taker_notional_usd >= ? AS is_large,
                    mi.title AS market_title,
                    mi.subtitle AS market_subtitle
                FROM aggregated a
                LEFT JOIN market_info_dedup mi ON mi.ticker = a.ticker
                WHERE a.taker_notional_usd >= ?
                ORDER BY a.traded_at DESC
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


def kalshi_tooltip_data() -> Dict[str, Dict[str, Any]]:
    """Return per-Kalshi-market hover content keyed by market_id.

    For every ``market_map`` row where ``book='kalshi'``, join to the latest
    ``kalshi_odds`` snapshot (by fetch_time) for the corresponding ticker. The
    join key is ``(book_label, book_subject)`` because market_map stores the
    ticker *prefix* (e.g. ``KXNFLDRAFTPICK-26-5``) while the full ticker has a
    candidate shortcode appended (``KXNFLDRAFTPICK-26-5-CTAT``).

    Returns:
        ``{market_id: {"ticker": str, "yes_bid": int, "yes_ask": int,
                       "last_price": int}}``
        Rows where no matching kalshi_odds snapshot exists are omitted.
        Returns ``QueryLocked`` on DuckDB lock contention.
    """
    def _query() -> Dict[str, Dict[str, Any]]:
        with read_connection() as con:
            rows = con.execute(
                f"""
                WITH latest_ko AS (
                  SELECT ticker, candidate, yes_bid, yes_ask, last_price,
                         ROW_NUMBER() OVER (PARTITION BY ticker ORDER BY fetch_time DESC) AS rn
                  FROM kalshi_odds
                  WHERE fetch_time > NOW() - INTERVAL '{MAX_AGE_MINUTES} minutes'
                )
                SELECT mm.market_id, ko.ticker, ko.yes_bid, ko.yes_ask, ko.last_price
                FROM market_map mm
                JOIN latest_ko ko
                  ON ko.ticker LIKE mm.book_label || '-%'
                 AND ko.candidate = mm.book_subject
                 AND ko.rn = 1
                WHERE mm.book = 'kalshi'
                """
            ).fetchall()
        return {
            market_id: {
                "ticker": ticker,
                "yes_bid": yes_bid,
                "yes_ask": yes_ask,
                "last_price": last_price,
            }
            for market_id, ticker, yes_bid, yes_ask, last_price in rows
        }

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
