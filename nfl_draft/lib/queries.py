"""Dashboard query functions — separated from Dash callbacks for testability.

These are pure(ish) functions that read from DuckDB and return plain Python
dicts/lists. No Dash imports here, so they can be unit-tested with a fresh
temp DuckDB file in each test.
"""

import statistics
from typing import List, Dict, Any, Optional
from datetime import datetime

from nfl_draft.lib.db import read_connection


def cross_book_grid(threshold_pp: float = 10.0) -> List[Dict[str, Any]]:
    """Return one row per market with all-venue prices, median, and outlier flags.

    For each market, take the most-recent devig_prob per book, compute the
    median across all posting books, and flag any book whose probability
    differs from that median by >= threshold_pp percentage points.

    A market with only one posting venue has no median to compare against,
    so it gets no flags (outlier_count = 0) — still listed for completeness.
    """
    with read_connection() as con:
        rows = con.execute(
            """
            WITH latest AS (
              SELECT market_id, book, devig_prob,
                     ROW_NUMBER() OVER (PARTITION BY market_id, book ORDER BY fetched_at DESC) AS rn
              FROM draft_odds
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
            # Single venue: no cross-book signal possible
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


def ev_candidates(threshold_pp: float = 10.0) -> List[Dict[str, Any]]:
    """Flat list of flagged (market, venue) outliers, sorted by |delta| desc.

    Each row represents a single book/market combo that sits more than
    threshold_pp points away from the cross-book median. Signed delta
    preserves direction: positive = book is higher than consensus (bet NO),
    negative = book is lower than consensus (bet YES).
    """
    grid = cross_book_grid(threshold_pp)
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


def trade_tape(limit: int = 200, large_threshold_usd: float = 500.0) -> List[Dict[str, Any]]:
    """Recent Kalshi trades. Computes is_large at read time."""
    with read_connection() as con:
        rows = con.execute(
            """
            SELECT trade_id, ticker, side, price_cents, count, notional_usd, traded_at,
                   notional_usd >= ? AS is_large
            FROM kalshi_trades
            ORDER BY traded_at DESC
            LIMIT ?
            """,
            [large_threshold_usd, limit],
        ).fetchall()
    cols = ["trade_id", "ticker", "side", "price_cents", "count", "notional_usd", "traded_at", "is_large"]
    return [dict(zip(cols, r)) for r in rows]


def bet_log_rows() -> List[Dict[str, Any]]:
    """All logged bets, most-recent first."""
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


def all_market_ids() -> List[str]:
    """Every known market_id from either the market catalog or the odds feed."""
    with read_connection() as con:
        rows = con.execute(
            """
            SELECT DISTINCT market_id FROM draft_markets
            UNION
            SELECT DISTINCT market_id FROM draft_odds
            """
        ).fetchall()
    return sorted({r[0] for r in rows if r[0]})


def latest_max_fetched_at(table: str) -> Optional[datetime]:
    """Return MAX(fetched_at) for a table. Used by the cheap-poll guard to
    decide whether the tab needs to re-render on the interval tick."""
    with read_connection() as con:
        result = con.execute(f"SELECT MAX(fetched_at) FROM {table}").fetchone()
    return result[0] if result and result[0] else None
