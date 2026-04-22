"""Write OddsRow batches to draft_odds with quarantine for unmapped rows."""

from typing import List, Tuple
from nfl_draft.lib.db import write_connection, read_connection
from nfl_draft.lib.devig import american_to_implied
from nfl_draft.scrapers._base import OddsRow


def write_or_quarantine(rows: List[OddsRow]) -> Tuple[int, int]:
    """Resolve each row's market_id; write mapped to draft_odds, unmapped to quarantine.

    Pre-fetches the entire market_map table once before the write loop so that
    per-row resolution is an O(1) in-memory dict lookup. Previous implementation
    opened a fresh read_connection() for each unique (book, book_label,
    book_subject) key -- with hundreds of unique keys per scrape, connection
    setup was the dominant cost.

    If ``row.implied_prob`` / ``row.devig_prob`` are set, they are written
    verbatim to draft_odds; otherwise both columns are derived from
    ``american_to_implied(row.american_odds)`` (legacy sportsbook behavior).

    Returns: (mapped_count, unmapped_count).
    """
    if not rows:
        return (0, 0)

    # Single read pass -- load every market_map entry into memory.
    with read_connection() as con:
        map_rows = con.execute(
            "SELECT book, book_label, book_subject, market_id FROM market_map"
        ).fetchall()
    lookup = {(b, l, s): m for b, l, s, m in map_rows}

    mapped = 0
    unmapped = 0
    with write_connection() as con:
        for row in rows:
            mid = lookup.get((row.book, row.book_label, row.book_subject))
            if mid is None:
                con.execute(
                    "INSERT INTO draft_odds_unmapped VALUES (?, ?, ?, ?, ?, ?)",
                    [row.book, row.book_label, row.book_subject, row.american_odds, row.fetched_at, "no market_map entry"],
                )
                unmapped += 1
                continue
            implied = (
                row.implied_prob
                if row.implied_prob is not None
                else american_to_implied(row.american_odds)
            )
            devig = (
                row.devig_prob
                if row.devig_prob is not None
                else implied
            )
            con.execute(
                "INSERT INTO draft_odds VALUES (?, ?, ?, ?, ?, ?)",
                [mid, row.book, row.american_odds, implied, devig, row.fetched_at],
            )
            mapped += 1
    return (mapped, unmapped)
