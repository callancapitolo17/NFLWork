"""Write OddsRow batches to draft_odds with quarantine for unmapped rows."""

from typing import List, Tuple
from nfl_draft.lib.db import write_connection, read_connection
from nfl_draft.lib.devig import american_to_implied
from nfl_draft.scrapers._base import OddsRow


def write_or_quarantine(rows: List[OddsRow]) -> Tuple[int, int]:
    """Resolve each row's market_id; write mapped to draft_odds, unmapped to quarantine.

    Returns: (mapped_count, unmapped_count).
    """
    # Pre-fetch all market_id mappings before opening write connection to avoid
    # mixing read + write connections in same context (DuckDB restriction).
    market_ids = {}
    for row in rows:
        key = (row.book, row.book_label, row.book_subject)
        if key not in market_ids:
            with read_connection() as con:
                result = con.execute(
                    "SELECT market_id FROM market_map WHERE book = ? AND book_label = ? AND book_subject = ?",
                    [row.book, row.book_label, row.book_subject],
                ).fetchone()
            market_ids[key] = result[0] if result else None

    mapped = 0
    unmapped = 0
    with write_connection() as con:
        for row in rows:
            key = (row.book, row.book_label, row.book_subject)
            mid = market_ids[key]
            if mid is None:
                con.execute(
                    "INSERT INTO draft_odds_unmapped VALUES (?, ?, ?, ?, ?, ?)",
                    [row.book, row.book_label, row.book_subject, row.american_odds, row.fetched_at, "no market_map entry"],
                )
                unmapped += 1
                continue
            implied = american_to_implied(row.american_odds)
            con.execute(
                "INSERT INTO draft_odds VALUES (?, ?, ?, ?, ?, ?)",
                [mid, row.book, row.american_odds, implied, implied, row.fetched_at],
            )
            mapped += 1
    return (mapped, unmapped)
