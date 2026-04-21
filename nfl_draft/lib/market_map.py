"""Per-book market label -> canonical market_id mapping + market_id construction."""

from typing import Optional
from nfl_draft.lib.db import read_connection


def slug(name: str) -> str:
    """Canonicalize a name to a URL-safe slug (used in market_id construction)."""
    return name.lower().replace(" ", "-").replace(".", "").replace("'", "")


def build_market_id(market_type: str, **kwargs) -> str:
    """Deterministic market_id constructor -- must be applied identically in every scraper.

    Tested by test_market_map.py.
    """
    if market_type == "first_at_position":
        return f"first_{kwargs['position'].lower()}_{slug(kwargs['player'])}"
    if market_type == "pick_outright":
        return f"pick_{kwargs['pick_number']}_overall_{slug(kwargs['player'])}"
    if market_type == "top_n_range":
        return f"top_{kwargs['range_high']}_{slug(kwargs['player'])}"
    if market_type == "team_first_pick":
        return f"team_{kwargs['team'].lower()}_first_pick_{slug(kwargs['player'])}"
    if market_type == "prop":
        return f"prop_{slug(kwargs['short_description'])}"
    raise ValueError(f"Unknown market_type: {market_type}")


def resolve_market_id(book: str, book_label: str, book_subject: str) -> Optional[str]:
    """Look up canonical market_id for a per-book (label, subject) pair."""
    with read_connection() as con:
        result = con.execute(
            "SELECT market_id FROM market_map WHERE book = ? AND book_label = ? AND book_subject = ?",
            [book, book_label, book_subject],
        ).fetchone()
    return result[0] if result else None
