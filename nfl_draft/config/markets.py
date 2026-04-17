"""Per-book market label -> canonical market_id mappings.

Populated incrementally during scraper reconnaissance (Tasks 16-20).
On first scraper run, unmapped rows land in draft_odds_unmapped -- add
entries here, then re-run `python -m nfl_draft.lib.seed`.
"""

# Schema: list of (book, book_label, book_subject, market_id) tuples.
# market_id must match what build_market_id() produces in lib/market_map.py.
MARKET_MAP: list[tuple[str, str, str, str]] = [
    # Kalshi: ticker prefix + candidate name
    # ("kalshi", "KXNFLDRAFT1", "Cam Ward", "pick_1_overall_cam-ward"),
    # ("kalshi", "KXNFLDRAFTQB", "Cam Ward", "first_qb_cam-ward"),
    # DK / FD / Bookmaker / Wagerzon entries added during reconnaissance.
]
