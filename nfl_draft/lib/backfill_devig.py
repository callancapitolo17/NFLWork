"""One-shot: re-compute draft_odds.devig_prob using the bucketed devig rule.

Usage (run once from main after merge, then delete this file):
    python -m nfl_draft.lib.backfill_devig

Groups rows by (date_trunc(minute, fetched_at), book, outright_group_key)
so concurrent scrape batches don't split. Skips Kalshi rows (their
devig_prob was always correct — written as mid by the scraper).
"""

from collections import defaultdict
from nfl_draft.lib.db import read_connection, write_connection
from nfl_draft.lib.devig import devig_n_way
from nfl_draft.lib.market_map import outright_group_key


def _load_rows():
    """Return (ts_bucket, book, market_id, american_odds, fetched_at) for every non-Kalshi draft_odds row.

    market_group is not persisted in draft_odds; we re-derive it from
    market_id structure via _derive_market_group below.
    """
    with read_connection() as con:
        return con.execute(
            """
            SELECT date_trunc('minute', o.fetched_at) AS ts_bucket,
                   o.book, o.market_id, o.american_odds, o.fetched_at
            FROM draft_odds o
            WHERE o.book != 'kalshi'
            ORDER BY o.fetched_at
            """
        ).fetchall()


def _derive_market_group(market_id: str) -> str:
    """Heuristic fallback: derive market_group from market_id prefix.

    We don't persist market_group in draft_odds, so for the backfill we
    reverse-engineer it from market_id structure. This is narrower than
    outright_group_key's input requirements but sufficient for backfill.
    """
    import re
    if market_id.startswith("pick_") and "_overall_" in market_id:
        return "pick_outright"
    if market_id.startswith("first_"):
        return "first_at_position"
    if market_id.startswith("top_"):
        m = re.match(r"^top_(\d+)_", market_id)
        return f"top_{m.group(1)}_range" if m else ""
    if market_id.startswith("team_") and "_first_pick_" in market_id:
        return "team_first_pick"
    if "_first_pick_pos_" in market_id:
        return "team_first_pick_position"
    if market_id.startswith("draft_position_ou_"):
        return "draft_position_over_under"
    if market_id.startswith("matchup_"):
        return "matchup_before"
    if market_id.startswith("mr_irrelevant_"):
        return "mr_irrelevant_position"
    if market_id.startswith("prop_"):
        return "prop_unknown"
    m = re.match(r"^(\d+)_[a-z]+_", market_id)
    if m:
        return f"nth_at_position_{m.group(1)}"
    return ""


def run():
    print("Loading draft_odds rows (excluding Kalshi)...")
    rows = _load_rows()
    print(f"  Loaded {len(rows)} rows.")

    buckets = defaultdict(list)
    ungrouped = []
    for ts_bucket, book, market_id, american_odds, fetched_at in rows:
        mg = _derive_market_group(market_id)
        group = outright_group_key(mg, market_id)
        if group is None:
            ungrouped.append((market_id, book, fetched_at))
            continue
        buckets[(ts_bucket, book, group)].append((market_id, american_odds, fetched_at))

    print(f"  Bucketed into {len(buckets)} groups; {len(ungrouped)} rows ungrouped (kept as-is).")

    updates = 0
    with write_connection() as con:
        for (ts_bucket, book, group), items in buckets.items():
            if len(items) < 2:
                continue
            devigged = devig_n_way([ao for _, ao, _ in items])
            for (market_id, _, fetched_at), dev in zip(items, devigged):
                con.execute(
                    """
                    UPDATE draft_odds
                    SET devig_prob = ?
                    WHERE market_id = ? AND book = ? AND fetched_at = ?
                    """,
                    [dev, market_id, book, fetched_at],
                )
                updates += 1
    print(f"  Updated devig_prob on {updates} rows.")
    print("Done. Delete this file after committing.")


if __name__ == "__main__":
    run()
