"""Write OddsRow batches to draft_odds with quarantine for unmapped rows.

Devig is applied at ingest time: rows are bucketed by (book,
outright_group_key(market_group, market_id)) and proportional_devig runs
per bucket. Rows whose scraper already pre-set ``devig_prob`` (Kalshi)
pass through untouched. Single-row buckets and rows whose market_group
maps to no group key (props, matchups, mr_irrelevant) get
``devig_prob = implied_prob`` — no devig, comparison against the
cross-venue median still works against raw implieds.
"""

from typing import List, Tuple, Dict
from nfl_draft.lib.db import write_connection, read_connection
from nfl_draft.lib.devig import american_to_implied, devig_n_way
from nfl_draft.lib.market_map import outright_group_key
from nfl_draft.scrapers._base import OddsRow


def write_or_quarantine(rows: List[OddsRow]) -> Tuple[int, int]:
    """Resolve each row's market_id; write mapped to draft_odds (with devig
    applied per outright group), unmapped to quarantine.

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

    # Resolve market_ids up front; collect mapped rows so we can bucket-devig them.
    mapped: List[Tuple[OddsRow, str]] = []  # (row, market_id)
    unmapped_rows: List[OddsRow] = []
    for row in rows:
        mid = lookup.get((row.book, row.book_label, row.book_subject))
        if mid is None:
            unmapped_rows.append(row)
        else:
            mapped.append((row, mid))

    # Bucket mapped rows by (book, outright_group_key). Rows whose group
    # resolves to None go into solo buckets keyed on identity so they remain
    # 1-row buckets (no devig).
    buckets: Dict[object, List[int]] = {}
    for idx, (row, mid) in enumerate(mapped):
        group = outright_group_key(row.market_group or "", mid)
        key: object = (row.book, group) if group is not None else ("__solo__", idx)
        buckets.setdefault(key, []).append(idx)

    # Compute devig_prob per mapped row index.
    devig_by_idx: Dict[int, float] = {}
    for key, idx_list in buckets.items():
        bucket_rows = [mapped[i] for i in idx_list]
        # If ANY row in the bucket arrived with devig_prob pre-set (Kalshi),
        # respect every row's scraper-supplied values verbatim. row.book is the
        # per-scraper isolation: within one (book, group) bucket every row came
        # from the same scraper, so either all have devig_prob pre-set (Kalshi)
        # or none do (sportsbooks). The fallback for a "mixed" bucket is
        # belt-and-suspenders and shouldn't fire in practice.
        if any(row.devig_prob is not None for row, _ in bucket_rows):
            for i, (row, _) in zip(idx_list, bucket_rows):
                if row.devig_prob is not None:
                    devig_by_idx[i] = row.devig_prob
                else:
                    # Shouldn't happen in practice; fall back to implied.
                    implied = (row.implied_prob
                               if row.implied_prob is not None
                               else american_to_implied(row.american_odds))
                    devig_by_idx[i] = implied
            continue
        if len(bucket_rows) >= 2:
            # n-way proportional devig over the bucket's American odds.
            devigged = devig_n_way([row.american_odds for row, _ in bucket_rows])
            for i, dev in zip(idx_list, devigged):
                devig_by_idx[i] = dev
        else:
            # 1-row bucket (or market_group with no group key): pass implied through.
            row, _ = bucket_rows[0]
            devig_by_idx[idx_list[0]] = american_to_implied(row.american_odds)

    # Single write pass: unmapped -> quarantine, mapped -> draft_odds.
    with write_connection() as con:
        for row in unmapped_rows:
            con.execute(
                "INSERT INTO draft_odds_unmapped VALUES (?, ?, ?, ?, ?, ?)",
                [row.book, row.book_label, row.book_subject,
                 row.american_odds, row.fetched_at, "no market_map entry"],
            )
        for idx, (row, mid) in enumerate(mapped):
            implied = (row.implied_prob
                       if row.implied_prob is not None
                       else american_to_implied(row.american_odds))
            devig = devig_by_idx[idx]
            con.execute(
                "INSERT INTO draft_odds VALUES (?, ?, ?, ?, ?, ?)",
                [mid, row.book, row.american_odds, implied, devig, row.fetched_at],
            )
    return (len(mapped), len(unmapped_rows))
