"""Write OddsRow batches to draft_odds with quarantine for unmapped rows.

Devig is applied at ingest time: rows are bucketed by (book,
outright_group_key(market_group, market_id)) and a devig helper runs per
bucket. Which helper depends on the bucket shape:

  * ``top_<N>`` bucket -- ``devig_pool(odds, N)``: N-winner pool, so
    posted YES implieds are normalized to sum to N. Sparsity guardrail:
    if ``sum(implieds) < TOP_N_POOL_MIN_COVERAGE * N`` the posted set is
    too thin to trust, so fall back to raw implieds.
  * Any other multi-row bucket -- ``devig_n_way(odds)``: classic 1-winner
    outright, normalized to sum to 1.
  * Single-row bucket / no group key -- raw ``implied_prob`` passes
    through as ``devig_prob``. Cross-venue median still compares take
    prices directly.

Rows whose scraper already pre-set ``devig_prob`` (Kalshi) pass through
untouched regardless of bucket shape.
"""

import re
from typing import List, Tuple, Dict
from nfl_draft.lib.db import write_connection, read_connection
from nfl_draft.lib.devig import american_to_implied, devig_n_way, devig_pool
from nfl_draft.lib.market_map import outright_group_key
from nfl_draft.scrapers._base import OddsRow

# Sparsity guardrail for top_N pool devig: only trust pool normalization
# when posted vigged implieds sum to at least this fraction of the pool
# size. Below that, the posted candidate set is too thin to contain ~all
# the top-N probability mass and normalizing would inflate every fair.
# 0.9 leaves headroom for typical hold (~5-8% on round-1 boards) and
# still catches thin top_5 boards like Bookmaker (which often post <10
# candidates and sum to well under 5).
TOP_N_POOL_MIN_COVERAGE = 0.9

# Group keys shaped ``top_<N>`` come from outright_group_key for
# top_N_range markets; the N parses out here so quarantine can size the
# pool. Tightly coupled to market_map.outright_group_key's output format.
_TOP_N_GROUP_RE = re.compile(r"^top_(\d+)$")


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
        if len(bucket_rows) < 2:
            # 1-row bucket (or market_group with no group key): pass implied through.
            row, _ = bucket_rows[0]
            devig_by_idx[idx_list[0]] = american_to_implied(row.american_odds)
            continue
        odds = [row.american_odds for row, _ in bucket_rows]
        # Detect top_N pool buckets. ``key`` is either ``(book, group)`` or
        # ``("__solo__", idx)``; only the former can be a top_N group.
        pool_n: int | None = None
        if isinstance(key, tuple) and isinstance(key[1], str):
            m = _TOP_N_GROUP_RE.match(key[1])
            if m:
                pool_n = int(m.group(1))
        if pool_n is not None:
            implieds = [american_to_implied(o) for o in odds]
            coverage = sum(implieds)
            if coverage >= TOP_N_POOL_MIN_COVERAGE * pool_n:
                devigged = devig_pool(odds, pool_n)
            else:
                # Posted candidate set too thin to plausibly contain ~all
                # of the pool_n probability mass. Pool-normalizing would
                # inflate every fair; pass raw implieds through instead so
                # the cross-venue median still has something meaningful to
                # compare against.
                devigged = implieds
        else:
            # Classic 1-winner outright: normalize to sum=1.
            devigged = devig_n_way(odds)
        for i, dev in zip(idx_list, devigged):
            devig_by_idx[i] = dev

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
