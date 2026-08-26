"""Structure-route ingest: one flight per (book, game), ~36 devigs (issue #96).

Inputs: an SGPService and the current slate. Outputs: SurfaceRows + counts.
Side effects: live HTTP to one book. No DB writes (the runner owns those).

Books on this route are FanDuel, BetMGM, Novig and Caesars (#95 coverage).
DraftKings carries no structure odds at all and ProphetX 403s at the events
stage on the first request of a session, so neither is here.

Cost shape, measured in #95: one cold structure fetch per game (median
0.31-0.42s) and then ~0.00s per additional leg on that game, with ZERO SGP
price calls — the whole ladder for a game comes out of one payload.
"""
from __future__ import annotations

import logging
import time

from kalshi_mlb_mm.leg_surface.devig import ExclusionCounts, devig_rung
from kalshi_mlb_mm.leg_surface.slate import rungs

log = logging.getLogger(__name__)

ROUTE = "structure"


def _decimals_for_rung(rung_legs, leg_index, odds) -> dict:
    """side -> decimal for one rung, from ``StructureLegOdds.odds``.

    Prefers a side's OWN resolved price; falls back to the other side's
    ``opposite_decimal`` when the book resolved only one leg of the pair but
    still published both prices.
    """
    decimals: dict = {}
    for leg in rung_legs:
        priced = odds.get(leg_index.get(leg))
        if priced is None:
            continue
        decimals[leg.side] = priced[0]
    for leg in rung_legs:
        priced = odds.get(leg_index.get(leg))
        if priced is None or priced[1] is None:
            continue
        other = next((l.side for l in rung_legs if l.side != leg.side), None)
        if other is not None and other not in decimals:
            decimals[other] = priced[1]
    return decimals


def price_game(book: str, service, game, *, skip_markets, band_min: float,
               band_max: float) -> tuple[list, ExclusionCounts, str]:
    """One (book, game): fetch structure once, devig every rung.

    ``skip_markets`` are the (market_type, period) pairs this book's SINGLES
    route owns — FanDuel's FG/F5 totals. Pricing them here too would write
    the same surface key from two routes: duplicate rows in the mirror, and
    in memory a coin flip between two slices. Skipped rungs are NOT counted
    as misses; this route was never asked for them.

    Returns (rows, counts, outcome). ``outcome`` is the StructureLegOdds
    outcome, so a dark book ('transport_error', 'no_event') is never
    mistaken for a book that simply offers none of these alt rungs.
    """
    counts = ExclusionCounts()
    legs = list(game.legs)
    # CanonicalLeg is a frozen dataclass, so it keys by VALUE; the
    # slate deduped legs, so the map is 1:1.
    leg_index = {leg: i for i, leg in enumerate(legs)}
    result = service.structure_leg_odds(book, game.game_ref(), legs)
    if result.outcome != "ok":
        # Whole-game miss. Charged once as game_unmatched, not once per rung:
        # inflating it by the ladder size would make one dark book look like
        # a coverage collapse across every line.
        counts.game_unmatched += 1
        return [], counts, result.outcome

    rows = []
    for (period, market_type, _line), rung_legs in rungs(legs).items():
        if (market_type, period) in skip_markets:
            continue
        decimals = _decimals_for_rung(rung_legs, leg_index, result.odds)
        outcome = devig_rung(book=book, route=ROUTE, game=game,
                             legs=rung_legs, decimals=decimals,
                             built_at=result.fetched_at, band_min=band_min,
                             band_max=band_max)
        if outcome.reason is not None:
            counts.bump(outcome.reason)
            continue
        rows.extend(outcome.rows)
    return rows, counts, "ok"


def run_pass(book: str, service, games, *, skip_markets=(),
             band_min: float, band_max: float, max_req_per_sec: float
             ) -> tuple[list, ExclusionCounts, int, int]:
    """One book's pass over the slate.

    Returns (rows, counts, games_priced, transport_failures).

    ``max_req_per_sec`` paces the per-game fetches. A whole-book pull is one
    fetch PER GAME — not the one request the epic's cost note reads it as —
    so a mistuned cadence is how this loop would talk itself into a 403.
    Pacing stretches a pass instead of firing it.
    """
    rows: list = []
    counts = ExclusionCounts()
    games_priced = 0
    transport_failures = 0
    skip_markets = set(skip_markets)
    min_gap_sec = (1.0 / max_req_per_sec) if max_req_per_sec > 0 else 0.0
    next_allowed = time.monotonic()
    for game in games:
        wait = next_allowed - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        next_allowed = time.monotonic() + min_gap_sec
        game_rows, game_counts, outcome = price_game(
            book, service, game, skip_markets=skip_markets,
            band_min=band_min, band_max=band_max)
        counts.add(game_counts)
        if game_rows:
            games_priced += 1
            rows.extend(game_rows)
        elif outcome in ("transport_error", "error"):
            transport_failures += 1
            log.warning("surface %s: %s failed on %s", book, outcome,
                        game.game_id)
    return rows, counts, games_priced, transport_failures
