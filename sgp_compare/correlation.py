"""How correlated a venue thinks a parlay's legs are. Pure: no network, no clock.

Definitions (the page shows each one):
  * parlay chance      = 1 / the venue's parlay decimal odds.
  * independent chance = the product of each leg's own chance, i.e. what the
                         parlay would pay if the legs had nothing to do with
                         each other.
  * Correlation factor = parlay chance / independent chance - 1. "+19%" means
                         the venue thinks the legs hit together 19% more often
                         than unrelated legs would; negative means they work
                         against each other.

It is computed twice:
  * with vig    — straight off the prices you can take: the parlay price over
                  the single prices multiplied together. The parlay's margin and
                  the singles' margins are both in it, so it mixes correlation
                  with how much more the venue charges on a parlay.
  * vig removed — from the venue's own parlay prices for EVERY combination of
                  the legs' outcomes (2^N prices, 2-3 legs only). Each
                  combination's chance is divided by their total (removing the
                  parlay margin), and the parlay's chance is compared with the
                  product of each leg's chance read off those same combinations.
                  Using one price list for both sides of the ratio means a
                  margin spread evenly over the combinations cancels exactly;
                  mixing a devigged parlay with separately devigged singles
                  does not (review 2026-11-18: independent 20%/20% legs at 30%
                  parlay vig read -31% that way). None, with the reason, when a
                  combination is missing, a leg has no other side, a line can
                  push, or there are 4+ legs.

The singles are also devigged (kalshi_common's probit two-way devig) for the
"legs as separate bets, no vig" price shown next to the parlay.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

from kalshi_common.fair_value import two_way_fair

MAX_PARTITION_LEGS = 3      # 2^3 = 8 parlay requests per venue is the ceiling
# A coherent partition's chances add to at most 100% + 25% per leg (the same
# envelope as kalshi_common.fair_value.devig_partition); past that a cell is
# stale or garbage and the factor would be noise.
PARTITION_OVERROUND_PER_LEG = 0.25
# A coherent two-sided single at an exchange: its two chances add to 100-125%.
# Wider than kalshi_common's sportsbook default (100.5-120%) on both ends:
# exchange spreads can be a cent wide, and prop ladders are often thin.
SINGLE_OVERROUND_MIN = 1.0
SINGLE_OVERROUND_MAX = 1.25
SINGLE_REJECT_REASONS = {
    "crossed": "its two sides add to under 100%",
    "overround": "its two sides are more than 25% apart",
    "bad_price": "a price is not usable",
}


def partition_flips(n_legs: int) -> list[tuple[bool, ...]]:
    """Every combination of the legs' outcomes, the chosen parlay first.
    Entry i says, per leg, whether that leg is flipped to its other side."""
    return [tuple(flips) for flips in itertools.product((False, True), repeat=n_legs)]


@dataclass(frozen=True)
class LegChance:
    leg_id: str
    raw_chance: float | None       # 1 / decimal price to take it; None = nobody offers it
    fair_chance: float | None      # devigged against its other side
    fair_unavailable_reason: str | None


@dataclass(frozen=True)
class CorrelationReport:
    parlay_raw_chance: float
    independent_raw_chance: float | None        # product of the singles as offered
    correlation_raw: float | None               # with vig
    independent_fair_chance: float | None       # product of the devigged singles (display only)
    parlay_fair_chance: float | None            # the parlay's share of the outcome combinations
    partition_independent_chance: float | None  # product of each leg's share of the combinations
    correlation_fair: float | None              # vig removed
    partition_overround: float | None           # combinations' chances added up, minus 1 (the parlay vig)
    raw_unavailable_reason: str | None
    fair_unavailable_reason: str | None


def leg_chance(leg_id: str, decimal_odds: float | None, opposite_decimal_odds: float | None) -> LegChance:
    if decimal_odds is None:
        return LegChance(leg_id, None, None, "nobody is offering this side alone")
    raw = 1.0 / decimal_odds
    if opposite_decimal_odds is None:
        return LegChance(leg_id, raw, None, "nobody is offering its other side")
    fair, reject = two_way_fair(decimal_odds, opposite_decimal_odds,
                                band_min=SINGLE_OVERROUND_MIN, band_max=SINGLE_OVERROUND_MAX)
    return LegChance(leg_id, raw, fair, SINGLE_REJECT_REASONS.get(reject, reject))


def correlation_report(parlay_decimal_odds: float, legs: list[LegChance],
                       partition_decimal_odds: list[float] | None,
                       partition_unavailable_reason: str | None) -> CorrelationReport:
    """`partition_decimal_odds` follows `partition_flips` order (index 0 = the
    parlay itself) and is complete, or is None with `partition_unavailable_reason`."""
    parlay_raw = 1.0 / parlay_decimal_odds
    unpriced = [leg for leg in legs if leg.raw_chance is None]
    independent_raw = None if unpriced else math.prod(leg.raw_chance for leg in legs)
    correlation_raw = parlay_raw / independent_raw - 1.0 if independent_raw else None
    raw_reason = f"{len(unpriced)} leg(s) not offered alone, so there is no separate-bets price" if unpriced else None

    fair_singles = [leg.fair_chance for leg in legs]
    independent_fair = math.prod(fair_singles) if all(c is not None for c in fair_singles) else None

    if partition_decimal_odds is None:
        if partition_unavailable_reason is None:
            raise ValueError("expected partition prices or a reason they are missing")
        parlay_fair = partition_independent = overround = None
        fair_reason = partition_unavailable_reason
    else:
        parlay_fair, partition_independent, overround, fair_reason = _partition_shares(
            partition_decimal_odds, len(legs))
    correlation_fair = parlay_fair / partition_independent - 1.0 if parlay_fair is not None else None
    return CorrelationReport(parlay_raw, independent_raw, correlation_raw, independent_fair, parlay_fair,
                             partition_independent, correlation_fair, overround, raw_reason, fair_reason)


def _partition_shares(partition_decimal_odds: list[float], n_legs: int
                      ) -> tuple[float | None, float | None, float | None, str | None]:
    """(parlay's share, product of the legs' shares, overround, reason)."""
    flips = partition_flips(n_legs)
    if len(partition_decimal_odds) != len(flips):
        raise ValueError(f"expected {len(flips)} partition prices for {n_legs} legs, "
                         f"got {len(partition_decimal_odds)}")
    raw = [1.0 / d for d in partition_decimal_odds]
    total = sum(raw)
    if total <= 1.0:
        return None, None, total - 1.0, (f"the outcome combinations add up to {total:.1%} (100% or less), "
                                         "so the prices are crossed")
    if total > 1.0 + PARTITION_OVERROUND_PER_LEG * n_legs:
        return None, None, total - 1.0, (f"the outcome combinations add up to {total:.0%}, more than "
                                         f"{1 + PARTITION_OVERROUND_PER_LEG * n_legs:.0%}, so a price is stale or wrong")
    shares = [r / total for r in raw]
    leg_shares = [sum(share for share, cell in zip(shares, flips) if not cell[leg]) for leg in range(n_legs)]
    return shares[0], math.prod(leg_shares), total - 1.0, None
