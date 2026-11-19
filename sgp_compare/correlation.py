"""How correlated a venue thinks a parlay's legs are. Pure: no network, no clock.

Definitions (the page shows each one):
  * parlay chance      = 1 / the venue's parlay decimal odds.
  * independent chance = the product of each leg's own chance at the venue,
                         i.e. what the parlay would pay if the legs had nothing
                         to do with each other.
  * % correlated       = parlay chance / independent chance - 1. "+19%" means the
                         venue thinks the legs hit together 19% more often than
                         unrelated legs would; negative means they work against
                         each other.

It is computed twice:
  * with vig    — straight off the prices you can take. The parlay's margin and
                  the singles' margins are both in it, so it mixes correlation
                  with how much more the venue charges on a parlay.
  * vig removed — each single devigged against its other side (probit), and the
                  parlay devigged across every combination of its legs'
                  outcomes (2^N parlay prices). Only for 2-3 legs whose every
                  combination the venue priced; otherwise None, never guessed.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

from scipy.optimize import brentq
from scipy.stats import norm

MAX_PARTITION_LEGS = 3      # 2^3 = 8 parlay requests per venue is the ceiling
PROBABILITY_CLIP = 1e-9
# A coherent two-sided single: its two chances add to between 100% and 125%.
TWO_WAY_OVERROUND_MIN = 1.0
TWO_WAY_OVERROUND_MAX = 1.25


def probit_devig(raw_probabilities: list[float]) -> list[float]:
    """Remove the vig from a set of mutually exclusive, exhaustive outcomes by
    shifting every outcome's z-score by one constant until they sum to 1.
    Same method as kalshi_common.fair_value (probit handles longshots better
    than dividing by the overround)."""
    if len(raw_probabilities) < 2:
        raise ValueError(f"expected >= 2 outcomes to devig, got {len(raw_probabilities)}")
    z_scores = [norm.ppf(min(max(p, PROBABILITY_CLIP), 1 - PROBABILITY_CLIP)) for p in raw_probabilities]
    if len(z_scores) == 2:
        shift = -(z_scores[0] + z_scores[1]) / 2
    else:
        shift = brentq(lambda c: sum(norm.cdf(z + c) for z in z_scores) - 1.0, -10, 10, xtol=1e-12)
    return [float(norm.cdf(z + shift)) for z in z_scores]


def partition_flips(n_legs: int) -> list[tuple[bool, ...]]:
    """Every combination of the legs' outcomes, the chosen parlay first.
    Entry i says, per leg, whether that leg is flipped to its other side."""
    return [tuple(flips) for flips in itertools.product((False, True), repeat=n_legs)]


@dataclass(frozen=True)
class LegChance:
    leg_id: str
    raw_chance: float | None      # 1 / decimal price to take it; None = nobody offers it
    fair_chance: float | None     # devigged against its other side; None = no other side


@dataclass(frozen=True)
class CorrelationReport:
    parlay_raw_chance: float
    independent_raw_chance: float | None
    correlation_raw: float | None           # with vig
    parlay_fair_chance: float | None
    independent_fair_chance: float | None
    correlation_fair: float | None          # vig removed
    partition_overround: float | None       # sum of every combination's chance - 1 (the parlay vig)
    fair_unavailable_reason: str | None


def leg_chance(leg_id: str, decimal_odds: float | None, opposite_decimal_odds: float | None) -> LegChance:
    """A side whose two prices add to under 100% (crossed) or over
    TWO_WAY_OVERROUND_MAX (a stale or lopsided quote, e.g. Collins 71.5 at
    -128 / -138 is 14% wide; -400 / -200 is 47%) gets no fair chance rather than a misleading one."""
    raw = 1.0 / decimal_odds if decimal_odds else None
    if raw is None or not opposite_decimal_odds:
        return LegChance(leg_id, raw, None)
    two_way_total = raw + 1.0 / opposite_decimal_odds
    if not TWO_WAY_OVERROUND_MIN <= two_way_total <= TWO_WAY_OVERROUND_MAX:
        return LegChance(leg_id, raw, None)
    fair, _ = probit_devig([raw, 1.0 / opposite_decimal_odds])
    return LegChance(leg_id, raw, fair)


def _product(values: list[float]) -> float:
    return math.prod(values)


def correlation_report(parlay_decimal_odds: float, legs: list[LegChance],
                       partition_decimal_odds: list[float | None] | None) -> CorrelationReport:
    """`partition_decimal_odds` follows `partition_flips` order (index 0 = the
    parlay itself) or is None when the partition was not priced."""
    parlay_raw = 1.0 / parlay_decimal_odds
    raw_singles = [leg.raw_chance for leg in legs]
    independent_raw = _product(raw_singles) if all(c is not None for c in raw_singles) else None
    correlation_raw = parlay_raw / independent_raw - 1.0 if independent_raw else None

    fair_singles = [leg.fair_chance for leg in legs]
    independent_fair = _product(fair_singles) if all(c is not None for c in fair_singles) else None

    parlay_fair, overround, reason = _devig_partition(partition_decimal_odds, len(legs))
    if independent_fair is None and reason is None:
        reason = "a leg has no usable price on its other side (none offered, crossed, or wider than 25%)"
    correlation_fair = (parlay_fair / independent_fair - 1.0
                        if parlay_fair is not None and independent_fair else None)
    return CorrelationReport(parlay_raw, independent_raw, correlation_raw, parlay_fair,
                             independent_fair, correlation_fair, overround, reason)


def _devig_partition(partition_decimal_odds: list[float | None] | None,
                     n_legs: int) -> tuple[float | None, float | None, str | None]:
    if n_legs > MAX_PARTITION_LEGS:
        return None, None, f"vig removal needs every outcome combination; only done up to {MAX_PARTITION_LEGS} legs"
    if partition_decimal_odds is None:
        return None, None, "a leg has no other side to flip"
    expected = 2 ** n_legs
    if len(partition_decimal_odds) != expected:
        raise ValueError(f"expected {expected} partition prices for {n_legs} legs, got {len(partition_decimal_odds)}")
    if any(d is None for d in partition_decimal_odds):
        missing = sum(d is None for d in partition_decimal_odds)
        return None, None, f"the venue would not price {missing} of the {expected} outcome combinations"
    raw = [1.0 / d for d in partition_decimal_odds]
    overround = sum(raw) - 1.0
    if overround <= 0:
        return None, overround, f"outcome combinations add up to {sum(raw):.3f} (< 1), prices are crossed"
    return probit_devig(raw)[0], overround, None
