"""Devig math — port of Answer Keys/Tools.R devigging functions."""

from typing import List, Tuple


def american_to_implied(american_odds: float) -> float:
    """Convert American odds to implied probability (vig included)."""
    o = float(american_odds)
    if o > 0:
        return 100.0 / (o + 100.0)
    return -o / (-o + 100.0)


def devig_two_way(odds_a: float, odds_b: float) -> Tuple[float, float]:
    """Remove vig from a two-way market (proportional method)."""
    a = american_to_implied(odds_a)
    b = american_to_implied(odds_b)
    total = a + b
    return (a / total, b / total)


def devig_n_way(odds_list: List[float]) -> List[float]:
    """Remove vig from an n-way market where the complementary set is fully posted."""
    implieds = [american_to_implied(o) for o in odds_list]
    total = sum(implieds)
    return [p / total for p in implieds]


def proportional_devig(odds_list: List[float]) -> List[float]:
    """Devig when the complementary set is INCOMPLETE (e.g., book only posts top 5).

    KNOWN LIMITATION: assumes posted outcomes contain ~100% of probability mass.
    Inflates the posted candidates' probs vs. ground truth. Spec annotates this
    in the dashboard.
    """
    return devig_n_way(odds_list)


def devig_pool(odds_list: List[float], pool_size: float) -> List[float]:
    """Devig an N-winner futures market where exactly ``pool_size`` outcomes win.

    Applies to markets like "will player X be drafted in top N" where up to N
    distinct players are YES winners and the YES-only sportsbook doesn't post
    a tradable NO side. Standard two-way devig is impossible, but the
    structural constraint is known: the true fair YES probabilities across
    ALL players sum to exactly ``pool_size``. If the posted candidates
    capture ~all of that probability mass and the book's vig is
    approximately proportional across outcomes, normalizing posted vigged
    implieds to sum to ``pool_size`` recovers the fair:

        fair_i = q_i * pool_size / sum(q_j)

    where ``q_i = american_to_implied(odds_i)``.

    Assumptions (caller must validate):
      * Posted outcomes contain ~all of the pool_size probability mass.
        If the book only posts a handful of candidates, sum(q) << pool_size
        and this formula will inflate each fair. Guard upstream (e.g. only
        call when ``sum(q) >= 0.9 * pool_size``).
      * Vig is proportional across outcomes. Real books steepen vig on
        longshots, so favorites are understated slightly, longshots
        overstated. Within ~1-2pp for round-1 draft candidates.
    """
    implieds = [american_to_implied(o) for o in odds_list]
    total = sum(implieds)
    if total == 0:
        return implieds
    return [p * pool_size / total for p in implieds]
