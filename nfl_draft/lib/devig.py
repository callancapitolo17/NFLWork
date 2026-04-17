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
