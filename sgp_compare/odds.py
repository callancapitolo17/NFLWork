"""American / decimal / probability conversions. Pure."""
from __future__ import annotations

import math


def american_to_decimal(american_odds: float) -> float:
    """+150 -> 2.5, -200 -> 1.5. Venues send fractional American odds (+217.3)."""
    if american_odds is None or not math.isfinite(american_odds) or -100 < american_odds < 100:
        raise ValueError(f"expected American odds <= -100 or >= +100, got {american_odds!r}")
    if american_odds >= 100:
        return 1.0 + american_odds / 100.0
    return 1.0 + 100.0 / -american_odds


def decimal_to_american(decimal_odds: float) -> float:
    """2.5 -> +150, 1.5 -> -200 (unrounded)."""
    if not math.isfinite(decimal_odds) or decimal_odds <= 1.0:
        raise ValueError(f"expected decimal odds > 1, got {decimal_odds!r}")
    if decimal_odds >= 2.0:
        return (decimal_odds - 1.0) * 100.0
    return -100.0 / (decimal_odds - 1.0)


def probability_to_american(probability: float) -> float:
    if not 0.0 < probability < 1.0:
        raise ValueError(f"expected a probability in (0, 1), got {probability!r}")
    return decimal_to_american(1.0 / probability)
