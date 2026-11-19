"""American -> decimal odds. Pure."""
from __future__ import annotations

import math


def american_to_decimal(american_odds: float) -> float:
    """+150 -> 2.5, -200 -> 1.5. Venues send fractional American odds (+217.3)."""
    if american_odds is None or not math.isfinite(american_odds) or -100 < american_odds < 100:
        raise ValueError(f"expected American odds <= -100 or >= +100, got {american_odds!r}")
    if american_odds >= 100:
        return 1.0 + american_odds / 100.0
    return 1.0 + 100.0 / -american_odds

