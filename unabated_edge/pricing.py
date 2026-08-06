import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root for kalshi_common
from kalshi_common.fair_value import _probit_devig_n

MIN_ABS_AMERICAN_ODDS = 100  # American odds are never in (-100, 100) — that gap is the
                             # moneyline's own "pick 'em" boundary, never a real price.


def american_to_prob(odds: float) -> float:
    """Devig input: American odds -> implied probability.

    line_american_price (feed.py) is the fail-closed choke point that keeps
    a malformed value (e.g. a stray decimal-odds price like 1.91) from ever
    reaching here — so an out-of-range value at this point means a bug
    upstream, not bad feed data. Raise loudly instead of silently producing
    a plausible-but-wrong probability (e.g. american_to_prob(1.91) would
    otherwise return ~0.98, a "fair" prob for garbage input)."""
    odds = float(odds)
    if -MIN_ABS_AMERICAN_ODDS < odds < MIN_ABS_AMERICAN_ODDS:
        raise ValueError(
            f"invalid American odds shape: {odds} (must be <= -{MIN_ABS_AMERICAN_ODDS} "
            f"or >= {MIN_ABS_AMERICAN_ODDS})"
        )
    return (-odds) / (-odds + 100) if odds < 0 else 100 / (odds + 100)

def devig(probs: list[float]) -> list[float]:
    return _probit_devig_n(list(probs))
