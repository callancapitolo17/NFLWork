import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root for kalshi_common
from kalshi_common.fair_value import _probit_devig_n

def american_to_prob(odds: float) -> float:
    odds = float(odds)
    return (-odds) / (-odds + 100) if odds < 0 else 100 / (odds + 100)

def devig(probs: list[float]) -> list[float]:
    return _probit_devig_n(list(probs))


def overround_reject(overround: float, *, band_min: float, band_max: float) -> str | None:
    """Feed-integrity gate on a two-way pair's raw implied sum, checked BEFORE
    devig (devig clips-and-solves any input, so a poisoned pair would come back
    as a plausible-but-wrong fair). Returns a reject reason or None to pass.

    "anchor_crossed": sum < 1 is arithmetically not a two-way market — one side
    refreshed mid-tick while the other is stale. Hard reject regardless of band.
    "anchor_overround": sum outside [band_min, band_max] is not a plausible
    single-book two-way quote (near-zero vig or blown vig)."""
    if overround < 1.0:
        return "anchor_crossed"
    if not (band_min <= overround <= band_max):
        return "anchor_overround"
    return None
