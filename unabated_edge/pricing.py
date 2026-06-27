import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root for kalshi_common
from kalshi_common.fair_value import _probit_devig_n

def american_to_prob(odds: float) -> float:
    odds = float(odds)
    return (-odds) / (-odds + 100) if odds < 0 else 100 / (odds + 100)

def devig(probs: list[float]) -> list[float]:
    return _probit_devig_n(list(probs))
