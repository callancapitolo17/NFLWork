import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from kalshi_common.ev_calc import fee_per_contract


def edge_for_yes(fair_prob: float, yes_ask: float) -> tuple[float, float]:
    fee = fee_per_contract(yes_ask)
    ev_d = fair_prob * (1 - yes_ask) - (1 - fair_prob) * yes_ask - fee
    return ev_d, (ev_d / yes_ask if yes_ask > 0 else 0.0)
