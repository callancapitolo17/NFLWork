"""Fixed-margin quote pricing: target a constant ROI per side.

ROI = (p - cost)/cost, where p = side win-prob and cost = bid + maker_fee.
Setting ROI = target gives cost = p/(1+target), so bid = p/(1+target) - fee.
The fee depends on the bid, so we iterate once.  Because the maker fee is a
stepped (ceil-to-cent) function, flooring the bid to the grid or crossing a
fee tier can leave bid+fee(bid) > raw, pushing realized ROI below the target.
A post-floor guard steps the bid down one tick at a time until cost ≤ raw,
guaranteeing realized ROI ≥ target_roi by construction.
"""
import math
from dataclasses import dataclass

from kalshi_common.ev_calc import maker_fee_per_contract

GRID_STEP = 0.001  # Kalshi MVE markets price in $0.001 (deci-cent) steps


@dataclass(frozen=True)
class Quote:
    yes_bid: float
    no_bid: float


def _round_down_to_grid(price: float) -> float:
    # Floor to the grid: a lower bid is conservative (cheaper for us / more margin).
    return math.floor(round(price / GRID_STEP, 6)) * GRID_STEP


def _price_for_side(p: float, target_roi: float) -> float:
    raw = p / (1.0 + target_roi)             # target all-in cost (bid + fee)
    bid = raw - maker_fee_per_contract(raw)
    bid = raw - maker_fee_per_contract(bid)  # one refinement (fee depends on bid)
    bid = _round_down_to_grid(max(bid, 0.0))
    # The maker fee is stepped (ceil-to-cent), so flooring or a fee-tier
    # crossing can leave bid+fee(bid) > raw, i.e. cost above target → ROI below
    # target. Step down one tick until the all-in cost is within target (or the
    # bid hits zero, in which case quote() declines via its <=0 guard). This
    # guarantees realized ROI >= target_roi by construction (cost <= raw).
    while bid > 0 and bid + maker_fee_per_contract(bid) > raw + 1e-12:
        bid = _round_down_to_grid(bid - GRID_STEP)
    return max(bid, 0.0)


def quote(fair: float, target_roi: float) -> "Quote | None":
    if not (0.0 < fair < 1.0):
        return None
    yes_bid = _price_for_side(fair, target_roi)
    no_bid = _price_for_side(1.0 - fair, target_roi)
    if yes_bid <= 0 or no_bid <= 0 or yes_bid + no_bid >= 1.0:
        return None
    return Quote(yes_bid=round(yes_bid, 4), no_bid=round(no_bid, 4))
