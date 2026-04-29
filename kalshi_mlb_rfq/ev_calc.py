"""Kalshi taker fee + post-fee EV math.

Fee formula (Feb 2026 schedule): ceil(0.07 * P * (1-P) * 100) / 100 per contract,
where P is contract price in dollars.
"""

import math


def fee_per_contract(price_dollars: float) -> float:
    """Quadratic fee, rounded UP to the nearest cent."""
    return math.ceil(0.07 * price_dollars * (1 - price_dollars) * 100) / 100


def post_fee_ev_buy_yes(blended_fair: float, no_bid: float) -> tuple[float, float]:
    """We BUY YES by accepting a maker's no_bid (selling NO to them = buying YES from them).

    Effective YES purchase price = 1 - no_bid.
    Returns (ev_dollars_per_contract, ev_pct_of_stake).
    """
    yes_ask = 1.0 - no_bid
    if yes_ask <= 0 or yes_ask >= 1:
        return 0.0, 0.0
    fee = fee_per_contract(yes_ask)
    ev = blended_fair * (1 - yes_ask) - (1 - blended_fair) * yes_ask - fee
    return ev, ev / yes_ask


def post_fee_ev_buy_no(blended_fair: float, yes_bid: float) -> tuple[float, float]:
    """We BUY NO by accepting a maker's yes_bid (selling YES to them = buying NO from them).

    Effective NO purchase price = 1 - yes_bid.
    Returns (ev_dollars_per_contract, ev_pct_of_stake).
    """
    no_ask = 1.0 - yes_bid
    if no_ask <= 0 or no_ask >= 1:
        return 0.0, 0.0
    fee = fee_per_contract(no_ask)
    fair_no = 1.0 - blended_fair
    ev = fair_no * (1 - no_ask) - blended_fair * no_ask - fee
    return ev, ev / no_ask
