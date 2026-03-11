"""
Quoting logic: converts answer key fair values into Kalshi bid/ask prices.

MVP: mechanical fair_value ± half_spread with inventory skew.
Phase 2: orderbook-aware quoting.
"""

import math
from config import (
    HALF_SPREAD_CENTS, SKEW_PER_CONTRACT, MIN_QUOTE_SPREAD_CENTS,
    MIN_FAIR_VALUE, MAX_FAIR_VALUE, CONTRACT_SIZE
)


def compute_quotes(fair_prob, net_position=0):
    """Compute bid/ask in Kalshi cents given a fair probability.

    Args:
        fair_prob: Model's fair probability (0-1) for the YES side
        net_position: Current net YES contracts (positive = long YES)

    Returns:
        dict with bid_yes, ask_yes, or None if no quote should be posted.
        Also returns the side labels for clarity.
    """
    fair_cents = fair_prob * 100

    # Don't quote at extremes where model is unreliable
    if fair_cents < MIN_FAIR_VALUE or fair_cents > MAX_FAIR_VALUE:
        return None

    half_spread = max(HALF_SPREAD_CENTS, MIN_QUOTE_SPREAD_CENTS // 2)

    # Inventory skew: shift both quotes to attract offsetting flow
    skew = net_position * SKEW_PER_CONTRACT

    bid_yes = math.floor(fair_cents - half_spread - skew)
    ask_yes = math.ceil(fair_cents + half_spread - skew)

    # Enforce: never bid above fair, never ask below fair
    bid_yes = min(bid_yes, math.floor(fair_cents) - 1)
    ask_yes = max(ask_yes, math.ceil(fair_cents) + 1)

    # Enforce minimum spread
    if (ask_yes - bid_yes) < MIN_QUOTE_SPREAD_CENTS:
        # Widen symmetrically
        mid = (bid_yes + ask_yes) / 2
        bid_yes = math.floor(mid - MIN_QUOTE_SPREAD_CENTS / 2)
        ask_yes = math.ceil(mid + MIN_QUOTE_SPREAD_CENTS / 2)

    # Clamp to valid Kalshi range
    bid_yes = max(1, min(99, bid_yes))
    ask_yes = max(1, min(99, ask_yes))

    # Final sanity: bid must be less than ask
    if bid_yes >= ask_yes:
        return None

    return {
        "bid_yes": bid_yes,
        "ask_yes": ask_yes,
        "fair_cents": round(fair_cents, 2),
        "spread": ask_yes - bid_yes,
        "skew": skew,
        "size": CONTRACT_SIZE,
    }


def should_amend(current_price, desired_price, threshold=1):
    """Check if a resting order should be amended (price moved enough)."""
    return abs(current_price - desired_price) >= threshold


def format_quote_summary(ticker, quote):
    """Format a quote for logging."""
    if quote is None:
        return f"  {ticker}: NO QUOTE (out of range or invalid)"
    return (
        f"  {ticker}: "
        f"fair={quote['fair_cents']:.1f}c  "
        f"bid={quote['bid_yes']}c  ask={quote['ask_yes']}c  "
        f"spread={quote['spread']}c  skew={quote['skew']}  "
        f"size={quote['size']}"
    )
