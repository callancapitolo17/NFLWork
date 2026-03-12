"""
Quoting logic: converts answer key fair values into Kalshi bid/ask prices.

Orderbook-aware: 5% minimum EV floor + penny-the-book for top-of-book priority.
"""

import math
from config import (
    MIN_EV_PCT, SKEW_PER_CONTRACT, MIN_QUOTE_SPREAD_CENTS,
    MIN_FAIR_VALUE, MAX_FAIR_VALUE, CONTRACT_SIZE
)


def compute_quotes(fair_prob, net_position=0, book_bid=0, book_ask=0):
    """Compute bid/ask in Kalshi cents given a fair probability and orderbook.

    Uses 5% minimum EV floor, then pennies the book to get top-of-book
    priority while maintaining edge.

    Args:
        fair_prob: Model's fair probability (0-1) for the YES side
        net_position: Current net YES contracts (positive = long YES)
        book_bid: Current best YES bid on Kalshi (cents, 0 if empty)
        book_ask: Current best YES ask on Kalshi (cents, 0 if empty)

    Returns:
        dict with bid_yes, ask_yes, etc., or None if no quote should be posted.
    """
    # Guard against NaN/invalid probabilities
    if fair_prob is None or not (0 <= fair_prob <= 1):
        return None

    fair_cents = fair_prob * 100

    # Don't quote at extremes where model is unreliable
    if fair_cents < MIN_FAIR_VALUE or fair_cents > MAX_FAIR_VALUE:
        return None

    # --- 5% EV floor ---
    # Bid: EV% = (fair - bid) / bid >= MIN_EV_PCT → bid <= fair / (1 + MIN_EV_PCT)
    max_bid = math.floor(fair_cents / (1 + MIN_EV_PCT))
    # Ask: EV% = (ask - fair) / (100 - ask) >= MIN_EV_PCT
    #   → ask >= (fair + 100 * MIN_EV_PCT) / (1 + MIN_EV_PCT)
    min_ask = math.ceil((fair_cents + 100 * MIN_EV_PCT) / (1 + MIN_EV_PCT))

    # Inventory skew: shift both quotes to attract offsetting flow
    skew = net_position * SKEW_PER_CONTRACT
    max_bid -= skew
    min_ask -= skew

    # --- Orderbook-aware pricing: penny the book ---
    if book_bid > 0 and book_bid < max_bid:
        # Improve best bid by 1c, but cap at our EV floor
        bid_yes = min(book_bid + 1, max_bid)
    else:
        # No book or book already at/above our limit
        bid_yes = max_bid

    if book_ask > 0 and book_ask > min_ask:
        # Improve best ask by 1c, but floor at our EV limit
        ask_yes = max(book_ask - 1, min_ask)
    else:
        # No book or book already at/below our limit
        ask_yes = min_ask

    # Safety: never bid above fair, never ask below fair
    bid_yes = min(bid_yes, math.floor(fair_cents) - 1)
    ask_yes = max(ask_yes, math.ceil(fair_cents) + 1)

    # Enforce minimum spread
    if (ask_yes - bid_yes) < MIN_QUOTE_SPREAD_CENTS:
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
        "book_bid": book_bid,
        "book_ask": book_ask,
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
        f"book=[{quote.get('book_bid', 0)}/{quote.get('book_ask', 0)}]  "
        f"size={quote['size']}"
    )
