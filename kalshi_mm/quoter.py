"""
Quoting logic: converts answer key fair values into Kalshi bid/ask prices.

Orderbook-aware: 5% minimum EV floor + penny-the-book for top-of-book priority.
Anti-penny-loop: detects when a counterparty is walking our price and stops chasing.
"""

import math
import time
from collections import defaultdict
from config import (
    MIN_EV_PCT, SKEW_PER_CONTRACT, MIN_QUOTE_SPREAD_CENTS,
    MIN_FAIR_VALUE, MAX_FAIR_VALUE, CONTRACT_SIZE
)

# Anti-penny tracking: {ticker: [(timestamp, book_bid, book_ask), ...]}
# If the book moves by exactly 1c more than PENNY_LOOP_MAX times in
# PENNY_LOOP_WINDOW seconds, stop pennying and quote at EV floor instead.
_book_history = defaultdict(list)
PENNY_LOOP_WINDOW = 120   # 2 minutes
PENNY_LOOP_MAX = 4        # 4 x 1c moves in 2 min = someone is walking us


def _detect_penny_loop(ticker, book_bid, book_ask):
    """Detect if someone is pennying us in a loop.

    Returns True if we should stop pennying and quote at EV floor.
    """
    now = time.time()
    history = _book_history[ticker]

    # Prune old entries
    history[:] = [(t, b, a) for t, b, a in history if now - t < PENNY_LOOP_WINDOW]

    # Count 1c bid increases (someone walking our bid up)
    # AND 1c ask decreases (someone walking our ask down)
    penny_count = 0
    for i in range(1, len(history)):
        prev_bid = history[i - 1][1]
        curr_bid = history[i][1]
        prev_ask = history[i - 1][2]
        curr_ask = history[i][2]
        if curr_bid - prev_bid == 1 and prev_bid > 0:
            penny_count += 1
        if prev_ask - curr_ask == 1 and curr_ask > 0:
            penny_count += 1

    # Record current observation
    history.append((now, book_bid, book_ask))

    return penny_count >= PENNY_LOOP_MAX


def compute_quotes(fair_prob, net_position=0, book_bid=0, book_ask=0,
                   bid_size=None, ask_size=None):
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
    # But first check for penny-loop (counterparty walking our price)
    pennied = False
    if book_bid > 0 and book_bid < max_bid:
        # Improve best bid by 1c, but cap at our EV floor
        bid_yes = min(book_bid + 1, max_bid)
        pennied = True
    else:
        # No book or book already at/above our limit
        bid_yes = max_bid

    if book_ask > 0 and book_ask > min_ask:
        # Improve best ask by 1c, but floor at our EV limit
        ask_yes = max(book_ask - 1, min_ask)
        pennied = True
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
        "bid_size": bid_size if bid_size is not None else CONTRACT_SIZE,
        "ask_size": ask_size if ask_size is not None else CONTRACT_SIZE,
        "book_bid": book_bid,
        "book_ask": book_ask,
        "pennied": pennied,
    }


def should_amend(current_price, desired_price, threshold=1):
    """Check if a resting order should be amended (price moved enough)."""
    return abs(current_price - desired_price) >= threshold


def format_quote_summary(ticker, quote):
    """Format a quote for logging."""
    if quote is None:
        return f"  {ticker}: NO QUOTE (out of range or invalid)"
    penny_flag = " [PENNY]" if quote.get("pennied") else ""
    bid_sz = quote.get("bid_size", quote.get("size", "?"))
    ask_sz = quote.get("ask_size", quote.get("size", "?"))
    size_str = f"size={bid_sz}" if bid_sz == ask_sz else f"bid_sz={bid_sz} ask_sz={ask_sz}"
    return (
        f"  {ticker}: "
        f"fair={quote['fair_cents']:.1f}c  "
        f"bid={quote['bid_yes']}c  ask={quote['ask_yes']}c  "
        f"spread={quote['spread']}c  skew={quote['skew']}  "
        f"book=[{quote.get('book_bid', 0)}/{quote.get('book_ask', 0)}]  "
        f"{size_str}{penny_flag}"
    )
