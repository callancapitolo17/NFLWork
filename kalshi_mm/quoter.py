"""
Quoting logic: Avellaneda-Stoikov inspired pricing for Kalshi binary contracts.

Reservation price shifts quotes against inventory to attract offsetting flow.
Dynamic half-spread widens near tipoff to protect against informed counterparties.
Penny-the-book for top-of-book priority when edge allows.
Anti-penny-loop: detects when a counterparty is walking our price and stops chasing.
"""

import math
import time
import datetime as _dt
from datetime import timezone
from collections import defaultdict
from config import (
    MIN_QUOTE_SPREAD_CENTS, MIN_FAIR_VALUE, MAX_FAIR_VALUE, CONTRACT_SIZE,
    AS_GAMMA, AS_BASE_HALF_SPREAD, AS_URGENCY_SCALE, AS_URGENCY_CLAMP_HRS,
    BANKROLL
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


def _hours_to_tipoff(commence_time):
    """Returns hours until tipoff (float, >=0), or None if unavailable."""
    if not commence_time:
        return None
    try:
        if isinstance(commence_time, str):
            # Kalshi API returns RFC 3339 strings like "2026-03-29T19:00:00Z"
            tip = _dt.datetime.fromisoformat(commence_time.replace("Z", "+00:00"))
        else:
            tip = commence_time
        if tip.tzinfo is None:
            tip = tip.replace(tzinfo=timezone.utc)
        delta = (tip - _dt.datetime.now(timezone.utc)).total_seconds()
        return max(0.0, delta / 3600)
    except Exception:
        return None


def compute_quotes(fair_prob, net_position=0, book_bid=0, book_ask=0,
                   bid_size=None, ask_size=None, commence_time=None):
    """Compute bid/ask using Avellaneda-Stoikov reservation price + dynamic spread.

    The reservation price shifts fair value against our inventory — if we're
    long, it drops below fair to make our ask cheaper (attracting sells to us).
    The half-spread widens near tipoff to protect against informed flow.

    Args:
        fair_prob: Model's fair probability (0-1) for the YES side
        net_position: Current net YES contracts (positive = long YES)
        book_bid: Current best YES bid on Kalshi (cents, 0 if empty)
        book_ask: Current best YES ask on Kalshi (cents, 0 if empty)
        bid_size: Contract count for bid (passed through, not used for pricing)
        ask_size: Contract count for ask (passed through, not used for pricing)
        commence_time: Game start time (ISO string or datetime) for urgency calc

    Returns:
        dict with bid_yes, ask_yes, reservation, etc., or None if no quote.
    """
    # Guard against NaN/invalid probabilities
    if fair_prob is None or not (0 <= fair_prob <= 1):
        return None

    fair_cents = fair_prob * 100

    # Don't quote at extremes where model is unreliable
    if fair_cents < MIN_FAIR_VALUE or fair_cents > MAX_FAIR_VALUE:
        return None

    # --- Reservation price (AS inventory adjustment) ---
    # σ² = p(1-p): binary variance — maximal at 50c, shrinks at extremes.
    # This naturally makes the bot more conservative at 50c (high uncertainty)
    # and more aggressive at 20c or 80c (low uncertainty).
    sigma_sq = fair_prob * (1 - fair_prob)

    # Convert contract position to dollar exposure at current fair value
    pos_dollars = net_position * (fair_cents / 100)
    pos_fraction = pos_dollars / BANKROLL

    # inventory_shift (cents): how far to move our midpoint away from fair.
    # Positive when long → reservation drops → cheaper ask → attracts sellers.
    # Negative when short → reservation rises → cheaper bid → attracts buyers.
    inventory_shift = pos_fraction * AS_GAMMA * sigma_sq * 100
    reservation = fair_cents - inventory_shift

    # --- Dynamic half-spread (time urgency) ---
    # Near tipoff: information asymmetry rises (lineups confirmed, sharp action).
    # Wider spread protects us from informed flow at the cost of fewer fills.
    hours = _hours_to_tipoff(commence_time)
    if hours is None:
        urgency = 1.0  # No tipoff info — use base spread only
    else:
        # normalized: 1.0 when hours=clamp (far away), 0.0 at tipoff
        normalized = min(hours, AS_URGENCY_CLAMP_HRS) / AS_URGENCY_CLAMP_HRS
        # urgency: 1.0 far from tipoff, up to (1 + AS_URGENCY_SCALE) at tipoff
        urgency = 1.0 + (1.0 - normalized) * AS_URGENCY_SCALE

    half_spread = round(AS_BASE_HALF_SPREAD * urgency)
    # Ensure half_spread is at least 1c (can't bid and ask at same price)
    half_spread = max(half_spread, 1)

    # --- Bid/ask from reservation ---
    bid_yes = round(reservation - half_spread)
    ask_yes = round(reservation + half_spread)

    # Safety: never cross fair value. If extreme inventory pushes the ask
    # below fair (e.g., long 1200 contracts → reservation = 46c, ask = 49c
    # but fair = 50c → ask is still above fair, OK). But if inventory_shift
    # exceeds half_spread, ask could dip below fair — clamp it.
    bid_yes = min(bid_yes, math.floor(fair_cents) - 1)
    ask_yes = max(ask_yes, math.ceil(fair_cents) + 1)

    # --- Orderbook-aware pricing: penny the book ---
    pennied = False
    if book_bid > 0 and book_bid < bid_yes:
        # Improve best bid by 1c, but cap at our AS limit
        bid_yes = min(book_bid + 1, bid_yes)
        pennied = True

    if book_ask > 0 and book_ask > ask_yes:
        # Improve best ask by 1c, but floor at our AS limit
        ask_yes = max(book_ask - 1, ask_yes)
        pennied = True

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
        "reservation": round(reservation, 2),
        "inventory_shift": round(inventory_shift, 2),
        "urgency": round(urgency, 3),
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
        f"spread={quote['spread']}c  "
        f"rsv={quote.get('reservation', '?')}  "
        f"inv_shift={quote.get('inventory_shift', '?')}  "
        f"urg={quote.get('urgency', '?')}  "
        f"book=[{quote.get('book_bid', 0)}/{quote.get('book_ask', 0)}]  "
        f"{size_str}{penny_flag}"
    )
