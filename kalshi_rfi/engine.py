"""Pure quote-decision logic for the one-sided YRFI maker.

Network-free and config-free (everything injected) so it is fully
unit-testable. The daemon (main.py) owns all I/O.
"""
import math
from dataclasses import dataclass

from kalshi_common.ev_calc import maker_fee_per_contract


@dataclass(frozen=True)
class QuoteDecision:
    action: str                 # "quote" | "no_quote"
    reason: str
    price_cents: int | None = None
    count: int | None = None


def desired_price_cents(fair_yes: float, yes_ask_cents: int | None,
                        margin_cents: int, min_edge_cents: float) -> tuple[int | None, str]:
    """The YES bid we want resting, or (None, reason).

    price = floor(fair − margin), clamped to stay a MAKER order (never at or
    through the ask — a cross would take liquidity and pay taker fees, a
    different trade than the one this bot is sized for).
    """
    if not (0.0 < fair_yes < 1.0):
        return None, "fair_out_of_range"
    price = math.floor(fair_yes * 100.0) - margin_cents
    if yes_ask_cents is not None and 0 < yes_ask_cents <= 99:
        price = min(price, yes_ask_cents - 1)
    if price < 1:
        return None, "price_below_1c"
    # Post-fee edge floor. maker_fee is dollars/contract -> cents.
    fee_cents = maker_fee_per_contract(price / 100.0) * 100.0
    edge_cents = fair_yes * 100.0 - price - fee_cents
    if edge_cents < min_edge_cents:
        return None, "edge_below_floor"
    return price, "ok"


def size_contracts(price_cents: int, per_game_cap_usd: float,
                   game_committed_usd: float,
                   daily_remaining_usd: float) -> int:
    """Contracts for a YES bid. Worst case of a YES buy = its cost, so the
    caps are enforced on cost: min(per-game headroom, daily headroom)."""
    if price_cents < 1:
        return 0
    headroom = min(per_game_cap_usd - game_committed_usd, daily_remaining_usd)
    if headroom <= 0:
        return 0
    return max(0, math.floor(headroom * 100.0 / price_cents))


def decide(fair_yes: float | None, yes_ask_cents: int | None, *,
           margin_cents: int, min_edge_cents: float,
           per_game_cap_usd: float, game_committed_usd: float,
           daily_remaining_usd: float,
           sec_to_start: float, pull_before_start_sec: float) -> QuoteDecision:
    """Full decision for one game this cycle.

    game_committed_usd / daily_remaining_usd must EXCLUDE the resting order
    this decision replaces (it is cancelled before the new one is placed) but
    INCLUDE fills and every other game's resting cost — otherwise a reprice
    could never size, or a burst of fills could breach the caps.

    "quote": the caller ensures an order rests at (price, count), replacing
    any resting order whose price differs (count drift alone is not worth the
    cancel/replace round trip). "no_quote": cancel anything resting — with
    reason "caps_exhausted" that is correct, not just safe: fills alone have
    consumed the cap, so a resting add-on order would breach it.
    """
    if sec_to_start <= pull_before_start_sec:
        return QuoteDecision("no_quote", "past_pull_deadline")
    if fair_yes is None:
        return QuoteDecision("no_quote", "no_consensus")
    price, reason = desired_price_cents(fair_yes, yes_ask_cents,
                                        margin_cents, min_edge_cents)
    if price is None:
        return QuoteDecision("no_quote", reason)
    count = size_contracts(price, per_game_cap_usd, game_committed_usd,
                           daily_remaining_usd)
    if count < 1:
        return QuoteDecision("no_quote", "caps_exhausted")
    return QuoteDecision("quote", "ok", price_cents=price, count=count)
