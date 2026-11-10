"""Kalshi cash pots ("exchange shards") for the Ticket's Place button and the bot.

Kalshi runs one matching engine per market group, and each keeps its own cash:
an order on a basketball market can only use the basketball pot's cash, even
with thousands in the main pot (docs.kalshi.com/getting_started/exchange_sharding,
2026-11-02: two college basketball Place orders were refused
`insufficient_balance` with $5,247 in pot 0 and $113.86 in pot 3). Kalshi's
own app moves cash between pots by itself; API orders do not. So before a
Place order this tops the market's pot up from the main pot (Cal, 2026-11-02:
option C, move exactly the shortfall, nothing moves back).

Inputs: the market's exchange_index (the runner's /kalshi_quote.json reads it
off GET /markets/<ticker>) and the order's worst-case cost.
Side effects: GET /portfolio/balance?exchange_index=N on Cal's main account;
when short, ONE POST /portfolio/intra_exchange_instance_transfer from pot 0
(money moves between Cal's own pots, never out of the account), then GETs of
that transfer until Kalshi says it is complete. Needs a key with transfer
permission. Writes no database.

The bot (taker_orders.py, Cal 2026-11-03: option A2) does the same inside its
own subaccount, but only after Kalshi refused an order `insufficient_balance`
(Cal 2026-11-06: no read before every order, which met Kalshi's rate limit):
the short pot is topped up by exactly the
shortfall from the bot's other pot (top_up_bot_pot), with the bot's own key,
which Kalshi lets move money only between that subaccount's own pots. No
money ever comes in from the primary account.
"""
from __future__ import annotations

import logging
import math
import time
import uuid
from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # Type hints only: taker_orders imports this module for the bot's top-up.
    from unabated_ticket.bets_service.taker_orders import KalshiOrderClient

log = logging.getLogger(__name__)

MAIN_SHARD = 0
# Kalshi's shard table (exchange_sharding, 2026-11-02); an index not named
# here still works, it only reads as "pot N".
SHARD_NAMES = {0: "main", 1: "combos", 2: "crypto/commodities", 3: "basketball/baseball/tennis"}
BALANCE_PATH = "/portfolio/balance"
TRANSFER_PATH = "/portfolio/intra_exchange_instance_transfer"
TRANSFER_STATUS_PATH = "/portfolio/intra_exchange_instance_transfers"
# Kalshi's transfer is accepted at once and completes moments later; Place
# waits this long for "complete" before giving up without ordering.
TRANSFER_WAIT_SEC = 5.0
TRANSFER_POLL_SEC = 0.25
# A subaccount-locked key gets HTTP 403 from the transfer status endpoint, so
# the bot pauses this long after Kalshi accepts its move, then resends.
BOT_TRANSFER_SETTLE_SEC = 1.0
# The main account's balance read meets Kalshi's rate limit now and then: one
# retry, then the order goes without the check (Kalshi still refuses a short pot).
RATE_LIMITED = 429
BALANCE_RETRY_SEC = 0.3
TRANSFER_COMPLETE = "complete"
# Kalshi's 403 code for a key without the scope a call needs.
INSUFFICIENT_SCOPE_CODE_PREFIX = "insufficient_scope"
CENTICENTS_PER_DOLLAR = 10_000
# Kalshi rounds an order's fee up to the cent, so the cost can top
# contracts x all-in price by under a cent.
FEE_ROUNDING_DOLLARS = 0.01


@dataclass(frozen=True)
class TopUp:
    """What topping up did. `refusal` set = do not send the order (in words);
    `note` = the money moved, for Cal's message (None = nothing moved)."""
    refusal: str | None = None
    note: str | None = None


def shard_name(shard: int) -> str:
    return SHARD_NAMES.get(shard, f"pot {shard}")


def order_cost_bound(count: int, all_in_price: float) -> float:
    """The most the order can cost in dollars: contracts x all-in price, up to
    the cent, plus Kalshi's fee rounding."""
    return math.ceil(count * all_in_price * 100 - 1e-9) / 100 + FEE_ROUNDING_DOLLARS


def top_up_shard(client: KalshiOrderClient, shard: int | None, cost: float,
                 sleep: Callable[[float], None] = time.sleep,
                 monotonic: Callable[[], float] = time.monotonic) -> TopUp:
    """Make sure the main account's pot `shard` holds `cost` dollars, moving the
    shortfall from the main pot when it does not. A main-pot market, an unknown
    pot or an unreadable balance sends the order as before (Kalshi itself still
    refuses a short pot); a top-up that fails or does not complete in time
    refuses the order."""
    if shard is None or shard == MAIN_SHARD:
        return TopUp()
    pot_cash = _shard_cash(client, shard, sleep)
    if pot_cash is None:
        return TopUp()  # _shard_cash logged why
    shortfall_cents = math.ceil((cost - pot_cash) * 100 - 1e-9)
    if shortfall_cents <= 0:
        return TopUp()
    return _transfer(client, MAIN_SHARD, shard, shortfall_cents, pot_cash, sleep, monotonic)


def top_up_bot_pot(client: KalshiOrderClient, subaccount: int, shard: int | None, cost: float,
                   bot_shards: tuple[int, ...],
                   sleep: Callable[[float], None] = time.sleep,
                   monotonic: Callable[[], float] = time.monotonic) -> TopUp:
    """Make sure the bot subaccount's pot `shard` holds `cost` dollars, moving
    the shortfall from its other pot in `bot_shards` when it does not. Unlike
    Place, every doubt refuses: a pot outside `bot_shards`, an unknown pot, an
    unreadable balance, or a move that fails or does not complete in time.
    The other pot is never read: Kalshi answers a second balance read right
    after the first with HTTP 429, so its refusal of the move is the check."""
    if shard is None:
        return TopUp(refusal="Kalshi did not say which cash pot the market draws on")
    if shard not in bot_shards:
        return TopUp(refusal=f"the market draws on the {shard_name(shard)} pot, which the bot's subaccount "
                             f"is not set up on (it has pots {list(bot_shards)})")
    pot_cash = _shard_cash(client, shard, sleep, subaccount)
    if pot_cash is None:
        return TopUp(refusal=f"the bot's {shard_name(shard)} pot balance is unreadable")
    shortfall_cents = math.ceil((cost - pot_cash) * 100 - 1e-9)
    if shortfall_cents <= 0:
        return TopUp()
    sources = [other for other in bot_shards if other != shard]
    if len(sources) != 1:
        return TopUp(refusal=f"the bot's {shard_name(shard)} pot is ${shortfall_cents / 100:.2f} short and it has "
                             f"{len(sources)} other pots to draw on, expected 1")
    return _transfer(client, sources[0], shard, shortfall_cents, pot_cash, sleep, monotonic, subaccount)


def _shard_cash(client: KalshiOrderClient, shard: int, sleep: Callable[[float], None],
                subaccount: int | None = None) -> float | None:
    """Dollars available in pot `shard` of the main account (or of `subaccount`),
    or None when Kalshi does not say (logged)."""
    path = f"{BALANCE_PATH}?exchange_index={shard}"
    if subaccount is not None:
        path += f"&subaccount={subaccount}"
    for attempt in range(2):
        try:
            status, reply = client.request("GET", path)
        except OSError as error:
            log.warning("kalshi pots: balance of pot %d got no answer: %s", shard, error)
            return None
        if status == RATE_LIMITED and attempt == 0:
            sleep(BALANCE_RETRY_SEC)
            continue
        cents = reply.get("balance") if isinstance(reply, dict) else None
        if status != 200 or not isinstance(cents, int) or isinstance(cents, bool):
            log.warning("kalshi pots: balance of pot %d: expected 200 with integer cents in balance, got HTTP %s %s",
                        shard, status, str(reply)[:300])
            return None
        return cents / 100
    return None


def _transfer(client: KalshiOrderClient, source: int, shard: int, shortfall_cents: int, pot_cash: float,
              sleep: Callable[[float], None], monotonic: Callable[[], float], subaccount: int | None = None) -> TopUp:
    """ONE transfer of `shortfall_cents` from pot `source` to pot `shard`, then
    a wait for Kalshi to finish it (for a subaccount, a fixed pause instead).
    Never re-sent: a subaccount's transfer takes no client_transfer_id (Kalshi
    rejects it there), so a resend could move the money twice."""
    name = shard_name(shard)
    amount = shortfall_cents / 100
    body = {"source": "event_contract", "destination": "event_contract",
            "source_exchange_shard": source, "destination_exchange_shard": shard,
            "amount": shortfall_cents * (CENTICENTS_PER_DOLLAR // 100)}
    if subaccount is None:
        body["client_transfer_id"] = str(uuid.uuid4())
    else:
        # A subaccount-locked key must name its own subaccount on both ends.
        body["source_subaccount"] = subaccount
        body["destination_subaccount"] = subaccount
    whose = "" if subaccount is None else f"subaccount {subaccount} "
    try:
        status, reply = client.request("POST", TRANSFER_PATH, body)
    except OSError as error:
        log.error("kalshi pots: %stop-up of $%.2f to pot %d got no answer: %s", whose, amount, shard, error)
        return TopUp(refusal=f"moving ${amount:.2f} into the {name} pot got no answer from Kalshi; "
                             f"check the pot before placing again.")
    transfer_id = reply.get("transfer_id") if isinstance(reply, dict) else None
    if status != 200 or not isinstance(transfer_id, str):
        log.warning("kalshi pots: %stop-up of $%.2f to pot %d refused: HTTP %s %s", whose, amount, shard, status,
                    str(reply)[:300])
        if _error_code(reply).startswith(INSUFFICIENT_SCOPE_CODE_PREFIX):
            if subaccount is not None:
                return TopUp(refusal=f"the {name} pot has ${pot_cash:.2f} and the bot's key may not move money "
                                     f"between its pots: re-run taker_subaccount_setup for a key with transfer permission.")
            return TopUp(refusal=f"the {name} pot has ${pot_cash:.2f} and the Place key may not move money between "
                                 f"pots: make a Place key with transfer permission (deploy/README.md step 12).")
        return TopUp(refusal=f"the {name} pot has ${pot_cash:.2f} and moving ${amount:.2f} into it from the "
                             f"{shard_name(source)} pot failed "
                             f"(HTTP {status}{': ' + _error_text(reply) if _error_text(reply) else ''}).")
    log.info("kalshi pots: %smoving $%.2f from pot %d to pot %d (transfer %s)", whose, amount, source, shard, transfer_id)
    if subaccount is not None:
        sleep(BOT_TRANSFER_SETTLE_SEC)
        return TopUp(note=f"sent ${amount:.2f} into the {name} pot first")
    deadline = monotonic() + TRANSFER_WAIT_SEC
    while True:
        if _transfer_complete(client, transfer_id):
            return TopUp(note=f"moved ${amount:.2f} into the {name} pot first")
        if monotonic() >= deadline:
            log.warning("kalshi pots: transfer %s to pot %d not complete after %.0f s", transfer_id, shard, TRANSFER_WAIT_SEC)
            return TopUp(refusal=f"${amount:.2f} is on its way into the {name} pot but Kalshi has not finished "
                                 f"moving it; place again in a few seconds.")
        sleep(TRANSFER_POLL_SEC)


def _transfer_complete(client: KalshiOrderClient, transfer_id: str) -> bool:
    try:
        status, reply = client.request("GET", f"{TRANSFER_STATUS_PATH}/{transfer_id}")
    except OSError:
        return False
    transfer = reply.get("transfer") if isinstance(reply, dict) else None
    return status == 200 and isinstance(transfer, dict) and transfer.get("status") == TRANSFER_COMPLETE


def _error_code(reply: object) -> str:
    error = reply.get("error") if isinstance(reply, dict) else None
    return str(error.get("code") or "") if isinstance(error, dict) else ""


def _error_text(reply: object) -> str:
    error = reply.get("error") if isinstance(reply, dict) else None
    if isinstance(error, dict):
        return str(error.get("message") or error.get("code") or "")
    return ""
