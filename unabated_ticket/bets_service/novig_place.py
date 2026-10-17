"""The Ticket's Place button on Novig (novig_place_plan.md, Cal 2026-10-16):
Cal clicks, and this desk sends ONE order from Place's own Novig sub-account
(never the bot's, never his app balance: no key trades that). Nothing here runs
on a timer: every order is one POST /place_novig.json from the panel or the
phone page.

Rules, in order, for every request (validate_novig_place_request's shape):
  1. the click was not already sent: 'nman-' + its clickId, held in
     bets.duckdb::manual_orders before the send (Novig does not dedupe its
     clientId, which is the bare clickId since Novig wants a UUID, so this
     row is the only guard);
  2. the server runner finds the board's Novig line, its outcome id (Unabated's
     rung sourceData) and the outcome on the other side, reads Novig's live
     price for it (1 minus the best resting order on the other outcome) and
     checks the market is open and pregame (GET /novig_quote.json);
  3. the game has not started, and the runner read the Novig this desk trades
     on (its base URL);
  4. the live price is at most PRICE_TOLERANCE (0.1c) worse than the Ticket's;
     a better one is taken (a pregame take pays no fee, so price = all-in);
  5. no open straight Novig bet is on the other outcome (sizing cannot price
     a hedge, as on Kalshi);
  6. the Novig reads the stake was sized on are fresh and ok;
  7. no earlier Place on this game may have filled after the read the stake
     was sized on: a Novig Place after `novigReadAt`, a Kalshi Place after
     `kalshiReadAt` (the stake counts bets at both; manual_orders.unseen_manual_fill);
  8. contracts = floor(stake x 100 / price) (novig_orders.qty_for_stake): a
     Novig contract pays 1c, not $1;
  9. the row is INSERTED (status 'sending', with the board line, so Place's
     sub-account read can name the bet) BEFORE the send;
 10. ONE immediate-or-cancel buy at the live price (NovigOrderClient, the
     Novig taker's), never re-sent; then the order is read until final and
     its fills are read, the bot desk's way (novig_taker_desk.py). An order
     Novig shows OPEN (an IOC never should) is cancelled. One that never
     reads final still has its fills read: all its contracts filled is a
     definite 'filled' (it cannot fill more); else it stays 'unknown' with
     what has filled so far, and Place's sub-account read settles it once
     Novig shows it final (sources/novig_subaccount.py).

Order statuses (manual_orders.status on a Novig row):
  sending   inserted, the POST not yet answered
  filled | partial | unfilled   what the fills read found
  rejected  Novig refused it before taking it (a 451: the phone's location
            check is over 3 days old); nothing placed
  unknown   no definite answer, the order never read final, it rested, or its
            fills could not be read: it may hold contracts, so the next Place
            on the game waits for a read that saw it. One Novig never answered
            has no order id: the novig_place read looks it up by clientId a
            minute on and gives it its id, or marks it rejected (never sent)

Side effects: a Novig order from Place's sub-account (money); INSERT / UPDATE
bets.duckdb::manual_orders (venue 'novig'). Network: GET
<runner>/novig_quote.json, then Novig REST signed with Place's trading key only.
"""
import logging
import math
import threading
import time
from collections.abc import Callable
from datetime import datetime, timedelta, timezone

from unabated_ticket.bets_service.manual_orders import (
    PRICE_TOLERANCE, STATUS_FILLED, STATUS_PARTIAL, STATUS_REFUSED, STATUS_REJECTED, STATUS_UNFILLED, STATUS_UNKNOWN,
    VENUE_NOVIG, _cents, _is_number, _iso_time, optional_read_at, unseen_manual_fill,
)
from unabated_ticket.bets_service.novig_orders import (FINAL_ORDER_STATUSES, GEO_STATUS, QTY_PER_DOLLAR, STATUS_OPEN,
                                                       NovigError, NovigOrderClient, error_code_of, is_on_price_grid,
                                                       order_totals, qty_for_stake)
from unabated_ticket.bets_service.novig_taker_desk import (DEFINITE_REJECTION_STATUSES, MAX_HOLDINGS_AGE_SEC,
                                                           UUID_PATTERN, line_context_of)
from unabated_ticket.bets_service.store import BetsStore

log = logging.getLogger(__name__)

NOVIG_ORDER_ID_PREFIX = "nman-"
RUNNER_NOVIG_QUOTE_PATH = "/novig_quote.json"
# The runner's Novig quote makes up to four Novig reads (10 s each at most), so
# it gets longer than a Kalshi quote: 2026-10-17 it took 7-18 s on the VM and
# every first click failed at the 10 s default. The phone waits 90 s for Place.
RUNNER_NOVIG_QUOTE_TIMEOUT_SEC = 30
# After a 201 the order is read this many times, this far apart, for its final
# status (about 9 s). The bot desk's 3 s was too short for Place on 2026-10-17:
# Novig kept the order unread past it while its fill had already landed.
PLACE_FINAL_READ_ATTEMPTS = 10
PLACE_FINAL_READ_SLEEP_SEC = 1.0
# ...and no new read starts once this long has passed since the first, so a
# Novig whose reads each hang to their 10 s timeout cannot hold the click past
# the phone's 90 s wait (runner quote 30 s + send 10 s + this + fills 10 s).
PLACE_FINAL_READ_WINDOW_SEC = 12.0
# The bets service's Novig reads the stake counted must be at most this old
# (the Novig bot's holdings gate).
NOVIG_READ_MAX_AGE = timedelta(seconds=MAX_HOLDINGS_AGE_SEC)
# Bet types an outcome id exists for: Unabated carries Novig ids on spread and
# total rungs only, never on a moneyline.
NOVIG_PLACE_BET_TYPES = (2, 3)
# Novig's sources whose reads the Ticket's stake counts: Cal's app account and
# Place's own sub-account (the bot's is added when it is read, service.py).
NOVIG_PLACE_READ_SOURCES = ("novig", "novig_place")


def validate_novig_place_request(body: object) -> dict | str:
    """A POST /place_novig.json body, or what was expected and what was found.
    {clickId, eventId, periodTypeId, betTypeId (2 spread | 3 total), sideIndex,
    points, stake (dollars), ticketAllIn (the Ticket's Novig price as a
    probability), novigReadAt (the Novig read the stake was sized on, ISO),
    kalshiReadAt (the Kalshi read, ISO, or null when Kalshi was never read)}"""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    click_id = body.get("clickId")
    if not isinstance(click_id, str) or not UUID_PATTERN.match(click_id):
        return f"clickId must be a UUID (Novig's clientId), got {click_id!r}"
    whole = {}
    for field, allowed in (("eventId", None), ("periodTypeId", None), ("betTypeId", NOVIG_PLACE_BET_TYPES),
                           ("sideIndex", (0, 1))):
        value = body.get(field)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0 or (allowed and value not in allowed):
            return f"{field} must be {'one of ' + str(allowed) if allowed else 'a whole number'}, got {value!r}"
        whole[field] = value
    points = body.get("points")
    if not _is_number(points) or not math.isfinite(points):
        return f"points must be a number (a spread or total), got {points!r}"
    stake = body.get("stake")
    if not _is_number(stake) or not math.isfinite(stake) or stake <= 0:
        return f"stake must be dollars above 0, got {stake!r}"
    ticket_all_in = body.get("ticketAllIn")
    if not _is_number(ticket_all_in) or not 0 < ticket_all_in < 1:
        return f"ticketAllIn must be a probability between 0 and 1, got {ticket_all_in!r}"
    novig_read_at = _iso_time(body.get("novigReadAt"))
    if novig_read_at is None:
        return f"novigReadAt must be the ISO time of the Novig read the stake was sized on, got {body.get('novigReadAt')!r}"
    kalshi_read_at = optional_read_at(body, "kalshiReadAt")
    if isinstance(kalshi_read_at, str):
        return kalshi_read_at
    return {"click_id": click_id, **whole, "points": float(points), "stake": float(stake),
            "ticket_all_in": float(ticket_all_in), "novig_read_at": novig_read_at, "kalshi_read_at": kalshi_read_at}


class ManualNovigDesk:
    """Sends Cal's Novig Place-button orders from Place's sub-account, one per click."""

    def __init__(self, client: NovigOrderClient, store: BetsStore, quote_reader: Callable[[dict], dict],
                 read_sources: tuple[str, ...] = NOVIG_PLACE_READ_SOURCES,
                 now: Callable[[], datetime] = lambda: datetime.now(timezone.utc),
                 sleep: Callable[[float], None] = time.sleep, clock: Callable[[], float] = time.monotonic):
        """`client` signs with Place's own trading key (never the bot's).
        `read_sources`: the bets service's source names whose reads the
        Ticket's stake counts; each must be ok and fresh."""
        self._client = client
        self._store = store
        self._quote_reader = quote_reader
        self._read_sources = read_sources
        self._now = now
        self._sleep = sleep
        self._clock = clock
        # One click at a time (panel and phone can click at once): the checks
        # read the orders already sent, so they must see each other's.
        self._lock = threading.Lock()

    def place(self, request: dict) -> dict:
        """One order. Returns {status, message, ...}: refused (nothing sent), or
        the sent order's filled | partial | unfilled | rejected | unknown."""
        with self._lock:
            return self._place(request)

    def _place(self, request: dict) -> dict:
        client_order_id = f"{NOVIG_ORDER_ID_PREFIX}{request['click_id']}"
        quote = self._quote_reader({"eventId": request["eventId"], "periodTypeId": request["periodTypeId"],
                                    "betTypeId": request["betTypeId"], "sideIndex": request["sideIndex"],
                                    "points": request["points"]})
        refusal = self._refusal(request, quote)
        if refusal:
            return {"status": STATUS_REFUSED, "message": f"Not placed: {refusal}"}
        context = line_context_of(quote.get("line"))
        if isinstance(context, str):
            return {"status": STATUS_REFUSED, "message": f"Not placed: {context}."}
        quantity = qty_for_stake(request["stake"], quote["price"])
        if quantity < 1:
            return {"status": STATUS_REFUSED,
                    "message": f"Not placed: ${request['stake']:.2f} buys no contract at {_cents(quote['price'])}."}
        inserted = self._store.insert_manual_order({
            "client_order_id": client_order_id, "placed_at": self._now(), "venue": VENUE_NOVIG, "event_id": request["eventId"],
            "ticker": quote["outcomeId"], "contract_side": "buy", "label": quote.get("label"),
            "ticket_all_in": request["ticket_all_in"], "ask_price": quote["price"], "all_in_price": quote["price"],
            "stake": request["stake"], "count_requested": quantity, "status": "sending", "updated_at": self._now(),
            "market_id": quote["marketId"], "event_start": _iso_time(quote["eventStart"]), **context,
        })
        if not inserted:
            return {"status": STATUS_REFUSED, "message": "Not placed: this click was already sent."}
        log.info("novig place %s: outcome %s x%d @ %s (Ticket %s)", client_order_id, quote["outcomeId"], quantity,
                 quote["price"], request["ticket_all_in"])
        return self._send(client_order_id, request["click_id"], quote, quantity)

    def _refusal(self, request: dict, quote: dict) -> str | None:
        """Why this click must not become an order, in words, or None."""
        if not quote.get("ok"):
            return f"{quote.get('detail') or quote.get('reason') or 'no Novig outcome for this line'}."
        start = _iso_time(quote.get("eventStart"))
        if start is None:
            return "the game's start time is unknown."
        if start <= self._now():
            return "the game has started."
        if quote.get("novigBaseUrl") != self._client.base_url:
            return (f"the runner read Novig at {quote.get('novigBaseUrl')!r} but this desk trades "
                    f"{self._client.base_url}.")
        price = quote.get("price")
        if not _is_number(price) or not is_on_price_grid(price):
            return f"Novig's price {price!r} is not a price on its grid."
        if not all(isinstance(quote.get(field), str) and quote[field]
                   for field in ("outcomeId", "otherOutcomeId", "marketId")):
            return "the runner did not name the market and both its outcomes."
        if price > request["ticket_all_in"] + PRICE_TOLERANCE:
            return (f"Novig is at {_cents(price)} now, more than 0.1¢ worse than the Ticket's "
                    f"{_cents(request['ticket_all_in'])}.")
        if self._store.novig_outcome_is_held(quote["otherOutcomeId"]):
            return "you hold the other side of this Novig market; the stake cannot size a hedge."
        stale = self._stale_read()
        if stale:
            return stale
        if unseen_manual_fill(self._store, request["eventId"], request["kalshi_read_at"], request["novig_read_at"]):
            return ("your last Place on this game is not on the Bets tab yet, so this stake does not count it. "
                    "Wait for the next read (up to a minute) and place again.")
        return None

    def _stale_read(self) -> str | None:
        """Words for a Novig read the stake counted that is failing or old, or None."""
        status = self._store.source_status()
        for source in self._read_sources:
            read = status.get(source)
            fetched_at = _iso_time(read.get("fetchedAt")) if read else None
            if not read or not read.get("ok") or fetched_at is None:
                return f"the {source} read is failing, so the stake may miss a bet. Try again once it reads."
            if self._now() - fetched_at > NOVIG_READ_MAX_AGE:
                return f"the {source} read is over {NOVIG_READ_MAX_AGE.seconds // 60} minutes old, so the stake may miss a bet."
        return None

    def _send(self, client_order_id: str, novig_client_id: str, quote: dict, quantity: int) -> dict:
        """ONE IOC buy, never re-sent; then the order read until final and its fills."""
        try:
            status, reply = self._client.place_ioc_order(quote["outcomeId"], quote["price"], quantity, novig_client_id)
        except (OSError, NovigError) as error:
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": f"no answer: {error}"})
            log.error("novig place %s: no answer from Novig: %s", client_order_id, error)
            return {"status": STATUS_UNKNOWN, "message": "Novig did not answer. Check Novig before placing again."}
        order_id = reply.get("orderId") if isinstance(reply, dict) else None
        if status in DEFINITE_REJECTION_STATUSES:
            code = error_code_of(reply)
            self._update(client_order_id, {"status": STATUS_REJECTED, "error": f"HTTP {status} {code or ''}: {str(reply)[:400]}"})
            log.warning("novig place %s refused: HTTP %s %s", client_order_id, status, reply)
            if status == GEO_STATUS:
                return {"status": STATUS_REJECTED, "message": "Novig refused it: your phone's Novig location check is over "
                                                              "3 days old. Open the Novig app on your phone, then place again."}
            return {"status": STATUS_REJECTED,
                    "message": f"Novig refused it (HTTP {status}{' ' + code if code else ''}). Nothing was placed."}
        if status != 201 or not isinstance(order_id, str) or not order_id:
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": f"no definite answer: HTTP {status}: {str(reply)[:300]}"})
            log.error("novig place %s: no definite answer: HTTP %s %s", client_order_id, status, reply)
            return {"status": STATUS_UNKNOWN, "message": "Novig did not say whether it took the order. Check Novig before placing again."}
        self._update(client_order_id, {"novig_order_id": order_id})
        return self._follow_order(client_order_id, order_id, quote, quantity)

    def _follow_order(self, client_order_id: str, order_id: str, quote: dict, quantity: int) -> dict:
        """Read the order until Novig shows it final (a few tries), then its
        fills. An OPEN order (an IOC never should be) is cancelled."""
        order_status = None
        deadline = self._clock() + PLACE_FINAL_READ_WINDOW_SEC
        for attempt in range(PLACE_FINAL_READ_ATTEMPTS):
            if attempt:
                if self._clock() >= deadline:
                    break
                self._sleep(PLACE_FINAL_READ_SLEEP_SEC)
            try:
                status, reply = self._client.get_order(order_id)
            except OSError as error:
                log.warning("novig place %s: order read failed: %s", client_order_id, error)
                continue
            order_status = reply.get("status") if status == 200 and isinstance(reply, dict) else None
            if order_status in FINAL_ORDER_STATUSES:
                break
            if order_status == STATUS_OPEN:
                return self._cancel_resting(client_order_id, order_id)
        if order_status not in FINAL_ORDER_STATUSES:
            return self._unconfirmed(client_order_id, order_id, quote, quantity)
        return self._read_fills(client_order_id, order_id, quote, quantity)

    def _unconfirmed(self, client_order_id: str, order_id: str, quote: dict, quantity: int) -> dict:
        """The order never read final: its fills decide what can be said. All
        contracts filled is final whatever the order read says; anything less
        may still grow, so the row stays 'unknown' with what has filled."""
        log.warning("novig place %s: order %s never read final; reading its fills", client_order_id, order_id)
        try:
            filled, cost, fee = order_totals(self._client.fills(order_id=order_id))
        except (NovigError, OSError) as error:
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": f"never read final; fills read failed: {error}"})
            log.error("novig place %s: order %s never read final and its fills could not be read: %s",
                      client_order_id, order_id, error)
            return {"status": STATUS_UNKNOWN, "message": "Novig took the order but has not said what filled. "
                                                         "Check Novig before placing again."}
        if filled >= quantity:
            return self._record_fills(client_order_id, quote, quantity, filled, cost, fee)
        paid = round(cost + fee, 5)
        self._update(client_order_id, {"status": STATUS_UNKNOWN, "count_filled": filled, "cost_filled": paid,
                                       "error": f"never read final; {filled} of {quantity} filled so far"})
        log.error("novig place %s: order %s never read final, %d of %d filled so far", client_order_id, order_id,
                  filled, quantity)
        if filled == 0:
            return {"status": STATUS_UNKNOWN, "message": "Novig took the order but has not said what filled. "
                                                         "Check Novig before placing again."}
        return {"status": STATUS_UNKNOWN, "outcomeId": quote["outcomeId"], "count": quantity, "filled": filled, "cost": paid,
                "message": f"Bet ${paid:,.2f} at {_cents(quote['price'])} on Novig so far; Novig has not said the order "
                           "is done. Check Novig before placing again."}

    def _cancel_resting(self, client_order_id: str, order_id: str) -> dict:
        log.error("novig place %s: order %s is OPEN on Novig after an IOC send; cancelling", client_order_id, order_id)
        try:
            status, reply = self._client.cancel_order(order_id)
        except OSError as error:
            status, reply = None, f"no answer: {error}"
        self._update(client_order_id, {"status": STATUS_UNKNOWN,
                                       "error": f"OPEN after IOC; cancel answered {status} {str(reply)[:200]}"})
        return {"status": STATUS_UNKNOWN, "message": "Novig left the order resting, so a cancel was sent. "
                                                     "Check Novig before placing again."}

    def _read_fills(self, client_order_id: str, order_id: str, quote: dict, quantity: int) -> dict:
        try:
            filled, cost, fee = order_totals(self._client.fills(order_id=order_id))
        except (NovigError, OSError) as error:
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": f"fills read failed: {error}"})
            log.error("novig place %s: fills read failed: %s", client_order_id, error)
            return {"status": STATUS_UNKNOWN, "message": "Novig took the order but its fills could not be read. "
                                                         "Check Novig before placing again."}
        return self._record_fills(client_order_id, quote, quantity, filled, cost, fee)

    def _record_fills(self, client_order_id: str, quote: dict, quantity: int, filled: int, cost: float,
                      fee: float) -> dict:
        paid = round(cost + fee, 5)
        outcome = STATUS_UNFILLED if filled == 0 else STATUS_FILLED if filled >= quantity else STATUS_PARTIAL
        self._update(client_order_id, {"status": outcome, "count_filled": filled, "cost_filled": paid, "error": None})
        log.info("novig place %s: %s, %d of %d for $%.2f", client_order_id, outcome, filled, quantity, paid)
        return {"status": outcome, "outcomeId": quote["outcomeId"], "count": quantity, "filled": filled, "cost": paid,
                "message": _outcome_message(outcome, quote, quantity, filled, paid)}

    def _update(self, client_order_id: str, fields: dict) -> None:
        self._store.update_manual_order(client_order_id, fields, self._now())


def _outcome_message(outcome: str, quote: dict, quantity: int, filled: int, paid: float) -> str:
    price = _cents(quote["price"])
    to_win = filled / QTY_PER_DOLLAR - paid
    if outcome == STATUS_FILLED:
        return f"Bet ${paid:,.2f} at {price} on Novig to win ${to_win:,.2f}."
    if outcome == STATUS_PARTIAL:
        asked = quantity * quote["price"] / QTY_PER_DOLLAR
        return (f"Bet ${paid:,.2f} of ${asked:,.2f} at {price} on Novig to win ${to_win:,.2f}; "
                "the rest was not there. Nothing rests.")
    return f"Nothing filled at {price}: the price was gone. Nothing rests."
