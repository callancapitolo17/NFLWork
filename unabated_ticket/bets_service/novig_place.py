"""The Ticket's Place button on Novig (novig_place_plan.md, Cal 2026-10-16):
Cal clicks, and this desk sends ONE order from Place's own Novig sub-account
(never the bot's, never his app balance: no key trades that). Nothing here runs
on a timer: every order is one POST /place_novig.json from the panel or the
phone page.

Rules, in order, for every request (validate_novig_place_request's shape):
  1. the click was not already sent: client id 'nman-' + its clickId, held in
     bets.duckdb::manual_orders before the send (Novig does not dedupe ids,
     so this row is the only guard);
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
     `kalshiReadAt` (the stake counts bets at both);
  8. contracts = floor(stake / (price x 1c)): a Novig contract pays 1c, not $1;
  9. one immediate-or-cancel limit buy at the live price through the order
     client, which looks a lost send up and never re-sends it.

Side effects: a Novig order from Place's sub-account (money); INSERT / UPDATE
bets.duckdb::manual_orders (venue 'novig'). Network: GET
<runner>/novig_quote.json, then the order client's Novig calls.
"""
import logging
import math
import threading
from collections.abc import Callable
from datetime import datetime, timedelta, timezone
from typing import Protocol

from unabated_ticket.bets_service.manual_orders import (
    CLICK_ID_PATTERN, PRICE_TOLERANCE, STATUS_FILLED, STATUS_PARTIAL, STATUS_REFUSED, STATUS_REJECTED,
    STATUS_UNFILLED, STATUS_UNKNOWN, _cents, _is_number, _iso_time,
)
from unabated_ticket.bets_service.store import BetsStore

log = logging.getLogger(__name__)

NOVIG_ORDER_ID_PREFIX = "nman-"
RUNNER_NOVIG_QUOTE_PATH = "/novig_quote.json"
# A Novig contract pays 1c (qty 100 = $1 payout), unlike Kalshi's $1.
NOVIG_CONTRACT_PAYOUT_DOLLARS = 0.01
# Novig's price grid runs 0.001 to 0.999.
MIN_NOVIG_PRICE = 0.001
MAX_NOVIG_PRICE = 0.999
# The bets service's Novig reads the stake counted must be at most this old
# (the Kalshi bot's poll gate).
NOVIG_READ_MAX_AGE = timedelta(minutes=5)
# Bet types an outcome id exists for: Unabated carries Novig ids on spread and
# total rungs only, never on a moneyline.
NOVIG_PLACE_BET_TYPES = (2, 3)
STATUS_GEO_BLOCKED = "geo_blocked"  # Novig's 451: nothing placed
VENUE_NOVIG = "novig"
VENUE_KALSHI = "kalshi"


class NovigPlaceClient(Protocol):
    """The shared Novig order client (bets_service/novig_orders.py, the
    Novig taker's), configured with Place's key id, key file and sub-account.

    place_ioc_buy sends ONE immediate-or-cancel limit buy and answers
    {outcome, orderId, filledQuantity, cost, detail}: outcome one of filled |
    partial | unfilled | rejected | geo_blocked | unknown, cost in dollars.
    A send with no answer is looked up by clientId on that outcome, never
    sent again; when the lookup is not definite either, outcome is unknown."""

    base_url: str

    def place_ioc_buy(self, outcome_id: str, price: float, quantity: int, client_id: str) -> dict: ...


def validate_novig_place_request(body: object) -> dict | str:
    """A POST /place_novig.json body, or what was expected and what was found.
    {clickId, eventId, periodTypeId, betTypeId (2 spread | 3 total), sideIndex,
    points, stake (dollars), ticketAllIn (the Ticket's Novig price as a
    probability), novigReadAt (the Novig read the stake was sized on, ISO),
    kalshiReadAt (the Kalshi read, ISO, or null when Kalshi was never read)}"""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    click_id = body.get("clickId")
    if not isinstance(click_id, str) or not CLICK_ID_PATTERN.match(click_id):
        return f"clickId must be 8-48 letters, digits or dashes, got {click_id!r}"
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
    if not _is_number(ticket_all_in) or not MIN_NOVIG_PRICE <= ticket_all_in <= MAX_NOVIG_PRICE:
        return f"ticketAllIn must be a probability from {MIN_NOVIG_PRICE} to {MAX_NOVIG_PRICE}, got {ticket_all_in!r}"
    novig_read_at = _iso_time(body.get("novigReadAt"))
    if novig_read_at is None:
        return f"novigReadAt must be the ISO time of the Novig read the stake was sized on, got {body.get('novigReadAt')!r}"
    kalshi_read_at = None
    if body.get("kalshiReadAt") is not None:
        kalshi_read_at = _iso_time(body.get("kalshiReadAt"))
        if kalshi_read_at is None:
            return f"kalshiReadAt must be an ISO time or null, got {body.get('kalshiReadAt')!r}"
    return {"click_id": click_id, **whole, "points": float(points), "stake": float(stake),
            "ticket_all_in": float(ticket_all_in), "novig_read_at": novig_read_at, "kalshi_read_at": kalshi_read_at}


def novig_contracts_for(stake: float, price: float) -> int:
    """Whole Novig contracts `stake` dollars buys at `price` (each pays 1c and
    costs price x 1c). Rounded to 6 places first so 50 / 0.5 is 10000, not 9999."""
    return math.floor(round(stake / (price * NOVIG_CONTRACT_PAYOUT_DOLLARS), 6))


class ManualNovigDesk:
    """Sends Cal's Novig Place-button orders from Place's sub-account, one per click."""

    def __init__(self, client: NovigPlaceClient, store: BetsStore, quote_reader: Callable[[dict], dict],
                 read_sources: tuple[str, ...] = (VENUE_NOVIG,),
                 now: Callable[[], datetime] = lambda: datetime.now(timezone.utc)):
        """`read_sources`: the bets service's source names whose reads the
        Ticket's stake counts (the app Portfolio read, plus Place's own
        sub-account read once that source exists); each must be ok and fresh."""
        self._client = client
        self._store = store
        self._quote_reader = quote_reader
        self._read_sources = read_sources
        self._now = now
        # One click at a time (panel and phone can click at once): the checks
        # read the orders already sent, so they must see each other's.
        self._lock = threading.Lock()

    def place(self, request: dict) -> dict:
        """One order. Returns {status, message, ...}: refused (nothing sent), or
        the sent order's filled | partial | unfilled | rejected | geo_blocked | unknown."""
        with self._lock:
            return self._place(request)

    def _place(self, request: dict) -> dict:
        client_id = f"{NOVIG_ORDER_ID_PREFIX}{request['click_id']}"
        quote = self._quote_reader({"eventId": request["eventId"], "periodTypeId": request["periodTypeId"],
                                    "betTypeId": request["betTypeId"], "sideIndex": request["sideIndex"],
                                    "points": request["points"]})
        refusal = self._refusal(request, quote)
        if refusal:
            return {"status": STATUS_REFUSED, "message": f"Not placed: {refusal}"}
        quantity = novig_contracts_for(request["stake"], quote["price"])
        if quantity < 1:
            return {"status": STATUS_REFUSED,
                    "message": f"Not placed: ${request['stake']:.2f} buys no contract at {_cents(quote['price'])}."}
        inserted = self._store.insert_manual_order({
            "client_order_id": client_id, "placed_at": self._now(), "venue": VENUE_NOVIG, "event_id": request["eventId"],
            "ticker": quote["outcomeId"], "contract_side": "buy", "label": quote.get("label"),
            "ticket_all_in": request["ticket_all_in"], "ask_price": quote["price"], "all_in_price": quote["price"],
            "stake": request["stake"], "count_requested": quantity, "status": "sending", "updated_at": self._now(),
        })
        if not inserted:
            return {"status": STATUS_REFUSED, "message": "Not placed: this click was already sent."}
        log.info("novig place %s: outcome %s x%d @ %s (Ticket %s)", client_id, quote["outcomeId"], quantity,
                 quote["price"], request["ticket_all_in"])
        return self._send(client_id, quote, quantity)

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
        if not _is_number(price) or not MIN_NOVIG_PRICE <= price <= MAX_NOVIG_PRICE:
            return f"Novig's price {price!r} is not a price."
        if not isinstance(quote.get("outcomeId"), str) or not isinstance(quote.get("otherOutcomeId"), str):
            return "the runner did not name both outcomes of this market."
        if price > request["ticket_all_in"] + PRICE_TOLERANCE:
            return (f"Novig is at {_cents(price)} now, more than 0.1¢ worse than the Ticket's "
                    f"{_cents(request['ticket_all_in'])}.")
        if self._store.novig_outcome_is_held(quote["otherOutcomeId"]):
            return "you hold the other side of this Novig market; the stake cannot size a hedge."
        stale = self._stale_read()
        if stale:
            return stale
        if self._maybe_unseen_fill(request):
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

    def _maybe_unseen_fill(self, request: dict) -> bool:
        """An earlier Place on this game that may have filled after the read
        of its venue the stake was sized on (no Kalshi read: any Kalshi Place)."""
        novig_orders = self._store.manual_orders_maybe_filled_since(request["eventId"], request["novig_read_at"],
                                                                    venue=VENUE_NOVIG)
        kalshi_since = request["kalshi_read_at"] or datetime(1970, 1, 1, tzinfo=timezone.utc)
        kalshi_orders = self._store.manual_orders_maybe_filled_since(request["eventId"], kalshi_since, venue=VENUE_KALSHI)
        return bool(novig_orders or kalshi_orders)

    def _send(self, client_id: str, quote: dict, quantity: int) -> dict:
        try:
            result = self._client.place_ioc_buy(quote["outcomeId"], quote["price"], quantity, client_id)
        except OSError as error:
            self._update(client_id, {"status": STATUS_UNKNOWN, "error": f"no answer: {error}"})
            log.error("novig place %s: no answer from Novig: %s", client_id, error)
            return {"status": STATUS_UNKNOWN, "message": "Novig did not answer. Check Novig before placing again."}
        outcome = result.get("outcome") if isinstance(result, dict) else None
        detail = str(result.get("detail") or "")[:500] if isinstance(result, dict) else f"not a result: {result!r}"[:500]
        if outcome == STATUS_GEO_BLOCKED:
            self._update(client_id, {"status": STATUS_REJECTED, "error": f"451: {detail}"})
            return {"status": STATUS_REJECTED, "message": "Novig refused it: your phone's Novig location check is over "
                                                          "3 days old. Open the Novig app on your phone, then place again."}
        if outcome == STATUS_REJECTED:
            self._update(client_id, {"status": STATUS_REJECTED, "error": detail})
            return {"status": STATUS_REJECTED,
                    "message": f"Novig refused it{': ' + detail if detail else ''}. Nothing was placed."}
        filled = result.get("filledQuantity")
        cost = result.get("cost")
        definite = outcome in (STATUS_FILLED, STATUS_PARTIAL, STATUS_UNFILLED) and _is_number(filled) and filled >= 0
        if not definite:
            self._update(client_id, {"status": STATUS_UNKNOWN, "novig_order_id": result.get("orderId"),
                                     "error": detail or f"no definite answer: {str(result)[:300]}"})
            log.error("novig place %s: no definite answer: %s", client_id, result)
            return {"status": STATUS_UNKNOWN, "message": "Novig did not say what filled. Check Novig before placing again."}
        self._update(client_id, {"status": outcome, "novig_order_id": result.get("orderId"), "count_filled": filled,
                                 "cost_filled": cost if _is_number(cost) else None, "error": None})
        return {"status": outcome, "outcomeId": quote["outcomeId"], "count": quantity, "filled": filled,
                "cost": cost if _is_number(cost) else None,
                "message": _outcome_message(outcome, quote, quantity, filled, cost)}

    def _update(self, client_id: str, fields: dict) -> None:
        self._store.update_manual_order(client_id, fields, self._now())


def _outcome_message(outcome: str, quote: dict, quantity: int, filled: float, cost: object) -> str:
    price = _cents(quote["price"])
    paid = cost if _is_number(cost) else filled * quote["price"] * NOVIG_CONTRACT_PAYOUT_DOLLARS
    to_win = filled * NOVIG_CONTRACT_PAYOUT_DOLLARS - paid
    if outcome == STATUS_FILLED:
        return f"Bet ${paid:,.2f} at {price} on Novig to win ${to_win:,.2f}."
    if outcome == STATUS_PARTIAL:
        asked = quantity * quote["price"] * NOVIG_CONTRACT_PAYOUT_DOLLARS
        return (f"Bet ${paid:,.2f} of ${asked:,.2f} at {price} on Novig to win ${to_win:,.2f}; "
                "the rest was not there. Nothing rests.")
    return f"Nothing filled at {price}: the price was gone. Nothing rests."
