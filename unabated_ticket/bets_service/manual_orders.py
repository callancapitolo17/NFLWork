"""The Ticket's Place button on Kalshi (manual_place_plan.md): Cal clicks, and
this desk sends ONE order on his own Kalshi key (KALSHI_API_KEY_ID, the key
the Kalshi source reads with). Nothing here runs on a timer: every order is
one POST /place_kalshi.json from the panel or the phone page.

Rules, in order, for every request (validate_place_request's shape):
  1. the click was not already sent: client_order_id 'uman-' + its clickId,
     held in bets.duckdb::manual_orders before the send;
  2. the server runner finds the board's Kalshi line, its contract and the
     live ask, and checks Kalshi's market is that line's (GET
     /kalshi_quote.json: the automated taker's own resolveKalshiContract,
     quoteOf and marketMismatch);
  3. the game has not started, and the runner read the Kalshi this desk
     trades on (KALSHI_BASE_URL);
  4. the live all-in price is at most PRICE_TOLERANCE (0.1c) worse than the
     Ticket's (Cal, 2026-10-16); a better one is taken;
  5. the other side of the market is not held: on Kalshi a buy of the other
     side sells it;
  6. no earlier Place on this game may have filled after the Kalshi read the
     Ticket's stake was sized on (`kalshiReadAt`): that stake would not know
     about it (the bot's unseen-fill rule, on time instead of ids);
  7. contracts = floor(stake / live all-in), at least one;
  8. one immediate-or-cancel limit buy at the live ask
     (taker_orders.order_body: nothing rests), then the order's fills read.

Side effects: a Kalshi order on Cal's own key, his primary account (no
subaccount named; the bot trades its own numbered one) (money); INSERT / UPDATE
bets.duckdb::manual_orders. Network: GET <runner>/kalshi_quote.json, then
Kalshi REST signed with Cal's key (the order, and its fills).
"""
import json
import logging
import math
import re
import threading
import urllib.error
import urllib.request
from collections.abc import Callable
from datetime import datetime, timezone
from urllib.parse import urlencode

from unabated_ticket.bets_service.sources.kalshi import PRIMARY_SUBACCOUNT
from unabated_ticket.bets_service.store import BetsStore
from unabated_ticket.bets_service.taker_orders import (
    DEFINITE_REJECTION_STATUSES, ORDERS_PATH, KalshiOrderClient, fill_totals, held_id_of, order_body,
)

log = logging.getLogger(__name__)

MANUAL_ORDER_ID_PREFIX = "uman-"
# The click's own id, made by the page: letters, digits and dashes.
CLICK_ID_PATTERN = re.compile(r"^[A-Za-z0-9-]{8,48}$")
# Unabated's Kalshi price sometimes lags the book by a tick: 0.1c worse still
# places (Cal, 2026-10-16, the bot's own tolerance); a better live price is taken.
PRICE_TOLERANCE = 0.001
MIN_PRICE = 0.01
MAX_PRICE = 0.99
RUNNER_QUOTE_PATH = "/kalshi_quote.json"
RUNNER_TIMEOUT_SEC = 10
FILLS_PATH = "/portfolio/fills"
FILLS_PAGE_LIMIT = 200
# The runner is local (or on the tailnet), never behind an HTTP proxy.
_RUNNER_OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))

STATUS_REFUSED = "refused"      # nothing sent
STATUS_FILLED = "filled"
STATUS_PARTIAL = "partial"
STATUS_UNFILLED = "unfilled"    # sent, nothing filled, nothing rests
STATUS_REJECTED = "rejected"    # Kalshi refused it: nothing placed
STATUS_UNKNOWN = "unknown"      # no definite answer: it may have filled


def validate_place_request(body: object) -> dict | str:
    """A POST /place_kalshi.json body, or what was expected and what was found.
    {clickId, eventId, periodTypeId, betTypeId, sideIndex, points (null on a
    moneyline), stake (dollars), ticketAllIn (the Ticket's price, all-in, as a
    probability), kalshiReadAt (the Kalshi read the stake was sized on, ISO)}"""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    click_id = body.get("clickId")
    if not isinstance(click_id, str) or not CLICK_ID_PATTERN.match(click_id):
        return f"clickId must be 8-48 letters, digits or dashes, got {click_id!r}"
    whole = {}
    for field, allowed in (("eventId", None), ("periodTypeId", None), ("betTypeId", (1, 2, 3)), ("sideIndex", (0, 1))):
        value = body.get(field)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0 or (allowed and value not in allowed):
            return f"{field} must be {'one of ' + str(allowed) if allowed else 'a whole number'}, got {value!r}"
        whole[field] = value
    points = body.get("points")
    if points is not None and (not _is_number(points) or not math.isfinite(points)):
        return f"points must be a number or null, got {points!r}"
    stake = body.get("stake")
    if not _is_number(stake) or not math.isfinite(stake) or stake <= 0:
        return f"stake must be dollars above 0, got {stake!r}"
    ticket_all_in = body.get("ticketAllIn")
    if not _is_number(ticket_all_in) or not MIN_PRICE <= ticket_all_in < 1:
        return f"ticketAllIn must be a probability from {MIN_PRICE} to under 1, got {ticket_all_in!r}"
    kalshi_read_at = _iso_time(body.get("kalshiReadAt"))
    if kalshi_read_at is None:
        return f"kalshiReadAt must be the ISO time of the Kalshi read the stake was sized on, got {body.get('kalshiReadAt')!r}"
    return {"click_id": click_id, **whole, "points": None if points is None else float(points),
            "stake": float(stake), "ticket_all_in": float(ticket_all_in), "kalshi_read_at": kalshi_read_at}


def _is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _iso_time(value: object) -> datetime | None:
    if not isinstance(value, str):
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return parsed if parsed.tzinfo is not None else None


def read_runner_quote(runner_url: str, request: dict, quote_path: str = RUNNER_QUOTE_PATH) -> dict:
    """GET <runner><quote_path> (/kalshi_quote.json, or Novig Place's
    /novig_quote.json) for the request's line: the runner's answer, or
    {ok: False, reason: 'runner', detail} when it could not be read."""
    query = {key: request[key] for key in ("eventId", "periodTypeId", "betTypeId", "sideIndex")}
    if request["points"] is not None:
        query["points"] = f"{request['points']:g}"
    url = f"{runner_url.rstrip('/')}{quote_path}?{urlencode(query)}"
    try:
        with _RUNNER_OPENER.open(url, timeout=RUNNER_TIMEOUT_SEC) as response:
            answer = json.loads(response.read())
    except urllib.error.HTTPError as error:
        return {"ok": False, "reason": "runner", "detail": f"the server runner answered HTTP {error.code}: {error.read()[:300]!r}"}
    except (urllib.error.URLError, OSError, ValueError) as error:
        return {"ok": False, "reason": "runner", "detail": f"the server runner at {runner_url} could not be read: {error}"}
    if not isinstance(answer, dict):
        return {"ok": False, "reason": "runner", "detail": f"the server runner answered {type(answer).__name__}, not an object"}
    return answer


def _cents(probability: float) -> str:
    return f"{probability * 100:.1f}¢"


class ManualKalshiDesk:
    """Sends Cal's Place-button orders on his own Kalshi key, one per click."""

    def __init__(self, client: KalshiOrderClient, store: BetsStore, quote_reader: Callable[[dict], dict],
                 now: Callable[[], datetime] = lambda: datetime.now(timezone.utc)):
        self._client = client
        self._store = store
        self._quote_reader = quote_reader
        self._now = now
        # One click at a time (panel and phone can click at once): the checks
        # read the orders already sent, so they must see each other's.
        self._lock = threading.Lock()

    def place(self, request: dict) -> dict:
        """One order. Returns {status, message, ...}: refused (nothing sent),
        or the sent order's filled | partial | unfilled | rejected | unknown."""
        with self._lock:
            return self._place(request)

    def _place(self, request: dict) -> dict:
        client_order_id = f"{MANUAL_ORDER_ID_PREFIX}{request['click_id']}"
        quote = self._quote_reader({"eventId": request["eventId"], "periodTypeId": request["periodTypeId"],
                                    "betTypeId": request["betTypeId"], "sideIndex": request["sideIndex"],
                                    "points": request["points"]})
        refusal = self._refusal(request, quote)
        if refusal:
            return {"status": STATUS_REFUSED, "message": f"Not placed: {refusal}"}
        count = math.floor(request["stake"] / quote["allIn"])
        if count < 1:
            return {"status": STATUS_REFUSED,
                    "message": f"Not placed: ${request['stake']:.2f} buys no contract at {_cents(quote['allIn'])}."}
        inserted = self._store.insert_manual_order({
            "client_order_id": client_order_id, "placed_at": self._now(), "venue": "kalshi", "event_id": request["eventId"],
            "ticker": quote["ticker"], "contract_side": quote["side"], "label": quote.get("label"),
            "ticket_all_in": request["ticket_all_in"], "ask_price": quote["ask"], "all_in_price": quote["allIn"],
            "stake": request["stake"], "count_requested": count, "status": "sending", "updated_at": self._now(),
        })
        if not inserted:
            return {"status": STATUS_REFUSED, "message": "Not placed: this click was already sent."}
        log.info("manual order %s: %s %s x%d @ %s (Ticket %s)", client_order_id, quote["ticker"], quote["side"],
                 count, quote["ask"], request["ticket_all_in"])
        return self._send(client_order_id, quote, count)

    def _refusal(self, request: dict, quote: dict) -> str | None:
        """Why this click must not become an order, in words, or None."""
        if not quote.get("ok"):
            return f"{quote.get('detail') or quote.get('reason') or 'no Kalshi contract for this line'}."
        start = _iso_time(quote.get("eventStart"))
        if start is None:
            return "the game's start time is unknown."
        if start <= self._now():
            return "the game has started."
        if quote.get("kalshiBaseUrl") != self._client.base_url:
            return (f"the runner read Kalshi at {quote.get('kalshiBaseUrl')!r} but this desk trades "
                    f"{self._client.base_url}.")
        ask = quote.get("ask")
        if not _is_number(ask) or not MIN_PRICE <= ask <= MAX_PRICE or not _is_number(quote.get("allIn")):
            return f"Kalshi's ask {quote.get('ask')!r} is not a price."
        if quote["allIn"] > request["ticket_all_in"] + PRICE_TOLERANCE:
            return (f"Kalshi is at {_cents(quote['allIn'])} all-in now, more than 0.1¢ worse than the "
                    f"Ticket's {_cents(request['ticket_all_in'])}.")
        other_side = "no" if quote["side"] == "yes" else "yes"
        if self._store.bet_is_open(held_id_of(quote["ticker"], other_side, PRIMARY_SUBACCOUNT)):
            return f"you hold the {other_side.upper()} side of {quote['ticker']}; buying {quote['side'].upper()} would sell it."
        recent = self._store.manual_orders_maybe_filled_since(request["eventId"], request["kalshi_read_at"])
        if recent:
            return ("your last order on this game is not on the Bets tab yet, so this stake does not count it. "
                    "Wait for the next Kalshi read (up to a minute) and place again.")
        return None

    def _send(self, client_order_id: str, quote: dict, count: int) -> dict:
        body = order_body(quote["ticker"], quote["side"], quote["ask"], count, client_order_id, PRIMARY_SUBACCOUNT)
        # Cal's MAIN account (his choice A): an order naming no subaccount is the
        # primary's, as every order was before the bot moved to its own.
        del body["subaccount"]
        try:
            status, reply = self._client.request("POST", ORDERS_PATH, body)
        except OSError as error:
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": f"no answer: {error}"})
            log.error("manual order %s: no answer from Kalshi: %s", client_order_id, error)
            return {"status": STATUS_UNKNOWN, "message": "Kalshi did not answer. Check Kalshi before placing again."}
        if status in DEFINITE_REJECTION_STATUSES:
            detail = _kalshi_error_text(reply)
            self._update(client_order_id, {"status": STATUS_REJECTED, "error": f"HTTP {status}: {str(reply)[:500]}"})
            log.warning("manual order %s refused by Kalshi: HTTP %s %s", client_order_id, status, reply)
            return {"status": STATUS_REJECTED, "message": f"Kalshi refused it (HTTP {status}{': ' + detail if detail else ''}). Nothing was placed."}
        order = reply.get("order") if isinstance(reply, dict) and isinstance(reply.get("order"), dict) else reply
        order_id = order.get("order_id") if isinstance(order, dict) else None
        if status not in (200, 201) or not isinstance(order_id, str):
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": f"HTTP {status}: {str(reply)[:500]}"})
            log.error("manual order %s: no definite answer from Kalshi: HTTP %s %s", client_order_id, status, reply)
            return {"status": STATUS_UNKNOWN, "message": f"Kalshi answered HTTP {status}. Check Kalshi before placing again."}
        filled = _count_of(order, "fill_count_fp", "fill_count")
        remaining = _count_of(order, "remaining_count_fp", "remaining_count")
        if filled is None or remaining is None or remaining > 0:
            # Immediate-or-cancel leaves nothing open; a reply that does not say so gets one cancel, then a check by hand.
            cancelled = self._cancel(order_id, quote["ticker"])
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "kalshi_order_id": order_id,
                                           "error": f"reply without a zero remaining count (cancel "
                                                    f"{'confirmed' if cancelled else 'NOT confirmed'}): {str(reply)[:500]}"})
            log.error("manual order %s %s: Kalshi's reply does not show it closed (cancel %s): %s", client_order_id,
                      order_id, "confirmed" if cancelled else "not confirmed", reply)
            rest = "Its rest was cancelled." if cancelled else "It may still rest: cancel it on Kalshi."
            return {"status": STATUS_UNKNOWN, "message": f"Kalshi took it but did not say what filled. {rest} Check Kalshi."}
        cost = self._fill_cost(order_id, quote["side"], filled)
        outcome = STATUS_UNFILLED if filled <= 0 else STATUS_FILLED if filled >= count else STATUS_PARTIAL
        self._update(client_order_id, {"status": outcome, "kalshi_order_id": order_id, "count_filled": filled,
                                       "cost_filled": cost, "error": None})
        return {"status": outcome, "ticker": quote["ticker"], "side": quote["side"], "count": count, "filled": filled,
                "cost": cost, "message": _outcome_message(outcome, quote, count, filled, cost)}

    def _cancel(self, order_id: str, ticker: str) -> bool:
        """One cancel of an order Kalshi did not report closed, routed by
        market (the bot's own route). True only on Kalshi's 200/204."""
        try:
            status, _reply = self._client.request("DELETE", f"{ORDERS_PATH}/{order_id}?market_ticker={ticker}")
        except OSError as error:
            log.error("manual order %s: cancel got no answer: %s", order_id, error)
            return False
        return status in (200, 204)

    def _fill_cost(self, order_id: str, contract_side: str, filled: float) -> float | None:
        """Dollars paid for the order's fills, fees included, or None when the
        fills feed does not yet show all `filled` contracts (it can lag)."""
        if filled <= 0:
            return 0.0
        try:
            status, reply = self._client.request("GET", f"{FILLS_PATH}?order_id={order_id}&limit={FILLS_PAGE_LIMIT}")
        except OSError:
            return None
        if status != 200 or not isinstance(reply, dict) or not isinstance(reply.get("fills"), list):
            return None
        totals = fill_totals([fill for fill in reply["fills"] if isinstance(fill, dict)], contract_side)
        if isinstance(totals, str) or totals[0] < filled:
            return None
        return totals[1]

    def _update(self, client_order_id: str, fields: dict) -> None:
        self._store.update_manual_order(client_order_id, fields, self._now())


def _count_of(order: object, *keys: str) -> float | None:
    if not isinstance(order, dict):
        return None
    for key in keys:
        if key in order:
            try:
                value = float(order[key])
            except (TypeError, ValueError):
                return None
            return value if math.isfinite(value) and value >= 0 else None
    return None


def _kalshi_error_text(reply: object) -> str:
    """Kalshi's own words for a refusal ({"error": {"message": ...}}), or ""."""
    if isinstance(reply, dict):
        error = reply.get("error")
        if isinstance(error, dict):
            return str(error.get("message") or error.get("code") or "")[:200]
        if isinstance(error, str):
            return error[:200]
    return ""


def _outcome_message(outcome: str, quote: dict, count: int, filled: float, cost: float | None) -> str:
    price = _cents(quote["allIn"])
    paid = f" for ${cost:,.2f}" if cost is not None else ""
    if outcome == STATUS_FILLED:
        return f"Bought {filled:,.0f} contracts at {price} all-in{paid}."
    if outcome == STATUS_PARTIAL:
        return f"Bought {filled:,.0f} of {count:,} contracts at {price} all-in{paid}; the rest was not there. Nothing rests."
    return f"Nothing filled at {price}: the price was gone. Nothing rests."
