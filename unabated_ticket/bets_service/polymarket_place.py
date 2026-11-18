"""The Ticket's Place button at Polymarket US (Cal 2026-10-21, option A): Cal
clicks, and this desk sends ONE order on his main Polymarket US account,
signed with the key the polymarket_us source reads with (POLYMARKET_US_*;
Cal 2026-10-22, a polymarket.us key trades too). Nothing here runs on a timer: every order is
one POST /place_polymarket.json from the panel or the phone page.

Rules, in order, for every request (validate_polymarket_place_request's shape):
  1. the click was not already sent: 'pman-' + its clickId, held in
     bets.duckdb::manual_orders before the send. Polymarket US takes no
     client order id, so this row is the only guard against a second order;
  2. the server runner finds the board's Polymarket US line and the one
     market and contract side in Polymarket US's own event list that is that
     line (the game by its two teams, the market by its own sides), and reads
     the live price (GET /polymarket_quote.json, server/polymarket.js);
  3. the game has not started (the earlier of Unabated's and Polymarket US's
     start), and the runner read the Polymarket US this desk trades on; a live
     click (`live`, a Ticket at a break) instead needs the line on Unabated's
     live board, Polymarket US's event live and the panel's in-game fair read
     within 90 s (manual_orders.start_refusal);
  4. the live all-in price (ask + the taker fee) is at most PRICE_TOLERANCE
     (0.1c) worse than the Ticket's; a better one is taken. The fee is the
     market's own `feeCoefficient` x p x (1 - p), read on this quote (no
     coefficient, no order) — unless Polymarket US charged more than its
     market showed on the newest Place fill of at least FEE_LEARN_MIN_CONTRACTS
     at that same coefficient: then the charged one prices the order
     (charged_fee_coefficient), so a fee the market misstates costs at most
     that one fill. The order's limit is the ask, so this all-in is the most a
     contract can cost;
  5. the other side of the market is not held: buying NO against a held YES
     sells it (one instrument per market);
  6. the polymarket_us read the stake counted is ok and fresh;
  7. no earlier Place on this game, at any venue, may have filled after the
     read of its venue the stake was sized on (manual_orders.unseen_manual_fill);
  8. contracts = floor(stake / live all-in), at least one;
  9. the row is INSERTED (status 'sending') BEFORE the send;
 10. ONE immediate-or-cancel limit buy at the live ask (polymarket_orders.py),
     never re-sent; then the order is read until final for what filled. An
     order Polymarket US shows resting (an IOC never should) gets one cancel.
     A fill that cost more than the Ticket's price answers `overpaid` (red on
     the Ticket) and is logged as an error.

Order statuses (manual_orders.status on a 'polymarket_us' row): as Kalshi
Place's (manual_orders.py): sending, filled, partial, unfilled, rejected,
unknown. An 'unknown' row blocks the next Place on the game until a
polymarket_us read after it (rule 7); that read shows any fill as a position.

Side effects: a Polymarket US order on Cal's main account (money); INSERT /
UPDATE bets.duckdb::manual_orders (venue 'polymarket_us', incl. the fee
coefficient priced at and the fill's average price and commission). Network: GET
<runner>/polymarket_quote.json, then Polymarket US REST signed with Place's key.
"""
import logging
import math
import threading
import time
from collections.abc import Callable
from datetime import datetime, timedelta, timezone

from unabated_ticket.bets_service.manual_orders import (
    CLICK_ID_PATTERN, PRICE_TOLERANCE, STATUS_FILLED, STATUS_PARTIAL, STATUS_REFUSED, STATUS_REJECTED,
    STATUS_UNFILLED, STATUS_UNKNOWN, VENUE_POLYMARKET, _cents, _is_number, _iso_time, contracts_detail, filled_bet_message,
    optional_read_at, start_refusal, unpriced_fill_message, unseen_manual_fill, validate_live_fields,
)
from unabated_ticket.bets_service.polymarket_orders import (
    FINAL_ORDER_STATES, MAX_PRICE, MIN_PRICE, RESTING_ORDER_STATES, PolymarketError, PolymarketOrderClient,
    order_fill_fee, order_fill_totals,
)
from unabated_ticket.bets_service.store import BetsStore
from unabated_ticket.bets_service.taker_orders import DEFINITE_REJECTION_STATUSES

log = logging.getLogger(__name__)

POLYMARKET_ORDER_ID_PREFIX = "pman-"
RUNNER_POLYMARKET_QUOTE_PATH = "/polymarket_quote.json"
# The quote reads Polymarket US's events and one book (10 s each at most).
RUNNER_POLYMARKET_QUOTE_TIMEOUT_SEC = 30
# After the order is taken it is read this many times, this far apart, for
# its final state; no new read starts once PLACE_FINAL_READ_WINDOW_SEC has
# passed since the first, so one click holds well under the phone's 90 s wait.
PLACE_FINAL_READ_ATTEMPTS = 6
PLACE_FINAL_READ_SLEEP_SEC = 1.0
PLACE_FINAL_READ_WINDOW_SEC = 8.0
# The polymarket_us read the stake counted must be at most this old.
POLYMARKET_READ_SOURCE = "polymarket_us"
POLYMARKET_READ_MAX_AGE = timedelta(minutes=5)
POLYMARKET_PLACE_BET_TYPES = (1, 2, 3)
BET_TYPE_MONEYLINE = 1
CONTRACT_SIDES = ("yes", "no")
# A market's feeCoefficient above this is taken as misread, not as a fee (server/polymarket.js MAX_FEE_COEFFICIENT).
MAX_FEE_COEFFICIENT = 0.25
# Commission is rounded to the cent per order, so a fill this small says little about the rate
# (fills under ~17 contracts read 0.0 to 0.12 on 2026-11-08); only larger ones correct the market's.
FEE_LEARN_MIN_CONTRACTS = 50
# A charged coefficient this far over the market's is a real difference, not that rounding
# (0.0948-0.0951 on 15 fills charged 0.095).
FEE_LEARN_SLACK = 0.002
# The commission's cent rounding a fill may show above the all-in it was priced at, in dollars per order.
COMMISSION_ROUNDING_DOLLARS = 0.01


def validate_polymarket_place_request(body: object) -> dict | str:
    """A POST /place_polymarket.json body, or what was expected and what was found.
    {clickId, eventId, periodTypeId, betTypeId (1|2|3), sideIndex, points (null on
    a moneyline), stake (dollars), ticketAllIn (the Ticket's Polymarket US price as
    a probability), polymarketReadAt (the polymarket_us read the stake was sized
    on, ISO), kalshiReadAt / novigReadAt / prophetxReadAt / bfaReadAt (ISO, or null when never read)}"""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    click_id = body.get("clickId")
    if not isinstance(click_id, str) or not CLICK_ID_PATTERN.match(click_id):
        return f"clickId must be 8-48 letters, digits or dashes, got {click_id!r}"
    whole = {}
    for field, allowed in (("eventId", None), ("periodTypeId", None), ("betTypeId", POLYMARKET_PLACE_BET_TYPES),
                           ("sideIndex", (0, 1))):
        value = body.get(field)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0 or (allowed and value not in allowed):
            return f"{field} must be {'one of ' + str(allowed) if allowed else 'a whole number'}, got {value!r}"
        whole[field] = value
    points = body.get("points")
    if whole["betTypeId"] == BET_TYPE_MONEYLINE:
        points = None
    elif not _is_number(points) or not math.isfinite(points):
        return f"points must be a number (a spread or total), got {points!r}"
    stake = body.get("stake")
    if not _is_number(stake) or not math.isfinite(stake) or stake <= 0:
        return f"stake must be dollars above 0, got {stake!r}"
    ticket_all_in = body.get("ticketAllIn")
    if not _is_number(ticket_all_in) or not MIN_PRICE <= ticket_all_in < 1:
        return f"ticketAllIn must be a probability from {MIN_PRICE} to under 1, got {ticket_all_in!r}"
    polymarket_read_at = _iso_time(body.get("polymarketReadAt"))
    if polymarket_read_at is None:
        return ("polymarketReadAt must be the ISO time of the Polymarket US read the stake was sized on, "
                f"got {body.get('polymarketReadAt')!r}")
    read_ats = {}
    for field in ("kalshiReadAt", "novigReadAt", "prophetxReadAt", "bfaReadAt"):
        read_at = optional_read_at(body, field)
        if isinstance(read_at, str):
            return read_at
        read_ats[field] = read_at
    live = validate_live_fields(body)
    if isinstance(live, str):
        return live
    return {"click_id": click_id, **whole, "points": None if points is None else float(points), **live,
            "stake": float(stake), "ticket_all_in": float(ticket_all_in), "polymarket_read_at": polymarket_read_at,
            "kalshi_read_at": read_ats["kalshiReadAt"], "novig_read_at": read_ats["novigReadAt"],
            "prophetx_read_at": read_ats["prophetxReadAt"], "bfa_read_at": read_ats["bfaReadAt"]}


def all_in_of(ask: float, fee_coefficient: float) -> float:
    """One contract's price with Polymarket US's taker fee (server/polymarket.js allInOf)."""
    return round(ask + fee_coefficient * ask * (1 - ask), 6)


def fee_coefficient_of_fill(fill_avg_price: float, fill_commission: float, count_filled: float) -> float | None:
    """The coefficient Polymarket US charged on one fill: commission / (contracts x p x (1 - p)), or None."""
    spread = count_filled * fill_avg_price * (1 - fill_avg_price)
    if spread <= 0:
        return None
    return fill_commission / spread


def charged_fee_coefficient(market_coefficient: float, latest_fill: dict | None) -> tuple[float, str | None]:
    """The coefficient to price this order at, and words when it is not the market's.

    The newest sizable Place fill (store.latest_polymarket_fee_fill) overrides the
    market's coefficient only when it was priced at that same coefficient and charged
    more than it: Polymarket US misstated its fee, so its word no longer holds. A
    market showing another coefficient since is a new word and is taken."""
    if latest_fill is None or abs(latest_fill["fee_coefficient"] - market_coefficient) > 1e-9:
        return market_coefficient, None
    charged = fee_coefficient_of_fill(latest_fill["fill_avg_price"], latest_fill["fill_commission"],
                                      latest_fill["count_filled"])
    if charged is None or charged <= market_coefficient + FEE_LEARN_SLACK:
        return market_coefficient, None
    return charged, (f"your last Polymarket US fill was charged a {charged:.4f} fee rate, above the "
                     f"{market_coefficient:.4f} its market shows, so this one is priced at {charged:.4f}")


def held_id_of(slug: str, contract_side: str) -> str:
    """The bets row id the polymarket_us source gives a held side (sources/polymarket_us.py)."""
    return f"{POLYMARKET_READ_SOURCE}:{slug}:{contract_side}"


class ManualPolymarketDesk:
    """Sends Cal's Polymarket US Place-button orders on his main account, one per click."""

    def __init__(self, client: PolymarketOrderClient, store: BetsStore, quote_reader: Callable[[dict], dict],
                 now: Callable[[], datetime] = lambda: datetime.now(timezone.utc),
                 sleep: Callable[[float], None] = time.sleep, clock: Callable[[], float] = time.monotonic):
        """`client` signs with the POLYMARKET_US_* key."""
        self._client = client
        self._store = store
        self._quote_reader = quote_reader
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
        client_order_id = f"{POLYMARKET_ORDER_ID_PREFIX}{request['click_id']}"
        quote = self._quote_reader({"eventId": request["eventId"], "periodTypeId": request["periodTypeId"],
                                    "betTypeId": request["betTypeId"], "sideIndex": request["sideIndex"],
                                    "points": request["points"], "live": request["live"]})
        refusal = self._refusal(request, quote)
        if refusal:
            return {"status": STATUS_REFUSED, "message": f"Not placed: {refusal}"}
        all_in = quote["chargedAllIn"]
        count = math.floor(request["stake"] / all_in)
        if count < 1:
            return {"status": STATUS_REFUSED,
                    "message": f"Not placed: ${request['stake']:.2f} buys no contract at {_cents(all_in)}."}
        line = quote.get("line") if isinstance(quote.get("line"), dict) else {}
        inserted = self._store.insert_manual_order({
            "client_order_id": client_order_id, "placed_at": self._now(), "venue": VENUE_POLYMARKET,
            "event_id": request["eventId"], "ticker": quote["slug"], "contract_side": quote["side"],
            "label": quote.get("label"), "ticket_all_in": request["ticket_all_in"], "ask_price": quote["ask"],
            "all_in_price": all_in, "fee_coefficient": quote["feeCoefficient"],
            "stake": request["stake"], "count_requested": count, "status": "sending",
            "updated_at": self._now(), "event_start": _iso_time(quote.get("eventStart")),
            "league": line.get("league"), "bet_type": line.get("betType"), "period": line.get("period"),
            "side": line.get("side"), "points": line.get("points"), "away_team": line.get("awayTeam"),
            "home_team": line.get("homeTeam"), "away_team_id": line.get("awayTeamId"),
            "home_team_id": line.get("homeTeamId"),
        })
        if not inserted:
            return {"status": STATUS_REFUSED, "message": "Not placed: this click was already sent."}
        log.info("polymarket place %s: %s %s x%d @ %s (all-in %s at fee %s, market fee %s, Ticket %s)", client_order_id,
                 quote["slug"], quote["side"], count, quote["ask"], all_in, quote["chargedFeeCoefficient"],
                 quote["feeCoefficient"], request["ticket_all_in"])
        return self._send(client_order_id, {**quote, "ticketAllIn": request["ticket_all_in"]}, count)

    def _refusal(self, request: dict, quote: dict) -> str | None:
        """Why this click must not become an order, in words, or None."""
        if not quote.get("ok"):
            return f"{quote.get('detail') or quote.get('reason') or 'no Polymarket US market for this line'}."
        clock = start_refusal(request, quote, self._now())
        if clock:
            return clock
        if quote.get("polymarketGatewayUrl") != self._client.gateway_url:
            return (f"the runner read Polymarket US at {quote.get('polymarketGatewayUrl')!r} but this desk trades "
                    f"{self._client.gateway_url}.")
        if not isinstance(quote.get("slug"), str) or not quote["slug"] or quote.get("side") not in CONTRACT_SIDES:
            return "the runner did not name the market and its side."
        ask, fee_coefficient = quote.get("ask"), quote.get("feeCoefficient")
        if not _is_number(ask) or not MIN_PRICE <= ask <= MAX_PRICE:
            return f"Polymarket US's ask {ask!r} is not a price."
        if not _is_number(fee_coefficient) or not 0 <= fee_coefficient <= MAX_FEE_COEFFICIENT:
            return f"the runner did not say Polymarket US's fee on this market (feeCoefficient {fee_coefficient!r})."
        if not _is_number(quote.get("allIn")) or abs(quote["allIn"] - all_in_of(ask, fee_coefficient)) > 1e-6:
            return (f"the runner's all-in {quote.get('allIn')!r} is not {_cents(ask)} plus the market's "
                    f"{fee_coefficient} fee.")
        charged, why = charged_fee_coefficient(fee_coefficient,
                                               self._store.latest_polymarket_fee_fill(FEE_LEARN_MIN_CONTRACTS))
        all_in = all_in_of(ask, charged)
        if all_in > request["ticket_all_in"] + PRICE_TOLERANCE:
            reason = f" ({why})" if why else ""
            return (f"Polymarket US is at {_cents(all_in)} all-in now ({_cents(ask)} + fee), more than 0.1¢ worse "
                    f"than the Ticket's {_cents(request['ticket_all_in'])}{reason}.")
        quote["chargedFeeCoefficient"], quote["chargedAllIn"] = charged, all_in
        other_side = "no" if quote["side"] == "yes" else "yes"
        if self._store.bet_is_open(held_id_of(quote["slug"], other_side)):
            return (f"you hold the {other_side.upper()} side of {quote['slug']}; buying {quote['side'].upper()} "
                    "would sell it.")
        stale = self._stale_read()
        if stale:
            return stale
        if unseen_manual_fill(self._store, request["eventId"], request["kalshi_read_at"], request["novig_read_at"],
                              request["polymarket_read_at"], request["prophetx_read_at"],
                              request["bfa_read_at"]):
            return ("your last Place on this game is not on the Bets tab yet, so this stake does not count it. "
                    "Wait for the next read (up to a minute) and place again.")
        return None

    def _stale_read(self) -> str | None:
        """Words for a polymarket_us read that is failing or old, or None."""
        read = self._store.source_status().get(POLYMARKET_READ_SOURCE)
        fetched_at = _iso_time(read.get("fetchedAt")) if read else None
        if not read or not read.get("ok") or fetched_at is None:
            return "the Polymarket US read is failing, so the stake may miss a bet. Try again once it reads."
        if self._now() - fetched_at > POLYMARKET_READ_MAX_AGE:
            return (f"the Polymarket US read is over {POLYMARKET_READ_MAX_AGE.seconds // 60} minutes old, "
                    "so the stake may miss a bet.")
        return None

    def _send(self, client_order_id: str, quote: dict, count: int) -> dict:
        """ONE IOC buy, never re-sent; then the order read until final."""
        try:
            status, reply = self._client.place_ioc_order(quote["slug"], quote["side"], quote["ask"], count)
        except OSError as error:
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": f"no answer: {error}"})
            log.error("polymarket place %s: no answer from Polymarket US: %s", client_order_id, error)
            return {"status": STATUS_UNKNOWN, "message": "Polymarket US did not answer. Check Polymarket US before placing again."}
        if status in DEFINITE_REJECTION_STATUSES:
            self._update(client_order_id, {"status": STATUS_REJECTED, "error": f"HTTP {status}: {str(reply)[:500]}"})
            log.warning("polymarket place %s refused: HTTP %s %s", client_order_id, status, reply)
            detail = _error_text(reply)
            return {"status": STATUS_REJECTED,
                    "message": f"Polymarket US refused it (HTTP {status}{': ' + detail if detail else ''}). Nothing was placed."}
        order_id = reply.get("id") if isinstance(reply, dict) else None
        if status not in (200, 201) or not isinstance(order_id, str) or not order_id:
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": f"HTTP {status}: {str(reply)[:500]}"})
            log.error("polymarket place %s: no definite answer: HTTP %s %s", client_order_id, status, reply)
            return {"status": STATUS_UNKNOWN,
                    "message": f"Polymarket US answered HTTP {status}. Check Polymarket US before placing again."}
        self._update(client_order_id, {"polymarket_order_id": order_id})
        return self._follow_order(client_order_id, order_id, quote, count)

    def _follow_order(self, client_order_id: str, order_id: str, quote: dict, count: int) -> dict:
        """Read the order until Polymarket US shows it final, then what filled.
        A resting order (an IOC never should be) gets one cancel."""
        order = None
        deadline = self._clock() + PLACE_FINAL_READ_WINDOW_SEC
        for attempt in range(PLACE_FINAL_READ_ATTEMPTS):
            if attempt:
                if self._clock() >= deadline:
                    break
                self._sleep(PLACE_FINAL_READ_SLEEP_SEC)
            try:
                status, reply = self._client.get_order(order_id)
            except OSError as error:
                log.warning("polymarket place %s: order read failed: %s", client_order_id, error)
                continue
            read = reply.get("order") if status == 200 and isinstance(reply, dict) else None
            if not isinstance(read, dict):
                continue
            order = read
            if order.get("state") in FINAL_ORDER_STATES:
                return self._record(client_order_id, order, quote, count)
        if order is not None and order.get("state") in RESTING_ORDER_STATES:
            return self._cancel_resting(client_order_id, order_id, order, quote)
        totals = order_fill_totals(order, quote["side"]) if order is not None else "the order never read"
        filled = totals[0] if isinstance(totals, tuple) else None
        self._update(client_order_id, {"status": STATUS_UNKNOWN, "count_filled": filled,
                                       "error": f"never read final (last state {order.get('state') if order else None})"})
        log.error("polymarket place %s: order %s never read final", client_order_id, order_id)
        return {"status": STATUS_UNKNOWN, "message": "Polymarket US took the order but has not said what filled. "
                                                     "Check Polymarket US before placing again."}

    def _cancel_resting(self, client_order_id: str, order_id: str, order: dict, quote: dict) -> dict:
        log.error("polymarket place %s: order %s rests (%s) after an IOC send; cancelling", client_order_id, order_id,
                  order.get("state"))
        try:
            status, reply = self._client.cancel_order(order_id, quote["slug"])
        except OSError as error:
            status, reply = None, f"no answer: {error}"
        self._update(client_order_id, {"status": STATUS_UNKNOWN,
                                       "error": f"rested after IOC; cancel answered {status} {str(reply)[:200]}"})
        rest = "A cancel was sent" if status in (200, 204) else "The cancel was NOT confirmed: cancel it in the app"
        return {"status": STATUS_UNKNOWN, "message": f"Polymarket US left the order resting. {rest}. "
                                                     "Check Polymarket US before placing again."}

    def _record(self, client_order_id: str, order: dict, quote: dict, count: int) -> dict:
        totals = order_fill_totals(order, quote["side"])
        if isinstance(totals, str):
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": f"final but unreadable: {totals}"})
            log.error("polymarket place %s: %s", client_order_id, totals)
            return {"status": STATUS_UNKNOWN, "message": "Polymarket US finished the order but did not say what filled. "
                                                         "Check Polymarket US before placing again."}
        filled, paid = totals
        if filled <= 0 and order.get("state") == "ORDER_STATE_REJECTED":
            outcome = STATUS_REJECTED
        else:
            outcome = STATUS_UNFILLED if filled <= 0 else STATUS_FILLED if filled >= count else STATUS_PARTIAL
        fee = order_fill_fee(order, quote["side"]) if filled > 0 else None
        overpaid = _overpaid_words(filled, paid, quote["ticketAllIn"])
        error = f"rejected: {str(order)[:300]}" if outcome == STATUS_REJECTED else overpaid
        self._update(client_order_id, {"status": outcome, "count_filled": filled, "cost_filled": paid, "error": error,
                                       "fill_avg_price": fee[0] if fee else None,
                                       "fill_commission": fee[1] if fee else None})
        log.info("polymarket place %s: %s, %s of %d for %s (fill fee %s)", client_order_id, outcome, filled, count, paid,
                 fee)
        if overpaid:
            log.error("polymarket place %s: %s", client_order_id, overpaid)
        message = _outcome_message(outcome, quote, count, filled, paid)
        if overpaid:
            message += f" {overpaid[0].upper()}{overpaid[1:]}."
        return {"status": outcome, "slug": quote["slug"], "side": quote["side"], "count": count, "filled": filled,
                "cost": paid, "overpaid": bool(overpaid), "message": message}

    def _update(self, client_order_id: str, fields: dict) -> None:
        self._store.update_manual_order(client_order_id, fields, self._now())


def _error_text(reply: object) -> str:
    """Polymarket US's own words for a refusal, or ""."""
    if isinstance(reply, dict):
        for key in ("message", "error", "code"):
            if isinstance(reply.get(key), str):
                return reply[key][:200]
    if isinstance(reply, str):
        return reply[:200]
    return ""


def _overpaid_words(filled: float, paid: float | None, ticket_all_in: float) -> str | None:
    """Words when a fill cost more a contract than the Ticket's price (past the 0.1c tolerance and the
    commission's cent rounding), else None."""
    if not filled or filled <= 0 or paid is None:
        return None
    each = paid / filled
    if each <= ticket_all_in + PRICE_TOLERANCE + COMMISSION_ROUNDING_DOLLARS / filled:
        return None
    return (f"this cost {_cents(each)} a contract with Polymarket US's fee, more than the Ticket's "
            f"{_cents(ticket_all_in)}: its fee was above what the market showed, so Place there now prices at the "
            f"fee it charged")


def _outcome_message(outcome: str, quote: dict, count: int, filled: float, paid: float | None) -> str:
    """A Polymarket US contract pays $1, so a fill pays back `filled` dollars.
    Without the cost it says what filled, at the quoted price."""
    price = _cents(quote["chargedAllIn"])
    if outcome == STATUS_REJECTED:
        return "Polymarket US rejected the order. Nothing was placed."
    if outcome not in (STATUS_FILLED, STATUS_PARTIAL):
        return f"Nothing filled at {price}: the price was gone. Nothing rests."
    rest = "." if outcome == STATUS_FILLED else "; the rest was not there. Nothing rests."
    if paid is None:
        return unpriced_fill_message("Polymarket US", quote.get("label"), filled, count, quote["chargedAllIn"], rest)
    return filled_bet_message(quote.get("label"), paid, filled, contracts_detail(filled, count, paid)) + rest
