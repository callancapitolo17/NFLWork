"""The automated taker's Novig order desk (automated_taker_plan.md § Novig):
the only place a bot order is sent to Novig. The server runner
(server/runner.js + server/taker.js) decides; this module checks the stop
files, the holdings' freshness and the Novig $500 ceiling, sends ONE
immediate-or-cancel buy at the ask the decision was checked at, and reads what
filled. Kalshi's desk (taker_orders.py) is untouched apart from asking this
module about Novig bot orders on the same game.

Armed only on the bot's own Novig subaccount: a trading key
(NOVIG_TAKER_KEY_ID + NOVIG_TAKER_PRIVATE_KEY_PATH) that reads its own
subaccount's balance AND is refused Novig's key listing (only a management key
lists keys, and only a management key moves money). The money in that
subaccount is then all the bot can lose. novig_subaccount_setup.py opens and
funds it. Without it, POST /taker/novig_order.json answers 503 and nothing is sent.

Novig facts this rests on (docs.novig.com, read 2026-10-16):
  * qty counts 1c contracts (100 = a $1 payout): cost = qty x price / 100;
  * an IOC that fills in part cancels the rest; one that fills nothing is
    REJECTED, so a bot order never rests;
  * 201 means queued, not placed; GET /v3/orders/{id} can 404 just after;
  * clientId is NOT deduplicated: an order is sent once and never retried; a
    send with no answer is looked up by outcome + clientId across every status;
  * pre-game game markets charge no taker fee (the runner refuses a market
    whose fee is charged ALWAYS); 451 GEOLOCATION_EXPIRED refuses placements
    once the account's phone has not checked its location for 3 days.

Rules, in order, for every request:
  1. a stop file exists (config.TAKER_STOP_PATH stops both bots,
     config.NOVIG_TAKER_STOP_PATH this one) -> refused, "stopped";
  2. the same clientId was already sent -> refused, "duplicate";
  3. the decision is fresh (MAX_DECISION_AGE_SEC), the game has not started,
     and the runner priced the Novig this desk trades -> else "stale_decision"
     / "started" / "env_mismatch";
  4. the holdings the stake was sized on are current: the bets service's
     reads of Cal's Novig app account and of the bot subaccount both ok and
     under MAX_HOLDINGS_AGE_SEC old -> else "stale_holdings";
  5. open orders are brought up to date when the last refresh is older than REFRESH_SEC;
  6. unseen fill: an earlier bot order on the SAME GAME, Novig or Kalshi,
     that may hold contracts must be among the held bets the runner sized
     this stake on (`heldIds`) -> else "unseen_fill";
  7. ceiling: qty cut so open Novig exposure + this order stays at or under
     OPEN_COST_CEILING; none left -> "ceiling";
  8. the row is INSERTED (status "sending", full cost) BEFORE the send;
  9. one IOC buy; then the order is read until final, and its fills.

Order statuses (taker.duckdb::novig_taker_orders.status):
  sending   inserted, the POST not yet answered
  unknown   no definite answer (it may have landed; looked up by clientId), or a read failed
  sent      Novig queued it (201); read until final
  resting   Novig shows it OPEN (never expected of an IOC): cancelled, full cost held
  filled | partial | unfilled   what the fills read found; `fills_final` once
            read FILLS_SETTLE_SEC after the order went final
  rejected  Novig refused it before taking it (a 451 answers the runner
            "geo_blocked"), or it never reached Novig

Side effects: Novig orders on the bot's key (money); INSERTs / UPDATEs
taker.duckdb::novig_taker_orders (taker_store.py). Network: Novig REST, signed
with the bot's trading key only; public catalog reads for settlement.
"""
import logging
import math
import threading
import time
from collections.abc import Callable
from datetime import datetime, timedelta, timezone
from pathlib import Path

from unabated_ticket.bets_service.novig_orders import (FINAL_ORDER_STATUSES, GEO_STATUS, QTY_PER_DOLLAR, STATUS_OPEN, NovigError, NovigOrderClient, cost_of,
                                                       error_code_of, is_on_price_grid, order_totals)
from unabated_ticket.bets_service.taker_orders import kalshi_unseen_fill
from unabated_ticket.bets_service.taker_store import TakerStore

log = logging.getLogger(__name__)

# Cal, 2026-10-16: Novig's own $500 cap on open cost, apart from Kalshi's. The
# Kelly bankroll that sizes the stake is the runner's (taker.js NOVIG_KELLY_BANKROLL).
OPEN_COST_CEILING = 500.0
# The bot subaccount's bets are recorded as "novig:<order id>:bot" (sources/novig_subaccount.py).
BOT_ACCOUNT_NAME = "bot"
REFRESH_SEC = 30
MAX_DECISION_AGE_SEC = 60
# The stake is conditional on the held Novig bets: the app poll is 60 s, the subaccount poll 30 s.
MAX_HOLDINGS_AGE_SEC = 5 * 60
HOLDINGS_SOURCES = ("novig", "novig_bot")
# A fills read this long after the order went final is taken as final.
FILLS_SETTLE_SEC = 60
# A send with no answer is looked up by clientId once it is this old.
LOST_ORDER_LOOKUP_AFTER_SEC = 60
# Novig answered the order and refused it outright: nothing was placed. 423
# (locked), 429 and 5xx may hide a queued order, so they are looked up instead.
DEFINITE_REJECTION_STATUSES = (400, 401, 403, 404, 413, 422, GEO_STATUS)
# After a 201 the order is read this many times, this far apart, for its final status.
FINAL_READ_ATTEMPTS = 6
FINAL_READ_SLEEP_SEC = 0.5
SETTLEMENT_CHECK_SEC = 5 * 60
MARKET_SETTLED = "SETTLED"
STATUSES_MAYBE_HOLDING = ("sending", "unknown", "sent", "resting", "filled", "partial")
MAX_HELD_IDS = 1000
MAX_HELD_ID_LENGTH = 200
BET_TYPES = ("spread", "total", "moneyline")
SIDES = ("away", "home", "over", "under")


def held_id_of(novig_order_id: str) -> str:
    """The bets service's record id for a bot-subaccount Novig order."""
    return f"novig:{novig_order_id}:{BOT_ACCOUNT_NAME}"


def _may_hold_contracts(order: dict) -> bool:
    if order["status"] in STATUSES_MAYBE_HOLDING:
        return True
    return order["status"] == "unfilled" and not order["fills_final"]


def novig_unseen_fill(store: TakerStore, event_id: int, held_ids: set[str]) -> str | None:
    """Why an add on Unabated event `event_id` would be sized blind, or None:
    an open Novig bot order there that may hold contracts and is not among
    `held_ids` (one with no Novig order id yet can never be). Store reads only;
    Kalshi's desk asks it too (taker_orders.TakerDesk other_venue_unseen)."""
    for order in store.open_novig_orders():
        if order["event_id"] != event_id or not _may_hold_contracts(order):
            continue
        if order["novig_order_id"] is None or held_id_of(order["novig_order_id"]) not in held_ids:
            shown = held_id_of(order["novig_order_id"]) if order["novig_order_id"] else "its Novig order id (none yet)"
            return (f"Novig bot order {order['client_id']} ({order['status']}) is not yet among the held bets "
                    f"the stake was sized on: no {shown}")
    return None


def _iso_time(value: object) -> datetime | None:
    if not isinstance(value, str):
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return parsed if parsed.tzinfo is not None else None


def _uuid_like(value: object) -> bool:
    return isinstance(value, str) and len(value) == 36 and value.count("-") == 4


def _optional_text(value: object) -> bool:
    return value is None or (isinstance(value, str) and len(value) <= 200)


def validate_order_request(body: object) -> dict | str:
    """A POST /taker/novig_order.json body, or what was expected and what was found."""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    for field in ("clientId", "outcomeId", "marketId"):
        if not _uuid_like(body.get(field)):
            return f"{field} must be a UUID, got {body.get(field)!r}"
    if not isinstance(body.get("rowKey"), str) or not body["rowKey"]:
        return f"rowKey must be a non-empty string, got {body.get('rowKey')!r}"
    price = body.get("price")
    if not isinstance(price, (int, float)) or isinstance(price, bool) or not is_on_price_grid(price):
        return f"price must be a price on Novig's grid, got {price!r}"
    qty = body.get("qty")
    if not isinstance(qty, int) or isinstance(qty, bool) or qty < 1:
        return f"qty must be a whole number of 1c contracts >= 1, got {qty!r}"
    event_id = body.get("eventId")
    if not isinstance(event_id, int) or isinstance(event_id, bool) or event_id < 1:
        return f"eventId must be Unabated's event id, got {event_id!r}"
    held_ids = body.get("heldIds")
    if (not isinstance(held_ids, list) or len(held_ids) > MAX_HELD_IDS
            or any(not isinstance(held, str) or len(held) > MAX_HELD_ID_LENGTH for held in held_ids)):
        return f"heldIds must be a list of at most {MAX_HELD_IDS} bet ids the stake was sized on, got {held_ids!r}"
    start = _iso_time(body.get("eventStart"))
    if start is None:
        return f"eventStart must be an ISO time, got {body.get('eventStart')!r}"
    decided_at = _iso_time(body.get("decidedAt"))
    if decided_at is None:
        return f"decidedAt must be an ISO time, got {body.get('decidedAt')!r}"
    novig_base_url = body.get("novigBaseUrl")
    if not isinstance(novig_base_url, str) or not novig_base_url.startswith("http"):
        return f"novigBaseUrl must be the Novig API the decision was priced on, got {novig_base_url!r}"
    line = body.get("line")
    if not isinstance(line, dict):
        return f"line must be the board line priced, got {line!r}"
    if not isinstance(line.get("league"), str) or not line["league"]:
        return f"line.league must be the feed's league path, got {line.get('league')!r}"
    if line.get("betType") not in BET_TYPES:
        return f"line.betType must be one of {BET_TYPES}, got {line.get('betType')!r}"
    if not isinstance(line.get("period"), str) or not line["period"]:
        return f"line.period must be a period name, got {line.get('period')!r}"
    if line.get("side") not in SIDES:
        return f"line.side must be one of {SIDES}, got {line.get('side')!r}"
    points = line.get("points")
    if points is not None and (not isinstance(points, (int, float)) or isinstance(points, bool)):
        return f"line.points must be a number or null, got {points!r}"
    for field in ("awayTeam", "homeTeam", "awayTeamId", "homeTeamId"):
        if not _optional_text(line.get(field)):
            return f"line.{field} must be a short string or null, got {line.get(field)!r}"
    return {"client_id": body["clientId"], "outcome_id": body["outcomeId"], "market_id": body["marketId"],
            "row_key": body["rowKey"], "price": float(price), "qty": qty, "event_id": event_id,
            "held_ids": set(held_ids), "event_start": start, "decided_at": decided_at,
            "novig_base_url": novig_base_url.rstrip("/"),
            "league": line["league"], "bet_type": line["betType"], "period": line["period"], "side": line["side"],
            "points": None if points is None else float(points), "away_team": line.get("awayTeam"),
            "home_team": line.get("homeTeam"), "away_team_id": line.get("awayTeamId"),
            "home_team_id": line.get("homeTeamId")}


def account_problem(client: NovigOrderClient) -> tuple[str | None, float | None]:
    """(None, the subaccount's balance) when the key is a trading key on its
    own subaccount: it reads that balance and Novig refuses it the key
    listing. Otherwise (why not, None). GETs only."""
    try:
        balance = client.balance()
    except (NovigError, OSError) as error:
        return f"the key cannot read its subaccount's balance: {error}", None
    try:
        management = client.is_management_key()
    except OSError as error:
        return f"could not check the key's scope: no answer ({error})", None
    if management is True:
        return "the key lists Novig's keys: it is a management key, which can move money; the bot needs a trading key", None
    if management is None:
        return "Novig's answer to the key listing says neither trading nor management", None
    return None, balance


class NovigTakerDesk:
    """Places the runner's Novig orders, one at a time, under the Novig ceiling."""

    def __init__(self, client: NovigOrderClient, store: TakerStore, stop_paths: tuple[Path, ...],
                 holdings_status: Callable[[], dict[str, dict]],
                 now: Callable[[], datetime] = lambda: datetime.now(timezone.utc),
                 sleep: Callable[[float], None] = time.sleep):
        self._client = client
        self._store = store
        self._stop_paths = tuple(Path(path) for path in stop_paths)
        # The bets service's source status (store.source_status): name -> {fetchedAt, ok, error}.
        self._holdings_status = holdings_status
        self._now = now
        self._sleep = sleep
        self._lock = threading.Lock()
        self._last_refresh_at: datetime | None = None
        self._settlement_checked_at: dict[str, datetime] = {}

    def run_refresh_loop(self, stop: threading.Event) -> None:
        """The desk's own thread: open orders brought up to date every REFRESH_SEC."""
        while not stop.wait(REFRESH_SEC):
            try:
                self.refresh()
            except Exception:  # noqa: BLE001 — the sweep must outlive one bad answer
                log.exception("novig taker desk refresh failed")

    def refresh(self) -> None:
        with self._lock:
            self._refresh_open_orders()

    def place(self, request: dict) -> dict:
        """One order (validate_order_request's shape). Returns {status, ...}:
        stopped | duplicate | stale_decision | started | env_mismatch |
        stale_holdings | unseen_fill | ceiling (nothing sent), or the sent
        order's filled | partial | unfilled | resting | rejected | unknown, or
        geo_blocked (451: Novig refused it for the phone's location; nothing placed)."""
        with self._lock:
            for stop_path in self._stop_paths:
                if stop_path.exists():
                    return {"status": "stopped", "detail": f"stop file {stop_path} is present"}
            if self._store.has_novig_order(request["client_id"]):
                return {"status": "duplicate", "detail": f"{request['client_id']} was already sent"}
            refused = self._refuse_unfit(request)
            if refused:
                return refused
            if self._last_refresh_at is None or self._now() - self._last_refresh_at >= timedelta(seconds=REFRESH_SEC):
                self._refresh_open_orders()
                refused = self._refuse_unfit(request)
                if refused:
                    return refused
            unseen = (novig_unseen_fill(self._store, request["event_id"], request["held_ids"])
                      or kalshi_unseen_fill(self._store, request["event_id"], request["held_ids"]))
            if unseen:
                return {"status": "unseen_fill", "detail": unseen}
            headroom = OPEN_COST_CEILING - self._store.novig_open_exposure()
            qty = min(request["qty"], math.floor(headroom * QTY_PER_DOLLAR / request["price"]))
            if qty < 1:
                return {"status": "ceiling", "detail": f"${headroom:.2f} of ${OPEN_COST_CEILING:.0f} left"}
            return self._send(request, qty)

    def _refuse_unfit(self, request: dict) -> dict | None:
        now = self._now()
        if request["novig_base_url"] != self._client.base_url:
            return {"status": "env_mismatch",
                    "detail": f"the runner priced {request['novig_base_url']}, this desk trades {self._client.base_url}"}
        age = now - request["decided_at"]
        if age > timedelta(seconds=MAX_DECISION_AGE_SEC):
            return {"status": "stale_decision", "detail": f"decided {age.total_seconds():.0f} s ago, over {MAX_DECISION_AGE_SEC} s"}
        if request["event_start"] <= now:
            return {"status": "started", "detail": f"the game started at {request['event_start'].isoformat()}"}
        stale = self._stale_holdings(now)
        if stale:
            return {"status": "stale_holdings", "detail": stale}
        return None

    def _stale_holdings(self, now: datetime) -> str | None:
        """Why the held Novig bets the stake was sized on may be wrong, or None."""
        status = self._holdings_status()
        for name in HOLDINGS_SOURCES:
            source = status.get(name)
            if not source:
                return f"the bets service does not read {name}"
            if source.get("ok") is not True:
                return f"the {name} read failed: {source.get('error') or 'no detail'}"
            fetched_at = _iso_time(source.get("fetchedAt"))
            if fetched_at is None or now - fetched_at > timedelta(seconds=MAX_HOLDINGS_AGE_SEC):
                return f"the {name} read is older than {MAX_HOLDINGS_AGE_SEC // 60} min ({source.get('fetchedAt') or 'never'})"
        return None

    def _send(self, request: dict, qty: int) -> dict:
        client_id = request["client_id"]
        self._store.insert_novig_order({
            **request, "placed_at": self._now(), "qty_requested": qty, "status": "sending",
            "open_cost": round(cost_of(qty, request["price"]), 5)})
        try:
            status, reply = self._client.place_ioc_order(request["outcome_id"], request["price"], qty, client_id)
        except (OSError, NovigError) as error:
            self._store.update_novig_order(client_id, {"status": "unknown", "error": f"no answer: {error}"})
            log.error("novig taker order %s %s: no answer from Novig: %s", client_id, request["outcome_id"], error)
            return {"status": "unknown", "qty": qty, "detail": f"no answer from Novig: {error}; full cost held"}
        order_id = reply.get("orderId") if isinstance(reply, dict) else None
        if status in DEFINITE_REJECTION_STATUSES:
            code = error_code_of(reply)
            self._store.update_novig_order(client_id, {"status": "rejected", "open_cost": 0.0, "settled_at": self._now(),
                                                       "error": f"HTTP {status} {code or ''}: {str(reply)[:400]}"})
            if status == GEO_STATUS:
                log.error("novig taker: Novig refused placement for location (%s): open the Novig app on Cal's phone", code)
            else:
                log.warning("novig taker order %s refused: HTTP %s %s", client_id, status, reply)
            if status == GEO_STATUS:
                return {"status": "geo_blocked", "qty": qty,
                        "detail": f"Novig HTTP 451 {code or ''}: open the Novig app on Cal's phone to check its location".strip()}
            return {"status": "rejected", "qty": qty, "detail": f"Novig HTTP {status} {code or ''}".strip()}
        if status != 201 or not isinstance(order_id, str) or not order_id:
            self._store.update_novig_order(client_id, {"status": "unknown",
                                                       "error": f"no definite answer: HTTP {status}: {str(reply)[:300]}"})
            log.error("novig taker order %s: no definite answer: HTTP %s %s", client_id, status, reply)
            return {"status": "unknown", "qty": qty, "detail": f"no definite answer from Novig: HTTP {status}; full cost held"}
        self._store.update_novig_order(client_id, {"status": "sent", "novig_order_id": order_id, "error": None})
        return self._follow_order(client_id, order_id, qty, wait=True)

    def _follow_order(self, client_id: str, order_id: str, qty: int, wait: bool) -> dict:
        """Read the order until Novig shows it final (a few tries when `wait`),
        then its fills. An OPEN order (an IOC never should be) is cancelled."""
        attempts = FINAL_READ_ATTEMPTS if wait else 1
        order = None
        for attempt in range(attempts):
            if attempt:
                self._sleep(FINAL_READ_SLEEP_SEC)
            try:
                status, reply = self._client.get_order(order_id)
            except OSError as error:
                log.warning("novig taker order %s: read failed: %s", order_id, error)
                continue
            if status == 200 and isinstance(reply, dict) and reply.get("status") in FINAL_ORDER_STATUSES:
                order = reply
                break
            if status == 200 and isinstance(reply, dict) and reply.get("status") == STATUS_OPEN:
                return self._cancel_resting(client_id, order_id, qty)
        if order is None:
            return {"status": "sent", "qty": qty, "orderId": order_id,
                    "detail": "Novig has not shown the order final yet; full cost held, the refresh reads it"}
        self._store.update_novig_order(client_id, {"final_at": self._now()})
        return self._read_fills(client_id, order_id, qty, final=False)

    def _cancel_resting(self, client_id: str, order_id: str, qty: int) -> dict:
        log.error("novig taker order %s is OPEN on Novig after an IOC send; cancelling", order_id)
        try:
            status, reply = self._client.cancel_order(order_id)
        except OSError as error:
            status, reply = None, f"no answer: {error}"
        self._store.update_novig_order(client_id, {"status": "resting",
                                                   "error": f"OPEN after IOC; cancel answered {status} {str(reply)[:200]}"})
        return {"status": "resting", "qty": qty, "orderId": order_id,
                "detail": "Novig showed the order resting; cancel sent, full cost held until it reads final"}

    def _read_fills(self, client_id: str, order_id: str, qty: int, final: bool) -> dict:
        """Read the order's fills into its row; `final` fixes open_cost to what filled."""
        try:
            filled, cost, fee = order_totals(self._client.fills(order_id=order_id))
        except (NovigError, OSError) as error:
            self._store.update_novig_order(client_id, {"status": "unknown", "error": f"fills read failed: {error}"})
            return {"status": "unknown", "qty": qty, "orderId": order_id, "detail": "fills read failed; full cost held"}
        paid = round(cost + fee, 5)
        outcome = "unfilled" if filled == 0 else "filled" if filled >= qty else "partial"
        fields = {"status": outcome, "qty_filled": filled, "cost_filled": paid, "fills_final": final, "error": None}
        if final:
            fields["open_cost"] = paid
        self._store.update_novig_order(client_id, fields)
        return {"status": outcome, "qty": qty, "filled": filled, "cost": paid, "orderId": order_id}

    def _refresh_open_orders(self) -> None:
        """Bring every open Novig order up to date, fail closed (a failed read
        holds the cost): a lost send is looked up by clientId, a sent or
        resting order read to final, fills read final FILLS_SETTLE_SEC after,
        and a started market checked for settlement, which frees its cost."""
        now = self._now()
        self._last_refresh_at = now
        open_orders = self._store.open_novig_orders()
        open_markets = {order["market_id"] for order in open_orders}
        self._settlement_checked_at = {market: at for market, at in self._settlement_checked_at.items() if market in open_markets}
        for order in open_orders:
            if order["novig_order_id"] is None:
                self._recover_lost_order(order, now)
                continue
            if order["status"] in ("sent", "resting") or order["final_at"] is None:
                self._follow_order(order["client_id"], order["novig_order_id"], order["qty_requested"], wait=False)
                continue
            if not order["fills_final"] and now - order["final_at"] >= timedelta(seconds=FILLS_SETTLE_SEC):
                self._read_fills(order["client_id"], order["novig_order_id"], order["qty_requested"], final=True)
            self._free_if_settled(order, now)

    def _recover_lost_order(self, order: dict, now: datetime) -> None:
        """A send with no definite answer: adopt it from Novig's listings by
        clientId once it is old enough; listed nowhere, it never arrived."""
        if order["status"] == "rejected" or now - order["placed_at"] < timedelta(seconds=LOST_ORDER_LOOKUP_AFTER_SEC):
            return
        try:
            match = self._client.find_order_by_client_id(order["outcome_id"], order["client_id"])
        except (NovigError, OSError) as error:
            log.error("novig taker order %s: lookup failed, holding its cost: %s", order["client_id"], error)
            return
        if match is None:
            self._store.update_novig_order(order["client_id"], {"status": "rejected", "open_cost": 0.0, "settled_at": now,
                                                                "error": "never reached Novig: in none of its order lists"})
            log.warning("novig taker order %s: not in Novig's order lists, counted as never sent", order["client_id"])
            return
        self._store.update_novig_order(order["client_id"], {"status": "sent", "novig_order_id": match["orderId"], "error": None})
        log.warning("novig taker order %s: found on Novig as %s (%s) after a lost answer",
                    order["client_id"], match["orderId"], match.get("status"))
        self._follow_order(order["client_id"], match["orderId"], order["qty_requested"], wait=False)

    def _free_if_settled(self, order: dict, now: datetime) -> None:
        if order["event_start"] > now:
            return
        last = self._settlement_checked_at.get(order["market_id"])
        if last is not None and now - last < timedelta(seconds=SETTLEMENT_CHECK_SEC):
            return
        self._settlement_checked_at[order["market_id"]] = now
        try:
            market = self._client.market(order["market_id"])
        except (NovigError, OSError):
            return
        if market.get("status") == MARKET_SETTLED:
            self._store.update_novig_order(order["client_id"], {"settled_at": now})
