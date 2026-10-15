"""The automated taker's order desk (automated_taker_plan.md): the only place a
bot order is signed and sent. The server runner (server/runner.js) decides;
this module checks the stop file and the $500 ceiling, sends one limit order
at the ask the decision was checked at, cancels whatever did not fill at
once, and reads what filled.

Armed only when a dedicated Kalshi key is configured
(KALSHI_TAKER_API_KEY_ID + KALSHI_TAKER_PRIVATE_KEY_PATH, a key made for the
bot alone so it can be revoked alone). Without it, POST /taker/order.json
answers 503 and nothing is sent.

Rules, in order, for every request:
  1. the stop file (config.TAKER_STOP_PATH) exists -> refused, "stopped";
  2. the same client_order_id was already sent -> refused, "duplicate";
  3. open orders are brought up to date when the last refresh is older than
     REFRESH_SEC (the desk's own thread does the same every REFRESH_SEC, so
     nothing below waits on the runner): see _refresh_open_orders;
  4. unseen fill: an earlier bot order on the SAME GAME that may still hold
     contracts must be among the held bets the runner sized this stake on
     (`heldIds`, the ids bets.js matched to the row: "kalshi:<ticker>:<side>"),
     else refused, "unseen_fill". Conditional Kelly is only right once it
     sees the earlier fill; until then an add would be sized as if alone;
  5. ceiling: contracts are cut so that open exposure + this order's cost
     stays at or under BOT_BANKROLL; none left -> refused, "ceiling";
  6. the order row is INSERTED (status "sending", full cost) BEFORE it is
     sent, so a crash mid-send still counts against the ceiling;
  7. one limit buy at the ask (never a market order), then a cancel of the
     remainder routed by market (a bare DELETE goes to exchange shard 0 and
     404s while the order keeps resting: kalshi_rfi, 2026-08-27), a 404
     verified against Kalshi's resting list, then a fills read.

Order statuses (taker.duckdb::taker_orders.status):
  sending   inserted, the POST not yet answered
  unknown   the POST got no answer (it may have landed; recovered from
            Kalshi's order list by client_order_id), or a fills read failed
  sent      Kalshi took it; the cancel and fills read follow at once
  resting   the cancel could not be confirmed, so it may still rest: the
            full cost stays held and every refresh retries the cancel
  filled | partial | unfilled   what the fills read found; `fills_final`
            once read FILLS_SETTLE_SEC after the confirmed cancel
  rejected  Kalshi refused it, or it never reached Kalshi: nothing held

Side effects: Kalshi orders on the bot's key (money); INSERTs / UPDATEs
taker.duckdb::taker_orders (taker_store.py). Network: Kalshi REST, signed
with the bot's key only.
"""
import base64
import json
import logging
import math
import threading
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from datetime import datetime, timedelta, timezone
from pathlib import Path

from unabated_ticket.bets_service.taker_store import TakerStore

log = logging.getLogger(__name__)

# Cal, 2026-10-14: $500, and the same $500 caps what the bot has open.
BOT_BANKROLL = 500.0
# Kalshi's fill feed can trail the match by a moment; a read this long after
# the confirmed cancel is taken as final. Until then the whole requested cost is held.
FILLS_SETTLE_SEC = 60
# Kalshi rounds each order's fee up to the cent: one cent of slack per order.
FEE_ROUNDING_SLACK = 0.01
SETTLED_MARKET_STATUSES = ("settled", "finalized")
REQUEST_TIMEOUT_SEC = 10
SIGNED_PATH_PREFIX = "/trade-api/v2"
ORDERS_PATH = "/portfolio/events/orders"
ORDER_LIST_PATH = "/portfolio/orders"
MIN_PRICE = 0.01
MAX_PRICE = 0.99
# Open orders are brought up to date this often: by the desk's thread, and
# before a place when the last refresh is older.
REFRESH_SEC = 30
# A started market is asked for settlement this often, per ticker.
SETTLEMENT_CHECK_SEC = 5 * 60
# A send with no answer and not resting: Kalshi's order list for the ticker
# is asked for it by client_order_id once it is this old; absent then, it
# never reached Kalshi and its cost is freed.
LOST_ORDER_LOOKUP_AFTER_SEC = 10 * 60
ORDER_LIST_PAGE_LIMIT = 1000
ORDER_LIST_MAX_PAGES = 10
MAX_HELD_IDS = 1000
MAX_HELD_ID_LENGTH = 200
BOT_ORDER_ID_PREFIX = "ubot-"
# Statuses under which the order may hold contracts (see _may_hold_contracts).
STATUSES_MAYBE_HOLDING = ("sending", "unknown", "sent", "resting", "filled", "partial")


class KalshiOrderClient:
    """Signed Kalshi REST on ONE key (RSA-PSS SHA-256, Kalshi's scheme). Its
    own key, not kalshi_common.auth_client's process-wide one, so the bets
    service's read key and the bot's trading key never mix."""

    def __init__(self, api_key_id: str, private_key_path: str, base_url: str):
        from cryptography.hazmat.primitives import serialization
        self._api_key_id = api_key_id
        self._base_url = base_url.rstrip("/")
        with open(private_key_path, "rb") as key_file:
            self._private_key = serialization.load_pem_private_key(key_file.read(), password=None)

    def _signature(self, timestamp_ms: str, method: str, path: str) -> str:
        from cryptography.hazmat.primitives import hashes
        from cryptography.hazmat.primitives.asymmetric import padding
        # Kalshi signs the path without its query (kalshi_common/auth_client.py does the same, live).
        message = f"{timestamp_ms}{method}{SIGNED_PATH_PREFIX}{path.split('?', 1)[0]}".encode()
        signature = self._private_key.sign(
            message, padding.PSS(mgf=padding.MGF1(hashes.SHA256()), salt_length=padding.PSS.DIGEST_LENGTH),
            hashes.SHA256())
        return base64.b64encode(signature).decode()

    def request(self, method: str, path: str, body: dict | None = None) -> tuple[int, object]:
        """(status, parsed body or text). Never retried: a POST that timed out
        may still have reached Kalshi. Raises OSError on no answer at all."""
        timestamp_ms = str(int(time.time() * 1000))
        data = json.dumps(body).encode() if body is not None else None
        request = urllib.request.Request(f"{self._base_url}{path}", data=data, method=method)
        request.add_header("KALSHI-ACCESS-KEY", self._api_key_id)
        request.add_header("KALSHI-ACCESS-SIGNATURE", self._signature(timestamp_ms, method, path))
        request.add_header("KALSHI-ACCESS-TIMESTAMP", timestamp_ms)
        request.add_header("Content-Type", "application/json")
        try:
            with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT_SEC) as response:
                return response.status, _parse(response.read())
        except urllib.error.HTTPError as error:
            return error.code, _parse(error.read())


def _parse(raw: bytes) -> object:
    text = raw.decode(errors="replace")
    if not text:
        return {}
    try:
        return json.loads(text)
    except ValueError:
        return text


def order_body(ticker: str, contract_side: str, ask: float, count: int, client_order_id: str) -> dict:
    """Kalshi v2's order body for BUYING `contract_side` at `ask`. v2 trades the
    YES leg only: buy YES = bid at the ask; buy NO = ask (sell YES) at 1 - the
    NO ask (unabated_edge/maker/gateway.py, the same API live)."""
    if contract_side == "yes":
        v2_side, price = "bid", ask
    else:
        v2_side, price = "ask", 1 - ask
    return {"ticker": ticker, "side": v2_side, "price": f"{price:.4f}", "count": f"{float(count):.2f}",
            "time_in_force": "good_till_canceled", "self_trade_prevention_type": "taker_at_cross",
            "client_order_id": client_order_id}


def fill_totals(fills: list[dict], contract_side: str) -> tuple[float, float]:
    """(contracts, dollars paid incl. fees) over one order's fills, priced on
    the bot's own contract side whichever leg Kalshi reports the fill on."""
    contracts = 0.0
    cost = 0.0
    for fill in fills:
        count = float(fill.get("count_fp") or fill.get("count") or 0)
        yes_price = float(fill.get("yes_price_dollars") or 0)
        no_price = float(fill.get("no_price_dollars") or (1 - yes_price))
        price = yes_price if contract_side == "yes" else no_price
        contracts += count
        cost += count * price + float(fill.get("fee_cost") or 0)
    return contracts, round(cost, 4)


def held_id_of(ticker: str, contract_side: str) -> str:
    """The bets service's record id for a Kalshi position (sources/kalshi.py)."""
    return f"kalshi:{ticker}:{contract_side}"


def _order_id_of(reply: object) -> str | None:
    if not isinstance(reply, dict):
        return None
    order = reply.get("order") if isinstance(reply.get("order"), dict) else reply
    order_id = order.get("order_id")
    return order_id if isinstance(order_id, str) and order_id else None


def _exchange_index_of(reply: object) -> int | None:
    if not isinstance(reply, dict):
        return None
    order = reply.get("order") if isinstance(reply.get("order"), dict) else reply
    index = order.get("exchange_index")
    return int(index) if isinstance(index, (int, float)) and not isinstance(index, bool) else None


def validate_order_request(body: object) -> dict | str:
    """A POST /taker/order.json body, or what was expected and what was found."""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    client_order_id = body.get("clientOrderId")
    if (not isinstance(client_order_id, str) or not client_order_id.startswith(BOT_ORDER_ID_PREFIX)
            or len(client_order_id) > 64):
        return f"clientOrderId must be a '{BOT_ORDER_ID_PREFIX}' id of at most 64 characters, got {client_order_id!r}"
    for field in ("ticker", "rowKey"):
        if not isinstance(body.get(field), str) or not body[field]:
            return f"{field} must be a non-empty string, got {body.get(field)!r}"
    if body.get("side") not in ("yes", "no"):
        return f"side must be 'yes' or 'no', got {body.get('side')!r}"
    ask = body.get("ask")
    if not isinstance(ask, (int, float)) or isinstance(ask, bool) or not MIN_PRICE <= ask <= MAX_PRICE:
        return f"ask must be a price from {MIN_PRICE} to {MAX_PRICE}, got {ask!r}"
    all_in = body.get("allIn")
    if not isinstance(all_in, (int, float)) or isinstance(all_in, bool) or not ask <= all_in < 1:
        return f"allIn must be at least the ask and under 1, got {all_in!r}"
    count = body.get("count")
    if not isinstance(count, int) or isinstance(count, bool) or count < 1:
        return f"count must be a whole number of contracts >= 1, got {count!r}"
    event_id = body.get("eventId")
    if not isinstance(event_id, int) or isinstance(event_id, bool) or event_id < 1:
        return f"eventId must be Unabated's event id, got {event_id!r}"
    held_ids = body.get("heldIds")
    if (not isinstance(held_ids, list) or len(held_ids) > MAX_HELD_IDS
            or any(not isinstance(held, str) or len(held) > MAX_HELD_ID_LENGTH for held in held_ids)):
        return f"heldIds must be a list of at most {MAX_HELD_IDS} bet ids the stake was sized on, got {held_ids!r}"
    event_start = body.get("eventStart")
    try:
        start = datetime.fromisoformat(event_start.replace("Z", "+00:00")) if isinstance(event_start, str) else None
    except ValueError:
        start = None
    if start is None:
        return f"eventStart must be an ISO time, got {event_start!r}"
    return {"client_order_id": client_order_id, "ticker": body["ticker"], "row_key": body["rowKey"],
            "contract_side": body["side"], "ask_price": float(ask), "all_in_price": float(all_in), "count": count,
            "event_id": event_id, "held_ids": set(held_ids), "event_start": start}


def _may_hold_contracts(order: dict) -> bool:
    """Could this open order hold contracts? Yes until a FINAL fills read says
    none did: a fill can trail the cancel in Kalshi's feed."""
    if order["status"] in STATUSES_MAYBE_HOLDING:
        return True
    return order["status"] == "unfilled" and not order["fills_final"]


class TakerDesk:
    """Places the runner's orders, one at a time, under the ceiling."""

    def __init__(self, client: KalshiOrderClient, store: TakerStore, stop_path: Path,
                 now: Callable[[], datetime] = lambda: datetime.now(timezone.utc)):
        self._client = client
        self._store = store
        self._stop_path = Path(stop_path)
        self._now = now
        self._lock = threading.Lock()
        self._last_refresh_at: datetime | None = None
        self._settlement_checked_at: dict[str, datetime] = {}

    def status(self) -> dict:
        return {"armed": True, "stopped": self._stop_path.exists(), "bankroll": BOT_BANKROLL,
                "openExposure": round(self._store.open_exposure(), 2),
                "lastRefreshAt": self._last_refresh_at.isoformat() if self._last_refresh_at else None}

    def run_refresh_loop(self, stop: threading.Event) -> None:
        """The desk's own thread (service.serve): open orders brought up to
        date every REFRESH_SEC whether or not the runner sends anything, so a
        resting order is cancelled and settled money freed on time."""
        while not stop.wait(REFRESH_SEC):
            try:
                self.refresh()
            except Exception:  # noqa: BLE001 — the sweep must outlive one bad answer
                log.exception("taker desk refresh failed")

    def refresh(self) -> None:
        with self._lock:
            self._refresh_open_orders()

    def place(self, request: dict) -> dict:
        """One order (validate_order_request's shape). Returns {status, ...}:
        stopped | duplicate | unseen_fill | ceiling (nothing sent), or the
        sent order's filled | partial | unfilled | resting | rejected | unknown."""
        with self._lock:
            if self._stop_path.exists():
                return {"status": "stopped", "detail": f"stop file {self._stop_path} is present"}
            if self._store.has_order(request["client_order_id"]):
                return {"status": "duplicate", "detail": f"{request['client_order_id']} was already sent"}
            if self._last_refresh_at is None or self._now() - self._last_refresh_at >= timedelta(seconds=REFRESH_SEC):
                self._refresh_open_orders()
            unseen = self._unseen_fill(request)
            if unseen:
                return {"status": "unseen_fill", "detail": unseen}
            headroom = BOT_BANKROLL - self._store.open_exposure()
            count = min(request["count"], math.floor((headroom - FEE_ROUNDING_SLACK) / request["all_in_price"]))
            if count < 1:
                return {"status": "ceiling", "detail": f"${headroom:.2f} of ${BOT_BANKROLL:.0f} left"}
            return self._send(request, count)

    def _unseen_fill(self, request: dict) -> str | None:
        """Why this order would add to a game blind, or None: an earlier bot
        order on the game that may hold contracts and is not in `held_ids`."""
        for order in self._store.open_orders():
            if order["event_id"] != request["event_id"] or not _may_hold_contracts(order):
                continue
            held_id = held_id_of(order["ticker"], order["contract_side"])
            if held_id not in request["held_ids"]:
                return (f"{order['client_order_id']} ({order['status']}) is not yet among the held bets the stake "
                        f"was sized on: no {held_id}")
        return None

    def _send(self, request: dict, count: int) -> dict:
        client_order_id = request["client_order_id"]
        reserved = round(count * request["all_in_price"] + FEE_ROUNDING_SLACK, 4)
        self._store.insert_order({**request, "placed_at": self._now(), "count_requested": count, "open_cost": reserved})
        body = order_body(request["ticker"], request["contract_side"], request["ask_price"], count, client_order_id)
        try:
            status, reply = self._client.request("POST", ORDERS_PATH, body)
        except OSError as error:
            # It may have reached Kalshi: keep the full cost held; the refresh looks it up by client_order_id.
            self._store.update_order(client_order_id, {"status": "unknown", "error": f"no answer: {error}"})
            log.error("taker order %s %s: no answer from Kalshi: %s", client_order_id, request["ticker"], error)
            return {"status": "unknown", "count": count, "detail": f"no answer from Kalshi: {error}"}
        order_id = _order_id_of(reply)
        if status not in (200, 201) or not order_id:
            # Nothing rests or filled: no cost, and no market to watch for settlement.
            self._store.update_order(client_order_id, {"status": "rejected", "open_cost": 0.0, "settled_at": self._now(),
                                                       "error": f"HTTP {status}: {str(reply)[:500]}"})
            log.warning("taker order %s %s refused by Kalshi: HTTP %s %s", client_order_id, request["ticker"], status, reply)
            return {"status": "rejected", "count": count, "detail": f"Kalshi HTTP {status}"}
        exchange_index = _exchange_index_of(reply)
        self._store.update_order(client_order_id, {"status": "sent", "kalshi_order_id": order_id, "exchange_index": exchange_index})
        order = {"client_order_id": client_order_id, "kalshi_order_id": order_id, "ticker": request["ticker"],
                 "contract_side": request["contract_side"], "exchange_index": exchange_index, "count_requested": count}
        return self._cancel_then_read(order, resting=None)

    def _cancel_then_read(self, order: dict, resting: list[dict] | None) -> dict:
        """Cancel the rest of a sent order and read its fills. An unconfirmed
        cancel leaves it "resting" with the full cost held; the refresh retries."""
        client_order_id = order["client_order_id"]
        if not self._cancel_rest(order, resting):
            self._store.update_order(client_order_id, {"status": "resting", "error": "cancel unconfirmed; the refresh retries it"})
            return {"status": "resting", "count": order["count_requested"], "orderId": order["kalshi_order_id"],
                    "detail": "the cancel could not be confirmed; full cost held, cancel retried"}
        self._store.update_order(client_order_id, {"cancelled_at": self._now(), "error": None})
        return self._read_fills(client_order_id, order["kalshi_order_id"], order["contract_side"],
                                order["count_requested"], final=False)

    def _cancel_rest(self, order: dict, resting: list[dict] | None) -> bool:
        """Cancel what did not fill, routed by market (and shard when Kalshi
        named one). True only when Kalshi confirmed it, or answered 404 and
        its resting list does not carry the order."""
        order_id = order["kalshi_order_id"]
        path = f"{ORDERS_PATH}/{order_id}?market_ticker={order['ticker']}"
        if order.get("exchange_index") is not None:
            path += f"&exchange_index={int(order['exchange_index'])}"
        try:
            status, reply = self._client.request("DELETE", path)
        except OSError as error:
            log.error("taker cancel %s: no answer: %s", order_id, error)
            return False
        if status in (200, 204):
            return True
        if status == 404:
            # Terminal (filled, expired) or misrouted: only Kalshi's resting list tells.
            listing = self._resting_orders() if resting is None else resting
            if listing is None:
                log.error("taker cancel %s: 404 and the resting list could not be read; holding", order_id)
                return False
            still_resting = any(resting_order.get("order_id") == order_id for resting_order in listing)
            if still_resting:
                log.error("taker cancel %s: 404 but still resting on Kalshi (exchange_index=%s); holding",
                          order_id, order.get("exchange_index"))
            return not still_resting
        log.error("taker cancel %s failed: HTTP %s %s", order_id, status, reply)
        return False

    def _read_fills(self, client_order_id: str, order_id: str, contract_side: str, count: int, final: bool) -> dict:
        try:
            status, reply = self._client.request("GET", f"/portfolio/fills?order_id={order_id}")
        except OSError as error:
            status, reply = None, str(error)
        if status != 200 or not isinstance(reply, dict) or not isinstance(reply.get("fills"), list):
            self._store.update_order(client_order_id, {"status": "unknown", "error": f"fills read failed: {status} {str(reply)[:300]}"})
            return {"status": "unknown", "count": count, "orderId": order_id, "detail": "fills read failed; full cost held"}
        filled, cost = fill_totals(reply["fills"], contract_side)
        outcome = "unfilled" if filled == 0 else "filled" if filled >= count else "partial"
        fields = {"status": outcome, "count_filled": filled, "cost_filled": cost, "fills_final": final, "error": None}
        if final:
            fields["open_cost"] = cost
        self._store.update_order(client_order_id, fields)
        return {"status": outcome, "count": count, "filled": filled, "cost": cost, "orderId": order_id}

    def _list_orders(self, query: str) -> list[dict] | None:
        """GET /portfolio/orders?<query>, every page; None when the read
        failed or the listing did not end (never a partial list)."""
        orders: list[dict] = []
        cursor = ""
        for _ in range(ORDER_LIST_MAX_PAGES):
            path = f"{ORDER_LIST_PATH}?{query}&limit={ORDER_LIST_PAGE_LIMIT}" + (f"&cursor={cursor}" if cursor else "")
            try:
                status, reply = self._client.request("GET", path)
            except OSError as error:
                log.error("taker order list %s: no answer: %s", query, error)
                return None
            if status != 200 or not isinstance(reply, dict) or not isinstance(reply.get("orders"), list):
                log.error("taker order list %s failed: HTTP %s %s", query, status, str(reply)[:300])
                return None
            orders.extend(order for order in reply["orders"] if isinstance(order, dict))
            cursor = reply.get("cursor") or ""
            if not cursor:
                return orders
        log.error("taker order list %s: more than %d pages; holding", query, ORDER_LIST_MAX_PAGES)
        return None

    def _resting_orders(self) -> list[dict] | None:
        return self._list_orders("status=resting")

    def _refresh_open_orders(self) -> None:
        """Bring every open order up to date, fail closed (a failed read holds
        the cost): a bot order resting on Kalshi that no row knows is
        cancelled; a send with no answer is looked up by client_order_id; an
        unconfirmed cancel is retried; fills are final FILLS_SETTLE_SEC after
        the cancel; a started market is checked for settlement, which frees it."""
        now = self._now()
        self._last_refresh_at = now
        open_orders = self._store.open_orders()
        open_tickers = {order["ticker"] for order in open_orders}
        self._settlement_checked_at = {ticker: at for ticker, at in self._settlement_checked_at.items() if ticker in open_tickers}
        if not open_orders:
            return
        resting = self._resting_orders()
        if resting is not None:
            self._cancel_untracked(resting, {order["client_order_id"] for order in open_orders})
        for order in open_orders:
            if not order["kalshi_order_id"]:
                self._recover_lost_order(order, resting, now)
                continue
            if order["status"] in ("sent", "resting"):
                self._cancel_then_read(order, resting)
                continue
            cancelled_at = order["cancelled_at"]
            if (not order["fills_final"] and cancelled_at is not None
                    and now - cancelled_at >= timedelta(seconds=FILLS_SETTLE_SEC)):
                self._read_fills(order["client_order_id"], order["kalshi_order_id"], order["contract_side"],
                                 order["count_requested"], final=True)
            if (order["event_start"] is not None and order["event_start"] <= now
                    and self._settlement_check_due(order["ticker"], now) and self._market_settled(order["ticker"])):
                self._store.update_order(order["client_order_id"], {"settled_at": now})

    def _cancel_untracked(self, resting: list[dict], known_ids: set[str]) -> None:
        """A bot order (client_order_id "ubot-…") resting on Kalshi that no
        open row names (a store lost or reset) is the bot's to cancel."""
        for resting_order in resting:
            client_order_id = resting_order.get("client_order_id")
            if not isinstance(client_order_id, str) or not client_order_id.startswith(BOT_ORDER_ID_PREFIX):
                continue
            if client_order_id in known_ids or not resting_order.get("order_id"):
                continue
            log.error("taker: untracked bot order %s resting on %s; cancelling", client_order_id, resting_order.get("ticker"))
            self._cancel_rest({"kalshi_order_id": resting_order["order_id"], "ticker": resting_order.get("ticker") or "",
                               "exchange_index": resting_order.get("exchange_index")}, resting)

    def _recover_lost_order(self, order: dict, resting: list[dict] | None, now: datetime) -> None:
        """An order whose POST got no answer: adopt it from Kalshi's resting
        list, else from the ticker's order list once it is old enough, else
        (listed nowhere) it never arrived and holds nothing."""
        client_order_id = order["client_order_id"]
        match = _by_client_order_id(resting, client_order_id)
        if match is None and now - order["placed_at"] >= timedelta(seconds=LOST_ORDER_LOOKUP_AFTER_SEC):
            listed = self._list_orders(f"ticker={order['ticker']}")
            if listed is None:
                return
            match = _by_client_order_id(listed, client_order_id)
            if match is None:
                self._store.update_order(client_order_id, {"status": "rejected", "open_cost": 0.0, "settled_at": now,
                                                           "error": "never reached Kalshi: not in its order list"})
                log.warning("taker order %s %s: not in Kalshi's order list, counted as never sent", client_order_id, order["ticker"])
                return
        if match is None:
            return
        adopted = {**order, "kalshi_order_id": match["order_id"], "exchange_index": _exchange_index_of(match)}
        self._store.update_order(client_order_id, {"status": "sent", "kalshi_order_id": match["order_id"],
                                                   "exchange_index": adopted["exchange_index"]})
        log.warning("taker order %s %s: found on Kalshi as %s after a lost answer", client_order_id, order["ticker"], match["order_id"])
        self._cancel_then_read(adopted, resting)

    def _settlement_check_due(self, ticker: str, now: datetime) -> bool:
        last = self._settlement_checked_at.get(ticker)
        if last is not None and now - last < timedelta(seconds=SETTLEMENT_CHECK_SEC):
            return False
        self._settlement_checked_at[ticker] = now
        return True

    def _market_settled(self, ticker: str) -> bool:
        try:
            status, reply = self._client.request("GET", f"/markets/{ticker}")
        except OSError:
            return False
        market = reply.get("market") if isinstance(reply, dict) else None
        return status == 200 and isinstance(market, dict) and market.get("status") in SETTLED_MARKET_STATUSES


def _by_client_order_id(orders: list[dict] | None, client_order_id: str) -> dict | None:
    for order in orders or []:
        if order.get("client_order_id") == client_order_id and order.get("order_id"):
            return order
    return None
