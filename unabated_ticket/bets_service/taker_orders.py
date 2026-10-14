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
  3. open orders are brought up to date: fills re-read once the order is
     FILLS_SETTLE_SEC old (until then the whole requested cost counts), and a
     market whose game has started is checked for settlement, which frees it;
  4. ceiling: contracts are cut so that open exposure + this order's cost
     stays at or under BOT_BANKROLL; none left -> refused, "ceiling";
  5. the order row is INSERTED (status "sending", full cost) BEFORE it is
     sent, so a crash mid-send still counts against the ceiling;
  6. one limit buy at the ask (never a market order), then an immediate
     cancel of the remainder, then a fills read.

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
# placing is taken as final. Until then the whole requested cost is held.
FILLS_SETTLE_SEC = 60
# Kalshi rounds each order's fee up to the cent: one cent of slack per order.
FEE_ROUNDING_SLACK = 0.01
SETTLED_MARKET_STATUSES = ("settled", "finalized")
REQUEST_TIMEOUT_SEC = 10
SIGNED_PATH_PREFIX = "/trade-api/v2"
ORDERS_PATH = "/portfolio/events/orders"
MIN_PRICE = 0.01
MAX_PRICE = 0.99


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


def validate_order_request(body: object) -> dict | str:
    """A POST /taker/order.json body, or what was expected and what was found."""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    client_order_id = body.get("clientOrderId")
    if not isinstance(client_order_id, str) or not client_order_id.startswith("ubot-") or len(client_order_id) > 64:
        return f"clientOrderId must be a 'ubot-' id of at most 64 characters, got {client_order_id!r}"
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
    event_start = body.get("eventStart")
    try:
        start = datetime.fromisoformat(event_start.replace("Z", "+00:00")) if isinstance(event_start, str) else None
    except ValueError:
        start = None
    if start is None:
        return f"eventStart must be an ISO time, got {event_start!r}"
    return {"client_order_id": client_order_id, "ticker": body["ticker"], "row_key": body["rowKey"],
            "contract_side": body["side"], "ask_price": float(ask), "all_in_price": float(all_in), "count": count,
            "event_start": start}


class TakerDesk:
    """Places the runner's orders, one at a time, under the ceiling."""

    def __init__(self, client: KalshiOrderClient, store: TakerStore, stop_path: Path,
                 now: Callable[[], datetime] = lambda: datetime.now(timezone.utc)):
        self._client = client
        self._store = store
        self._stop_path = Path(stop_path)
        self._now = now
        self._lock = threading.Lock()

    def status(self) -> dict:
        return {"armed": True, "stopped": self._stop_path.exists(), "bankroll": BOT_BANKROLL,
                "openExposure": round(self._store.open_exposure(), 2)}

    def place(self, request: dict) -> dict:
        """One order (validate_order_request's shape). Returns {status, ...}:
        stopped | duplicate | ceiling (nothing sent), or the sent order's
        filled | partial | unfilled | rejected | unknown."""
        with self._lock:
            if self._stop_path.exists():
                return {"status": "stopped", "detail": f"stop file {self._stop_path} is present"}
            if self._store.has_order(request["client_order_id"]):
                return {"status": "duplicate", "detail": f"{request['client_order_id']} was already sent"}
            self._refresh_open_orders()
            headroom = BOT_BANKROLL - self._store.open_exposure()
            count = min(request["count"], math.floor((headroom - FEE_ROUNDING_SLACK) / request["all_in_price"]))
            if count < 1:
                return {"status": "ceiling", "detail": f"${headroom:.2f} of ${BOT_BANKROLL:.0f} left"}
            return self._send(request, count)

    def _send(self, request: dict, count: int) -> dict:
        client_order_id = request["client_order_id"]
        reserved = round(count * request["all_in_price"] + FEE_ROUNDING_SLACK, 4)
        self._store.insert_order({**request, "placed_at": self._now(), "count_requested": count, "open_cost": reserved})
        body = order_body(request["ticker"], request["contract_side"], request["ask_price"], count, client_order_id)
        try:
            status, reply = self._client.request("POST", ORDERS_PATH, body)
        except OSError as error:
            # It may have reached Kalshi: keep the full cost held, the fills read settles it later.
            self._store.update_order(client_order_id, {"status": "unknown", "error": f"no answer: {error}"})
            log.error("taker order %s %s: no answer from Kalshi: %s", client_order_id, request["ticker"], error)
            return {"status": "unknown", "count": count, "detail": f"no answer from Kalshi: {error}"}
        order_id = (reply.get("order_id") or (reply.get("order") or {}).get("order_id")) if isinstance(reply, dict) else None
        if status not in (200, 201) or not order_id:
            # Nothing rests or filled: no cost, and no market to watch for settlement.
            self._store.update_order(client_order_id, {"status": "rejected", "open_cost": 0.0, "settled_at": self._now(),
                                                       "error": f"HTTP {status}: {str(reply)[:500]}"})
            log.warning("taker order %s %s refused by Kalshi: HTTP %s %s", client_order_id, request["ticker"], status, reply)
            return {"status": "rejected", "count": count, "detail": f"Kalshi HTTP {status}"}
        self._store.update_order(client_order_id, {"status": "sent", "kalshi_order_id": order_id})
        self._cancel_rest(order_id)
        return self._read_fills(client_order_id, order_id, request["contract_side"], count, final=False)

    def _cancel_rest(self, order_id: str) -> None:
        """Cancel what did not fill at once; 404 means nothing was left."""
        try:
            status, reply = self._client.request("DELETE", f"{ORDERS_PATH}/{order_id}")
        except OSError as error:
            log.error("taker cancel %s: no answer: %s", order_id, error)
            return
        if status not in (200, 204, 404):
            log.error("taker cancel %s failed: HTTP %s %s", order_id, status, reply)

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

    def _refresh_open_orders(self) -> None:
        """Final fills for orders old enough; settlement for games under way.
        A failed read changes nothing, so the cost stays held (fail closed)."""
        now = self._now()
        for order in self._store.open_orders():
            placed = order["placed_at"]
            if (not order["fills_final"] and order["kalshi_order_id"]
                    and now - placed >= timedelta(seconds=FILLS_SETTLE_SEC)):
                self._read_fills(order["client_order_id"], order["kalshi_order_id"], order["contract_side"],
                                 order["count_requested"], final=True)
            if order["event_start"] is not None and order["event_start"] <= now and self._market_settled(order["ticker"]):
                self._store.update_order(order["client_order_id"], {"settled_at": now})

    def _market_settled(self, ticker: str) -> bool:
        try:
            status, reply = self._client.request("GET", f"/markets/{ticker}")
        except OSError:
            return False
        market = reply.get("market") if isinstance(reply, dict) else None
        return status == 200 and isinstance(market, dict) and market.get("status") in SETTLED_MARKET_STATUSES
