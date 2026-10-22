"""Signed Polymarket US REST on ONE trading key: the order client the Ticket's
Polymarket US Place button sends through (polymarket_place.py). Pure transport
plus Polymarket US's own order rules; it decides nothing and keeps no state.

Polymarket US facts it rests on (the official SDK, github.com/Polymarket/
polymarket-us-python, read 2026-10-21; the docs site is unreachable from the
cloud container, so nothing here was read off docs.polymarket.us directly):
  * Signing, as sources/polymarket_us.py reads the account: headers
    X-PM-Access-Key, X-PM-Timestamp (unix ms, +-30 s) and X-PM-Signature =
    standard base64 of the Ed25519 signature over timestamp + METHOD + path,
    the path WITHOUT its query string (a signed query string is refused 401,
    live 2026-09-23). The body is not signed.
  * The secret key is base64 of the 32-byte Ed25519 seed followed by the
    public key (the "Secret Key" polymarket.us/developer shows once); it is
    read from a FILE here (the Place key moves to the VM by scp), never .env.
  * One instrument per market: the YES (long) contract. Buying NO is
    ORDER_INTENT_BUY_SHORT, and `price.value` is ALWAYS the YES price, on a
    short order too: buying NO at 83c sends 0.17 (yes_price_of).
  * POST /v1/orders {marketSlug, intent, type, price, quantity, tif,
    synchronousExecution} answers {id, executions}. There is NO client order
    id: a re-sent order is a second order. Nothing here ever re-sends; the
    desk's bets.duckdb::manual_orders row is the only guard.
  * GET /v1/order/<id> answers {order: {state, cumQuantity, leavesQuantity,
    avgPx, commissionNotionalTotalCollected, ...}}; POST
    /v1/order/<id>/cancel {marketSlug} cancels one.
  * Public reads (no key) live on the gateway: GET /v1/markets/<slug>/book
    {marketData: {bids, offers, state}}, bids/offers of the YES contract.

Side effects: Polymarket US orders and cancels on the configured key (money)
through place_ioc_order / cancel_order; everything else is a read. Nothing on disk.
"""
import base64
import json
import math
import time
from collections.abc import Callable
from pathlib import Path

import requests
from cryptography.hazmat.primitives.asymmetric import ed25519

from unabated_ticket.bets_service.sources.polymarket_us import API_BASE_URL, ED25519_SEED_BYTES, GATEWAY_BASE_URL

ACCESS_KEY_HEADER = "X-PM-Access-Key"
TIMESTAMP_HEADER = "X-PM-Timestamp"
SIGNATURE_HEADER = "X-PM-Signature"
REQUEST_TIMEOUT_SEC = 10

ORDERS_PATH = "/v1/orders"
ORDER_PATH = "/v1/order"
BALANCES_PATH = "/v1/account/balances"
INTENT_BUY_LONG = "ORDER_INTENT_BUY_LONG"
INTENT_BUY_SHORT = "ORDER_INTENT_BUY_SHORT"
# The contract side Place buys -> the order intent that buys it.
INTENT_OF_SIDE = {"yes": INTENT_BUY_LONG, "no": INTENT_BUY_SHORT}
ORDER_TYPE_LIMIT = "ORDER_TYPE_LIMIT"
TIF_IOC = "TIME_IN_FORCE_IMMEDIATE_OR_CANCEL"
CURRENCY = "USD"
# Order states after which nothing more fills (types/orders.py OrderState).
FINAL_ORDER_STATES = {"ORDER_STATE_FILLED", "ORDER_STATE_CANCELED", "ORDER_STATE_REJECTED",
                      "ORDER_STATE_EXPIRED", "ORDER_STATE_REPLACED"}
# Resting: an IOC never should be.
RESTING_ORDER_STATES = {"ORDER_STATE_NEW", "ORDER_STATE_PARTIALLY_FILLED"}
MARKET_STATE_OPEN = "MARKET_STATE_OPEN"
# A price is a whole number of tenths of a cent on Polymarket US's widest
# documented grid; the desk sends the book's own level, so this only formats it.
PRICE_DECIMALS = 4
MIN_PRICE = 0.01
MAX_PRICE = 0.99

Transport = Callable[[str, str, dict, bytes | None], tuple[int, bytes]]


class PolymarketError(RuntimeError):
    """A Polymarket US answer that is not what the call expects."""


def signer_of_key_file(path: str | Path) -> ed25519.Ed25519PrivateKey:
    """The Ed25519 key in a Place key file: the base64 "Secret Key" text."""
    text = Path(path).read_text().strip()
    try:
        seed = base64.b64decode(text, validate=True)[:ED25519_SEED_BYTES]
        return ed25519.Ed25519PrivateKey.from_private_bytes(seed)
    except ValueError as error:
        raise PolymarketError(f"{path} does not hold a base64 Ed25519 Secret Key from polymarket.us/developer: "
                              f"{error}") from None


def yes_price_of(contract_side: str, side_price: float) -> float:
    """The YES price an order on `contract_side` at `side_price` sends: a NO
    buy at 0.83 is a short at YES 0.17."""
    if contract_side == "yes":
        return side_price
    if contract_side == "no":
        return round(1 - side_price, PRICE_DECIMALS)
    raise ValueError(f"contract_side must be 'yes' or 'no', got {contract_side!r}")


def price_string(price: float) -> str:
    return f"{price:.{PRICE_DECIMALS}f}"


def ioc_order_body(slug: str, contract_side: str, side_price: float, quantity: int) -> dict:
    """ONE immediate-or-cancel limit buy of `quantity` contracts on
    `contract_side` at no worse than `side_price` (that side's own price)."""
    if not isinstance(quantity, int) or quantity < 1:
        raise ValueError(f"quantity must be a whole number of contracts >= 1, got {quantity!r}")
    if not MIN_PRICE <= side_price <= MAX_PRICE:
        raise ValueError(f"price must be {MIN_PRICE}-{MAX_PRICE}, got {side_price!r}")
    return {
        "marketSlug": slug, "intent": INTENT_OF_SIDE[contract_side], "type": ORDER_TYPE_LIMIT,
        "price": {"value": price_string(yes_price_of(contract_side, side_price)), "currency": CURRENCY},
        "quantity": quantity, "tif": TIF_IOC, "synchronousExecution": True,
    }


def amount_of(amount: object) -> float | None:
    """{value: "0.17", currency: "USD"} (or a bare number/string) as a float, or None."""
    value = amount.get("value") if isinstance(amount, dict) else amount
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def quantity_of(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) and number >= 0 else None


def order_fill_totals(order: dict, contract_side: str) -> tuple[float, float | None] | str:
    """(contracts filled, dollars paid incl. commission or None when the
    order does not say) off one order object, or what was missing. avgPx is
    the YES price like every price here, so a NO fill costs 1 - avgPx."""
    filled = quantity_of(order.get("cumQuantity"))
    if filled is None:
        return f"no cumQuantity in the order: {str(order)[:200]}"
    if filled == 0:
        return 0.0, 0.0
    average = amount_of(order.get("avgPx"))
    commission = amount_of(order.get("commissionNotionalTotalCollected"))
    if average is None or commission is None:
        return filled, None
    side_price = average if contract_side == "yes" else 1 - average
    return filled, round(filled * side_price + commission, 5)


def requests_transport() -> Transport:
    session = requests.Session()

    def send(method: str, url: str, headers: dict, body: bytes | None) -> tuple[int, bytes]:
        response = session.request(method, url, headers=headers, data=body, timeout=REQUEST_TIMEOUT_SEC)
        return response.status_code, response.content

    return send


def _parse(raw: bytes) -> object:
    if not raw:
        return None
    try:
        return json.loads(raw)
    except ValueError:
        return raw[:300].decode("utf-8", "replace")


class PolymarketOrderClient:
    """Signed order calls on one Polymarket US key, plus the public book."""

    def __init__(self, key_id: str, secret_key_path: str | Path, base_url: str = API_BASE_URL,
                 gateway_url: str = GATEWAY_BASE_URL, transport: Transport | None = None,
                 clock: Callable[[], float] = time.time):
        """Reads the key file now, so a bad file fails at start, not on a click."""
        self.key_id = key_id
        self.base_url = base_url.rstrip("/")
        self.gateway_url = gateway_url.rstrip("/")
        self._signer = signer_of_key_file(secret_key_path)
        self._transport = transport or requests_transport()
        self._clock = clock

    def _headers(self, method: str, path: str) -> dict:
        timestamp = str(int(self._clock() * 1000))
        signature = self._signer.sign(f"{timestamp}{method}{path}".encode())
        return {ACCESS_KEY_HEADER: self.key_id, TIMESTAMP_HEADER: timestamp,
                SIGNATURE_HEADER: base64.b64encode(signature).decode(), "Accept": "application/json",
                "Content-Type": "application/json"}

    def request(self, method: str, path: str, body: dict | None = None) -> tuple[int, object]:
        """One signed call: (HTTP status, parsed JSON or a text snippet).
        Raises OSError / requests' errors when nothing came back."""
        raw = json.dumps(body).encode() if body is not None else None
        try:
            status, content = self._transport(method, f"{self.base_url}{path}", self._headers(method, path), raw)
        except requests.RequestException as error:
            raise OSError(str(error)) from error
        return status, _parse(content)

    def place_ioc_order(self, slug: str, contract_side: str, side_price: float, quantity: int) -> tuple[int, object]:
        """ONE order, never re-sent (no client id makes a resend safe)."""
        return self.request("POST", ORDERS_PATH, ioc_order_body(slug, contract_side, side_price, quantity))

    def get_order(self, order_id: str) -> tuple[int, object]:
        return self.request("GET", f"{ORDER_PATH}/{order_id}")

    def cancel_order(self, order_id: str, slug: str) -> tuple[int, object]:
        return self.request("POST", f"{ORDER_PATH}/{order_id}/cancel", {"marketSlug": slug})

    def balance(self) -> tuple[int, object]:
        """GET /v1/account/balances: proves the key signs (the desk's start check)."""
        return self.request("GET", BALANCES_PATH)
