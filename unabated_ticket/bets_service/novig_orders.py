"""Signed Novig v3 REST on ONE trading key: the order client every Novig order
goes through (the automated taker's Novig desk, novig_taker_desk.py, and the
Ticket's Novig Place button after it). Pure transport plus Novig's own rules;
it decides nothing and keeps no state.

Novig facts it rests on (docs.novig.com, read 2026-10-16; plan:
automated_taker_plan.md § Novig):
  * Signing (NOVIG-V3): headers Novig-Key-Id, Novig-Timestamp (unix ms,
    +-30 s) and Novig-Signature = standard base64 of the signature over six
    LF-joined lines: "NOVIG-V3", the timestamp, the METHOD, the path verbatim,
    the canonical query (canonical_query), the lowercase hex SHA-256 of the
    exact body bytes. Ed25519 signs the string; P-256 signs it with SHA-256,
    DER-encoded. The algorithm is read from the key file.
  * A trading key reaches ONE subaccount: it places and cancels there and can
    never move money. Orders always buy an outcome (`outcomeId`).
  * `qty` counts 1c contracts: 100 = $1 paid if the outcome wins. A price is a
    probability on a three-band grid (is_on_price_grid); off-grid is refused.
  * POST /v3/orders answers 201 = queued, not done; the order then moves
    PENDING -> OPEN / FILLED / CANCELED / REJECTED. An IOC order never rests:
    the unfilled rest is cancelled, and one that fills nothing is REJECTED.
  * `clientId` is never checked for uniqueness: a re-sent order is a second
    order. Nothing here ever re-sends; a lost answer is looked up by clientId
    (find_order_by_client_id).
  * Every placement needs the key holder's phone app to have checked its
    location within 3 days, else 451 GEOLOCATION_EXPIRED; reads and cancels
    do not. A data-center address is fine; a VPN or proxy is not.

Side effects: Novig orders and cancels on the configured key (money) through
place_ioc_order / cancel_order; everything else is a read. Nothing on disk.
"""
import base64
import hashlib
import json
import time
import urllib.parse
from collections.abc import Callable
from pathlib import Path

SCHEME = "NOVIG-V3"
KEY_ID_HEADER = "Novig-Key-Id"
TIMESTAMP_HEADER = "Novig-Timestamp"
SIGNATURE_HEADER = "Novig-Signature"
REQUEST_TIMEOUT_SEC = 10
# CloudFront refused Python's default urllib agent on 2026-10-10; name ours.
USER_AGENT = "unabated-ticket-bets-service/1.0"

ORDERS_PATH = "/v3/orders"
FILLS_PATH = "/v3/portfolio/fills"
POSITIONS_PATH = "/v3/portfolio/positions"
TRANSACTION_SETTLEMENT = "SETTLEMENT"
KEYS_PATH = "/v3/keys"
PUBLIC_MARKET_PATH = "/v3/public/catalog/markets"
TIME_IN_FORCE_IOC = "IOC"
# Order statuses (OrderStatus). There is no partly-filled status: track `remaining`.
STATUS_PENDING = "PENDING"
STATUS_OPEN = "OPEN"
# A public market's status once its outcomes are graded; nothing about it changes after.
MARKET_SETTLED = "SETTLED"
# The error code of a market Novig's public catalog no longer lists.
CODE_MARKET_NOT_FOUND = "MARKET_NOT_FOUND"
STATUS_FILLED = "FILLED"
STATUS_CANCELED = "CANCELED"
STATUS_REJECTED = "REJECTED"
ORDER_STATUSES = (STATUS_PENDING, STATUS_OPEN, STATUS_FILLED, STATUS_CANCELED, STATUS_REJECTED)
FINAL_ORDER_STATUSES = (STATUS_FILLED, STATUS_CANCELED, STATUS_REJECTED)
# One contract pays 1c: `qty` units per dollar of payout.
QTY_PER_DOLLAR = 100
PAGE_LIMIT = 500
MAX_PAGES = 20
# Price grid (PlaceOrder.price), in thousandths: 1-50 by 1, 55-945 by 5, 950-999 by 1.
GRID_LOW_MAX = 50
GRID_MID_MIN, GRID_MID_MAX, GRID_MID_STEP = 55, 945, 5
GRID_HIGH_MIN, GRID_HIGH_MAX = 950, 999
GEO_STATUS = 451

# Transport: (method, url, headers, body bytes or None) -> (status, body bytes).
# Raises OSError when no answer came back at all.
Transport = Callable[[str, str, dict, bytes | None], tuple[int, bytes]]


class NovigError(RuntimeError):
    """A Novig answer this module could not use, with what was expected and found."""


# ---- signing --------------------------------------------------------------------------------

def _unreserved_encode(text: str) -> str:
    # RFC 3986 unreserved stays; everything else is %XX in uppercase hex (quote does uppercase).
    return urllib.parse.quote(text, safe="-._~")


def canonical_query(raw_query: str) -> str:
    """Line 5 of the string to sign, from the raw query exactly as sent: split
    on '&', split each pair at the FIRST '=', percent-decode as UTF-8 (a bare
    '+' stays '+'), re-encode keeping only ALPHA DIGIT - . _ ~, sort bytewise
    by name then value (repeats kept), join with '&'. Empty when there is none."""
    if not raw_query:
        return ""
    pairs = []
    for token in raw_query.split("&"):
        name, _, value = token.partition("=")
        pairs.append((_unreserved_encode(urllib.parse.unquote(name)), _unreserved_encode(urllib.parse.unquote(value))))
    pairs.sort(key=lambda pair: (pair[0].encode(), pair[1].encode()))
    return "&".join(f"{name}={value}" for name, value in pairs)


def string_to_sign(timestamp_ms: int, method: str, path: str, raw_query: str, body: bytes) -> str:
    """The six NOVIG-V3 lines, LF-joined, no trailing newline."""
    body_hash = hashlib.sha256(body).hexdigest()
    return "\n".join([SCHEME, str(timestamp_ms), method.upper(), path, canonical_query(raw_query), body_hash])


class NovigSigner:
    """One private key (PKCS#8 PEM, Ed25519 or P-256). Signs canonical strings."""

    def __init__(self, private_key_pem: bytes):
        from cryptography.hazmat.primitives import serialization
        from cryptography.hazmat.primitives.asymmetric import ec, ed25519
        key = serialization.load_pem_private_key(private_key_pem, password=None)
        if isinstance(key, ed25519.Ed25519PrivateKey):
            self.algorithm = "Ed25519"
        elif isinstance(key, ec.EllipticCurvePrivateKey) and key.curve.name == "secp256r1":
            self.algorithm = "P-256"
        else:
            raise NovigError(f"expected an Ed25519 or P-256 private key, got {type(key).__name__}")
        self._key = key

    def sign(self, message: str) -> str:
        from cryptography.hazmat.primitives import hashes
        from cryptography.hazmat.primitives.asymmetric import ec
        data = message.encode()
        if self.algorithm == "Ed25519":
            signature = self._key.sign(data)
        else:
            # cryptography's ECDSA signature is DER-encoded SEQUENCE {r, s}, as Novig requires.
            signature = self._key.sign(data, ec.ECDSA(hashes.SHA256()))
        return base64.b64encode(signature).decode()


def public_key_pem_of(private_key_pem: bytes) -> str:
    """The SPKI PEM of a private key: what Novig stores when a key is registered."""
    from cryptography.hazmat.primitives import serialization
    key = serialization.load_pem_private_key(private_key_pem, password=None)
    return key.public_key().public_bytes(serialization.Encoding.PEM,
                                         serialization.PublicFormat.SubjectPublicKeyInfo).decode()


# ---- price and size ---------------------------------------------------------------------------

def is_on_price_grid(price: float) -> bool:
    """Is `price` one of Novig's submittable prices?"""
    thousandths = round(price * 1000)
    if abs(price * 1000 - thousandths) > 1e-6:
        return False
    if 1 <= thousandths <= GRID_LOW_MAX:
        return True
    if GRID_MID_MIN <= thousandths <= GRID_MID_MAX:
        return thousandths % GRID_MID_STEP == 0
    return GRID_HIGH_MIN <= thousandths <= GRID_HIGH_MAX


def price_string(price: float) -> str:
    """The wire form of a grid price ("0.465"); raises off the grid."""
    if not is_on_price_grid(price):
        raise NovigError(f"price {price!r} is not on Novig's price grid")
    return f"{round(price * 1000) / 1000:.3f}"


def qty_for_stake(stake_dollars: float, price: float) -> int:
    """Whole 1c contracts a stake buys at `price`: floor(stake x 100 / price).
    A Novig contract pays 1c, never Kalshi's $1."""
    if not (0 < price < 1) or stake_dollars < 0:
        raise NovigError(f"qty_for_stake: expected a stake >= 0 and a price in (0, 1), got {stake_dollars!r}, {price!r}")
    return int((stake_dollars * QTY_PER_DOLLAR) // price + 1e-9)


def cost_of(qty: int, price: float) -> float:
    """Dollars `qty` contracts cost at `price`."""
    return qty * price / QTY_PER_DOLLAR


def ioc_order_body(outcome_id: str, price: float, qty: int, client_id: str) -> dict:
    """POST /v3/orders body: buy `qty` of `outcome_id` at no worse than `price`,
    immediate-or-cancel (nothing rests)."""
    if isinstance(qty, bool) or not isinstance(qty, int) or qty < 1:
        raise NovigError(f"qty must be a whole number of contracts >= 1, got {qty!r}")
    return {"outcomeId": outcome_id, "price": price_string(price), "qty": qty, "tif": TIME_IN_FORCE_IOC,
            "clientId": client_id}


def error_code_of(reply: object) -> str | None:
    """Novig's stable error `code` (the only field to branch on), or None."""
    code = reply.get("code") if isinstance(reply, dict) else None
    return code if isinstance(code, str) else None


# ---- transport ---------------------------------------------------------------------------------

def requests_transport() -> Transport:
    """The default transport: one requests.Session, never retried (a POST that
    timed out may still have reached Novig)."""
    import requests
    session = requests.Session()

    def send(method: str, url: str, headers: dict, body: bytes | None) -> tuple[int, bytes]:
        try:
            response = session.request(method, url, headers=headers, data=body, timeout=REQUEST_TIMEOUT_SEC)
        except requests.RequestException as error:
            raise OSError(f"{type(error).__name__}: {error}") from error
        return response.status_code, response.content

    return send


def _parse(raw: bytes) -> object:
    text = raw.decode(errors="replace")
    if not text:
        return {}
    try:
        return json.loads(text)
    except ValueError:
        return text


class NovigOrderClient:
    """Signed Novig v3 REST on ONE key. Every argument comes in (no module-level
    key), so the bot's desk and the Place button each run their own instance on
    their own subaccount's key."""

    def __init__(self, key_id: str, private_key_path: str | Path, base_url: str,
                 transport: Transport | None = None, clock_ms: Callable[[], int] | None = None):
        if not key_id:
            raise NovigError("a Novig key id is required")
        self.key_id = key_id
        self.base_url = base_url.rstrip("/")
        self._signer = NovigSigner(Path(private_key_path).read_bytes())
        self._transport = transport or requests_transport()
        self._clock_ms = clock_ms or (lambda: int(time.time() * 1000))

    def request(self, method: str, path: str, query: list[tuple[str, str]] | None = None,
                body: dict | None = None) -> tuple[int, object]:
        """(status, parsed body or text). Raises OSError on no answer at all.
        The body is serialised once and those exact bytes are hashed and sent."""
        raw_query = urllib.parse.urlencode(query or [], quote_via=urllib.parse.quote, safe="-._~")
        data = json.dumps(body, separators=(",", ":")).encode() if body is not None else b""
        timestamp_ms = self._clock_ms()
        headers = {
            KEY_ID_HEADER: self.key_id,
            TIMESTAMP_HEADER: str(timestamp_ms),
            SIGNATURE_HEADER: self._signer.sign(string_to_sign(timestamp_ms, method, path, raw_query, data)),
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        }
        if body is not None:
            # Without a content type Novig hashes zero bytes instead of the body.
            headers["Content-Type"] = "application/json"
        url = f"{self.base_url}{path}" + (f"?{raw_query}" if raw_query else "")
        status, raw = self._transport(method.upper(), url, headers, data if body is not None else None)
        return status, _parse(raw)

    # ---- orders (money) ----

    def place_ioc_order(self, outcome_id: str, price: float, qty: int, client_id: str) -> tuple[int, object]:
        """POST /v3/orders, IOC. Sent once, never retried. 201 = queued."""
        return self.request("POST", ORDERS_PATH, body=ioc_order_body(outcome_id, price, qty, client_id))

    def cancel_order(self, order_id: str) -> tuple[int, object]:
        """DELETE /v3/orders/{id}. Idempotent at Novig."""
        return self.request("DELETE", f"{ORDERS_PATH}/{order_id}")

    # ---- reads ----

    def get_order(self, order_id: str) -> tuple[int, object]:
        """GET /v3/orders/{id}: may 404 right after a 201 (still queued)."""
        return self.request("GET", f"{ORDERS_PATH}/{order_id}")

    def _paged(self, path: str, query: list[tuple[str, str]]) -> list[dict]:
        """Every item of a cursor-paged listing; raises NovigError on any
        failed page or a listing that does not end (never a partial list)."""
        items: list[dict] = []
        after = None
        for _ in range(MAX_PAGES):
            page_query = [*query, ("limit", str(PAGE_LIMIT))] + ([("after", after)] if after else [])
            status, reply = self.request("GET", path, page_query)
            if status != 200 or not isinstance(reply, dict) or not isinstance(reply.get("items"), list):
                raise NovigError(f"GET {path}: expected 200 with an items list, got HTTP {status} {str(reply)[:300]}")
            items.extend(item for item in reply["items"] if isinstance(item, dict))
            after = reply.get("next")
            if not after:
                return items
        raise NovigError(f"GET {path}: more than {MAX_PAGES} pages; refusing a partial list")

    def orders(self, outcome_id: str, status: str) -> list[dict]:
        """This subaccount's orders on one outcome in one status (the listing defaults to OPEN only)."""
        if status not in ORDER_STATUSES:
            raise NovigError(f"order status must be one of {ORDER_STATUSES}, got {status!r}")
        return self._paged(ORDERS_PATH, [("outcome", outcome_id), ("status", status)])

    def fills(self, order_id: str | None = None, outcome_id: str | None = None,
              starts_after_ms: int | None = None) -> list[dict]:
        """This subaccount's fills, filtered by order, outcome or event start."""
        query = []
        if order_id:
            query.append(("order", order_id))
        if outcome_id:
            query.append(("outcome", outcome_id))
        if starts_after_ms is not None:
            query.append(("startsAfter", str(starts_after_ms)))
        return self._paged(FILLS_PATH, query)

    def open_positions(self) -> list[dict]:
        """This subaccount's nonzero positions ([{marketId, outcomeId, qty, cost}]).
        A settled market's position leaves this list (Novig closes it)."""
        status, reply = self.request("GET", POSITIONS_PATH)
        if status != 200 or not isinstance(reply, list):
            raise NovigError(f"GET {POSITIONS_PATH}: expected 200 with a list, got HTTP {status} {str(reply)[:300]}")
        return [item for item in reply if isinstance(item, dict)]

    def settlements(self, starts_after_ms: int) -> list[dict]:
        """This subaccount's SETTLEMENT ledger rows ([{amount, ref, ts}]) on
        events that start after `starts_after_ms`: one per settled market that
        paid it something (ref = market id; a loss pays nothing, so writes none)."""
        return self._paged(f"/v3/account/subaccounts/{self.key_id}/transactions",
                           [("kind", TRANSACTION_SETTLEMENT), ("startsAfter", str(starts_after_ms))])

    def balance(self) -> float:
        """This subaccount's balance in dollars (GET .../subaccounts/{own key}/balance)."""
        status, reply = self.request("GET", f"/v3/account/subaccounts/{self.key_id}/balance")
        balance = reply.get("balance") if isinstance(reply, dict) else None
        if status != 200 or not isinstance(balance, str):
            raise NovigError(f"balance: expected 200 with a balance string, got HTTP {status} {str(reply)[:300]}")
        try:
            return float(balance)
        except ValueError as error:
            raise NovigError(f"balance: {balance!r} is not a decimal") from error

    def find_order_by_client_id(self, outcome_id: str, client_id: str) -> dict | None:
        """The order on `outcome_id` carrying `client_id`, in any status, or
        None when Novig lists it nowhere. Raises NovigError when a listing
        cannot be read: "not found" is only ever said off complete listings."""
        for status in ORDER_STATUSES:
            for order in self.orders(outcome_id, status):
                if order.get("clientId") == client_id and order.get("orderId"):
                    return order
        return None

    def market(self, market_id: str) -> dict:
        """The public catalog's market (status, strike, outcomes, fee, eventId)."""
        market = self.market_if_listed(market_id)
        if market is None:
            raise NovigError(f"market {market_id}: Novig's catalog does not list it")
        return market

    def market_if_listed(self, market_id: str) -> dict | None:
        """The public catalog's market, or None when Novig answers 404
        MARKET_NOT_FOUND (it drops some markets from the catalog)."""
        status, reply = self.request("GET", f"{PUBLIC_MARKET_PATH}/{market_id}")
        if status == 404 and error_code_of(reply) == CODE_MARKET_NOT_FOUND:
            return None
        if status != 200 or not isinstance(reply, dict) or not reply.get("marketId"):
            raise NovigError(f"market {market_id}: expected 200 with a market, got HTTP {status} {str(reply)[:300]}")
        return reply

    def is_management_key(self) -> bool | None:
        """True when this key lists keys (only a management key can), False
        when Novig refuses it the listing, None when the answer says neither."""
        status, _reply = self.request("GET", KEYS_PATH)
        if status == 200:
            return True
        if status in (401, 403):
            return False
        return None


def fill_outcome(filled: int, requested: int) -> str:
    """An order's row status from its contracts filled: 'unfilled' | 'filled' |
    'partial' (the bot desk's, Place's and Place's later confirmation)."""
    return "unfilled" if filled == 0 else "filled" if filled >= requested else "partial"


def order_totals(fills: list[dict]) -> tuple[int, float, float]:
    """(contracts, dollars paid, fees) over one order's fills. Raises NovigError
    on a fill whose qty or cost cannot be read: a missing field must never read
    as a free fill."""
    qty_total = 0
    cost_total = 0.0
    fee_total = 0.0
    for fill in fills:
        qty = fill.get("qty")
        try:
            cost = float(fill.get("cost"))
            fee = float(fill["fee"]) if fill.get("fee") not in (None, "") else 0.0
        except (TypeError, ValueError) as error:
            raise NovigError(f"fill {fill.get('fillId') or '?'} unreadable: {error}") from error
        if isinstance(qty, bool) or not isinstance(qty, int) or qty < 1 or cost < 0:
            raise NovigError(f"fill {fill.get('fillId') or '?'} unreadable: qty {qty!r}, cost {fill.get('cost')!r}")
        qty_total += qty
        cost_total += cost
        fee_total += fee
    return qty_total, round(cost_total, 5), round(fee_total, 5)
