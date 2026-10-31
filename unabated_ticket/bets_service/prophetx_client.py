"""ProphetX Trading API on ONE access/secret key pair: the transport both the
ProphetX bet source (sources/prophetx.py) and the Ticket's ProphetX Place desk
(prophetx_place.py) read and trade through. Pure transport plus ProphetX's own
session rules; it decides nothing and keeps nothing on disk.

ProphetX facts it rests on (docs.prophetx.co, read 2026-10-31; the API itself
was not reachable from the cloud container, so nothing here is a live capture):
  * Base URL https://cash.api.prophetx.co/partner (sandbox
    https://api.sandbox.prophetx.dev/partner).
  * POST /auth/login {access_key, secret_key} -> {data: {access_token,
    access_expire_time, refresh_token, refresh_expire_time}}; every other call
    sends `Authorization: Bearer <access_token>`. An access token lasts 10
    minutes, so the client logs in again a minute before it runs out (a fresh
    login, never /auth/refresh, whose body the docs leave unstated).
  * Prices are AMERICAN odds and must sit on ProphetX's price ladder
    (GET /v4/mm/get_price_ladder). An order's `quantity` is the stake in
    dollars ("random $1 orders" send quantity 1.0).
  * Orders only buy: there is no side flag, an order backs its strike.
  * POST /v4/mm/submit_order {strike_id, price, quantity, external_id,
    order_strategy}: `fillOrKill` fills what it can at the price or better and
    cancels the rest about 2 s after the last activity (it rests that long).
  * GET /v4/mm/get_order/<order_id> -> {data: order}; order status open /
    inactive / invalid / canceled / void / wiped / settled / closed /
    manually_settled, matching_status unmatched / partially_matched /
    fully_matched, filled_quantity / open_quantity in dollars, `fee`.

Side effects: ProphetX orders and cancels on the configured key (money) through
submit_fill_or_kill / cancel_order; everything else is a read.
"""
import json
import logging
import math
import threading
import time
from collections.abc import Callable
from urllib.parse import urlencode

import requests

log = logging.getLogger(__name__)

DEFAULT_BASE_URL = "https://cash.api.prophetx.co/partner"
REQUEST_TIMEOUT_SEC = 15
LOGIN_PATH = "/auth/login"
SPORT_EVENTS_PATH = "/mm/get_sport_events"
MARKETS_PATH = "/v4/mm/get_markets"
STRIKES_PATH = "/v4/mm/get_strikes"
ORDER_HISTORY_PATH = "/v4/mm/get_order_history"
TRADES_PATH = "/v4/mm/get_trades"
ORDER_PATH = "/v4/mm/get_order"
SUBMIT_ORDER_PATH = "/v4/mm/submit_order"
CANCEL_ORDER_PATH = "/v4/mm/cancel_order"
BALANCE_PATH = "/v4/mm/get_balance"
# An access token lasts 10 minutes; a new login goes out this long before that.
ACCESS_TOKEN_LIFETIME_SEC = 600
ACCESS_TOKEN_MARGIN_SEC = 60
# get_strikes takes at most this many ids per read.
MAX_STRIKES_PER_READ = 20
ORDER_STRATEGY_FILL_OR_KILL = "fillOrKill"
# Order statuses after which nothing more fills, and the ones still on the book.
FINAL_ORDER_STATUSES = {"canceled", "void", "wiped", "invalid", "settled", "closed", "manually_settled"}
RESTING_ORDER_STATUSES = {"open", "inactive"}
# The smallest American price magnitude ProphetX can quote (even money is +-100).
MIN_AMERICAN_MAGNITUDE = 100

Transport = Callable[[str, str, dict, bytes | None], tuple[int, bytes]]


class ProphetXError(RuntimeError):
    """A ProphetX answer that is not what the call expects."""


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


def is_american_price(price: object) -> bool:
    """A whole American price ProphetX can carry: +100 or more, -100 or less."""
    return (isinstance(price, (int, float)) and not isinstance(price, bool) and math.isfinite(price)
            and float(price).is_integer() and abs(price) >= MIN_AMERICAN_MAGNITUDE)


def fill_or_kill_body(strike_id: str, american_price: int, stake_dollars: float, external_id: str) -> dict:
    """ONE fill-or-kill buy of `strike_id` for `stake_dollars` at
    `american_price` or better."""
    if not isinstance(strike_id, str) or not strike_id:
        raise ValueError(f"strike_id must be ProphetX's strike id, got {strike_id!r}")
    if not is_american_price(american_price):
        raise ValueError(f"price must be a whole American price (|price| >= 100), got {american_price!r}")
    if not isinstance(stake_dollars, (int, float)) or not math.isfinite(stake_dollars) or stake_dollars <= 0:
        raise ValueError(f"stake must be dollars above 0, got {stake_dollars!r}")
    return {"strike_id": strike_id, "price": int(american_price), "quantity": round(float(stake_dollars), 2),
            "external_id": external_id, "order_strategy": ORDER_STRATEGY_FILL_OR_KILL}


def data_of(reply: object, what: str) -> object:
    """The `data` of a ProphetX reply, or ProphetXError naming what came back."""
    if not isinstance(reply, dict) or "data" not in reply:
        raise ProphetXError(f"{what}: expected {{data: ...}}, got {str(reply)[:200]}")
    return reply["data"]


class ProphetXClient:
    """Calls on one ProphetX key pair. Thread-safe: the bet source's poll and a
    Place click share one client and one session token."""

    def __init__(self, access_key: str, secret_key: str, base_url: str = DEFAULT_BASE_URL,
                 transport: Transport | None = None, clock: Callable[[], float] = time.time):
        if not access_key or not secret_key:
            raise ProphetXError("a ProphetX access key and secret key are both needed")
        self.access_key = access_key
        self._secret_key = secret_key
        self.base_url = base_url.rstrip("/")
        self._transport = transport or requests_transport()
        self._clock = clock
        self._token: str | None = None
        self._token_good_until = 0.0
        self._token_lock = threading.Lock()

    # ---- session -----------------------------------------------------------------

    def _send(self, method: str, path: str, headers: dict, body: dict | None) -> tuple[int, object]:
        raw = json.dumps(body).encode() if body is not None else None
        try:
            status, content = self._transport(method, f"{self.base_url}{path}", headers, raw)
        except requests.RequestException as error:
            raise OSError(str(error)) from error
        return status, _parse(content)

    def _login(self) -> None:
        status, reply = self._send("POST", LOGIN_PATH, {"Content-Type": "application/json", "Accept": "application/json"},
                                   {"access_key": self.access_key, "secret_key": self._secret_key})
        if status != 200:
            raise ProphetXError(f"ProphetX login refused: HTTP {status} {str(reply)[:200]}")
        data = data_of(reply, "ProphetX login")
        token = data.get("access_token") if isinstance(data, dict) else None
        if not isinstance(token, str) or not token:
            raise ProphetXError(f"ProphetX login answered no access_token: {str(reply)[:200]}")
        self._token = token
        self._token_good_until = self._clock() + ACCESS_TOKEN_LIFETIME_SEC - ACCESS_TOKEN_MARGIN_SEC

    def _bearer(self, force_login: bool = False) -> str:
        with self._token_lock:
            if force_login or self._token is None or self._clock() >= self._token_good_until:
                self._login()
            return self._token

    def request(self, method: str, path: str, params: dict | None = None, body: dict | None = None) -> tuple[int, object]:
        """One call with a live token: (HTTP status, parsed JSON or a text
        snippet). A 401 logs in again and sends once more: a 401 is a request
        ProphetX did not take. Raises OSError when nothing came back."""
        query = f"?{urlencode(params, doseq=False)}" if params else ""
        for attempt in range(2):
            headers = {"Authorization": f"Bearer {self._bearer(force_login=attempt > 0)}",
                       "Accept": "application/json", "Content-Type": "application/json"}
            status, reply = self._send(method, f"{path}{query}", headers, body)
            if status != 401:
                return status, reply
        return status, reply

    def get_data(self, path: str, params: dict | None = None) -> object:
        """GET's `data`, or ProphetXError on anything but a 200 with data."""
        status, reply = self.request("GET", path, params)
        if status != 200:
            raise ProphetXError(f"GET {path}: HTTP {status} {str(reply)[:200]}")
        return data_of(reply, f"GET {path}")

    # ---- reads ---------------------------------------------------------------------

    def sport_events(self, event_ids: list[int] | None = None) -> list[dict]:
        """GET /mm/get_sport_events: the listed (not ended) games, all of them
        or the ones named. ProphetX drops a game from this list once it ends."""
        params = {"event_ids": ",".join(str(event_id) for event_id in event_ids)} if event_ids else None
        data = self.get_data(SPORT_EVENTS_PATH, params)
        events = data.get("sport_events") if isinstance(data, dict) else None
        if not isinstance(events, list):
            raise ProphetXError(f"get_sport_events: no sport_events list in {str(data)[:200]}")
        return [event for event in events if isinstance(event, dict)]

    def markets(self, event_id: int) -> list[dict]:
        """GET /v4/mm/get_markets?event_id=: every market of one game, alt
        lines nested under `market_strikes`, each side's offers in `selections`."""
        data = self.get_data(MARKETS_PATH, {"event_id": event_id})
        markets = data.get("markets") if isinstance(data, dict) else None
        if not isinstance(markets, list):
            raise ProphetXError(f"get_markets {event_id}: no markets list in {str(data)[:200]}")
        return [market for market in markets if isinstance(market, dict)]

    def strikes(self, strike_ids: list[str]) -> list[dict]:
        """GET /v4/mm/get_strikes: {event_id, market_id, strike_id, type} per id
        (at most MAX_STRIKES_PER_READ)."""
        data = self.get_data(STRIKES_PATH, {"strike_ids": ",".join(strike_ids[:MAX_STRIKES_PER_READ])})
        if not isinstance(data, list):
            raise ProphetXError(f"get_strikes: expected a list, got {str(data)[:200]}")
        return [strike for strike in data if isinstance(strike, dict)]

    def _pages(self, path: str, params: dict, list_key: str, max_pages: int) -> list[dict]:
        rows: list[dict] = []
        cursor = None
        for _page in range(max_pages):
            data = self.get_data(path, {**params, **({"next_cursor": cursor} if cursor else {})})
            page = data.get(list_key) if isinstance(data, dict) else None
            if not isinstance(page, list):
                raise ProphetXError(f"GET {path}: no {list_key} list in {str(data)[:200]}")
            rows.extend(row for row in page if isinstance(row, dict))
            cursor = data.get("next_cursor")
            if not cursor or not page:
                return rows
        raise ProphetXError(f"GET {path}: still paging after {max_pages} pages")

    def order_history(self, created_from: int, created_to: int, max_pages: int = 50) -> list[dict]:
        """Every order created in [created_from, created_to] (unix seconds, at
        most 60 days apart), every status, newest pages first."""
        return self._pages(ORDER_HISTORY_PATH, {"from": created_from, "to": created_to, "limit": 1000},
                           "orders", max_pages)

    def trades(self, matched_from: int, matched_to: int, max_pages: int = 100) -> list[dict]:
        """Every fill matched in [matched_from, matched_to] (unix seconds):
        {order_id, price, quantity, profit, fee, matched_at, ...}."""
        return self._pages(TRADES_PATH, {"from": matched_from, "to": matched_to, "limit": 100}, "trades", max_pages)

    def get_order(self, order_id: str) -> tuple[int, object]:
        return self.request("GET", f"{ORDER_PATH}/{order_id}")

    def balance(self) -> tuple[int, object]:
        """GET /v4/mm/get_balance: proves the key logs in (the desk's start check)."""
        return self.request("GET", BALANCE_PATH)

    # ---- orders (money) --------------------------------------------------------------

    def submit_fill_or_kill(self, strike_id: str, american_price: int, stake_dollars: float,
                            external_id: str) -> tuple[int, object]:
        """ONE order. Never re-sent by the caller: only a 401 (not taken) is
        sent again, inside request()."""
        return self.request("POST", SUBMIT_ORDER_PATH,
                            body=fill_or_kill_body(strike_id, american_price, stake_dollars, external_id))

    def cancel_order(self, order_id: str, external_id: str) -> tuple[int, object]:
        return self.request("POST", CANCEL_ORDER_PATH, body={"order_id": order_id, "external_id": external_id})
