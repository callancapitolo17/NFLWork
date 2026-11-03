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
    minutes; a minute before it runs out the client calls POST /auth/refresh
    {refresh_token} -> {data: {access_token, access_expire_time}}, and logs in
    afresh only when it holds no refresh token or the refresh is refused.
    Every login opens a SESSION, and a key may hold only 20 at once
    (`session_num_exceed`, 403), freed only as they expire, so logging in
    every 10 minutes would lock the key out within hours.
  * Rate limit: the docs say 50 requests/s, but the live API answers 429
    "please send only 1 requests/second" (2026-11-03), per KEY, so the source's
    poll and Place share one budget. Every send on the client is spaced
    MIN_SEND_INTERVAL_SEC apart; a read that still gets a 429 waits and is
    sent once more (never an order: a 429 order was not taken, but the click
    is the one decision to send it).
  * ProphetX's CloudFront refuses some IPs outright: a 403 HTML page ("Request
    blocked") before the API. This VM's IP got one from the old MLB maker's
    SGP scraping (2026-09); retrying once a minute keeps such a block alive,
    so a refusal stops EVERY call on the client for BLOCKED_BACKOFF_START_SEC,
    doubling to BLOCKED_BACKOFF_MAX_SEC, until a call gets through.
  * Prices are AMERICAN odds and must sit on ProphetX's price ladder
    (GET /v4/mm/get_price_ladder). An order's `quantity` is the stake in
    dollars ("random $1 orders" send quantity 1.0).
  * Orders only buy: there is no side flag, an order backs its strike.
  * POST /v4/mm/submit_order {strike_id, price, quantity, external_id,
    order_strategy}: `fillOrKill` fills what it can at the price or better and
    cancels the rest about 2 s after the last activity (it rests that long).
  * GET /v4/mm/get_order_history rows (the order shape): order status open /
    inactive / invalid / canceled / void / wiped / settled / closed /
    manually_settled, matching_status unmatched / partially_matched /
    fully_matched, filled_quantity / open_quantity in dollars, `fee`.
    GET /v4/mm/get_order/<order_id> is not used: live, it answers 404
    data_not_found for 8 s+ after a send (2026-11-03) while the history
    already lists the order.

Side effects: ProphetX orders on the configured key (money) through
submit_fill_or_kill; everything else is a read. Nothing here cancels: a
fill-or-kill kills its own remainder.
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
REFRESH_PATH = "/auth/refresh"
NOT_FOUND_ERROR = "data_not_found"
SPORT_EVENTS_PATH = "/mm/get_sport_events"
MARKETS_PATH = "/v4/mm/get_markets"
STRIKES_PATH = "/v4/mm/get_strikes"
ORDER_HISTORY_PATH = "/v4/mm/get_order_history"
TRADES_PATH = "/v4/mm/get_trades"
SUBMIT_ORDER_PATH = "/v4/mm/submit_order"
BALANCE_PATH = "/v4/mm/get_balance"
# An access token lasts 10 minutes; a new login goes out this long before that.
ACCESS_TOKEN_LIFETIME_SEC = 600
ACCESS_TOKEN_MARGIN_SEC = 60
# After a refusal (CloudFront block, or a refused login: a 403 session cap or
# a bad key) nothing is sent for this long, doubling per refusal up to the max.
BLOCKED_BACKOFF_START_SEC = 3600.0
BLOCKED_BACKOFF_MAX_SEC = 86400.0
# ProphetX takes 1 request per second per key; a little over that leaves room for clock drift.
MIN_SEND_INTERVAL_SEC = 1.1
# A read that still gets a 429 is sent once more after this long.
RATE_LIMITED_WAIT_SEC = 2.0
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


class ProphetXBlocked(ProphetXError):
    """ProphetX refused this server (or the key's login) and the client is
    sending nothing until `until` (unix seconds)."""

    def __init__(self, message: str, until: float, server_refused: bool):
        super().__init__(message)
        self.until = until
        # True for CloudFront refusing this server (says nothing about the key);
        # False for a refused login (the key, or its session cap).
        self.server_refused = server_refused


def is_cdn_block(status: int, reply: object) -> bool:
    """A 403 that is CloudFront's HTML page, not the API's JSON error."""
    return status == 403 and not isinstance(reply, dict)


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
                 transport: Transport | None = None, clock: Callable[[], float] = time.time,
                 min_send_interval_sec: float = MIN_SEND_INTERVAL_SEC,
                 monotonic: Callable[[], float] = time.monotonic, sleep: Callable[[float], None] = time.sleep):
        if not access_key or not secret_key:
            raise ProphetXError("a ProphetX access key and secret key are both needed")
        self.access_key = access_key
        self._secret_key = secret_key
        self.base_url = base_url.rstrip("/")
        self._transport = transport or requests_transport()
        self._clock = clock
        self._token: str | None = None
        self._token_good_until = 0.0
        self._refresh_token: str | None = None
        self._token_lock = threading.Lock()
        self._min_send_interval_sec = min_send_interval_sec
        self._monotonic = monotonic
        self._sleep = sleep
        self._send_lock = threading.Lock()
        self._last_send_at = -math.inf
        self._blocked_until = 0.0
        self._blocked_reason = ""
        self._blocked_by_server = False
        self._next_backoff_sec = BLOCKED_BACKOFF_START_SEC

    # ---- refusals --------------------------------------------------------------

    @property
    def blocked_until(self) -> float:
        """Unix seconds until which nothing is sent (0 when not blocked)."""
        return self._blocked_until if self._clock() < self._blocked_until else 0.0

    def _raise_if_blocked(self) -> None:
        if self._clock() < self._blocked_until:
            raise ProphetXBlocked(f"{self._blocked_reason}; not retrying until "
                                  f"{time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime(self._blocked_until))}",
                                  self._blocked_until, self._blocked_by_server)

    def _block(self, reason: str, server_refused: bool) -> None:
        """Stop every call for the current backoff, and double the next one."""
        self._blocked_until = self._clock() + self._next_backoff_sec
        self._blocked_reason = reason
        self._blocked_by_server = server_refused
        log.error("prophetx: %s; sending nothing for %.0f min", reason, self._next_backoff_sec / 60)
        self._next_backoff_sec = min(self._next_backoff_sec * 2, BLOCKED_BACKOFF_MAX_SEC)
        self._raise_if_blocked()

    def _got_through(self) -> None:
        self._next_backoff_sec = BLOCKED_BACKOFF_START_SEC

    # ---- session -----------------------------------------------------------------

    def _send(self, method: str, path: str, headers: dict, body: dict | None) -> tuple[int, object]:
        raw = json.dumps(body).encode() if body is not None else None
        with self._send_lock:
            # One key, one budget: the poll thread and a Place click take turns.
            wait = self._last_send_at + self._min_send_interval_sec - self._monotonic()
            if wait > 0:
                self._sleep(wait)
            self._last_send_at = self._monotonic()
            try:
                status, content = self._transport(method, f"{self.base_url}{path}", headers, raw)
            except requests.RequestException as error:
                raise OSError(str(error)) from error
        reply = _parse(content)
        if is_cdn_block(status, reply):
            self._block(f"ProphetX's CloudFront refused this server (HTTP 403 on {path})", server_refused=True)
        return status, reply

    def _take_token(self, data: object, what: str) -> None:
        token = data.get("access_token") if isinstance(data, dict) else None
        if not isinstance(token, str) or not token:
            raise ProphetXError(f"{what} answered no access_token: {str(data)[:200]}")
        self._token = token
        self._token_good_until = self._clock() + ACCESS_TOKEN_LIFETIME_SEC - ACCESS_TOKEN_MARGIN_SEC
        refresh = data.get("refresh_token")
        if isinstance(refresh, str) and refresh:
            self._refresh_token = refresh

    def _login(self) -> None:
        status, reply = self._send("POST", LOGIN_PATH, {"Content-Type": "application/json", "Accept": "application/json"},
                                   {"access_key": self.access_key, "secret_key": self._secret_key})
        if status != 200:
            # A refused login (403 session cap, 401 bad key) only gets worse retried each minute.
            self._block(f"ProphetX login refused: HTTP {status} {str(reply)[:200]}", server_refused=False)
        self._take_token(data_of(reply, "ProphetX login"), "ProphetX login")
        log.info("prophetx: logged in (a new session)")

    def _refresh(self) -> bool:
        """A new access token on the held refresh token, without a new session.
        False when ProphetX refuses it (then the caller logs in)."""
        status, reply = self._send("POST", REFRESH_PATH, {"Content-Type": "application/json", "Accept": "application/json"},
                                   {"refresh_token": self._refresh_token})
        if status != 200:
            log.warning("prophetx: token refresh refused (HTTP %s %s); logging in again", status, str(reply)[:200])
            self._refresh_token = None
            return False
        self._take_token(data_of(reply, "ProphetX refresh"), "ProphetX refresh")
        return True

    def _bearer(self, force_login: bool = False) -> str:
        with self._token_lock:
            self._raise_if_blocked()
            expired = self._token is None or self._clock() >= self._token_good_until
            if force_login or expired:
                if not (self._refresh_token and self._refresh()):
                    self._login()
            return self._token

    def request(self, method: str, path: str, params: dict | None = None, body: dict | None = None) -> tuple[int, object]:
        """One call with a live token: (HTTP status, parsed JSON or a text
        snippet). A 401 gets a new token and sends once more: a 401 is a
        request ProphetX did not take. Raises OSError when nothing came back,
        ProphetXBlocked while the client is backing off a refusal."""
        query = f"?{urlencode(params, doseq=False)}" if params else ""
        for attempt in range(2):
            headers = {"Authorization": f"Bearer {self._bearer(force_login=attempt > 0)}",
                       "Accept": "application/json", "Content-Type": "application/json"}
            status, reply = self._send(method, f"{path}{query}", headers, body)
            if status == 429 and method == "GET":
                log.warning("prophetx: 429 on GET %s; sending once more after %.1f s", path, RATE_LIMITED_WAIT_SEC)
                self._sleep(RATE_LIMITED_WAIT_SEC)
                status, reply = self._send(method, f"{path}{query}", headers, body)
            if status != 401:
                self._got_through()
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
        status, reply = self.request("GET", SPORT_EVENTS_PATH, params)
        if status == 404 and isinstance(reply, dict) and reply.get("error") == NOT_FOUND_ERROR:
            # ProphetX answers 404 data_not_found when none of the named games is listed
            # (seen live 2026-11-03 on ended games), not an empty list.
            return []
        if status != 200:
            raise ProphetXError(f"GET {SPORT_EVENTS_PATH}: HTTP {status} {str(reply)[:200]}")
        data = data_of(reply, f"GET {SPORT_EVENTS_PATH}")
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
