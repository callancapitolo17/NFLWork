"""Novig catalog reads for the runner, signed with the bot's Novig key.

The runner checks every Novig order and Place quote against Novig's catalog:
the open pregame market list, one market, its event, and (with no stream
book) its book. It used to read them from Novig's unkeyed `/v3/public`
routes, which share one small per-IP throttle (about 10 at once, then 2 a
second, measured by a third party). On 2026-10-28 and 10-30 Place quotes were
refused with HTTP 429 on the market and event reads. Novig serves the same
catalog on keyed routes (`/v3/catalog/...`), which spend the key's own `read`
bucket (64 at once, 16 a second, docs.novig.com/throttling). This module reads
those routes for the runner (POST /novig_read.json, service.py).

Only the four catalog routes the runner needs are allowed, each path and query
checked before anything is signed: the key is the bot's trading key, and this
module only ever sends GETs, so a bad body can neither place an order nor read
an account.

Side effects: GET requests to Novig on the bot's key. Nothing on disk, no DuckDB.
"""
import queue
import re
from collections.abc import Callable

from unabated_ticket.bets_service.novig_orders import (
    NovigOrderClient,
    Transport,
    requests_transport,
)

CATALOG_MARKETS_PATH = "/v3/catalog/markets"
# Novig's ids are UUIDs; anything else never reaches the signer.
_ID = r"[0-9A-Za-z-]{1,64}"
_MARKET_PATH = re.compile(rf"^/v3/catalog/markets/{_ID}$")
_BOOK_PATH = re.compile(rf"^/v3/catalog/markets/{_ID}/book$")
_EVENT_PATH = re.compile(rf"^/v3/catalog/events/{_ID}$")
# The market list's filters the runner sends (runner.js readNovigMarkets).
MARKET_LIST_QUERY_KEYS = frozenset({"league", "eventStatus", "limit", "after"})
MAX_QUERY_VALUE_CHARS = 200
# Concurrent reads (a Place quote's market and event go at once, while a taker
# pass may be reading the list): each takes its own HTTP session.
SESSIONS = 4


def validate_read_request(body: object) -> tuple[str, list[tuple[str, str]]] | str:
    """(path, query pairs) for one allowed catalog GET, or why the body is refused.
    body {path, query?: {name: value}}: query only on the market list."""
    if not isinstance(body, dict):
        return "expected a JSON object body"
    path = body.get("path")
    if not isinstance(path, str):
        return f"path must be a string, got {path!r}"
    query = body.get("query") or {}
    if not isinstance(query, dict):
        return f"query must be an object, got {query!r}"
    if path != CATALOG_MARKETS_PATH:
        if not (_MARKET_PATH.match(path) or _BOOK_PATH.match(path) or _EVENT_PATH.match(path)):
            return f"path must be {CATALOG_MARKETS_PATH}, a market, its book or an event under /v3/catalog, got {path!r}"
        if query:
            return f"only {CATALOG_MARKETS_PATH} takes a query, got {sorted(query)} on {path}"
        return path, []
    pairs = []
    for name, value in query.items():
        if name not in MARKET_LIST_QUERY_KEYS:
            return f"query {name!r} is not one of {sorted(MARKET_LIST_QUERY_KEYS)}"
        if not isinstance(value, str) or not value or len(value) > MAX_QUERY_VALUE_CHARS:
            return f"query {name} must be a non-empty string of at most {MAX_QUERY_VALUE_CHARS} characters, got {value!r}"
        pairs.append((name, value))
    return path, pairs


def pooled_transport(make_transport: Callable[[], Transport] = requests_transport, size: int = SESSIONS) -> Transport:
    """A transport that runs each request on one of `size` sessions (requests'
    Session is not made for concurrent use), waiting for a free one."""
    free: queue.Queue = queue.Queue()
    for _ in range(size):
        free.put(make_transport())

    def send(method: str, url: str, headers: dict, body: bytes | None) -> tuple[int, bytes]:
        transport = free.get()
        try:
            return transport(method, url, headers, body)
        finally:
            free.put(transport)

    return send


class NovigCatalogReads:
    """GETs of Novig's keyed catalog on one key's client. read() raises OSError
    when Novig did not answer; any answer, a 429 included, is returned as is."""

    def __init__(self, client: NovigOrderClient):
        self.client = client

    def read(self, path: str, query: list[tuple[str, str]]) -> dict:
        """{status, body}: Novig's status and parsed body for one allowed GET (validate_read_request first)."""
        status, body = self.client.request("GET", path, query)
        return {"status": status, "body": body}
