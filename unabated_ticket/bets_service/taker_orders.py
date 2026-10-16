"""The automated taker's order desk (automated_taker_plan.md): the only place a
bot order is signed and sent. The server runner (server/runner.js) decides;
this module checks the stop file and the $500 ceiling, sends one limit order
at the ask the decision was checked at, cancels whatever did not fill at
once, and reads what filled.

Armed only when a dedicated Kalshi key is configured
(KALSHI_TAKER_API_KEY_ID + KALSHI_TAKER_PRIVATE_KEY_PATH, a key made for the
bot alone so it can be revoked alone) together with the bot's own subaccount
(KALSHI_TAKER_SUBACCOUNT, 1-63; Cal, 2026-10-16), and only once Kalshi shows
the key is locked to that subaccount (key_lock_problem): it reads its balance
there and is refused the primary account's. The money in the subaccount is
then all the bot can lose. Without all of it, POST /taker/order.json answers
503 and nothing is sent. taker_subaccount_setup.py makes the subaccount, funds
it and mints the locked key.

Rules, in order, for every request:
  1. the stop file (config.TAKER_STOP_PATH) exists -> refused, "stopped";
  2. the same client_order_id was already sent -> refused, "duplicate";
  3. open orders are brought up to date when the last refresh is older than
     REFRESH_SEC (the desk's own thread does the same every REFRESH_SEC, so
     nothing below waits on the runner): see _refresh_open_orders;
  4. the decision is fresh (decidedAt within MAX_DECISION_AGE_SEC: a request
     that waited behind the desk lock is not the live price any more) and the
     game has not started, and the runner priced the same Kalshi (its
     kalshiBaseUrl is this desk's); else refused, "stale_decision" /
     "started" / "env_mismatch";
  5. unseen fill: an earlier bot order on the SAME GAME, on Kalshi or on
     Novig (novig_taker_desk.py), that may still hold
     contracts must be among the held bets the runner sized this stake on
     (`heldIds`, the ids bets.js matched to the row: "kalshi:<ticker>:<side>:s<N>"
     for the bot's subaccount N, sources/kalshi.py record_id_of), else
     refused, "unseen_fill". Conditional Kelly is only right once it
     sees the earlier fill; until then an add would be sized as if alone;
  6. ceiling: contracts are cut so that open exposure + this order's cost
     stays at or under BOT_BANKROLL; none left -> refused, "ceiling";
  7. the order row is INSERTED (status "sending", full cost) BEFORE it is
     sent, so a crash mid-send still counts against the ceiling;
  8. one limit buy at the ask (never a market order), then a cancel of the
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

Every order, cancel and listing names the subaccount, which a locked key must
match. An open order from another subaccount (the bot's first orders, placed
from the primary account on 2026-10-16) cannot be read or cancelled with the
locked key: its cost stays held until its market settles.

Side effects: Kalshi orders on the bot's key (money); INSERTs / UPDATEs
taker.duckdb::taker_orders (taker_store.py). Network: Kalshi REST, signed
with the bot's key only.
"""
import base64
import http.client
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

from unabated_ticket.bets_service.sources.kalshi import MAX_SUBACCOUNT, PRIMARY_SUBACCOUNT, record_id_of
from unabated_ticket.bets_service.taker_store import TakerStore

log = logging.getLogger(__name__)

# The ceiling on the bot's open cost and what the setup funds its subaccount with; the Kelly bankroll is taker.js KELLY_BANKROLL.
BOT_BANKROLL = 500.0
# Kalshi's fill feed can trail the match by a moment; a read this long after
# the confirmed cancel is taken as final. Until then the whole requested cost is held.
FILLS_SETTLE_SEC = 60
# Kalshi rounds each order's fee up to the cent: one cent of slack per order.
FEE_ROUNDING_SLACK = 0.01
# Kalshi's market status enum ends initialized/inactive/active/closed/determined/.../finalized.
SETTLED_MARKET_STATUSES = ("finalized",)
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
# An immediate-or-cancel order is final within a second; its listing can lag a little.
LOST_ORDER_LOOKUP_AFTER_SEC = 60
# Kalshi answered the order itself and refused it: nothing was placed. Any other
# answer (409 conflict, 429, a 5xx, a gateway timeout) may have reached the
# matching engine after all: the order is looked up by client_order_id instead.
DEFINITE_REJECTION_STATUSES = (400, 401, 403, 404, 422)
ORDER_LIST_STATUSES = ("resting", "executed", "canceled")
ORDER_LIST_PAGE_LIMIT = 1000
ORDER_LIST_MAX_PAGES = 10
MAX_HELD_IDS = 1000
MAX_HELD_ID_LENGTH = 200
# A decision older than this was priced on a book that may have moved (the
# runner's pass is 30 s; a request can wait behind the desk lock).
MAX_DECISION_AGE_SEC = 60
FILLS_PAGE_LIMIT = 1000
FILLS_MAX_PAGES = 10
BOT_ORDER_ID_PREFIX = "ubot-"
# Statuses under which the order may hold contracts (see _may_hold_contracts).
STATUSES_MAYBE_HOLDING = ("sending", "unknown", "sent", "resting", "filled", "partial")
BALANCE_PATH = "/portfolio/balance"
# The key-lock check at arming: a no-answer, 429 or 5xx is retried this many times, this far apart.
LOCK_CHECK_ATTEMPTS = 3
LOCK_CHECK_RETRY_SEC = 2.0


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
        except http.client.HTTPException as error:
            # A cut-off or malformed answer (IncompleteRead, BadStatusLine) is "no answer" to the callers.
            raise OSError(f"{type(error).__name__}: {error}") from error

    @property
    def base_url(self) -> str:
        return self._base_url


def _parse(raw: bytes) -> object:
    text = raw.decode(errors="replace")
    if not text:
        return {}
    try:
        return json.loads(text)
    except ValueError:
        return text


def order_body(ticker: str, contract_side: str, ask: float, count: int, client_order_id: str,
               subaccount: int) -> dict:
    """Kalshi v2's order body for BUYING `contract_side` at `ask`, immediate-
    or-cancel: what does not fill at once Kalshi cancels in the same request,
    so a bot order never rests (docs.kalshi.com, Create Order V2). v2 trades
    the YES leg only: buy YES = bid at the ask; buy NO = ask (sell YES) at
    1 - the NO ask (unabated_edge/maker/gateway.py, the same API live).
    `subaccount` is the bot's: a key locked to it may name only it."""
    if contract_side == "yes":
        v2_side, price = "bid", ask
    else:
        v2_side, price = "ask", 1 - ask
    return {"ticker": ticker, "side": v2_side, "price": f"{price:.4f}", "count": f"{float(count):.2f}",
            "time_in_force": "immediate_or_cancel", "self_trade_prevention_type": "taker_at_cross",
            "client_order_id": client_order_id, "subaccount": subaccount}


def fill_totals(fills: list[dict], contract_side: str) -> tuple[float, float] | str:
    """(contracts, dollars paid incl. fees) over one order's fills, priced on
    the bot's own contract side whichever leg Kalshi reports the fill on; or
    why a fill could not be read (a renamed or missing field must never read
    as a free fill)."""
    contracts = 0.0
    cost = 0.0
    for fill in fills:
        count = _positive_number(fill.get("count_fp", fill.get("count")))
        yes_price = _positive_number(fill.get("yes_price_dollars"))
        no_price = _positive_number(fill.get("no_price_dollars"))
        if no_price is None and yes_price is not None:
            no_price = 1 - yes_price
        fee = 0.0 if fill.get("fee_cost") in (None, "") else _number(fill.get("fee_cost"))
        if count is None or yes_price is None or no_price is None or fee is None:
            return f"fill {fill.get('fill_id') or fill.get('trade_id') or '?'} unreadable, fields {sorted(fill)}"
        price = yes_price if contract_side == "yes" else no_price
        contracts += count
        cost += count * price + fee
    return contracts, round(cost, 4)


def _number(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def _positive_number(value: object) -> float | None:
    number = _number(value)
    return number if number is not None and number > 0 else None


def held_id_of(ticker: str, contract_side: str, subaccount: int) -> str:
    """The bets service's record id for a Kalshi position (sources/kalshi.py)."""
    return record_id_of(ticker, contract_side, subaccount)


def _balance_dollars(reply: object) -> float | None:
    """GET /portfolio/balance's available balance (`balance`, in cents), in dollars."""
    cents = reply.get("balance") if isinstance(reply, dict) else None
    if isinstance(cents, bool) or not isinstance(cents, (int, float)):
        return None
    return cents / 100


def _retryable(status: int | None) -> bool:
    return status is None or status == 429 or status >= 500


def _get_status(client: "KalshiOrderClient", path: str) -> tuple[int | None, object]:
    """(status, reply), or (None, why) when Kalshi did not answer."""
    try:
        return client.request("GET", path)
    except OSError as error:
        return None, f"no answer: {error}"


def key_lock_problem(client: "KalshiOrderClient", subaccount: int,
                     sleep: Callable[[float], None] = time.sleep) -> tuple[str | None, float | None]:
    """(None, the subaccount's available dollars) when Kalshi shows the key is
    locked to `subaccount`, two ways that must agree: its balance reply there
    carries no `balance_breakdown` (docs.kalshi.com, Get Balance: "omitted only
    when using a subaccount-restricted API key"), and the primary account's
    balance is refused. Otherwise (why not, None). A no-answer, 429 or 5xx is
    retried LOCK_CHECK_ATTEMPTS times; never proven means not armed. GETs only."""
    last = "not checked"
    for attempt in range(LOCK_CHECK_ATTEMPTS):
        if attempt:
            sleep(LOCK_CHECK_RETRY_SEC)
        own_status, own_reply = _get_status(client, f"{BALANCE_PATH}?subaccount={subaccount}")
        if _retryable(own_status):
            last = f"its balance in subaccount {subaccount}: {own_status or own_reply}"
            continue
        if own_status != 200:
            return (f"the key cannot read subaccount {subaccount}: HTTP {own_status} {str(own_reply)[:200]}", None)
        if isinstance(own_reply, dict) and "balance_breakdown" in own_reply:
            return ("the key is not restricted: its balance reply carries balance_breakdown, which Kalshi omits "
                    "only for a subaccount-restricted key", None)
        primary_status, primary_reply = _get_status(client, f"{BALANCE_PATH}?subaccount={PRIMARY_SUBACCOUNT}")
        if _retryable(primary_status):
            last = f"the primary account's balance: {primary_status or primary_reply}"
            continue
        if primary_status == 200:
            return (f"the key reads the primary account too: it is not locked to subaccount {subaccount}", None)
        if 400 <= primary_status < 500:
            return None, _balance_dollars(own_reply)
        return f"the primary account's balance answered HTTP {primary_status}, expected a refusal", None
    return f"could not check after {LOCK_CHECK_ATTEMPTS} tries ({last})", None


def _order_id_of(reply: object) -> str | None:
    if not isinstance(reply, dict):
        return None
    order = reply.get("order") if isinstance(reply.get("order"), dict) else reply
    order_id = order.get("order_id")
    return order_id if isinstance(order_id, str) and order_id else None


def _remaining_count_of(reply: object) -> float | None:
    """Contracts Kalshi's create reply says are still open ("0.00" after an
    immediate-or-cancel), or None when the reply does not say."""
    return _reply_count(reply, "remaining_count_fp", "remaining_count")


def _fill_count_of(reply: object) -> float | None:
    """Contracts Kalshi's create reply says filled at once, or None."""
    return _reply_count(reply, "fill_count_fp", "fill_count")


def _reply_count(reply: object, *keys: str) -> float | None:
    if not isinstance(reply, dict):
        return None
    order = reply.get("order") if isinstance(reply.get("order"), dict) else reply
    for key in keys:
        if key in order:
            number = _number(order[key])
            return number if number is not None and number >= 0 else None
    return None


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
    start = _iso_time(body.get("eventStart"))
    if start is None:
        return f"eventStart must be an ISO time, got {body.get('eventStart')!r}"
    decided_at = _iso_time(body.get("decidedAt"))
    if decided_at is None:
        return f"decidedAt must be an ISO time, got {body.get('decidedAt')!r}"
    kalshi_base_url = body.get("kalshiBaseUrl")
    if not isinstance(kalshi_base_url, str) or not kalshi_base_url.startswith("http"):
        return f"kalshiBaseUrl must be the Kalshi API the decision was priced on, got {kalshi_base_url!r}"
    return {"client_order_id": client_order_id, "ticker": body["ticker"], "row_key": body["rowKey"],
            "contract_side": body["side"], "ask_price": float(ask), "all_in_price": float(all_in), "count": count,
            "event_id": event_id, "held_ids": set(held_ids), "event_start": start, "decided_at": decided_at,
            "kalshi_base_url": kalshi_base_url.rstrip("/")}


def _iso_time(value: object) -> datetime | None:
    if not isinstance(value, str):
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return parsed if parsed.tzinfo is not None else None


def _may_hold_contracts(order: dict) -> bool:
    """Could this open order hold contracts? Yes until a FINAL fills read says
    none did: a fill can trail the cancel in Kalshi's feed."""
    if order["status"] in STATUSES_MAYBE_HOLDING:
        return True
    return order["status"] == "unfilled" and not order["fills_final"]


def kalshi_unseen_fill(store: TakerStore, event_id: int, held_ids: set[str], ticker: str | None = None) -> str | None:
    """Why an add on Unabated event `event_id` (or on `ticker`) would be sized
    blind, or None: an open Kalshi bot order there that may hold contracts and
    whose position is not among `held_ids`. Store reads only, no desk lock, so
    the Novig desk asks it too."""
    for order in store.open_orders():
        same_game = order["event_id"] == event_id or (ticker is not None and order["ticker"] == ticker)
        if not same_game or not _may_hold_contracts(order):
            continue
        held_id = held_id_of(order["ticker"], order["contract_side"], order["subaccount"])
        if held_id not in held_ids:
            return (f"{order['client_order_id']} ({order['status']}) is not yet among the held bets the stake "
                    f"was sized on: no {held_id}")
    return None


class TakerDesk:
    """Places the runner's orders, one at a time, under the ceiling."""

    def __init__(self, client: KalshiOrderClient, store: TakerStore, stop_path: Path, subaccount: int,
                 now: Callable[[], datetime] = lambda: datetime.now(timezone.utc),
                 other_venue_unseen: Callable[[int, set[str]], str | None] | None = None):
        if isinstance(subaccount, bool) or not isinstance(subaccount, int) or not 1 <= subaccount <= MAX_SUBACCOUNT:
            raise ValueError(f"the bot trades only in a numbered subaccount 1-{MAX_SUBACCOUNT}, got {subaccount!r}")
        self._client = client
        self._subaccount = subaccount
        self._store = store
        self._stop_path = Path(stop_path)
        self._now = now
        # The Novig bot's orders on the same game (novig_taker_desk.novig_unseen_fill); None = Kalshi's alone.
        self._other_venue_unseen = other_venue_unseen
        self._lock = threading.Lock()
        self._last_refresh_at: datetime | None = None
        self._settlement_checked_at: dict[str, datetime] = {}

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
            refused = self._refuse_unfit(request)
            if refused:
                return refused
            if self._last_refresh_at is None or self._now() - self._last_refresh_at >= timedelta(seconds=REFRESH_SEC):
                self._refresh_open_orders()
                # The refresh waits on Kalshi, which can be slow: the decision may have aged out meanwhile.
                refused = self._refuse_unfit(request)
                if refused:
                    return refused
            unseen = self._unseen_fill(request)
            if unseen:
                return {"status": "unseen_fill", "detail": unseen}
            headroom = BOT_BANKROLL - self._store.open_exposure()
            count = min(request["count"], math.floor((headroom - FEE_ROUNDING_SLACK) / request["all_in_price"]))
            if count < 1:
                return {"status": "ceiling", "detail": f"${headroom:.2f} of ${BOT_BANKROLL:.0f} left"}
            return self._send(request, count)

    def _refuse_unfit(self, request: dict) -> dict | None:
        """A request the desk must not act on, whatever the ceiling says, or None."""
        now = self._now()
        if request["kalshi_base_url"] != self._client.base_url:
            return {"status": "env_mismatch",
                    "detail": f"the runner priced {request['kalshi_base_url']}, this desk trades {self._client.base_url}"}
        age = now - request["decided_at"]
        if age > timedelta(seconds=MAX_DECISION_AGE_SEC):
            return {"status": "stale_decision", "detail": f"decided {age.total_seconds():.0f} s ago, over {MAX_DECISION_AGE_SEC} s"}
        if request["event_start"] <= now:
            return {"status": "started", "detail": f"the game started at {request['event_start'].isoformat()}"}
        return None

    def _unseen_fill(self, request: dict) -> str | None:
        """Why this order would add to a game blind, or None: an earlier bot
        order on the game, on Kalshi or (other_venue_unseen) on Novig, that
        may hold contracts and is not in `held_ids`."""
        unseen = kalshi_unseen_fill(self._store, request["event_id"], request["held_ids"], request["ticker"])
        if unseen is None and self._other_venue_unseen is not None:
            unseen = self._other_venue_unseen(request["event_id"], request["held_ids"])
        return unseen

    def _send(self, request: dict, count: int) -> dict:
        client_order_id = request["client_order_id"]
        reserved = round(count * request["all_in_price"] + FEE_ROUNDING_SLACK, 4)
        self._store.insert_order({**request, "placed_at": self._now(), "count_requested": count, "open_cost": reserved,
                                  "subaccount": self._subaccount})
        body = order_body(request["ticker"], request["contract_side"], request["ask_price"], count, client_order_id,
                          self._subaccount)
        try:
            status, reply = self._client.request("POST", ORDERS_PATH, body)
        except OSError as error:
            # It may have reached Kalshi: keep the full cost held; the refresh looks it up by client_order_id.
            self._store.update_order(client_order_id, {"status": "unknown", "error": f"no answer: {error}"})
            log.error("taker order %s %s: no answer from Kalshi: %s", client_order_id, request["ticker"], error)
            return {"status": "unknown", "count": count, "detail": f"no answer from Kalshi: {error}"}
        order_id = _order_id_of(reply)
        if status in DEFINITE_REJECTION_STATUSES:
            # Nothing rests or filled: no cost, and no market to watch for settlement.
            self._store.update_order(client_order_id, {"status": "rejected", "open_cost": 0.0, "settled_at": self._now(),
                                                       "error": f"HTTP {status}: {str(reply)[:500]}"})
            log.warning("taker order %s %s refused by Kalshi: HTTP %s %s", client_order_id, request["ticker"], status, reply)
            return {"status": "rejected", "count": count, "detail": f"Kalshi HTTP {status}"}
        if status not in (200, 201) or not order_id:
            # A gateway error or a conflict may hide a placed order: full cost held, looked up by client_order_id.
            self._store.update_order(client_order_id, {"status": "unknown",
                                                       "error": f"no definite answer: HTTP {status}: {str(reply)[:300]}"})
            log.error("taker order %s %s: no definite answer from Kalshi: HTTP %s %s", client_order_id, request["ticker"], status, reply)
            return {"status": "unknown", "count": count, "detail": f"no definite answer from Kalshi: HTTP {status}; full cost held"}
        exchange_index = _exchange_index_of(reply)
        self._store.update_order(client_order_id, {"status": "sent", "kalshi_order_id": order_id, "exchange_index": exchange_index,
                                                   "count_reported": _fill_count_of(reply)})
        order = {"client_order_id": client_order_id, "kalshi_order_id": order_id, "ticker": request["ticker"],
                 "contract_side": request["contract_side"], "exchange_index": exchange_index, "count_requested": count}
        remaining = _remaining_count_of(reply)
        if remaining is not None and remaining <= 0:
            # Immediate-or-cancel: Kalshi cancelled the rest in the same request, nothing rests.
            self._store.update_order(client_order_id, {"cancelled_at": self._now(), "error": None})
            return self._read_fills(client_order_id, order_id, request["contract_side"], count, final=False)
        # No count in the reply, or contracts left open against the order's terms: cancel as if it rested.
        log.warning("taker order %s %s: Kalshi reports %s contracts still open after an IOC send; cancelling",
                    client_order_id, request["ticker"], remaining)
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
        path = f"{ORDERS_PATH}/{order_id}?market_ticker={order['ticker']}&subaccount={self._subaccount}"
        if order.get("exchange_index") is not None:
            path += f"&exchange_index={int(order['exchange_index'])}"
        try:
            status, reply = self._client.request("DELETE", path)
        except OSError as error:
            log.error("taker cancel %s: no answer: %s", order_id, error)
            return False
        if status in (200, 204):
            return True
        if 400 <= status < 500:
            # Terminal (filled, expired), misrouted, or refused: only Kalshi's resting list tells.
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

    def _order_fills(self, order_id: str) -> list[dict] | str:
        """Every fill of one order (paged), or why the read failed."""
        fills: list[dict] = []
        cursor = ""
        for _ in range(FILLS_MAX_PAGES):
            path = (f"/portfolio/fills?order_id={order_id}&subaccount={self._subaccount}&limit={FILLS_PAGE_LIMIT}"
                    + (f"&cursor={cursor}" if cursor else ""))
            try:
                status, reply = self._client.request("GET", path)
            except OSError as error:
                return f"no answer: {error}"
            if status != 200 or not isinstance(reply, dict) or not isinstance(reply.get("fills"), list):
                return f"HTTP {status} {str(reply)[:300]}"
            fills.extend(fill for fill in reply["fills"] if isinstance(fill, dict))
            cursor = reply.get("cursor") or ""
            if not cursor:
                return fills
        return f"more than {FILLS_MAX_PAGES} pages of fills"

    def _read_fills(self, client_order_id: str, order_id: str, contract_side: str, count: int, final: bool,
                    count_reported: float | None = None) -> dict:
        """Read the order's fills into its row. `final` fixes open_cost to what
        filled; refused while the fills feed shows fewer contracts than Kalshi's
        create reply reported (`count_reported`), since the feed can lag."""
        fills = self._order_fills(order_id)
        totals = fills if isinstance(fills, str) else fill_totals(fills, contract_side)
        if isinstance(totals, str):
            self._store.update_order(client_order_id, {"status": "unknown", "error": f"fills read failed: {totals}"})
            return {"status": "unknown", "count": count, "orderId": order_id, "detail": "fills read failed; full cost held"}
        filled, cost = totals
        if final and count_reported is not None and filled < count_reported:
            log.warning("taker order %s: the fills feed shows %s of the %s contracts Kalshi reported filled; not final yet",
                        client_order_id, filled, count_reported)
            final = False
        outcome = "unfilled" if filled == 0 else "filled" if filled >= count else "partial"
        fields = {"status": outcome, "count_filled": filled, "cost_filled": cost, "fills_final": final, "error": None}
        if final:
            fields["open_cost"] = cost
        self._store.update_order(client_order_id, fields)
        return {"status": outcome, "count": count, "filled": filled, "cost": cost, "orderId": order_id}

    def _list_orders(self, query: str) -> list[dict] | None:
        """GET /portfolio/orders?<query> in the bot's subaccount, every page;
        None when the read failed or the listing did not end (never a partial list)."""
        orders: list[dict] = []
        cursor = ""
        for _ in range(ORDER_LIST_MAX_PAGES):
            path = (f"{ORDER_LIST_PATH}?{query}&subaccount={self._subaccount}&limit={ORDER_LIST_PAGE_LIMIT}"
                    + (f"&cursor={cursor}" if cursor else ""))
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

    def _orders_for_ticker(self, ticker: str) -> list[dict] | None:
        """Every order on the ticker in every status Kalshi lists (what an
        unfiltered listing returns is not documented), or None when any read failed."""
        orders: list[dict] = []
        for status in ORDER_LIST_STATUSES:
            listed = self._list_orders(f"ticker={ticker}&status={status}")
            if listed is None:
                return None
            orders.extend(listed)
        return orders

    def _refresh_open_orders(self) -> None:
        """Bring every open order up to date, fail closed (a failed read holds
        the cost): a bot order resting on Kalshi that no row knows is
        cancelled; a send with no answer is looked up by client_order_id; an
        unconfirmed cancel is retried; fills are final FILLS_SETTLE_SEC after
        the cancel; a started market is checked for settlement, which frees it.
        An order from another subaccount is only checked for settlement: the
        locked key can neither read nor cancel it, so its cost stays held."""
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
            if order["subaccount"] != self._subaccount:
                self._free_if_settled(order, now)
                continue
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
                                 order["count_requested"], final=True, count_reported=order["count_reported"])
            self._free_if_settled(order, now)

    def _free_if_settled(self, order: dict, now: datetime) -> None:
        """Marks a started order's market settled once Kalshi says so (a public
        market read, any subaccount): that frees its open_cost."""
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
            listed = self._orders_for_ticker(order["ticker"])
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
                                                   "exchange_index": adopted["exchange_index"],
                                                   "count_reported": _fill_count_of(match)})
        log.warning("taker order %s %s: found on Kalshi as %s (%s) after a lost answer",
                    client_order_id, order["ticker"], match["order_id"], match.get("status"))
        if match.get("status") in ("executed", "canceled"):
            # Done on Kalshi's side (an immediate-or-cancel always is): nothing to cancel, only fills to read.
            self._store.update_order(client_order_id, {"cancelled_at": now, "error": None})
            self._read_fills(client_order_id, match["order_id"], order["contract_side"], order["count_requested"], final=False)
            return
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
