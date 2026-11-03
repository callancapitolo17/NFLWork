"""ProphetX bet source: the account's orders and fills -> normalised bet records,
read on the account's own Trading API key (PROPHETX_ACCESS_KEY / _SECRET_KEY).

Replaces the content-script route of the unmerged #117 branch
(feature/unabated-prophetx-bets, 2026-09-16, which read the web app's wager
list while a ProphetX tab was open); its money and grading rules are kept
here: only the MATCHED dollars of an order are a bet, a resting remainder is a
separate unmatchable record, and a settled order grades off winning_status.

Two halves, like sources/polymarket_us.py:
  normalize_prophetx(...)  pure: orders + their fills + the games and markets
                           they name -> records. Pinned by
                           tests/test_prophetx.py.
  ProphetXSource           the network half (Source protocol): every poll reads
                           the order history and fills back PROPHETX_HISTORY_DAYS
                           (at most 60 days, ProphetX's window), older orders once
                           per boot in 60-day slices, and once per game its sport
                           event and markets (cached in memory).

Wire grammar (docs.prophetx.co v4 schemas, 2026-10-31; not a live capture):
  order   {order_id, external_id, sport_event_id, market_id, outcome_id,
          strike_id, strike, price (American), quantity / original_quantity
          (stake dollars), filled_quantity, open_quantity, status,
          matching_status, winning_status, created_at, updated_at, settled_at}
  trade   {order_id, strike_id, price, quantity (dollars matched), profit,
          fee, matched_at, status, winning_status}: one per fill. A fill can be
          at a better price than the order's limit, so the record's price is
          the fills' own: stake = the matched dollars plus fees, toWin = what
          they win at each fill's price less fees, price = the American price
          of that ratio. A voided or invalid fill matched nothing.
  status  void (voided by ProphetX) and invalid (never accepted) void the order;
          canceled / wiped end a resting remainder but keep what matched;
          settled / closed / manually_settled grade off winning_status
          (profit won, loss lost, push / draw push, no_result void, the manual_
          variants alike, tbd still open).
Games and markets come from GET /mm/get_sport_events and /v4/mm/get_markets
(prophetx_markets.py reads them). ProphetX drops a game from the event list once
it ends, so an order whose game is no longer listed keeps the game fields its
stored record already had (`known_records`, read from bets.duckdb); one never
resolved is unmatchable with the reason.

Outputs: records (id "prophetx:<order_id>", resting part
"prophetx:<order_id>:order"). Side effects: none on disk; network reads only.
Raises on a failed order-history or trades read (the service keeps the previous
records); a failed game or market read only leaves that game unresolved for one poll.
"""
import logging
import time
from collections.abc import Callable
from datetime import datetime, timezone

from unabated_ticket.bets_service import config
from unabated_ticket.bets_service.normalize import EASTERN, json_clean, round_cents, utc_now_iso
from unabated_ticket.bets_service.prophetx_client import ProphetXClient, ProphetXError
from unabated_ticket.bets_service.prophetx_markets import league_of, lines_of_event, number_of, profit_at
from unabated_ticket.bets_service.sources.betonline import american_from_payout

log = logging.getLogger(__name__)

SOURCE = "prophetx_api"
VENUE = "prophetx"
APPROX_ORDER = "prophetx_unmatched_order"
VOID_STATUSES = {"void", "invalid"}
RESTING_STATUSES = {"open", "inactive"}
SETTLED_STATUSES = {"settled", "closed", "manually_settled"}
GRADES = {"profit": "won", "manually_won": "won", "loss": "lost", "manually_lost": "lost", "push": "push",
          "draw": "push", "no_result": "void"}
# Game fields a record carries; a stored copy of them is reused once ProphetX stops listing the game.
GAME_FIELDS = ("league", "eventStart", "eventDate", "awayTeam", "homeTeam", "betType", "period", "side", "points")
# A game's markets are read again after this long (its alt numbers change; the ones an order bought do not).
MARKETS_TTL_SEC = 600.0
# A game ProphetX's event list no longer returns (it drops ended games) is asked for again after this long.
UNLISTED_RETRY_SEC = 3600.0
# Once per boot, orders older than the polled window are read in 60-day slices back
# until this many slices in a row hold no order (the account's start), capped at
# BACKFILL_MAX_WINDOWS slices (~3 years).
BACKFILL_EMPTY_WINDOWS_TO_STOP = 2
BACKFILL_MAX_WINDOWS = 18


def iso_utc_of(value: object) -> str | None:
    """ProphetX's timestamp (ISO text or unix seconds) as ISO UTC, or None."""
    if isinstance(value, (int, float)) and not isinstance(value, bool) and value > 0:
        return datetime.fromtimestamp(value / 1000 if value > 1e12 else value, timezone.utc).isoformat().replace("+00:00", "Z")
    if isinstance(value, str) and value.strip():
        try:
            parsed = datetime.fromisoformat(value.strip().replace("Z", "+00:00"))
        except ValueError:
            return None
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=timezone.utc)
        return parsed.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")
    return None


def eastern_date_of(iso: str | None) -> str | None:
    if not iso:
        return None
    return datetime.fromisoformat(iso.replace("Z", "+00:00")).astimezone(EASTERN).date().isoformat()


def matched_of(order: dict, fills: list[dict]) -> tuple[float, float | None, float]:
    """(matched dollars, what they win or None, fees) of one order: off its
    fills when it has any, else off filled_quantity at the order's price."""
    # A fill ProphetX voided or never accepted matched nothing.
    fills = [fill for fill in fills if str(fill.get("status") or "").lower() not in VOID_STATUSES]
    fees = sum(number_of(fill.get("fee")) or 0.0 for fill in fills)
    if fills:
        stake = 0.0
        to_win = 0.0
        for fill in fills:
            quantity = number_of(fill.get("quantity")) or 0.0
            price = number_of(fill.get("price"))
            if quantity <= 0:
                continue
            if price is None or abs(price) < 100:
                return round_cents(sum(number_of(f.get("quantity")) or 0.0 for f in fills)), None, fees
            stake += quantity
            to_win += profit_at(quantity, price)
        return round_cents(stake), round_cents(to_win), fees
    filled = number_of(order.get("filled_quantity")) or 0.0
    price = number_of(order.get("price"))
    if filled <= 0:
        return 0.0, 0.0, fees
    return round_cents(filled), None if price is None or abs(price) < 100 else round_cents(profit_at(filled, price)), fees


def status_of(order: dict, matched: float) -> str:
    """The status of an order's MATCHED part."""
    status = str(order.get("status") or "").lower()
    if matched <= 0 or status in VOID_STATUSES:
        return "void"
    grade = GRADES.get(str(order.get("winning_status") or "").lower())
    if grade is not None:
        return grade
    if status in SETTLED_STATUSES:
        return "unknown"
    return "open"


def resting_of(order: dict) -> float:
    """Dollars still resting on the book, while the order is."""
    if str(order.get("status") or "").lower() not in RESTING_STATUSES:
        return 0.0
    return round_cents(max(number_of(order.get("open_quantity")) or 0.0, 0.0))


def game_fields_of(order: dict, event: dict | None, lines: list[dict] | str | None) -> dict | str:
    """The order's board line, or why it has none."""
    if event is None or lines is None:
        return f"ProphetX does not list game {order.get('sport_event_id')} any more"
    if isinstance(lines, str):
        return lines
    found = [line for line in lines if line["strike_id"] == order.get("strike_id")]
    if len(found) != 1:
        return (f"strike {order.get('strike_id')} reads as {len(found)} lines on ProphetX game "
                f"{order.get('sport_event_id')} (not a game market, or a market that does not read)")
    line = found[0]
    start = iso_utc_of(event.get("scheduled"))
    return {"league": line["league"], "eventStart": start, "eventDate": eastern_date_of(start),
            "awayTeam": line["awayTeam"], "homeTeam": line["homeTeam"], "betType": line["betType"],
            "period": line["period"], "side": line["side"], "points": line["points"]}


def _unmatchable(record: dict, reason: str) -> dict:
    record.update({field: None for field in GAME_FIELDS})
    record["betType"] = "other"
    record["unmatchable"] = reason
    return record


def _stored_game(stored: dict | None) -> dict | None:
    """A stored record's game fields, when it had a line."""
    if not stored or stored.get("unmatchable") or stored.get("betType") not in ("moneyline", "spread", "total"):
        return None
    return {field: stored.get(field) for field in GAME_FIELDS}


def normalize_prophetx(orders: list[dict], trades: list[dict], events: dict[int, dict],
                       lines_by_event: dict[int, list[dict] | str], known_records: dict[str, dict],
                       fetched_at: str | None) -> list[dict]:
    """The records of every order. `events` / `lines_by_event` are keyed by
    ProphetX's sport_event_id (lines_of_event's answer); `known_records` are
    the stored records by id. Pure."""
    fills_by_order: dict[str, list[dict]] = {}
    for trade in trades:
        if trade.get("order_id"):
            fills_by_order.setdefault(str(trade["order_id"]), []).append(trade)
    records = []
    for order in orders:
        order_id = order.get("order_id")
        if not order_id:
            continue
        fills = fills_by_order.get(str(order_id), [])
        matched, to_win, fees = matched_of(order, fills)
        resting = resting_of(order)
        event_id = order.get("sport_event_id")
        game = game_fields_of(order, events.get(event_id), lines_by_event.get(event_id))
        placed_at = iso_utc_of(order.get("created_at"))
        event = events.get(event_id)
        raw = {"orderId": order_id, "externalId": order.get("external_id"), "sportEventId": event_id,
               "league": league_of(event) if event else None,
               "marketId": order.get("market_id"), "strikeId": order.get("strike_id"), "strike": order.get("strike"),
               "limitPrice": order.get("price"), "requestedStake": number_of(order.get("original_quantity") or order.get("quantity")),
               "matchedStake": matched, "restingStake": resting, "fees": round(fees, 4), "fillCount": len(fills),
               "orderStatus": order.get("status"), "matchingStatus": order.get("matching_status"),
               "winningStatus": order.get("winning_status")}
        base = {"source": SOURCE, "venue": VENUE, "awayKey": None, "homeKey": None, "rotation": None, "contracts": None,
                "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None, "sourceFetchedAt": fetched_at,
                "placedAt": placed_at, "venueIds": {"strikeId": order.get("strike_id"), "sportEventId": event_id},
                "raw": raw}
        if matched > 0:
            status = status_of(order, matched)
            # Fees count in the money like Polymarket US's: stake is what the bet cost, toWin what it nets.
            stake = round_cents(matched + fees)
            net_to_win = round_cents(to_win - fees) if to_win is not None else None
            record = {"id": f"{VENUE}:{order_id}", **base,
                      "price": american_from_payout(stake, net_to_win) if net_to_win else None,
                      "stake": stake, "toWin": net_to_win, "status": status, "approx": [], "unmatchable": None,
                      "closedAt": None if status == "open" else iso_utc_of(order.get("settled_at")) or iso_utc_of(order.get("updated_at"))}
            records.append(json_clean(_with_game(record, game, known_records)))
        if resting > 0:
            price = number_of(order.get("price"))
            record = {"id": f"{VENUE}:{order_id}:order", **base, "price": int(price) if price else None,
                      "stake": resting, "toWin": round_cents(profit_at(resting, price)) if price and abs(price) >= 100 else None,
                      "status": "open", "closedAt": None, "approx": [APPROX_ORDER]}
            record = _with_game(record, game, known_records)
            label = f" at {'+' if price and price > 0 else ''}{int(price)}" if price else ""
            record["unmatchable"] = record.get("unmatchable") or f"unmatched order (resting ${resting}{label})"
            records.append(json_clean(record))
        elif _was_resting(known_records.get(f"{VENUE}:{order_id}:order")):
            # The remainder stopped resting (matched, cancelled or wiped): its
            # stored open copy is closed as void, never left open.
            ended = dict(known_records[f"{VENUE}:{order_id}:order"])
            ended.update({"status": "void", "stake": 0.0, "toWin": 0.0, "sourceFetchedAt": fetched_at,
                          "closedAt": iso_utc_of(order.get("updated_at")) or fetched_at, "raw": raw})
            records.append(json_clean(ended))
    return records


def _was_resting(stored: dict | None) -> bool:
    return isinstance(stored, dict) and stored.get("status") == "open"


def _with_game(record: dict, game: dict | str, known_records: dict[str, dict]) -> dict:
    if isinstance(game, dict):
        record.update(game)
        record["unmatchable"] = None
        return record
    stored = _stored_game(known_records.get(record["id"]))
    if stored is not None:
        record.update(stored)
        record["unmatchable"] = None
        return record
    return _unmatchable(record, game)


def _has_strikes(lines: list[dict] | str, strikes: set[str]) -> bool:
    """Every strike an order bought is among a game's cached lines."""
    return isinstance(lines, list) and strikes <= {line.get("strike_id") for line in lines}


# ---- network half -----------------------------------------------------------------


class ProphetXSource:
    """Polls the account's ProphetX orders on its Trading API key."""

    name = VENUE

    def __init__(self, client: ProphetXClient, known_records: Callable[[], dict[str, dict]],
                 poll_sec: float = config.PROPHETX_POLL_SEC, history_days: int = config.PROPHETX_HISTORY_DAYS,
                 clock: Callable[[], float] = time.time):
        """`known_records`: the stored ProphetX records by id (store.load_records_of_venue)."""
        self.poll_sec = poll_sec
        self._client = client
        self._known_records = known_records
        self._history_days = history_days
        self._clock = clock
        self._events: dict[int, dict] = {}
        self._lines: dict[int, tuple[float, list[dict] | str]] = {}
        self._unlisted: dict[int, float] = {}
        # Orders older than the polled window, read once per boot (backfill_older).
        self._older: tuple[list[dict], list[dict]] | None = None

    @property
    def client(self) -> ProphetXClient:
        """The source's client: the Place desk trades on the same key pair and session."""
        return self._client

    def _read_games(self, strikes_by_event: dict[int, set[str]]) -> None:
        """Read the games not yet cached (one event-list read), and a game's
        markets when none are cached or an order's strike is not among them
        (a new alt number), at most once per MARKETS_TTL_SEC. A game ProphetX
        no longer lists is not asked for again for UNLISTED_RETRY_SEC. Caches
        keep only this poll's games. A failed read leaves a game unresolved."""
        now = self._clock()
        event_ids = set(strikes_by_event)
        for cache in (self._events, self._lines, self._unlisted):
            for event_id in [key for key in cache if key not in event_ids]:
                del cache[event_id]
        missing = sorted(event_id for event_id in event_ids if event_id not in self._events
                         and now - self._unlisted.get(event_id, -UNLISTED_RETRY_SEC) >= UNLISTED_RETRY_SEC)
        if missing:
            try:
                events = self._client.sport_events(missing)
                if not events and len(missing) > 1:
                    # Whether one ended game makes the whole list 404 is unverified:
                    # ask each alone (each answer is then remembered for an hour).
                    events = [event for event_id in missing for event in self._client.sport_events([event_id])]
                for event in events:
                    if event.get("event_id") is not None:
                        self._events[event["event_id"]] = event
                for event_id in missing:
                    if event_id not in self._events:
                        self._unlisted[event_id] = now
            except (OSError, ProphetXError) as error:
                log.warning("prophetx: sport events %s unread: %s", missing[:5], error)
        for event_id, strikes in strikes_by_event.items():
            event = self._events.get(event_id)
            cached = self._lines.get(event_id)
            if event is None or (cached and (now - cached[0] < MARKETS_TTL_SEC or _has_strikes(cached[1], strikes))):
                continue
            try:
                self._lines[event_id] = (now, lines_of_event(event, self._client.markets(event_id)))
            except (OSError, ProphetXError) as error:
                log.warning("prophetx: markets of game %s unread: %s", event_id, error)

    def fetch(self) -> list[dict]:
        now = int(self._clock())
        start = now - self._history_days * 86400
        if self._older is None:
            self._older = self._backfill_older(start)
        orders = self._client.order_history(start, now) + self._older[0]
        trades = self._client.trades(start, now) + self._older[1]
        strikes_by_event: dict[int, set[str]] = {}
        for order in orders:
            if isinstance(order.get("sport_event_id"), int):
                strikes = strikes_by_event.setdefault(order["sport_event_id"], set())
                if order.get("strike_id"):
                    strikes.add(str(order["strike_id"]))
        self._read_games(strikes_by_event)
        lines = {event_id: entry[1] for event_id, entry in self._lines.items()}
        return normalize_prophetx(orders, trades, self._events, lines, self._known_records(), utc_now_iso())

    def _backfill_older(self, window_end: int) -> tuple[list[dict], list[dict]]:
        """Every order (and its fills) before the polled window, read once per boot in
        ProphetX's 60-day slices back until BACKFILL_EMPTY_WINDOWS_TO_STOP slices in a
        row hold no order, at most BACKFILL_MAX_WINDOWS. These are settled bets the
        store already holds after the first boot; re-reading them each boot keeps
        their stored copy current."""
        orders: list[dict] = []
        trades: list[dict] = []
        empty_in_a_row = 0
        for _ in range(BACKFILL_MAX_WINDOWS):
            window_start = window_end - self._history_days * 86400
            window_orders = self._client.order_history(window_start, window_end - 1)
            # Trades are read by match time, and a fill can match a slice after its
            # order was placed, so every slice scanned has its trades read too.
            trades += self._client.trades(window_start, window_end - 1)
            if window_orders:
                orders += window_orders
                empty_in_a_row = 0
            else:
                empty_in_a_row += 1
                if empty_in_a_row >= BACKFILL_EMPTY_WINDOWS_TO_STOP:
                    break
            window_end = window_start
        log.info("prophetx: backfilled %d order(s) older than %d days", len(orders), self._history_days)
        return orders, trades

    def balance(self) -> dict:
        """The balances.py reading off GET /v4/mm/get_balance: `balance` is the
        cash free to bet, `matched_order_balance` the money in matched bets
        (field meanings read off their names; unverified until a live read)."""
        status, reply = self._client.balance()
        data = reply.get("data") if status == 200 and isinstance(reply, dict) else None
        cash = data.get("balance") if isinstance(data, dict) else None
        if isinstance(cash, bool) or not isinstance(cash, (int, float)):
            raise ProphetXError(f"get_balance: HTTP {status}, no balance number in {str(reply)[:200]}")
        in_bets = data.get("matched_order_balance")
        in_bets = float(in_bets) if isinstance(in_bets, (int, float)) and not isinstance(in_bets, bool) else None
        return {"cash": float(cash), "inBets": in_bets, "available": float(cash), "creditLimit": None}


def source_if_configured(known_records: Callable[[], dict[str, dict]] = dict) -> ProphetXSource | None:
    """The source when the key pair is set (the Polymarket US pattern). The
    service passes the store's stored ProphetX records; check_sources none."""
    client = client_if_configured()
    return ProphetXSource(client, known_records) if client is not None else None


def client_if_configured() -> ProphetXClient | None:
    """The account's ProphetX client when its key pair is set, else None (logged once)."""
    if not (config.PROPHETX_ACCESS_KEY and config.PROPHETX_SECRET_KEY):
        log.warning("prophetx: no PROPHETX_ACCESS_KEY / PROPHETX_SECRET_KEY — ProphetX source and Place not registered")
        return None
    return ProphetXClient(config.PROPHETX_ACCESS_KEY, config.PROPHETX_SECRET_KEY, config.PROPHETX_BASE_URL)
