"""Novig subaccount bet source (automated_taker_plan.md § Novig, 2026-10-16):
the bets of ONE Novig subaccount, read on its own trading key, kept apart from
Cal's app account (sources/novig.py) and from the other subaccount.

  "bot"    the automated taker's subaccount (NOVIG_TAKER_KEY_ID), source "novig_bot"
  "place"  the Ticket's Novig Place button's subaccount (NOVIG_PLACE_KEY_ID), source "novig_place"

Why a source of its own: the app's Portfolio feed reads Cal's cash wallet on
his Auth0 login; a subaccount is a separate wallet only its trading key reads.
Every bot bet must be held (sized against, and the unseen-fill rule's "seen")
the moment its fills land, so the bot subaccount is read every POLL_SEC.

    GET /v3/portfolio/fills?startsAfter=<ms>   (signed, NOVIG-V3, paged)
        -> one record per order: its fills summed
    GET /v3/public/catalog/markets/<id>        (public, once a game has started)
        -> settled outcome status: WIN / LOSS / PUSH / a decimal (fair-value settlement)

Novig's catalog names neither the board's teams nor its periods, so a bet's
game fields (league, teams with Unabated's team ids, market, period, side,
number) come from `contexts()`: the board line the order was placed on, by
Novig order id. For "bot" that is taker.duckdb::novig_taker_orders
(TakerStore.novig_order_contexts); an order with no context (placed outside the
desk) is recorded unmatchable, so it shows but sizes nothing.

Ids: "novig:<order id>:<name>", venue "novig" (bets.js joins a spread or total
on venueIds.outcomeId and the side on its own fields). A Novig contract pays
1c: qty 100 = $1; `contracts` here is dollars of payout like sources/novig.py.

Side effects: none on disk. Raises on any failed read, so the service records
a failed run and keeps the previous records.
"""
import logging
from collections.abc import Callable
from datetime import datetime, timedelta, timezone
from zoneinfo import ZoneInfo

from unabated_ticket.bets_service.normalize import js_round, json_clean, utc_now_iso
from unabated_ticket.bets_service.novig_orders import QTY_PER_DOLLAR, NovigError, NovigOrderClient, order_totals
from unabated_ticket.bets_service.sources.novig import probability_to_american

log = logging.getLogger(__name__)

NAMES = ("bot", "place")
SOURCE_PREFIX = "novig_"
VENUE = "novig"
EASTERN = ZoneInfo("America/New_York")
POLL_SEC = 30.0
# Fills on events that started within this many days are read every poll (the store keeps older ones).
HISTORY_DAYS = 14
MARKET_SETTLED = "SETTLED"
OUTCOME_STATUS = {"WIN": "won", "LOSS": "lost", "PUSH": "push"}
OUTCOME_OPEN = "TBD"
REASON_NO_CONTEXT = "no board line on record for this Novig order (placed outside the bot's desk)"

Contexts = Callable[[], dict[str, dict]]


def _iso_ms(value: datetime | None) -> str | None:
    if value is None:
        return None
    value = value.astimezone(timezone.utc)
    return value.strftime("%Y-%m-%dT%H:%M:%S.") + f"{value.microsecond // 1000:03d}Z"


def _iso_of_epoch_ms(epoch_ms: object) -> str | None:
    if isinstance(epoch_ms, bool) or not isinstance(epoch_ms, int):
        return None
    return _iso_ms(datetime.fromtimestamp(epoch_ms / 1000, tz=timezone.utc))


def _round_cents(dollars: float) -> float:
    return js_round(dollars * 100) / 100


def settled_status(outcome_status: object, stake: float, qty: int) -> str | None:
    """won / lost / push for a settled outcome status, or None while it is
    open (TBD) or unreadable. A decimal is a fair-value settlement, graded on
    what it paid against the stake, like sources/novig.py's "Settled"."""
    if outcome_status in OUTCOME_STATUS:
        return OUTCOME_STATUS[outcome_status]
    if not isinstance(outcome_status, str) or outcome_status == OUTCOME_OPEN:
        return None
    try:
        fraction = float(outcome_status)
    except ValueError:
        return None
    paid = qty / QTY_PER_DOLLAR * fraction
    if abs(paid - stake) < 0.005:
        return "push"
    return "won" if paid > stake else "lost"


def _venue_team(name: object, unabated_id: object) -> dict:
    return {"id": None, "name": name if isinstance(name, str) else None, "shortName": None, "symbol": None,
            "unabatedId": unabated_id if isinstance(unabated_id, str) else None}


def _game_fields(context: dict | None) -> dict:
    if context is None:
        return {"league": None, "eventStart": None, "eventDate": None, "awayTeam": None, "homeTeam": None,
                "awayKey": None, "homeKey": None, "betType": "other", "period": None, "side": None, "points": None,
                "unmatchable": REASON_NO_CONTEXT, "awayTeamVenue": None, "homeTeamVenue": None}
    start = context.get("event_start")
    return {
        "league": context["league"], "eventStart": _iso_ms(start),
        "eventDate": start.astimezone(EASTERN).strftime("%Y-%m-%d") if start else None,
        "awayTeam": context.get("away_team"), "homeTeam": context.get("home_team"), "awayKey": None, "homeKey": None,
        "betType": context["bet_type"], "period": context["period"], "side": context["side"],
        "points": context.get("points"), "unmatchable": None,
        "awayTeamVenue": _venue_team(context.get("away_team"), context.get("away_team_id")),
        "homeTeamVenue": _venue_team(context.get("home_team"), context.get("home_team_id")),
    }


def normalize_order(name: str, order_id: str, fills: list[dict], context: dict | None, market: dict | None,
                    fetched_at: str) -> dict:
    """One order's fills -> one record. `market` is the public market once its game has started (else None)."""
    qty, cost, fee = order_totals(fills)
    stake = round(cost + fee, 5)
    first = fills[0]
    probability = cost / qty * QTY_PER_DOLLAR if qty else None
    payout = qty / QTY_PER_DOLLAR
    status = "open"
    outcome_status = None
    if market is not None and market.get("status") == MARKET_SETTLED:
        outcome = next((item for item in market.get("outcomes") or [] if item.get("outcomeId") == first.get("outcomeId")), None)
        outcome_status = outcome.get("status") if outcome else None
        graded = settled_status(outcome_status, stake, qty)
        if graded is None:
            raise NovigError(f"market {market.get('marketId')} is SETTLED but outcome {first.get('outcomeId')} reads "
                             f"{outcome_status!r}")
        status = graded
    placed_ms = min(fill.get("ts") for fill in fills if isinstance(fill.get("ts"), int)) if any(
        isinstance(fill.get("ts"), int) for fill in fills) else None
    record = {
        "id": f"novig:{order_id}:{name}",
        "source": f"{SOURCE_PREFIX}{name}",
        "venue": VENUE,
        "rotation": None,
        "price": probability_to_american(probability),
        "stake": _round_cents(stake),
        "toWin": _round_cents(payout - stake),
        "contracts": _round_cents(payout),
        "placedAt": _iso_of_epoch_ms(placed_ms),
        "status": status,
        "closedAt": None if status == "open" else fetched_at,
        "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
        "approx": [],
        "sourceFetchedAt": fetched_at,
        "venueIds": {"marketId": first.get("marketId"), "outcomeId": first.get("outcomeId"),
                     "eventId": market.get("eventId") if market else None, "gameId": None},
        "raw": {"orderId": order_id, "subaccount": name, "qty": qty, "cost": cost, "fee": fee,
                "probability": probability, "qtyUnit": QTY_PER_DOLLAR, "outcomeStatus": outcome_status,
                "fillIds": [fill.get("fillId") for fill in fills]},
    }
    record.update(_game_fields(context))
    return json_clean(record)


class NovigSubaccountSource:
    """Source protocol for one Novig subaccount (signed GETs on its trading key)."""

    def __init__(self, name: str, client: NovigOrderClient, contexts: Contexts,
                 now: Callable[[], datetime] = lambda: datetime.now(timezone.utc), poll_sec: float = POLL_SEC):
        if name not in NAMES:
            raise ValueError(f"name must be one of {NAMES}, got {name!r}")
        self.name = f"{SOURCE_PREFIX}{name}"
        self.poll_sec = poll_sec
        self._account = name
        self._client = client
        self._contexts = contexts
        self._now = now
        # A market seen SETTLED never changes again: read once.
        self._settled_markets: dict[str, dict] = {}

    def fetch(self) -> list[dict]:
        now = self._now()
        starts_after_ms = int((now - timedelta(days=HISTORY_DAYS)).timestamp() * 1000)
        fills = self._client.fills(starts_after_ms=starts_after_ms)
        by_order: dict[str, list[dict]] = {}
        for fill in fills:
            order_id = fill.get("orderId")
            if not isinstance(order_id, str) or not order_id:
                raise NovigError(f"fill {fill.get('fillId') or '?'} names no order")
            by_order.setdefault(order_id, []).append(fill)
        contexts = self._contexts()
        fetched_at = utc_now_iso()
        records = []
        for order_id, order_fills in by_order.items():
            context = contexts.get(order_id)
            market = self._market_if_started(order_fills[0].get("marketId"), context, now)
            records.append(normalize_order(self._account, order_id, order_fills, context, market, fetched_at))
        log.info("%s: %d fills -> %d bets", self.name, len(fills), len(records))
        return records

    def _market_if_started(self, market_id: object, context: dict | None, now: datetime) -> dict | None:
        """The public market once the game may be over (started, or no context
        to say when), else None: before the start nothing can have settled."""
        if not isinstance(market_id, str):
            return None
        if market_id in self._settled_markets:
            return self._settled_markets[market_id]
        start = context.get("event_start") if context else None
        if start is not None and start > now:
            return None
        market = self._client.market(market_id)
        if market.get("status") == MARKET_SETTLED:
            self._settled_markets[market_id] = market
        return market
