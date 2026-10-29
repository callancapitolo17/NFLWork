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

Novig's catalog names neither the board's teams nor its periods, so a bet's
game fields (league, teams with Unabated's team ids, market, period, side,
number) come from `contexts()`: the board line the order was placed on, by
Novig order id. For "bot" that is taker.duckdb::novig_taker_orders
(TakerStore.novig_order_contexts), for "place" bets.duckdb::manual_orders
(BetsStore.novig_place_order_contexts); an order with no context (placed
outside the desk) is recorded unmatchable, so it shows but sizes nothing.

A context's `settled_at` is when the order's position was first seen closed
(its market settled), kept for good. The bot desk keeps it itself; Place has
no desk sweep, so its read keeps the first sighting through `note_settled` (BetsStore.note_novig_place_settled, an UPDATE of
manual_orders.settled_at) and a restart never moves it. For the same reason
this read also confirms Place's UNCONFIRMED orders (`unconfirmed_orders` /
`confirm_order`): one Novig named but never showed final while Cal waited is
read again here (by clientId when its own read fails) until Novig shows it
final, then given what its fills say. It
also adopts Place's LOST orders (`lost_orders` / `adopt_order`): a
Place send Novig never answered has no order id, so no fill could ever be
named; once it is LOST_ORDER_LOOKUP_AFTER_SEC old the order is looked up by
clientId across Novig's order lists (never re-sent) and given its id, or,
listed nowhere, marked as never having reached Novig.

Grading is off the subaccount's own money, on its own key, every
LEDGER_RECHECK, never off Novig's public catalog: Novig drops most markets
from it around when they settle (2026-10-18: none of the first 15 Place bets
ever read SETTLED there).

    GET /v3/portfolio/positions                                     (signed; open positions only)
    GET /v3/account/subaccounts/<key>/transactions?kind=SETTLEMENT  (signed, paged)

A SETTLEMENT row (ref = market id, amount = what the market paid this
subaccount) grades every order on that market by what it paid against its
stake: more is won, the stake back a push, less lost (a fair-value price pays
part). A loss pays nothing and writes no row, so a started order whose market
has left the open positions with no SETTLEMENT is lost once its settle time
(`settled_at`: the first sighting of the closed position, kept by note_settled
for Place and by the bot's desk) is LOSS_CONFIRM old, dated at that settle time but no later than
LOSS_CLOSED_AT_MAX_AFTER_START after its game's start; a SETTLEMENT lagging the
close turns it won. On 2026-10-18 this read 4 SETTLEMENT rows and 10 closed
positions for the 14 settled Place bets: the 4 won / 10 lost Cal's app shows.
An order with no context grades only off a SETTLEMENT row (no start time to
call a loss by). A failed ledger read keeps the last one and does not fail
the run; an order with a kept settle time it cannot grade (no read since a
restart) is left out of the poll, so its stored row stands.

Ids: "novig:<order id>:<name>", venue "novig" (bets.js joins a spread or total
on venueIds.outcomeId and the side on its own fields). A Novig contract pays
1c: qty 100 = $1; `contracts` here is dollars of payout like sources/novig.py.

balance() reads the subaccount's own balance for the Balances tab:

    GET /v3/account/subaccounts/<key>/balance   (signed) -> {"balance": "<dollars>"}

Side effects: none on disk but `note_settled`'s, `adopt_order`'s and `confirm_order`'s, when given. Raises on any other failed read, so the service
records a failed run and keeps the previous records.
"""
import logging
from collections.abc import Callable
from datetime import datetime, timedelta, timezone

from unabated_ticket.bets_service.normalize import EASTERN, json_clean, round_cents, utc_now_iso
from unabated_ticket.bets_service.novig_orders import FINAL_ORDER_STATUSES, QTY_PER_DOLLAR
from unabated_ticket.bets_service.novig_orders import STATUS_FILLED as NOVIG_STATUS_FILLED
from unabated_ticket.bets_service.novig_orders import NovigError, NovigOrderClient, fill_outcome, order_totals
from unabated_ticket.bets_service.sources.novig import probability_to_american

log = logging.getLogger(__name__)

NAMES = ("bot", "place")
SOURCE_PREFIX = "novig_"
VENUE = "novig"
POLL_SEC = 10.0  # the app source's NOVIG_POLL_SEC cadence (same poll thread)
# Fills on events that started within this many days are read every poll (the store keeps older ones).
HISTORY_DAYS = 14
REASON_NO_CONTEXT = "no board line on record for this Novig order (placed outside the bot's desk)"

Contexts = Callable[[], dict[str, dict]]
# (Novig order id, first SETTLED sighting) -> the settle time kept for good.
NoteSettled = Callable[[str, datetime], datetime]
# placed_before -> [{client_order_id, outcome_id}]: Place orders Novig never answered or named.
LostOrders = Callable[[datetime], list[dict]]
# (client_order_id, Novig order id or None when Novig lists it nowhere) -> None.
AdoptOrder = Callable[[str, str | None], None]
# placed_before -> [{client_order_id, novig_order_id, outcome_id, count_requested}]: Place orders Novig named but never showed final.
UnconfirmedOrders = Callable[[datetime], list[dict]]
# (client_order_id, 'filled' | 'partial' | 'unfilled', contracts filled, dollars paid) -> None.
ConfirmOrder = Callable[[str, str, int, float], None]
# A lost send is looked up only once it is this old: Novig's listings can lag a 201.
LOST_ORDER_LOOKUP_AFTER_SEC = 60
# The order row's client id is 'nman-' + Novig's clientId (novig_place.py).
PLACE_ORDER_ID_PREFIX = "nman-"
# How often the open positions and SETTLEMENT ledger are re-read.
LEDGER_RECHECK = timedelta(minutes=2)
# A closed position with no SETTLEMENT row is a loss only this long after it was first seen closed.
LOSS_CONFIRM = timedelta(minutes=15)
# A no-payout loss is dated no later than this after its game's start. Novig gives a loss no
# settle time, so its closedAt is when the read first saw the position closed; a read that was
# down or not yet grading (2026-10-19: #199's deploy dated 11 Place losses from Oct 17-18 to
# the deploy) would move the loss to another day. Long enough for any game to end, short
# enough that a late Pacific kickoff stays on its own day.
LOSS_CLOSED_AT_MAX_AFTER_START = timedelta(hours=4)


def _iso_ms(value: datetime | None) -> str | None:
    if value is None:
        return None
    value = value.astimezone(timezone.utc)
    return value.strftime("%Y-%m-%dT%H:%M:%S.") + f"{value.microsecond // 1000:03d}Z"


def _iso_of_epoch_ms(epoch_ms: object) -> str | None:
    if isinstance(epoch_ms, bool) or not isinstance(epoch_ms, int):
        return None
    return _iso_ms(datetime.fromtimestamp(epoch_ms / 1000, tz=timezone.utc))


def record_id_of(name: str, novig_order_id: str) -> str:
    """The bet record id of a subaccount's Novig order; the desk's unseen-fill rule looks for it."""
    return f"novig:{novig_order_id}:{name}"


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


def normalize_order(name: str, order_id: str, fills: list[dict], context: dict | None, fetched_at: str,
                    grade: tuple[str, str, float] | None = None) -> dict:
    """One order's fills -> one record, open unless `grade` is given: (status,
    closedAt, dollars paid) off the subaccount's ledger (_ledger_grade). The
    closedAt is fixed (a SETTLEMENT's time or a kept settle time), never the
    poll's: the record would change, and its settle day move, every poll."""
    qty, cost, fee = order_totals(fills)
    stake = round(cost + fee, 5)
    first = fills[0]
    probability = cost / qty * QTY_PER_DOLLAR if qty else None
    payout = qty / QTY_PER_DOLLAR
    status, closed_at, to_win = "open", None, round_cents(payout - stake)
    if grade is not None:
        status, closed_at, paid = grade
        if status == "won":
            # A fair-value settlement pays part of the contract value.
            to_win = round_cents(paid - stake)
    placed_ms = min(fill.get("ts") for fill in fills if isinstance(fill.get("ts"), int)) if any(
        isinstance(fill.get("ts"), int) for fill in fills) else None
    record = {
        "id": record_id_of(name, order_id),
        "source": f"{SOURCE_PREFIX}{name}",
        "venue": VENUE,
        "rotation": None,
        "price": probability_to_american(probability),
        "stake": round_cents(stake),
        "toWin": to_win,
        "contracts": round_cents(payout),
        "placedAt": _iso_of_epoch_ms(placed_ms),
        "status": status,
        "closedAt": closed_at,
        "isParlayLeg": False, "parlayId": None, "legIndex": None, "legCount": None,
        "approx": [],
        "sourceFetchedAt": fetched_at,
        "venueIds": {"marketId": first.get("marketId"), "outcomeId": first.get("outcomeId"),
                     "eventId": None, "gameId": None},
        "raw": {"orderId": order_id, "subaccount": name, "qty": qty, "cost": cost, "fee": fee,
                "probability": probability, "qtyUnit": QTY_PER_DOLLAR,
                "fillIds": [fill.get("fillId") for fill in fills]},
    }
    record.update(_game_fields(context))
    return json_clean(record)


class NovigSubaccountSource:
    """Source protocol for one Novig subaccount (signed GETs on its trading key)."""

    def __init__(self, name: str, client: NovigOrderClient, contexts: Contexts,
                 now: Callable[[], datetime] = lambda: datetime.now(timezone.utc), poll_sec: float = POLL_SEC,
                 note_settled: NoteSettled | None = None, lost_orders: LostOrders | None = None,
                 adopt_order: AdoptOrder | None = None, unconfirmed_orders: UnconfirmedOrders | None = None,
                 confirm_order: ConfirmOrder | None = None):
        if name not in NAMES:
            raise ValueError(f"name must be one of {NAMES}, got {name!r}")
        self.name = f"{SOURCE_PREFIX}{name}"
        self.poll_sec = poll_sec
        self._account = name
        self._client = client
        self._contexts = contexts
        self._note_settled = note_settled
        self._lost_orders = lost_orders
        self._adopt_order = adopt_order
        self._unconfirmed_orders = unconfirmed_orders
        self._confirm_order = confirm_order
        self._now = now
        # (open position market ids, market id -> (dollars paid, last SETTLEMENT ms)) as of _ledger_read_at.
        self._ledger: tuple[set[str], dict[str, tuple[float, int]]] | None = None
        self._ledger_read_at: datetime | None = None
        self._ledger_tried_at: datetime | None = None

    def fetch(self) -> list[dict]:
        now = self._now()
        self._adopt_lost_orders(now)
        starts_after_ms = int((now - timedelta(days=HISTORY_DAYS)).timestamp() * 1000)
        fills = self._client.fills(starts_after_ms=starts_after_ms)
        by_order: dict[str, list[dict]] = {}
        for fill in fills:
            order_id = fill.get("orderId")
            if not isinstance(order_id, str) or not order_id:
                raise NovigError(f"fill {fill.get('fillId') or '?'} names no order")
            by_order.setdefault(order_id, []).append(fill)
        self._confirm_unconfirmed_orders(now)
        contexts = self._contexts()
        self._read_ledger_if_due(now, starts_after_ms)
        held_qty: dict[str, int] = {}
        held_outcomes: dict[str, set] = {}
        for order_fills in by_order.values():
            market_id = order_fills[0].get("marketId")
            held_qty[market_id] = held_qty.get(market_id, 0) + order_totals(order_fills)[0]
            held_outcomes.setdefault(market_id, set()).update(fill.get("outcomeId") for fill in order_fills)
        fetched_at = utc_now_iso()
        records = []
        for order_id, order_fills in by_order.items():
            context = contexts.get(order_id)
            grade = self._ledger_grade(order_id, order_fills, context, now, held_qty, held_outcomes)
            if grade is None and self._may_be_graded(context, now):
                # It may already be stored graded, and no ledger read can grade it now: the stored row stands.
                continue
            records.append(normalize_order(self._account, order_id, order_fills, context, fetched_at, grade))
        log.info("%s: %d fills -> %d bets", self.name, len(fills), len(records))
        return records

    def balance(self) -> dict:
        """This subaccount's balance for the tracker's Balances tab (balances.py).
        A subaccount is a wallet of its own, outside the app account's summary
        (2026-10-23: the app showed the main wallet's $11,893 alone, Place's
        $2,446 apart). Novig gives one figure ("balance") and no money-in-bets
        split, so inBets is None. Raises NovigError on a failed read."""
        cash = self._client.balance()
        return {"cash": cash, "inBets": None, "available": cash, "creditLimit": None}

    def _may_be_graded(self, context: dict | None, now: datetime) -> bool:
        """True for an order an earlier run may have graded: its position was seen
        closed, or (before this run's first ledger read) it has started or has no context."""
        if context is not None and context.get("settled_at") is not None:
            return True
        if self._ledger is not None:
            return False
        start = context.get("event_start") if context else None
        return start is None or start <= now

    def _read_ledger_if_due(self, now: datetime, starts_after_ms: int) -> None:
        """Re-read the open positions, then the SETTLEMENT ledger (in that
        order: a market that settles between the two reads shows its payout),
        once per LEDGER_RECHECK, a failed one too (a 429 is not retried every
        poll). A failed read keeps the last one."""
        if self._ledger_tried_at is not None and now - self._ledger_tried_at < LEDGER_RECHECK:
            return
        self._ledger_tried_at = now
        try:
            open_markets = {position.get("marketId") for position in self._client.open_positions()}
            settlements: dict[str, tuple[float, int]] = {}
            for row in self._client.settlements(starts_after_ms):
                market_id, amount, ts_ms = row.get("ref"), float(row.get("amount")), row.get("ts")
                if not isinstance(market_id, str) or not isinstance(ts_ms, int):
                    raise ValueError(f"SETTLEMENT row {row.get('transactionId')} has ref {market_id!r}, ts {ts_ms!r}")
                paid, last_ms = settlements.get(market_id, (0.0, 0))
                settlements[market_id] = (paid + amount, max(last_ms, ts_ms))
        except (NovigError, OSError, TypeError, ValueError) as error:
            log.error("%s: positions / SETTLEMENT ledger read failed, grading off the last read: %s", self.name, error)
            return
        self._ledger = (open_markets, settlements)
        self._ledger_read_at = now

    def _ledger_grade(self, order_id: str, fills: list[dict], context: dict | None, now: datetime,
                      held_qty: dict[str, int], held_outcomes: dict[str, set]) -> tuple[str, str, float] | None:
        """(status, closedAt, dollars paid) for an order its market's SETTLEMENT
        row grades (shared out by contracts across this subaccount's orders
        on the market), or "lost" for a started order whose position closed
        with none LOSS_CONFIRM ago; None while the ledger cannot say."""
        market_id = fills[0].get("marketId")
        if self._ledger is None or len(held_outcomes.get(market_id, ())) != 1:
            return None
        open_markets, settlements = self._ledger
        qty, cost, fee = order_totals(fills)
        stake = cost + fee
        if market_id in settlements and held_qty[market_id] > 0:
            if context is not None and context.get("settled_at") is None and self._note_settled is not None:
                self._note_settled(order_id, now)  # so a restart's failed first read leaves it as stored
            market_paid, settled_ms = settlements[market_id]
            paid = market_paid * qty / held_qty[market_id]
            # A payout under the stake (a fair-value settle) is "won" with a negative toWin, since
            # the tracker scores every "lost" bet at the full stake.
            if abs(paid - stake) < 0.005:
                status = "push"
            else:
                status = "won" if paid >= 0.005 else "lost"
            return status, _iso_of_epoch_ms(settled_ms), paid
        start = context.get("event_start") if context else None
        # A positions read from before the start can predate the order's fills.
        if start is None or start > now or self._ledger_read_at <= start or market_id in open_markets:
            return None
        settled_at = context.get("settled_at")
        if settled_at is None and self._note_settled is not None:
            settled_at = self._note_settled(order_id, now)
        if settled_at is None or now - settled_at < LOSS_CONFIRM:
            return None
        return "lost", _iso_ms(min(settled_at, start + LOSS_CLOSED_AT_MAX_AFTER_START)), 0.0

    def _confirm_unconfirmed_orders(self, now: datetime) -> None:
        """Give each Place order Novig never showed final while Cal waited its
        fills, once Novig shows it final: its own read, or by clientId on its
        outcome's order lists when that read fails (Novig can 404 an order).
        The fills are read after the final status, so none can land between.
        A FILLED order whose fills have not all landed waits for the next
        poll; a failed read too."""
        if self._unconfirmed_orders is None or self._confirm_order is None:
            return
        for order in self._unconfirmed_orders(now - timedelta(seconds=LOST_ORDER_LOOKUP_AFTER_SEC)):
            client_order_id = order["client_order_id"]
            order_id = order["novig_order_id"]
            try:
                order_status = self._final_status_of(order)
                if order_status not in FINAL_ORDER_STATUSES:
                    continue
                filled, cost, fee = order_totals(self._client.fills(order_id=order_id))
            except (NovigError, OSError) as error:
                log.error("%s: unconfirmed order %s: read failed, still unconfirmed: %s", self.name, client_order_id, error)
                continue
            requested = order["count_requested"]
            if order_status == NOVIG_STATUS_FILLED and filled < requested:
                continue
            outcome = fill_outcome(filled, requested)
            self._confirm_order(client_order_id, outcome, filled, round(cost + fee, 5))
            log.warning("%s: unconfirmed order %s is %s on Novig: %s, %d of %d", self.name, client_order_id,
                        order_status, outcome, filled, requested)

    def _final_status_of(self, order: dict) -> str | None:
        """Novig's status for an unconfirmed Place order: its own read, else
        its listing by clientId; None when neither names it."""
        status, reply = self._client.get_order(order["novig_order_id"])
        if status == 200 and isinstance(reply, dict):
            return reply.get("status")
        novig_client_id = order["client_order_id"].removeprefix(PLACE_ORDER_ID_PREFIX)
        listed = self._client.find_order_by_client_id(order["outcome_id"], novig_client_id)
        if listed and listed.get("orderId") == order["novig_order_id"]:
            return listed.get("status")
        return None

    def _adopt_lost_orders(self, now: datetime) -> None:
        """Give each Place order Novig never answered its Novig id off Novig's
        order lists by clientId, or mark it as never sent when it is listed
        nowhere. A failed listing read leaves the order lost for the next poll."""
        if self._lost_orders is None or self._adopt_order is None:
            return
        for lost in self._lost_orders(now - timedelta(seconds=LOST_ORDER_LOOKUP_AFTER_SEC)):
            client_order_id = lost["client_order_id"]
            novig_client_id = client_order_id.removeprefix(PLACE_ORDER_ID_PREFIX)
            try:
                match = self._client.find_order_by_client_id(lost["outcome_id"], novig_client_id)
            except (NovigError, OSError) as error:
                log.error("%s: lost order %s: lookup failed, still lost: %s", self.name, client_order_id, error)
                continue
            order_id = match.get("orderId") if match else None
            self._adopt_order(client_order_id, order_id)
            log.warning("%s: lost order %s %s", self.name, client_order_id,
                        f"found on Novig as {order_id} ({match.get('status')})" if match else "is in none of Novig's order lists: never sent")
