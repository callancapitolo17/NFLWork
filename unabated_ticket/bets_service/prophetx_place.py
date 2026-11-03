"""The Ticket's Place button at ProphetX (Cal 2026-10-31, option A): Cal clicks,
and this desk sends ONE order on his main ProphetX account, on the key pair the
ProphetX bet source reads with (PROPHETX_ACCESS_KEY / _SECRET_KEY). Nothing here
runs on a timer: every order is one POST /place_prophetx.json from the panel or
the phone page.

Rules, in order, for every request (validate_prophetx_place_request's shape):
  1. the click was not already sent: 'xman-' + its clickId, held in
     bets.duckdb::manual_orders before the send, and sent as the order's
     external_id;
  2. the server runner finds the board's one ProphetX line and the ProphetX
     strike ids Unabated's feed carries for that game (GET
     /prophetx_quote.json); this desk finds the ProphetX game: by those strike
     ids (GET /v4/mm/get_strikes names their game), else by the two teams and
     the start among ProphetX's listed games, refusing on zero or two;
  3. the game is listed and not started (the earlier of Unabated's and
     ProphetX's start), and its away and home teams are the board's: ProphetX's
     two names resolve to Unabated team ids through the runner (GET
     /team_keys.json, teams.js + the crosswalk, as bet matching does) and equal
     the board line's ids (a swapped or unresolved team refuses: the side
     would be wrong);
  4. the one strike among the game's markets that reads as the board line
     (prophetx_markets.lines_of_event; the line's own strike id when Unabated
     names one), and its best offer now;
  5. that offer's implied probability is at most PRICE_TOLERANCE (0.1c) worse
     than the Ticket's; a better one is taken. ProphetX's fee is set per
     account and unpublished, so it is not counted (Cal, 2026-10-31: use
     Unabated's price as is); each fill's own fee is stored;
  6. the other side of the line is not held (as at every venue: Place adds to a
     side, never hedges it);
  7. the prophetx read the stake counted is ok and fresh;
  8. no earlier Place on this game, at any venue, may have filled after the
     read of its venue the stake was sized on (manual_orders.unseen_manual_fill);
  9. the row is INSERTED (status 'sending') BEFORE the send;
 10. ONE fill-or-kill buy at that best price for the stake (prophetx_client),
     never re-sent and never cancelled: ProphetX kills a fill-or-kill's
     remainder itself about 2 s after its last fill. The order is read off the
     account's order history (the last HISTORY_READ_LOOKBACK_SEC of orders;
     the single-order read lags 8 s+ after a send) for up to
     PLACE_FINAL_READ_WINDOW_SEC; one not yet read finished is left 'unknown'
     and Place says ProphetX is still matching it.
 11. The account's order history is the record (Cal 2026-11-03, option A): every
     ProphetX source poll hands its orders and fills to
     settle_unconfirmed_place_orders, which settles each 'unknown' ProphetX
     Place row (by its order id, else its external_id = client order id) once
     ProphetX shows the order finished. On 2026-11-03 two orders read
     "inactive, money open" for the whole window while ProphetX had already
     matched them in full, and the cancel that followed answered
     order_already_filled: one look is not the record.

Order statuses (manual_orders.status on a 'prophetx' row): as Kalshi Place's.
Side effects: a ProphetX order on Cal's main account (money); INSERT / UPDATE
bets.duckdb::manual_orders (venue 'prophetx'); settle_unconfirmed_place_orders
UPDATEs 'unknown' prophetx rows from the order history (no network). Network: GET
<runner>/prophetx_quote.json and /team_keys.json, then ProphetX's Trading API on the key pair.
"""
import logging
import math
import threading
import time
from collections.abc import Callable
from datetime import datetime, timedelta, timezone

from unabated_ticket.bets_service.manual_orders import (
    CLICK_ID_PATTERN, PRICE_TOLERANCE, STATUS_FILLED, STATUS_PARTIAL, STATUS_REFUSED, STATUS_REJECTED,
    STATUS_UNFILLED, STATUS_UNKNOWN, VENUE_PROPHETX, _cents, _is_number, _iso_time, filled_bet_message, optional_read_at,
    unseen_manual_fill,
)
from unabated_ticket.bets_service.prophetx_client import (
    FINAL_ORDER_STATUSES, MAX_STRIKES_PER_READ, ProphetXClient, ProphetXError,
)
from unabated_ticket.bets_service.prophetx_markets import (
    american_to_probability, best_offer_of, league_of, lines_of_event, number_of, profit_at, team_name_of, teams_of,
)
from unabated_ticket.bets_service.sources.prophetx import VENUE, iso_utc_of, matched_of
from unabated_ticket.bets_service.store import BetsStore
from unabated_ticket.bets_service.taker_orders import DEFINITE_REJECTION_STATUSES

log = logging.getLogger(__name__)

PROPHETX_ORDER_ID_PREFIX = "xman-"
RUNNER_PROPHETX_QUOTE_PATH = "/prophetx_quote.json"
RUNNER_PROPHETX_QUOTE_TIMEOUT_SEC = 10
PLACE_FINAL_READ_ATTEMPTS = 6
PLACE_FINAL_READ_SLEEP_SEC = 1.0
PLACE_FINAL_READ_WINDOW_SEC = 8.0
# The order-history window Place reads its order from: orders created this
# long before now, to this long after (room for ProphetX's clock running ahead).
HISTORY_READ_LOOKBACK_SEC = 120
HISTORY_READ_LOOKAHEAD_SEC = 60
HISTORY_READ_MAX_PAGES = 2
PROPHETX_READ_SOURCE = VENUE
PROPHETX_READ_MAX_AGE = timedelta(minutes=5)
BET_TYPE_MONEYLINE = 1
# ProphetX's game must start within this of the board's to be the same game
# when it is found by its teams (bets.js's own game rule).
START_TOLERANCE = timedelta(minutes=30)
# A stake ProphetX takes: whole cents, at least a dollar (its minimum is unpublished).
MIN_STAKE_DOLLARS = 1.0
PROPHETX_EVENT_NOT_STARTED = "not_started"


def validate_prophetx_place_request(body: object) -> dict | str:
    """A POST /place_prophetx.json body, or what was expected and what was found.
    {clickId, eventId, periodTypeId, betTypeId (1|2|3), sideIndex, points (null on
    a moneyline), stake (dollars), ticketAllIn (the Ticket's ProphetX price as a
    probability), prophetxReadAt (the prophetx read the stake was sized on, ISO),
    kalshiReadAt / novigReadAt / polymarketReadAt / bfaReadAt (ISO, or null when never read)}"""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    click_id = body.get("clickId")
    if not isinstance(click_id, str) or not CLICK_ID_PATTERN.match(click_id):
        return f"clickId must be 8-48 letters, digits or dashes, got {click_id!r}"
    whole = {}
    for field, allowed in (("eventId", None), ("periodTypeId", None), ("betTypeId", (1, 2, 3)), ("sideIndex", (0, 1))):
        value = body.get(field)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0 or (allowed and value not in allowed):
            return f"{field} must be {'one of ' + str(allowed) if allowed else 'a whole number'}, got {value!r}"
        whole[field] = value
    points = body.get("points")
    if whole["betTypeId"] == BET_TYPE_MONEYLINE:
        points = None
    elif not _is_number(points) or not math.isfinite(points):
        return f"points must be a number (a spread or total), got {points!r}"
    stake = body.get("stake")
    if not _is_number(stake) or not math.isfinite(stake) or stake < MIN_STAKE_DOLLARS:
        return f"stake must be at least ${MIN_STAKE_DOLLARS:.0f}, got {stake!r}"
    ticket_all_in = body.get("ticketAllIn")
    if not _is_number(ticket_all_in) or not 0 < ticket_all_in < 1:
        return f"ticketAllIn must be a probability between 0 and 1, got {ticket_all_in!r}"
    prophetx_read_at = _iso_time(body.get("prophetxReadAt"))
    if prophetx_read_at is None:
        return ("prophetxReadAt must be the ISO time of the ProphetX read the stake was sized on, "
                f"got {body.get('prophetxReadAt')!r}")
    read_ats = {}
    for field in ("kalshiReadAt", "novigReadAt", "polymarketReadAt", "bfaReadAt"):
        read_at = optional_read_at(body, field)
        if isinstance(read_at, str):
            return read_at
        read_ats[field] = read_at
    return {"click_id": click_id, **whole, "points": None if points is None else float(points),
            "stake": round(float(stake), 2), "ticket_all_in": float(ticket_all_in), "prophetx_read_at": prophetx_read_at,
            "kalshi_read_at": read_ats["kalshiReadAt"], "novig_read_at": read_ats["novigReadAt"],
            "polymarket_read_at": read_ats["polymarketReadAt"], "bfa_read_at": read_ats["bfaReadAt"]}


# A fill within a cent of the stake is all of it (ProphetX rounds dollars to the cent).
FULL_FILL_TOLERANCE_DOLLARS = 0.01
REVERSED_TEAMS_WORDS = "home and away the other way round from the board"


class TeamKeysUnread(Exception):
    """The runner could not resolve ProphetX's team names: the Place refuses, naming why."""


# (league, away name, home name) -> the runner's /team_keys.json answer
# {ok, awayKey, homeKey} (manual_orders.read_runner_team_keys).
TeamKeysReader = Callable[[str, str, str], dict]


def board_team_keys(line: dict) -> tuple[str, str] | str:
    """The board line's away and home as Unabated team keys ("<league>:<id>",
    teams.js keyOf), or why it names none."""
    league, away_id, home_id = line.get("league"), line.get("awayTeamId"), line.get("homeTeamId")
    if not league or away_id in (None, "") or home_id in (None, ""):
        return f"the board line names no team ids ({league} {away_id} @ {home_id})"
    return f"{league}:{away_id}", f"{league}:{home_id}"


def teams_mismatch(event: dict, line: dict, team_keys_of: TeamKeysReader) -> str | None:
    """Why ProphetX's away / home are not the board's, or None. ProphetX's
    names resolve to Unabated team ids through the runner (teams.js + the
    crosswalk, as every bet record does), compared with the board line's ids:
    "Sam Houston State" and the board's "Sam Houston St" are one team, and a
    name nothing resolves refuses. Raises TeamKeysUnread when the runner could
    not answer."""
    teams = teams_of(event)
    if isinstance(teams, str):
        return teams
    board = board_team_keys(line)
    if isinstance(board, str):
        return board
    px_away, px_home = team_name_of(teams["away"]), team_name_of(teams["home"])
    answer = team_keys_of(line["league"], px_away, px_home)
    if not answer.get("ok"):
        raise TeamKeysUnread(answer.get("detail") or answer.get("error") or "the runner named no team keys")
    keys = (answer.get("awayKey"), answer.get("homeKey"))
    unresolved = [name for name, key in zip((px_away, px_home), keys) if not key]
    if unresolved:
        return f"ProphetX's {' and '.join(unresolved)} resolves to no {line['league']} team on Unabated's board"
    if keys == board:
        return None
    if keys == board[::-1]:
        return f"ProphetX lists {px_away} @ {px_home}, {REVERSED_TEAMS_WORDS}"
    return f"ProphetX's {px_away} @ {px_home} does not read as the board's {line.get('awayTeam')} @ {line.get('homeTeam')}"


def board_line_of(line: dict) -> dict:
    """The board line as prophetx_markets names one: {betType, period, side, points}."""
    return {"betType": line.get("betType"), "period": line.get("period"), "side": line.get("side"),
            "points": None if line.get("betType") == "moneyline" else line.get("points")}


def matching_lines(lines: list[dict], wanted: dict) -> list[dict]:
    return [line for line in lines if all(line[key] == wanted[key] for key in ("betType", "period", "side", "points"))]


class ManualProphetXDesk:
    """Sends Cal's ProphetX Place-button orders on his main account, one per click."""

    def __init__(self, client: ProphetXClient, store: BetsStore, quote_reader: Callable[[dict], dict],
                 team_keys_reader: TeamKeysReader, now: Callable[[], datetime] = lambda: datetime.now(timezone.utc),
                 sleep: Callable[[float], None] = time.sleep, clock: Callable[[], float] = time.monotonic):
        self._client = client
        self._store = store
        self._quote_reader = quote_reader
        self._team_keys_reader = team_keys_reader
        self._now = now
        self._sleep = sleep
        self._clock = clock
        self._lock = threading.Lock()

    def place(self, request: dict) -> dict:
        """One order. Returns {status, message, ...}: refused (nothing sent), or
        the sent order's filled | partial | unfilled | rejected | unknown."""
        with self._lock:
            return self._place(request)

    def _place(self, request: dict) -> dict:
        client_order_id = f"{PROPHETX_ORDER_ID_PREFIX}{request['click_id']}"
        quote = self._quote_reader({"eventId": request["eventId"], "periodTypeId": request["periodTypeId"],
                                    "betTypeId": request["betTypeId"], "sideIndex": request["sideIndex"],
                                    "points": request["points"]})
        target = self._target(quote)
        if isinstance(target, str):
            return {"status": STATUS_REFUSED, "message": f"Not placed: {target}"}
        refusal = self._refusal(request, target)
        if refusal:
            return {"status": STATUS_REFUSED, "message": f"Not placed: {refusal}"}
        line = quote["line"]
        target["label"] = quote.get("label")
        probability = target["probability"]
        inserted = self._store.insert_manual_order({
            "client_order_id": client_order_id, "placed_at": self._now(), "venue": VENUE_PROPHETX,
            "event_id": request["eventId"], "ticker": target["strike_id"], "contract_side": "buy",
            "label": quote.get("label"), "ticket_all_in": request["ticket_all_in"], "ask_price": probability,
            "all_in_price": probability, "stake": request["stake"], "count_requested": int(round(request["stake"] * 100)),
            "status": "sending", "updated_at": self._now(), "event_start": target["eventStart"],
            "market_id": str(target["marketId"]), "league": line.get("league"), "bet_type": line.get("betType"),
            "period": line.get("period"), "side": line.get("side"), "points": line.get("points"),
            "away_team": line.get("awayTeam"), "home_team": line.get("homeTeam"),
            "away_team_id": line.get("awayTeamId"), "home_team_id": line.get("homeTeamId"),
        })
        if not inserted:
            return {"status": STATUS_REFUSED, "message": "Not placed: this click was already sent."}
        log.info("prophetx place %s: strike %s $%.2f @ %+d (%s, Ticket %s)", client_order_id, target["strike_id"],
                 request["stake"], target["price"], _cents(probability), _cents(request["ticket_all_in"]))
        return self._send(client_order_id, target, request["stake"])

    # ---- finding the strike ----------------------------------------------------------

    def _target(self, quote: dict) -> dict | str:
        """{strike_id, price, probability, depth, marketId, eventStart, otherStrikeId}
        for the quote's board line, or why there is none."""
        if not quote.get("ok"):
            return f"{quote.get('detail') or quote.get('reason') or 'no ProphetX line for this'}."
        line = quote.get("line") if isinstance(quote.get("line"), dict) else None
        if line is None or line.get("betType") not in ("moneyline", "spread", "total"):
            return "the runner did not name the board line."
        if line.get("betType") == "moneyline" and line.get("period") != "FG":
            return f"a {line.get('period')} moneyline: whether ProphetX voids a tied period is unchecked."
        try:
            event = self._event_of(quote, line)
            if isinstance(event, str):
                return event
            mismatch = teams_mismatch(event, line, self._team_keys_reader)
            # A total is the same bet with home and away swapped; any other mismatch refuses.
            if mismatch and not (line.get("betType") == "total" and mismatch.endswith(REVERSED_TEAMS_WORDS)):
                return f"{mismatch}."
            lines = lines_of_event(event, self._client.markets(event["event_id"]))
        except TeamKeysUnread as error:
            return f"ProphetX's teams could not be checked against the board's ({error})."
        except (OSError, ProphetXError) as error:
            return f"ProphetX could not be read ({error})."
        if isinstance(lines, str):
            return f"{lines}."
        found = matching_lines(lines, board_line_of(line))
        own_strike = quote.get("lineStrikeId")
        label = f"{line.get('side')} {line.get('points') if line.get('points') is not None else ''} {line.get('period')}".strip()
        if len(found) != 1:
            return f"{len(found)} ProphetX strikes read as {label} on game {event['event_id']} ({len(lines)} sides listed)."
        target = found[0]
        if own_strike and own_strike != target["strike_id"]:
            return f"Unabated's ProphetX strike {own_strike} is not the one that reads as {label} ({target['strike_id']})."
        offer = best_offer_of(target["ladder"])
        if offer is None:
            return f"nothing rests on ProphetX at {label}."
        px_start = _prophetx_time(event.get("scheduled"))
        board_start = _iso_time(quote.get("eventStart"))
        starts = [start for start in (px_start, board_start) if start is not None]
        return {"strike_id": target["strike_id"], "price": offer["price"], "depth": offer["depth"],
                "probability": round(american_to_probability(offer["price"]), 6), "marketId": target["marketId"],
                "eventStart": min(starts) if starts else None, "eventStatus": event.get("status"),
                "otherStrikeId": target["otherStrikeId"], "prophetxEventId": event["event_id"]}

    def _event_of(self, quote: dict, line: dict) -> dict | str:
        """ProphetX's one listed game the board line is on, or why there is none."""
        strike_ids = [strike for strike in (quote.get("strikeIds") or []) if isinstance(strike, str) and strike]
        if strike_ids:
            event_ids = {strike.get("event_id") for strike in self._client.strikes(strike_ids[:MAX_STRIKES_PER_READ])}
            event_ids.discard(None)
            if len(event_ids) != 1:
                return f"Unabated's ProphetX strikes for this game name {len(event_ids)} ProphetX games."
            listed = self._client.sport_events([event_ids.pop()])
            if len(listed) != 1:
                return "ProphetX no longer lists this game (started or ended)."
            return listed[0]
        board_start = _iso_time(quote.get("eventStart"))
        if board_start is None:
            return "the board has no start for this game."
        found = []
        for event in self._client.sport_events():
            start = _prophetx_time(event.get("scheduled"))
            if league_of(event) != line.get("league") or start is None or abs(start - board_start) > START_TOLERANCE:
                continue
            if teams_mismatch(event, line, self._team_keys_reader) is None:
                found.append(event)
        if len(found) != 1:
            return (f"{len(found)} of ProphetX's listed games are {line.get('awayTeam')} @ {line.get('homeTeam')} "
                    "at the board's start.")
        return found[0]

    # ---- refusals ------------------------------------------------------------------

    def _refusal(self, request: dict, target: dict) -> str | None:
        """Why this click must not become an order, in words, or None."""
        start = target["eventStart"]
        if start is None:
            return "the game's start time is unknown."
        if start <= self._now() or target.get("eventStatus") not in (None, PROPHETX_EVENT_NOT_STARTED):
            return "the game has started."
        if target["probability"] > request["ticket_all_in"] + PRICE_TOLERANCE:
            return (f"ProphetX is at {target['price']:+d} now ({_cents(target['probability'])}), more than 0.1¢ worse "
                    f"than the Ticket's {_cents(request['ticket_all_in'])}.")
        if self._holds_strike(target["otherStrikeId"]):
            return "you hold the other side of this line on ProphetX."
        stale = self._stale_read()
        if stale:
            return stale
        if unseen_manual_fill(self._store, request["eventId"], request["kalshi_read_at"], request["novig_read_at"],
                              request["polymarket_read_at"], request["prophetx_read_at"],
                              request["bfa_read_at"]):
            return ("your last Place on this game is not on the Bets tab yet, so this stake does not count it. "
                    "Wait for the next read (up to a minute) and place again.")
        return None

    def _holds_strike(self, strike_id: str) -> bool:
        return any(record.get("status") == "open" and not record.get("unmatchable")
                   and (record.get("venueIds") or {}).get("strikeId") == strike_id
                   for record in self._store.load_records_of_venue(VENUE).values())

    def _stale_read(self) -> str | None:
        read = self._store.source_status().get(PROPHETX_READ_SOURCE)
        fetched_at = _iso_time(read.get("fetchedAt")) if read else None
        if not read or not read.get("ok") or fetched_at is None:
            return "the ProphetX read is failing, so the stake may miss a bet. Try again once it reads."
        if self._now() - fetched_at > PROPHETX_READ_MAX_AGE:
            return (f"the ProphetX read is over {PROPHETX_READ_MAX_AGE.seconds // 60} minutes old, "
                    "so the stake may miss a bet.")
        return None

    # ---- the order ---------------------------------------------------------------------

    def _send(self, client_order_id: str, target: dict, stake: float) -> dict:
        """ONE fill-or-kill buy, never re-sent; then the order read until nothing rests."""
        try:
            status, reply = self._client.submit_fill_or_kill(target["strike_id"], target["price"], stake, client_order_id)
        except (OSError, ProphetXError) as error:
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": f"no answer: {error}"})
            log.error("prophetx place %s: no answer from ProphetX: %s", client_order_id, error)
            return {"status": STATUS_UNKNOWN, "message": "ProphetX did not answer. Check ProphetX before placing again."}
        if status in DEFINITE_REJECTION_STATUSES:
            self._update(client_order_id, {"status": STATUS_REJECTED, "error": f"HTTP {status}: {str(reply)[:500]}"})
            log.warning("prophetx place %s refused: HTTP %s %s", client_order_id, status, reply)
            detail = _error_text(reply)
            return {"status": STATUS_REJECTED,
                    "message": f"ProphetX refused it (HTTP {status}{': ' + detail if detail else ''}). Nothing was placed."}
        data = reply.get("data") if isinstance(reply, dict) else None
        order = data.get("order") if isinstance(data, dict) else None
        order_id = order.get("order_id") if isinstance(order, dict) else None
        if status not in (200, 201) or not isinstance(order_id, str) or not order_id:
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": f"HTTP {status}: {str(reply)[:500]}"})
            log.error("prophetx place %s: no definite answer: HTTP %s %s", client_order_id, status, reply)
            return {"status": STATUS_UNKNOWN, "message": f"ProphetX answered HTTP {status}. Check ProphetX before placing again."}
        self._update(client_order_id, {"prophetx_order_id": order_id})
        return self._follow_order(client_order_id, order_id, order, target, stake)

    def _follow_order(self, client_order_id: str, order_id: str, order: dict, target: dict, stake: float) -> dict:
        """Read the order off the account's order history until nothing of it
        can still fill (final, or fully matched), then what filled. ProphetX's
        single-order read (GET get_order/<id>) answers 404 data_not_found for
        8 s and more after a send (2026-11-03, while the order had filled in
        about 1 s) and the history held it within 3 s, so the history is read.
        One not read finished within the window is left 'unknown' for the
        source's next order-history poll to settle; it is never cancelled."""
        deadline = self._clock() + PLACE_FINAL_READ_WINDOW_SEC
        for attempt in range(PLACE_FINAL_READ_ATTEMPTS):
            if _done(order, stake):
                return self._record(client_order_id, order, target, stake)
            if self._clock() >= deadline:
                break
            self._sleep(PLACE_FINAL_READ_SLEEP_SEC)
            read = self._history_read(client_order_id, order_id, attempt + 1)
            if read is not None:
                order = read
        if _done(order, stake):
            return self._record(client_order_id, order, target, stake)
        return self._still_matching(client_order_id, order_id, order)

    def _history_read(self, client_order_id: str, order_id: str, attempt: int) -> dict | None:
        """The order as the order history shows it now, or None when the read
        fails or does not hold it yet."""
        now = int(self._now().timestamp())
        try:
            orders = self._client.order_history(now - HISTORY_READ_LOOKBACK_SEC, now + HISTORY_READ_LOOKAHEAD_SEC,
                                                max_pages=HISTORY_READ_MAX_PAGES)
        except (OSError, ProphetXError) as error:
            log.warning("prophetx place %s: order history read %d failed: %s", client_order_id, attempt, error)
            return None
        order = next((row for row in orders if str(row.get("order_id")) == order_id), None)
        if order is None:
            log.warning("prophetx place %s: order history read %d: order %s not listed yet (%d orders)",
                        client_order_id, attempt, order_id, len(orders))
            return None
        log.info("prophetx place %s: order history read %d: %s", client_order_id, attempt, _order_state(order))
        return order

    def _still_matching(self, client_order_id: str, order_id: str, order: dict) -> dict:
        log.warning("prophetx place %s: order %s not read finished in %.0f s (%s); the order history settles it",
                    client_order_id, order_id, PLACE_FINAL_READ_WINDOW_SEC, _order_state(order))
        self._update(client_order_id, {"status": STATUS_UNKNOWN,
                                       "error": f"not read finished in the window: {_order_state(order)}"})
        return {"status": STATUS_UNKNOWN,
                "message": "Sent. ProphetX was still matching it when Place stopped reading; what filled shows "
                           "in Bets in a few seconds. Check Bets before placing again."}

    def _record(self, client_order_id: str, order: dict, target: dict, stake: float) -> dict:
        filled = number_of(order.get("filled_quantity"))
        if filled is None:
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": f"final but unreadable: {str(order)[:300]}"})
            return {"status": STATUS_UNKNOWN, "message": "ProphetX finished the order but did not say what filled. "
                                                         "Check ProphetX before placing again."}
        fee = number_of(order.get("fee")) or 0.0
        outcome = _outcome_of(order, filled, stake)
        self._update(client_order_id, {"status": outcome, "count_filled": round(filled * 100),
                                       "cost_filled": round(filled + fee, 4),
                                       "error": None if outcome != STATUS_REJECTED else f"invalid: {str(order)[:300]}"})
        log.info("prophetx place %s: %s, $%.2f of $%.2f (fee %s)", client_order_id, outcome, filled, stake, fee)
        return {"status": outcome, "strikeId": target["strike_id"], "price": target["price"], "stake": stake,
                "filled": filled, "fee": fee, "message": _outcome_message(outcome, target, stake, filled, fee)}

    def _update(self, client_order_id: str, fields: dict) -> None:
        self._store.update_manual_order(client_order_id, fields, self._now())


def _prophetx_time(value: object) -> datetime | None:
    """A ProphetX time (ISO text, or unix seconds / ms) as an aware datetime."""
    return _iso_time(iso_utc_of(value))


def _done(order: dict, stake: float) -> bool:
    """Nothing of the order can still fill: ProphetX finished it, or all of
    the stake matched. An open order reading 0 open and less than the stake
    filled may be one matching has not reached yet (a fill-or-kill rests ~2 s),
    so it is not done."""
    status = str(order.get("status") or "").lower()
    if status in FINAL_ORDER_STATUSES or str(order.get("matching_status") or "").lower() == "fully_matched":
        return True
    filled = number_of(order.get("filled_quantity"))
    return filled is not None and filled >= stake - FULL_FILL_TOLERANCE_DOLLARS


def _outcome_of(order: dict, filled: float, stake: float) -> str:
    """A finished order's manual_orders status from the dollars it matched."""
    if filled <= 0 and str(order.get("status")).lower() == "invalid":
        return STATUS_REJECTED
    if filled <= 0:
        return STATUS_UNFILLED
    return STATUS_FILLED if filled >= stake - FULL_FILL_TOLERANCE_DOLLARS else STATUS_PARTIAL


def _order_state(order: dict) -> str:
    """An order read in a log line: what decides whether it is done."""
    return (f"status {order.get('status')}, matching {order.get('matching_status')}, "
            f"filled {order.get('filled_quantity')}, open {order.get('open_quantity')}")


def settle_unconfirmed_place_orders(store: BetsStore, orders: list[dict], trades: list[dict], now: datetime) -> None:
    """Settle every 'unknown' ProphetX Place row from the order history the
    ProphetX source just read (no network): found by its ProphetX order id,
    else by its external_id (a click whose submit answer was lost), and
    settled once ProphetX shows the order finished, at the dollars its fills
    matched plus their fees. An order the history does not hold, or still
    matching, stays 'unknown' and is looked for again next poll, for
    store.UNCONFIRMED_PROPHETX_PLACE_WINDOW. Side effect: UPDATE bets.duckdb::manual_orders."""
    unconfirmed = store.unconfirmed_prophetx_place_orders(now)
    if not unconfirmed:
        return
    by_order_id = {str(order["order_id"]): order for order in orders if order.get("order_id")}
    by_external_id = {str(order["external_id"]): order for order in orders if order.get("external_id")}
    fills_by_order: dict[str, list[dict]] = {}
    for trade in trades:
        if trade.get("order_id"):
            fills_by_order.setdefault(str(trade["order_id"]), []).append(trade)
    for row in unconfirmed:
        client_order_id = row["client_order_id"]
        order = by_order_id.get(row["prophetx_order_id"] or "") or by_external_id.get(client_order_id)
        stake = row["count_requested"] / 100
        if order is None or not _done(order, stake):
            continue
        order_id = str(order["order_id"])
        filled, _to_win, fee = matched_of(order, fills_by_order.get(order_id, []))
        outcome = _outcome_of(order, filled, stake)
        store.update_manual_order(client_order_id, {
            "status": outcome, "prophetx_order_id": order_id, "count_filled": round(filled * 100),
            "cost_filled": round(filled + fee, 4), "error": f"settled from the order history: {_order_state(order)}"}, now)
        log.info("prophetx place %s: settled from the order history: %s, $%.2f of $%.2f (fee %s)",
                 client_order_id, outcome, filled, stake, fee)


def _error_text(reply: object) -> str:
    if isinstance(reply, dict):
        for key in ("message", "error"):
            if isinstance(reply.get(key), str) and reply[key]:
                return reply[key][:200]
    if isinstance(reply, str):
        return reply[:200]
    return ""


def _outcome_message(outcome: str, target: dict, stake: float, filled: float, fee: float) -> str:
    """A fill in the bet's own terms, like the other venues' Place messages. The
    to-win is at the order's limit price: ProphetX's order read names no fill
    price in the docs, and a fill-or-kill buy fills at that price or better."""
    price = f"{target['price']:+d}"
    if outcome == STATUS_REJECTED:
        return "ProphetX rejected the order (invalid, often the balance). Nothing was placed."
    if outcome not in (STATUS_FILLED, STATUS_PARTIAL):
        return f"Nothing filled at {price}: the price was gone. Nothing rests."
    rest = "." if outcome == STATUS_FILLED else "; the rest was not there. Nothing rests."
    payout = filled + profit_at(filled, target["price"])
    matched = "" if outcome == STATUS_FILLED else f"${filled:,.2f} of ${stake:,.2f} matched, "
    detail = f"{matched}limit {price} on ProphetX" + (f", fee ${fee:,.2f}" if fee else "")
    return filled_bet_message(target.get("label"), filled + fee, payout, detail) + rest
