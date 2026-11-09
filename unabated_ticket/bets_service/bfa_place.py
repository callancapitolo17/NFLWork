"""The Ticket's Place button at BFA (Cal 2026-11-02): ONE straight bet on Cal's BFA account, on
the BFA source's own login (sources/bfa.py), the same way the Teasers tab's Place sends a teaser
(bfa_teaser.py). Nothing here runs on a timer: every bet is one POST /place_bfa.json from the
panel or the phone page.

The wager is the teaser's POST (bfa_teaser.py documents the recorded shape) with one pick and
the straight type: WagerType 0 and IdWagerType the account's own straight type for the sport,
both read off the account metadata, whose betType rows give a teaser wagerTypes 2 (the recorded
teaser's WagerType) and a straight wagerTypes 0. The pick sends TeaserPoints 0 and Amount = the
stake, as the risk (RiskOrWinType / RiskWin RISK_OR_WIN_RISK, at any price). A straight bet has
not been recorded off the site, so the placed ticket is read back and checked against what was
sent (rotation, market, side, number, price, risk and to-win); a difference is said loudly.

Rules, in order, for every request (validate_bfa_place_request's shape):
  1. the click was not already sent: 'bman-' + its clickId, held in bets.duckdb::manual_orders
     before the send (venue 'bfa'); one placement at a time;
  2. the server runner names the board's one Buckeye line (GET /bfa_quote.json): its side's
     rotation, its start and the line as a bet record names it. A main full-game spread, total
     or moneyline in a league BFA's board is read for; an alt or a period line refuses (BFA's
     board lists only main full-game lines);
  3. the account's straight type for the sport, and the stake within its limit;
  4. BFA's board now: the game by rotation within START_WINDOW_SEC of the board's start, not
     started, the spread's or moneyline's team on the side the board has it, the market open
     at the SAME number (bfa_teaser.match_leg's checks);
  5. BFA's price is at most PRICE_TOLERANCE (0.1c) worse than the Ticket's; a better one is
     taken (AcceptChanges 0: BFA refuses a price or number that moves after this read);
  5b. a free play request (freePlay true, Cal 2026-11-09): the straight type offers free play,
     the price is within its odds limit, the balance read now covers the stake and the to-win
     is within its cap (bfa_free_play.py); the pick then sends UseFreePlay true and the ticket
     it confirms on is a "STRAIGHT BET (FP)" one;
  6. the other side of the line is not held at BFA, cash or free play (Place adds to a side,
     never hedges it);
  7. the BFA read the stake counted is ok and fresh, and no earlier Place on this game, at
     any venue, may have filled after the read of its venue the stake was sized on
     (manual_orders.unseen_manual_fill);
  8. the row is INSERTED (status 'sending') BEFORE the POST, which is sent once and never
     retried. Only a 4xx counts as refused; anything else may have booked it, so the open bets
     decide: a new straight on the game within CONFIRM_TIMEOUT_SEC is 'filled' (BFA books a
     straight whole or not at all), none is 'unknown' (check BFA before placing again).

Order statuses (manual_orders.status on a 'bfa' row): sending, filled, rejected, unknown.
Inputs:  the BFA source (its in-memory login), the store, a runner quote reader.
Outputs: {status: filled | rejected | unknown | refused, message, ticket?, legsMatch?,
         legCheck?, openWagers?}; openWagers is the open-bets list read at confirmation, for the
         caller to store.
Side effects: one straight wager at BFA when every check passes (real money); INSERT / UPDATE
bets.duckdb::manual_orders (venue 'bfa'). Network: GET <runner>/bfa_quote.json, then BFA's
metadata, board, open bets and wager endpoints on the BFA source's login.
"""
import logging
import math
import threading
import time
import uuid
from collections.abc import Callable
from datetime import datetime, timezone

from unabated_ticket.bets_service.bfa_teaser import (
    CONTESTANT_SIDES, FULL_GAME_LABEL, FULL_GAME_PERIOD, HTTP_CLIENT_ERRORS, HTTP_TIMEOUT_SEC, MAIN_LINE_INDEX,
    MARKET_OPEN, MONEY_TOLERANCE, ODDS_OPEN, RISK_OR_WIN_RISK, START_WINDOW_SEC, TOTAL_SIDES, WAGER_URL,
    games_with_rotation, is_free_play,
    is_number, main_market, one_line, parse_iso, read_account_metadata, read_board, risk_matches, short_text,
    spread_side_problem, ticket_of,
)
from unabated_ticket.bets_service.bfa_free_play import (
    FreePlayReader, free_play_problem, free_play_rules_of, read_free_play_balance,
)
from unabated_ticket.bets_service.manual_orders import (
    CLICK_ID_PATTERN, PRICE_TOLERANCE, STATUS_FILLED, STATUS_REFUSED, STATUS_REJECTED, STATUS_UNKNOWN, VENUE_BFA,
    _cents, _iso_time, optional_read_at, unseen_manual_fill,
)
from unabated_ticket.bets_service.prophetx_markets import american_to_probability, profit_at
from unabated_ticket.bets_service.sources.bfa import BFASource, new_session, parse_leg
from unabated_ticket.bets_service.store import BetsStore

log = logging.getLogger(__name__)

BFA_ORDER_ID_PREFIX = "bman-"
RUNNER_BFA_QUOTE_PATH = "/bfa_quote.json"
RUNNER_BFA_QUOTE_TIMEOUT_SEC = 10

# The board's league -> BFA's board slug and the sport code its metadata and odds rows use
# (bfa_odds/scraper.py reads these slugs). A league outside these refuses; a slug or code BFA
# does not know fails closed on the board read or the straight type.
BOARD_SLUGS = {"nfl": "nfl", "cfb": "ncaa_f", "nba": "nba", "cbb": "ncaa_b", "nhl": "nhl"}
SPORT_CODES = {"nfl": "NFL", "cfb": "CFB", "nba": "NBA", "cbb": "CBB", "nhl": "NHL"}

# metadata betType.wagerTypes and the wager's WagerType: 0 = straight (2 = teaser).
WAGER_TYPE_STRAIGHT = 0
STRAIGHT_TEAMS = 1
MARKET_MONEYLINE = 1
BET_TYPE_MARKETS = {"moneyline": MARKET_MONEYLINE, "spread": 2, "total": 3}
BET_TYPE_MONEYLINE = 1
STRAIGHT_HEADER_PREFIX = "STRAIGHT BET"

BFA_READ_SOURCE = "bfa"
BFA_READ_MAX_AGE_SEC = 5 * 60
CONFIRM_POLL_SEC = 1.0
CONFIRM_TIMEOUT_SEC = 20.0
METADATA_TTL_SEC = 3600.0
MIN_STAKE_DOLLARS = 1
OPPOSITE_SIDES = {"away": "home", "home": "away", "over": "under", "under": "over"}


# ---- the panel's request ------------------------------------------------------------

def validate_bfa_place_request(body: object) -> dict | str:
    """A POST /place_bfa.json body, or what was expected and what was found.
    {clickId, eventId, periodTypeId, betTypeId (1|2|3), sideIndex, points (null on a
    moneyline), stake (whole dollars: BFA takes no cents), freePlay (optional, default false:
    bet the account's free play, not its cash), ticketAllIn (the Ticket's BFA price
    as a probability), bfaReadAt (the BFA read the stake was sized on, ISO), kalshiReadAt /
    novigReadAt / polymarketReadAt / prophetxReadAt (ISO, or null when never read)}"""
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
    elif not is_number(points) or not math.isfinite(points):
        return f"points must be a number (a spread or total), got {points!r}"
    stake = body.get("stake")
    if not isinstance(stake, int) or isinstance(stake, bool) or stake < MIN_STAKE_DOLLARS:
        return f"stake must be whole dollars, at least ${MIN_STAKE_DOLLARS} (BFA takes no cents), got {stake!r}"
    free_play = body.get("freePlay", False)
    if not isinstance(free_play, bool):
        return f"freePlay must be true or false, got {free_play!r}"
    ticket_all_in = body.get("ticketAllIn")
    if not is_number(ticket_all_in) or not 0 < ticket_all_in < 1:
        return f"ticketAllIn must be a probability between 0 and 1, got {ticket_all_in!r}"
    bfa_read_at = _iso_time(body.get("bfaReadAt"))
    if bfa_read_at is None:
        return f"bfaReadAt must be the ISO time of the BFA read the stake was sized on, got {body.get('bfaReadAt')!r}"
    read_ats = {}
    for field in ("kalshiReadAt", "novigReadAt", "polymarketReadAt", "prophetxReadAt"):
        read_at = optional_read_at(body, field)
        if isinstance(read_at, str):
            return read_at
        read_ats[field] = read_at
    return {"click_id": click_id, **whole, "points": None if points is None else float(points), "stake": stake,
            "free_play": free_play, "ticket_all_in": float(ticket_all_in), "bfa_read_at": bfa_read_at,
            "kalshi_read_at": read_ats["kalshiReadAt"], "novig_read_at": read_ats["novigReadAt"],
            "polymarket_read_at": read_ats["polymarketReadAt"], "prophetx_read_at": read_ats["prophetxReadAt"]}


# ---- the account's straight type ----------------------------------------------------

def straight_type_of(metadata: dict, sport_code: str) -> dict | str:
    """{wagerTypeId, maxRisk} of the account's straight bet for `sport_code`, or why there is
    none: its metadata betType rows with wagerTypes 0 and one team, for that sport, must name
    exactly one type."""
    rows = metadata.get("betType")
    if not isinstance(rows, list):
        return f"BFA's account metadata carries no betType list (keys: {sorted(metadata)[:12]})"
    straights = [row for row in rows if isinstance(row, dict) and row.get("wagerTypes") == WAGER_TYPE_STRAIGHT
                 and row.get("numTeams") == STRAIGHT_TEAMS]
    found = [row for row in straights if row.get("sportId") == sport_code]
    type_ids = sorted({row.get("wagerTypeId") for row in found})
    if len(type_ids) != 1 or not isinstance(type_ids[0], int):
        listed = sorted({str(row.get("sportId")) for row in straights})
        return (f"expected one straight-bet type for {sport_code} in BFA's account metadata, found "
                f"{type_ids or 'none'} (straight types listed for {listed or 'no sport'})")
    limits = [row.get("maxRisk") for row in found if is_number(row.get("maxRisk")) and row.get("maxRisk") > 0]
    return {"wagerTypeId": type_ids[0], "maxRisk": min(limits) if limits else None}


# ---- the board -------------------------------------------------------------------------

def _number_text(bet_type: str, points: float | None) -> str:
    if points is None:
        return "no number"
    text = f"{points:g}"
    return f"+{text}" if bet_type == "spread" and points > 0 else text


def _odds_row(market: dict, line: dict, contestant: dict) -> dict | None:
    for odds in market.get("odds") or []:
        if odds.get("index") != MAIN_LINE_INDEX:
            continue
        if line["betType"] == "total":
            if odds.get("side") == TOTAL_SIDES[line["side"]]:
                return odds
        elif odds.get("contestantId") == contestant.get("id"):
            return odds
    return None


def match_straight(games: list[dict], line: dict, rotation: int, event_start: datetime, now: datetime,
                   label: str) -> dict | str:
    """BFA's game, market and odds row for one straight bet, or why it cannot be bet as the
    Ticket shows it: {game, fixture, market, odds}. `line` is the runner's {betType, side,
    points}; `rotation` the side's own (a total's: either team's, it finds the game)."""
    found = [item for item in games_with_rotation(games, rotation)
             if (start := parse_iso(item[1].get("date"))) is not None
             and abs((start - event_start).total_seconds()) <= START_WINDOW_SEC]
    if not found:
        return f"{label}: BFA's board lists no game with rotation {rotation}"
    if len(found) > 1:
        return f"{label}: BFA's board lists {len(found)} games with rotation {rotation}"
    game, fixture, contestant = found[0]
    start = parse_iso(fixture.get("date"))
    if start is None or start <= now or fixture.get("isLive") or game.get("isLive"):
        return f"{label}: the game has started"
    # spread_side_problem reads a spread leg; a moneyline's team sits on its side the same way.
    side_problem = None if line["betType"] == "total" else spread_side_problem(
        {"betType": "spread", "side": line["side"], "label": label, "rotation": rotation}, contestant)
    if side_problem:
        return side_problem
    market = main_market(game, fixture, BET_TYPE_MARKETS[line["betType"]])
    odds = _odds_row(market, line, contestant) if market else None
    if market is None or odds is None:
        return f"{label}: BFA lists no full-game {line['betType']} for {game.get('name')}"
    if line["betType"] != "total" and odds.get("side") != CONTESTANT_SIDES[line["side"]]:
        return (f"{label}: BFA's odds row for {contestant.get('name')} carries side {odds.get('side')!r}, "
                f"the {line['side']} team's is {CONTESTANT_SIDES[line['side']]}")
    if market.get("status") != MARKET_OPEN or odds.get("status") != ODDS_OPEN:
        return f"{label}: BFA has the {line['betType']} closed"
    if line["betType"] != "moneyline" and (not is_number(odds.get("line")) or odds["line"] != line["points"]):
        return (f"{label}: BFA has {_number_text(line['betType'], odds.get('line'))} now, "
                f"the Ticket has {_number_text(line['betType'], line['points'])}")
    if (not is_number(odds.get("price")) or odds["price"] == 0 or not odds.get("dSpId") or odds.get("dGtId") is None
            or odds.get("dGmId") is None):
        return f"{label}: BFA's odds row carries no price, sport, game type or game id"
    return {"game": game, "fixture": fixture, "market": market, "odds": odds}


def straight_pick_of(match: dict, stake: int, free_play: bool = False) -> dict:
    """The one Picks[] entry: the teaser's pick (bfa_teaser.pick_of) with no points bought,
    risking `stake` (of the account's free play when `free_play`)."""
    odds = match["odds"]
    return {
        "EventId": match["game"]["id"],
        "FixtureId": match["fixture"]["id"],
        "MarketType": match["market"]["type"],
        "PeriodNumber": FULL_GAME_PERIOD,
        "Side": odds["side"],
        "Index": MAIN_LINE_INDEX,
        "ContestantId": odds.get("contestantId"),
        "Line": odds.get("line"),
        "TeaserPoints": 0.0,
        "Price": odds["price"],
        "Amount": float(stake),
        "PointsPurchased": 0.0,
        "PitcherAction": 0,
        "RoundRobinCombinations": 0,
        "UseFreePlay": free_play,
        "RiskOrWinType": RISK_OR_WIN_RISK,
        "IdGameType": odds["dGtId"],
        "IdSport": odds["dSpId"],
    }


def straight_wager_body(transaction_id: str, player_id: str, wager_type_id: int, pick: dict) -> list[dict]:
    """The POST body: one straight wager on the pick (its Amount the risk), no line changes accepted."""
    return [{
        "IdTransaction": transaction_id,
        "IdPlayer": int(player_id),
        "FillIdWager": -1,
        "WagerType": WAGER_TYPE_STRAIGHT,
        "OpenSpots": 0,
        "Picks": [pick],
        "RiskWin": RISK_OR_WIN_RISK,
        "AcceptChanges": 0,
        "IdWagerType": wager_type_id,
        "IsLive": False,
        "PhoneLine": None,
        "UserId": None,
    }]


# ---- the open bets ---------------------------------------------------------------------

def is_straight(wager: dict) -> bool:
    """A straight bet, cash or free play ("STRAIGHT BET (FP)")."""
    return str(wager.get("headerDescription") or "").upper().startswith(STRAIGHT_HEADER_PREFIX)


def _straight_leg(wager: dict) -> dict | None:
    """The one leg of a straight wager as parse_leg reads it, with its BFA game id, or None."""
    details = wager.get("betDetails") or []
    if not is_straight(wager) or len(details) != 1:
        return None
    leg = parse_leg(str(details[0].get("detailDescription") or ""))
    return {**leg, "gameId": details[0].get("idGame")} if isinstance(leg, dict) else None


def holds_other_side(open_wagers: list[dict], game_id: object, line: dict) -> bool:
    """Whether an open straight on this BFA game is the other side of the same full-game market."""
    opposite = OPPOSITE_SIDES[line["side"]]
    for wager in open_wagers:
        leg = _straight_leg(wager)
        if (leg and leg["gameId"] == game_id and leg["betType"] == line["betType"]
                and leg["period"] == FULL_GAME_LABEL and leg["side"] == opposite):
            return True
    return False


def new_straights_on_game(open_wagers: list[dict], game_id: object, known_ids: set,
                          free_play: bool = False) -> list[dict]:
    """The straights not open before the POST whose one leg is on this BFA game, bet with free
    play when `free_play` and with cash when not."""
    return [wager for wager in open_wagers if wager.get("idWager") not in known_ids
            and is_straight(wager) and is_free_play(wager) == free_play
            and [detail.get("idGame") for detail in wager.get("betDetails") or []] == [game_id]]


def straight_mismatches(wager: dict, want: dict, stake: int, free_play: bool = False) -> list[str]:
    """Every way a placed straight differs from what was sent, in words; [] when it is exactly
    it. `want` {rotations (the ones its leg may read as written), betType, side, points, price}.
    A free play's risk may read 0 (bfa_teaser.risk_matches)."""
    details = wager.get("betDetails") or []
    text = str(details[0].get("detailDescription") or "") if len(details) == 1 else ""
    got = parse_leg(text) if text else "no leg"
    problems = []
    if isinstance(got, str):
        problems.append(f"BFA's ticket could not be read ({got}): {one_line(text)}")
    elif (got["rotationAsWritten"] not in want["rotations"] or got["betType"] != want["betType"]
          or got["period"] != FULL_GAME_LABEL or got["side"] != want["side"] or got["points"] != want["points"]
          or got["price"] != want["price"]):
        problems.append(f"BFA has {one_line(text)}")
    risk, to_win = wager.get("riskAmount"), wager.get("winAmount")
    expected_win = profit_at(stake, want["price"])
    # BFA rounds the to-win to the cent, so it may land a cent off profit_at's.
    if (not risk_matches(risk, stake, free_play) or not is_number(to_win)
            or abs(to_win - expected_win) > 0.01 + MONEY_TOLERANCE):
        problems.append(f"BFA's ticket risks {risk} to win {to_win}, ${stake} to win ${expected_win:.2f} was sent")
    return problems


# ---- the placer ------------------------------------------------------------------------

class BFAStraightPlacer:
    """Places Cal's BFA Place-button straight bets, one at a time, on the BFA source's login.

    `quote_reader(request)` reads the runner's /bfa_quote.json; `session_factory()` builds an
    HTTP session with .get()/.post() (tests inject a fake); `clock()` / `sleep()` pace the
    confirmation reads; `new_transaction_id()` mints the wager's IdTransaction. The account
    metadata is held in memory for METADATA_TTL_SEC.
    """

    def __init__(self, bfa_source: BFASource, store: BetsStore, quote_reader: Callable[[dict], dict],
                 session_factory: Callable[[], object] = new_session,
                 now: Callable[[], datetime] = lambda: datetime.now(timezone.utc),
                 clock: Callable[[], float] = time.monotonic, sleep: Callable[[float], None] = time.sleep,
                 new_transaction_id: Callable[[], str] = lambda: str(uuid.uuid4())):
        self._bfa = bfa_source
        self._store = store
        self._quote_reader = quote_reader
        self._session_factory = session_factory
        self._now = now
        self._clock = clock
        self._sleep = sleep
        self._new_transaction_id = new_transaction_id
        self._placing = threading.Lock()
        self._metadata: dict | None = None
        self._metadata_at = 0.0
        self._free_play_reader = FreePlayReader(bfa_source, session_factory, read_account_metadata, clock)

    def free_play(self) -> dict:
        """GET /bfa_free_play.json's body (bfa_free_play.FreePlayReader.read). Never raises."""
        return self._free_play_reader.read()

    def place(self, request: dict) -> dict:
        """One validated request (validate_bfa_place_request) -> {status, message, ...}."""
        if not self._placing.acquire(blocking=False):
            return _refused("another bet is being placed at BFA right now")
        try:
            return self._place(request)
        finally:
            self._placing.release()

    def _place(self, request: dict) -> dict:
        client_order_id = f"{BFA_ORDER_ID_PREFIX}{request['click_id']}"
        quote = self._quote_reader({"eventId": request["eventId"], "periodTypeId": request["periodTypeId"],
                                    "betTypeId": request["betTypeId"], "sideIndex": request["sideIndex"],
                                    "points": request["points"]})
        board_line = _board_line_problem(quote)
        if isinstance(board_line, str):
            return _refused(board_line)
        line, label = board_line["line"], board_line["label"]
        # BFA's own start is checked on its board; the board's is the second clock.
        if board_line["eventStart"] <= self._now():
            return _refused(f"{label}: the game has started")
        headers, player_id = self._bfa.authorized()
        metadata = self._account_metadata(headers, player_id)
        straight_type = straight_type_of(metadata, SPORT_CODES[line["league"]])
        if isinstance(straight_type, str):
            return _refused(straight_type)
        if straight_type["maxRisk"] is not None and request["stake"] > straight_type["maxRisk"]:
            return _refused(f"${request['stake']} is over BFA's ${straight_type['maxRisk']:g} limit on this bet")
        games = read_board(self._session_factory(), BOARD_SLUGS[line["league"]], line["league"], headers, player_id,
                           metadata.get("agentId"))
        match = match_straight(games, line, board_line["rotation"], board_line["eventStart"], self._now(), label)
        if isinstance(match, str):
            return _refused(match)
        odds = match["odds"]
        probability = round(american_to_probability(odds["price"]), 6)
        if probability > request["ticket_all_in"] + PRICE_TOLERANCE:
            return _refused(f"BFA is at {odds['price']:+g} now ({_cents(probability)}), more than 0.1¢ worse than the "
                            f"Ticket's {_cents(request['ticket_all_in'])}")
        if request["free_play"]:
            refusal = self._free_play_refusal(metadata, straight_type["wagerTypeId"], SPORT_CODES[line["league"]],
                                              request["stake"], odds["price"], headers, player_id)
            if refusal:
                return _refused(refusal)
        open_before = self._bfa.open_wagers()
        if holds_other_side(open_before, odds["dGmId"], line):
            return _refused("you hold the other side of this line at BFA")
        refusal = self._stale_or_unseen(request)
        if refusal:
            return _refused(refusal)
        transaction_id = self._new_transaction_id()
        body = straight_wager_body(transaction_id, player_id, straight_type["wagerTypeId"],
                                   straight_pick_of(match, request["stake"], request["free_play"]))
        inserted = self._store.insert_manual_order({
            "client_order_id": client_order_id, "placed_at": self._now(), "venue": VENUE_BFA,
            "event_id": request["eventId"], "ticker": str(odds["dGmId"]), "contract_side": line["side"],
            "label": label, "ticket_all_in": request["ticket_all_in"], "ask_price": probability,
            "all_in_price": probability, "stake": float(request["stake"]), "count_requested": request["stake"],
            "status": "sending", "updated_at": self._now(), "event_start": board_line["eventStart"],
            "market_id": str(match["market"].get("id")), "league": line["league"], "bet_type": line["betType"],
            "period": line["period"], "side": line["side"], "points": line.get("points"),
            "away_team": line.get("awayTeam"), "home_team": line.get("homeTeam"),
            "away_team_id": line.get("awayTeamId"), "home_team_id": line.get("homeTeamId"),
            "free_play": request["free_play"],
        })
        if not inserted:
            return _refused("this click was already sent")
        log.info("bfa place %s: %s $%d%s at %+g (%s, Ticket %s), transaction %s", client_order_id, label,
                 request["stake"], " free play" if request["free_play"] else "", odds["price"], _cents(probability),
                 _cents(request["ticket_all_in"]), transaction_id)
        try:
            posted = self._post_wager(body, headers, player_id)
        except Exception as error:  # noqa: BLE001 — a lost reply is not a refusal: the wager may stand
            log.warning("bfa place %s: POST failed: %s: %s", client_order_id, type(error).__name__, error)
            posted = {"refused": None, "note": f"the request to BFA failed ({type(error).__name__})"}
        if posted["refused"]:
            self._update(client_order_id, {"status": STATUS_REJECTED, "error": posted["refused"]})
            log.warning("bfa place %s refused: %s", client_order_id, posted["refused"])
            return {"status": STATUS_REJECTED, "message": f"{posted['refused']}. Nothing was placed."}
        want = {"rotations": _rotations_of(board_line), "betType": line["betType"], "side": line["side"],
                "points": line.get("points") if line["betType"] != "moneyline" else None, "price": odds["price"]}
        try:
            return self._confirm(client_order_id, odds["dGmId"], want, request["stake"], request["free_play"], label,
                                 {wager.get("idWager") for wager in open_before}, posted["note"])
        except Exception as error:  # noqa: BLE001 — the wager went out: never "Not placed" from here
            log.exception("bfa place %s: confirming failed", client_order_id)
            self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": f"confirm failed: {error}"})
            return {"status": STATUS_UNKNOWN,
                    "message": (f"Sent to BFA, but reading the open bets failed ({type(error).__name__}). "
                                "Check BFA's open bets before placing again.")}

    def _account_metadata(self, headers: dict, player_id: str) -> dict:
        if self._metadata is not None and self._clock() - self._metadata_at < METADATA_TTL_SEC:
            return self._metadata
        body = read_account_metadata(self._session_factory(), headers, player_id)
        self._metadata, self._metadata_at = body, self._clock()
        return body

    def _free_play_refusal(self, metadata: dict, wager_type_id: int, sport_code: str, stake: int, price: float,
                           headers: dict, player_id: str) -> str | None:
        """Why this straight cannot go out on free play (bfa_free_play), or None. The balance
        is read now: free play spent since the Ticket last read it must not go out twice."""
        rules = free_play_rules_of(metadata, wager_type_id, {sport_code})
        if isinstance(rules, str):
            return rules
        balance = read_free_play_balance(self._session_factory(), headers, player_id)
        return free_play_problem(rules, balance, stake, profit_at(stake, price), price)

    def _stale_or_unseen(self, request: dict) -> str | None:
        read = self._store.source_status().get(BFA_READ_SOURCE)
        fetched_at = _iso_time(read.get("fetchedAt")) if read else None
        if not read or not read.get("ok") or fetched_at is None:
            return "the BFA read is failing, so the stake may miss a bet. Try again once it reads"
        if (self._now() - fetched_at).total_seconds() > BFA_READ_MAX_AGE_SEC:
            return f"the BFA read is over {BFA_READ_MAX_AGE_SEC // 60} minutes old, so the stake may miss a bet"
        if unseen_manual_fill(self._store, request["eventId"], request["kalshi_read_at"], request["novig_read_at"],
                              request["polymarket_read_at"], request["prophetx_read_at"], request["bfa_read_at"]):
            return ("your last Place on this game is not on the Bets tab yet, so this stake does not count it. "
                    "Wait for the next read (up to a minute) and place again")
        return None

    def _post_wager(self, body: list[dict], headers: dict, player_id: str) -> dict:
        """{refused, note}, as bfa_teaser's: only a 4xx is a refusal; a network error is raised."""
        response = self._session_factory().post(WAGER_URL, params={"playerId": player_id}, json=body,
                                                headers={**headers, "Content-Type": "application/json; charset=utf-8"},
                                                timeout=HTTP_TIMEOUT_SEC)
        reply = None
        try:
            reply = response.json()
        except ValueError:
            pass
        shown = short_text(str(reply if reply is not None else response.text))
        if HTTP_CLIENT_ERRORS[0] <= response.status_code <= HTTP_CLIENT_ERRORS[1]:
            return {"refused": f"BFA refused the bet: HTTP {response.status_code} {shown}", "note": None}
        if response.status_code == 200 and isinstance(reply, dict) and reply.get("state") == 0:
            return {"refused": None, "note": None}
        return {"refused": None, "note": f"BFA answered HTTP {response.status_code} {shown}"}

    def _confirm(self, client_order_id: str, game_id: object, want: dict, stake: int, free_play: bool, label: str,
                 known_ids: set, note: str | None) -> dict:
        deadline = self._clock() + CONFIRM_TIMEOUT_SEC
        while True:
            self._sleep(CONFIRM_POLL_SEC)
            try:
                open_wagers = self._bfa.open_wagers()
            except Exception as error:  # noqa: BLE001 — a failed read is one more try, not a failed wager
                log.warning("bfa place %s: open bets read failed while confirming: %s", client_order_id, error)
                open_wagers = []
            candidates = new_straights_on_game(open_wagers, game_id, known_ids, free_play)
            if candidates:
                return self._placed(client_order_id, candidates, want, stake, free_play, label, open_wagers)
            if self._clock() >= deadline:
                log.warning("bfa place %s: sent, no ticket after %.0f s (%s)", client_order_id, CONFIRM_TIMEOUT_SEC,
                            note or "BFA took it for processing")
                self._update(client_order_id, {"status": STATUS_UNKNOWN, "error": note or "no ticket showed"})
                reply = f" ({note})" if note else ""
                return {"status": STATUS_UNKNOWN,
                        "message": (f"Sent to BFA{reply}, but no ticket showed in {CONFIRM_TIMEOUT_SEC:.0f} s. "
                                    "Check BFA's open bets before placing again."),
                        "openWagers": open_wagers}

    def _placed(self, client_order_id: str, candidates: list[dict], want: dict, stake: int, free_play: bool,
                label: str, open_wagers: list[dict]) -> dict:
        """The filled result for the new ticket that reads exactly as sent, else the first new one
        on the game with every difference named (a booked bet cannot be undone). A free play
        costs no cash: its row's cost_filled is 0 and its count_filled the free play risked."""
        checked = [(wager, straight_mismatches(wager, want, stake, free_play)) for wager in candidates]
        wager, problems = next(((wager, found) for wager, found in checked if not found), checked[0])
        ticket = ticket_of(wager)
        risk = wager.get("riskAmount") if is_number(wager.get("riskAmount")) else stake
        to_win = wager.get("winAmount") if is_number(wager.get("winAmount")) else profit_at(stake, want["price"])
        risked = float(stake) if free_play else risk
        self._update(client_order_id, {"status": STATUS_FILLED, "count_filled": risked,
                                       "cost_filled": 0.0 if free_play else risk, "error": "; ".join(problems) or None})
        result = {"status": STATUS_FILLED, "ticket": ticket, "legsMatch": not problems, "legCheck": problems,
                  "openWagers": open_wagers}
        paid_with = " of free play" if free_play else ""
        bet = f"Bet ${risked:,.2f}{paid_with} on {label} at {want['price']:+g} to win ${to_win:,.2f}"
        if problems:
            log.error("bfa place %s: placed as ticket %s BUT it differs from what was sent: %s", client_order_id,
                      ticket["ticketNumber"], " | ".join(problems))
            result["message"] = (f"{bet} (BFA ticket {ticket['ticketNumber']}), but BFA's ticket differs from the "
                                 f"one sent: {'; '.join(problems)}. Check it at BFA.")
            return result
        log.info("bfa place %s: placed as ticket %s, as sent", client_order_id, ticket["ticketNumber"])
        result["message"] = f"{bet} (BFA ticket {ticket['ticketNumber']})."
        return result

    def _update(self, client_order_id: str, fields: dict) -> None:
        self._store.update_manual_order(client_order_id, fields, self._now())


def _board_line_problem(quote: dict) -> dict | str:
    """{line, label, rotation, eventStart} of the runner's quote when BFA Place can bet it, or why not."""
    if not quote.get("ok"):
        return f"{quote.get('detail') or quote.get('reason') or 'no BFA line for this'}"
    line = quote.get("line") if isinstance(quote.get("line"), dict) else None
    if line is None or line.get("betType") not in BET_TYPE_MARKETS or line.get("side") not in OPPOSITE_SIDES:
        return "the runner did not name the board line"
    label = quote.get("label") or f"{line.get('side')} {line.get('betType')}"
    if quote.get("isAlt"):
        return f"{label} is an alt line: BFA Place bets main lines only (BFA's board lists no alts)"
    if line.get("period") != FULL_GAME_LABEL:
        return f"{label} is a {line.get('period')} line: BFA Place bets full-game lines only"
    if line.get("league") not in BOARD_SLUGS:
        return f"BFA Place does not read BFA's {str(line.get('league')).upper()} board"
    rotation = quote.get("rotation")
    if not isinstance(rotation, int) or isinstance(rotation, bool) or rotation <= 0:
        return f"{label}: the board has no rotation number for this side"
    event_start = _iso_time(quote.get("eventStart"))
    if event_start is None:
        return f"{label}: the board has no start for this game"
    return {"line": line, "label": label, "rotation": rotation, "eventStart": event_start,
            "awayRotation": quote.get("awayRotation"), "homeRotation": quote.get("homeRotation")}


def _rotations_of(board_line: dict) -> set:
    """The rotations the placed leg may read as written: the side's own; a total's either team's."""
    if board_line["line"]["betType"] != "total":
        return {board_line["rotation"]}
    return {rotation for rotation in (board_line["rotation"], board_line["awayRotation"], board_line["homeRotation"])
            if isinstance(rotation, int)}


def _refused(message: str) -> dict:
    return {"status": STATUS_REFUSED, "message": f"Not placed: {message}."}
