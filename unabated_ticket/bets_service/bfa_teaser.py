"""Place a Buckeye 4-team 6-point teaser at BFA from the Teasers tab (POST /place_teaser.json).

BFA IS Buckeye. One wager is one JSON POST, recorded 2026-10-03 off a ticket Cal placed by
hand (ticket 356323496) — POST api.bfagaming.com/wagering/api/v1/wager?playerId=<id> with the
account's Keycloak bearer token and a ONE-element array:
  [{IdTransaction <uuid: the idempotency key>, IdPlayer, FillIdWager -1, WagerType 2 (teaser),
    OpenSpots 0, RiskWin 1 (the amount is the risk), AcceptChanges 0 (BFA refuses a moved number
    rather than taking the new one), IdWagerType <the account's own teaser type>, IsLive false,
    PhoneLine null, UserId null,
    Picks: [{EventId, FixtureId, MarketType 2 spread | 3 total, PeriodNumber 0, Side, Index 0,
             ContestantId (null on a total), Line <the UNTEASED number>, TeaserPoints 6,
             Price <BFA's juice>, Amount <the ticket stake, on every pick>, PointsPurchased 0,
             PitcherAction 0, RoundRobinCombinations 0, UseFreePlay false, RiskOrWinType 1,
             IdGameType, IdSport}]}]
No preflight call and no captcha. The reply is only {"state": 0} — taken for processing; the
outcome reaches the site ~1 s later on a SignalR hub this module does not speak, so it reads the
account's open bets until the new ticket shows. A total leg's pick is built the same way from its
odds row (side 4 over / 5 under, no contestant); the recording held spreads only.

place(request), in order — every check before the POST refuses with nothing sent:
  1. validate_place_request: 4 legs on 4 different games, NFL or CFB, a whole-dollar stake.
  2. one placement at a time, and never the same legs twice within RECENT_TICKET_SEC of a
     placement that went out (placed or unconfirmed).
  3. the account's teaser type from its metadata (cached METADATA_TTL_SEC): the one 4-team type
     at TEASER_POINTS for every leg's sport, paying +300 — what teaser.js prices.
  4. BFA's board now (the leagues' popular lists, as the site reads them): each leg's game by
     rotation, its main full-game market, open, and the SAME number the list was built on.
  5. the open bets before: an open teaser on the same four games and sides refuses the ticket.
  6. POST once, never retried. Only a 4xx reply counts as refused; a lost reply or any other
     answer may still have booked the wager, so the open bets decide.
  7. the open bets every CONFIRM_POLL_SEC until a new teaser on those four games shows
     (placed) or CONFIRM_TIMEOUT_SEC passes (unconfirmed: check BFA before trying again).

Inputs:  the BFA source (sources/bfa.py: its in-memory Keycloak session, one login per
         process) and the panel's request (validate_place_request).
Outputs: {status: placed | refused | unconfirmed, message, ticket?, openWagers?}; openWagers is
         the open-bets list read at confirmation, for the caller to store.
Side effects: one wager at BFA when every check passes — real money. Nothing on disk here; the
service's handler UPSERTs the returned open bets into bets.duckdb::bets so the panel sees the
ticket on its next poll.
"""
import logging
import re
import threading
import time
import uuid
from collections.abc import Callable
from datetime import datetime, timezone

from unabated_ticket.bets_service.sources.bfa import BFASource, new_session

log = logging.getLogger(__name__)

API_BASE = "https://api.bfagaming.com"
WAGER_URL = f"{API_BASE}/wagering/api/v1/wager"
# Lower case on purpose: the client's own GetPlayerMetadataByPlayerId spelling 404s here.
PLAYER_METADATA_URL = f"{API_BASE}/metadata/api/getplayermetadatabyplayerid"
BOARD_URL = f"{API_BASE}/oddsservice/events/popular/{{slug}}"
HTTP_TIMEOUT_SEC = 20

# The panel's league -> BFA's board slug and sport code.
BOARD_SLUGS = {"nfl": "nfl", "cfb": "ncaa_f"}
SPORT_CODES = {"nfl": "NFL", "cfb": "CFB"}
# A pregame board, as the site asks for it.
BOARD_FIXTURE_TYPE = 1

LEGS_PER_TICKET = 4
TEASER_POINTS = 6
# teaser.js prices a 4-team ticket at +300; any other payout makes its math wrong.
TICKET_PAYOUT_AMERICAN = 300
# metadata betType.wagerTypes and the wager's WagerType: 2 = teaser.
WAGER_TYPE_TEASER = 2
# The teaser points a type gives each sport's sides and totals.
TYPE_POINTS_FIELDS = {"NFL": ("nflSide", "nflTotal"), "CFB": ("cfbSide", "cfbTotal")}

MARKET_SPREAD = 2
MARKET_TOTAL = 3
FULL_GAME_PERIOD = 0
MAIN_LINE_INDEX = 0
# An odds row's side: a contestant's 1 home / 2 away; a total's 4 over / 5 under.
SIDE_OVER = 4
SIDE_UNDER = 5
MARKET_OPEN = 1
ODDS_OPEN = 0
BET_TYPE_MARKETS = {"spread": MARKET_SPREAD, "total": MARKET_TOTAL}
SPREAD_SIDES = ("away", "home")
TOTAL_SIDES = {"over": SIDE_OVER, "under": SIDE_UNDER}
# bets.js's own "start time differs" window: BFA's start against the list's.
START_WINDOW_SEC = 12 * 3600
# A rotation is a positive whole number; BFA's run past six digits for college.
MAX_ROTATION = 10_000_000
MAX_LABEL_CHARS = 120

# A 4xx books nothing; any other reply than {"state": 0} is left to the open bets to decide.
HTTP_CLIENT_ERRORS = (400, 499)
CONFIRM_POLL_SEC = 1.0
CONFIRM_TIMEOUT_SEC = 20.0
METADATA_TTL_SEC = 3600.0
# How long a ticket that went out blocks the same legs (the panel drops it once BFA lists it).
RECENT_TICKET_SEC = 15 * 60
TEASER_HEADER_RE = re.compile(r"TEASER", re.IGNORECASE)
DETAIL_ROTATION_RE = re.compile(r"\[(\d+)\]")
MAX_ERROR_TEXT = 300

STATUS_PLACED = "placed"
STATUS_REFUSED = "refused"
STATUS_UNCONFIRMED = "unconfirmed"


# ---- the panel's request ------------------------------------------------------------

def _is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and value == value


def _leg_error(index: int, raw: object) -> str | None:
    """Why `legs[index]` is not a leg, or None."""
    if not isinstance(raw, dict):
        return f"legs[{index}] must be an object, got {type(raw).__name__}"
    if raw.get("league") not in BOARD_SLUGS:
        return f"legs[{index}].league must be one of {sorted(BOARD_SLUGS)}, got {raw.get('league')!r}"
    bet_type = raw.get("betType")
    if bet_type not in BET_TYPE_MARKETS:
        return f"legs[{index}].betType must be one of {sorted(BET_TYPE_MARKETS)}, got {bet_type!r}"
    sides = SPREAD_SIDES if bet_type == "spread" else tuple(TOTAL_SIDES)
    if raw.get("side") not in sides:
        return f"legs[{index}].side must be one of {list(sides)} on a {bet_type}, got {raw.get('side')!r}"
    rotation = raw.get("rotation")
    if not isinstance(rotation, int) or isinstance(rotation, bool) or not 0 < rotation < MAX_ROTATION:
        return f"legs[{index}].rotation must be a positive whole number, got {rotation!r}"
    if not _is_number(raw.get("points")):
        return f"legs[{index}].points must be a number (Buckeye's number before the teaser), got {raw.get('points')!r}"
    start = raw.get("eventStart")
    if start is not None and _parse_iso(start) is None:
        return f"legs[{index}].eventStart must be an ISO time or null, got {start!r}"
    label = raw.get("label")
    if not isinstance(label, str) or not label or len(label) > MAX_LABEL_CHARS:
        return f"legs[{index}].label must be a non-empty string of at most {MAX_LABEL_CHARS} characters"
    return None


def _parse_iso(value: object) -> datetime | None:
    if not isinstance(value, str):
        return None
    try:
        moment = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return moment if moment.tzinfo else None


def validate_place_request(body: object) -> dict | str:
    """{stake, legs} of a POST /place_teaser.json body, or an error naming the first problem.

    body {stake: whole dollars, legs: [{league "nfl"|"cfb", betType "spread"|"total",
          side "away"|"home"|"over"|"under", rotation (the side's own, a total's either team's),
          points (Buckeye's number BEFORE the teaser), eventStart (ISO) or null, label}]}
    Four legs on four different games (teaser.js: one leg per game)."""
    if not isinstance(body, dict):
        return f"body must be an object, got {type(body).__name__}"
    stake = body.get("stake")
    if not isinstance(stake, int) or isinstance(stake, bool) or stake <= 0:
        return f"stake must be a whole number of dollars above 0, got {stake!r}"
    legs = body.get("legs")
    if not isinstance(legs, list) or len(legs) != LEGS_PER_TICKET:
        count = len(legs) if isinstance(legs, list) else type(legs).__name__
        return f"legs must be a list of {LEGS_PER_TICKET} legs, got {count}"
    for index, raw in enumerate(legs):
        error = _leg_error(index, raw)
        if error:
            return error
    clean = [{"league": raw["league"], "betType": raw["betType"], "side": raw["side"], "rotation": raw["rotation"],
              "points": float(raw["points"]), "eventStart": raw.get("eventStart"), "label": raw["label"]}
             for raw in legs]
    return {"stake": stake, "legs": clean}


def ticket_key(legs: list[dict]) -> tuple:
    """The same four sides at the same numbers, whatever order the panel sent them in."""
    return tuple(sorted((leg["league"], leg["rotation"], leg["betType"], leg["side"], leg["points"]) for leg in legs))


# ---- the account's teaser type ------------------------------------------------------

def teaser_type_of(metadata: dict, sport_codes: set[str]) -> dict | str:
    """{wagerTypeId, maxRisk} of the account's 4-team TEASER_POINTS teaser covering every
    sport in `sport_codes`, or why there is none. The metadata's betType lists one row per
    (type, sport); a type counts for a sport only when its sides and totals both tease
    TEASER_POINTS there, and it must pay +300 (what the Teasers tab prices)."""
    rows = metadata.get("betType")
    if not isinstance(rows, list):
        return f"BFA's account metadata carries no betType list (keys: {sorted(metadata)[:12]})"
    by_type: dict[int, dict] = {}
    for row in rows:
        if not isinstance(row, dict) or row.get("wagerTypes") != WAGER_TYPE_TEASER:
            continue
        if row.get("numTeams") != LEGS_PER_TICKET or row.get("sportId") not in TYPE_POINTS_FIELDS:
            continue
        side_field, total_field = TYPE_POINTS_FIELDS[row["sportId"]]
        if row.get(side_field) != TEASER_POINTS or row.get(total_field) != TEASER_POINTS:
            continue
        entry = by_type.setdefault(row["wagerTypeId"], {"sports": set(), "payouts": set(), "maxRisk": []})
        entry["sports"].add(row["sportId"])
        entry["payouts"].add(row.get("payOuts"))
        entry["maxRisk"].append(row.get("maxRisk"))
    covering = {type_id: entry for type_id, entry in by_type.items() if sport_codes <= entry["sports"]}
    if len(covering) != 1:
        return (f"expected one {LEGS_PER_TICKET}-team {TEASER_POINTS}-point teaser type for "
                f"{sorted(sport_codes)} in BFA's account metadata, found {sorted(covering) or 'none'}")
    type_id, entry = next(iter(covering.items()))
    if entry["payouts"] != {TICKET_PAYOUT_AMERICAN}:
        return (f"BFA's {LEGS_PER_TICKET}-team teaser (type {type_id}) pays {sorted(entry['payouts'])}, "
                f"the Teasers tab prices +{TICKET_PAYOUT_AMERICAN}")
    limits = [limit for limit in entry["maxRisk"] if _is_number(limit) and limit > 0]
    return {"wagerTypeId": type_id, "maxRisk": min(limits) if limits else None}


# ---- the board -------------------------------------------------------------------------

def _signed(points: float) -> str:
    text = f"{points:g}"
    return f"+{text}" if points > 0 else text


def _number_text(bet_type: str, points: float) -> str:
    return _signed(points) if bet_type == "spread" else f"{points:g}"


def _contestant_with_rotation(fixture: dict, rotation: int) -> dict | None:
    for contestant in fixture.get("contestants") or []:
        if contestant.get("rotNum") == rotation:
            return contestant
    return None


def _games_with_rotation(games: list[dict], rotation: int) -> list[tuple[dict, dict, dict]]:
    """(game, main fixture, contestant) for every game listing a team at `rotation`."""
    found = []
    for game in games:
        for fixture in game.get("fixtures") or []:
            if not fixture.get("isMain"):
                continue
            contestant = _contestant_with_rotation(fixture, rotation)
            if contestant is not None:
                found.append((game, fixture, contestant))
    return found


def _main_market(game: dict, fixture: dict, market_type: int) -> dict | None:
    for market in game.get("markets") or []:
        if (market.get("fixtureId") == fixture.get("id") and market.get("type") == market_type
                and market.get("periodNumber") == FULL_GAME_PERIOD):
            return market
    return None


def _odds_row(market: dict, leg: dict, contestant: dict) -> dict | None:
    for odds in market.get("odds") or []:
        if odds.get("index") != MAIN_LINE_INDEX:
            continue
        if leg["betType"] == "spread" and odds.get("contestantId") == contestant.get("id"):
            return odds
        if leg["betType"] == "total" and odds.get("side") == TOTAL_SIDES[leg["side"]]:
            return odds
    return None


def match_leg(games: list[dict], leg: dict, now: datetime) -> dict | str:
    """BFA's game, market and odds row for one leg, or why it cannot be bet as listed:
    {game, fixture, market, odds}. The game is the one listing the leg's rotation, starting
    within START_WINDOW_SEC of the list's start and not yet started; the odds row is the main
    full-game line, open, at the very number the list was built on (AcceptChanges 0 would
    refuse any other, so the check here only says why first)."""
    found = _games_with_rotation(games, leg["rotation"])
    expected_start = _parse_iso(leg["eventStart"])
    if expected_start is not None:
        found = [item for item in found
                 if (start := _parse_iso(item[1].get("date"))) is not None
                 and abs((start - expected_start).total_seconds()) <= START_WINDOW_SEC]
    if not found:
        return f"{leg['label']}: BFA's board lists no game with rotation {leg['rotation']}"
    if len(found) > 1:
        return f"{leg['label']}: BFA's board lists {len(found)} games with rotation {leg['rotation']}"
    game, fixture, contestant = found[0]
    start = _parse_iso(fixture.get("date"))
    if start is None or start <= now or fixture.get("isLive") or game.get("isLive"):
        return f"{leg['label']}: the game has started"
    market = _main_market(game, fixture, BET_TYPE_MARKETS[leg["betType"]])
    odds = _odds_row(market, leg, contestant) if market else None
    if market is None or odds is None:
        return f"{leg['label']}: BFA lists no full-game {leg['betType']} for {game.get('name')}"
    if market.get("status") != MARKET_OPEN or odds.get("status") != ODDS_OPEN:
        return f"{leg['label']}: BFA has the {leg['betType']} closed"
    if not _is_number(odds.get("line")) or odds["line"] != leg["points"]:
        return (f"{leg['label']}: BFA has {_number_text(leg['betType'], odds.get('line') or 0)} now, "
                f"the list has {_number_text(leg['betType'], leg['points'])}")
    if not _is_number(odds.get("price")) or not odds.get("dSpId") or odds.get("dGtId") is None:
        return f"{leg['label']}: BFA's odds row carries no price, sport or game type"
    return {"game": game, "fixture": fixture, "market": market, "odds": odds}


def pick_of(match: dict, stake: int) -> dict:
    """One Picks[] entry, field for field as the site sends it."""
    odds = match["odds"]
    return {
        "EventId": match["game"]["id"],
        "FixtureId": match["fixture"]["id"],
        "MarketType": match["market"]["type"],
        "PeriodNumber": FULL_GAME_PERIOD,
        "Side": odds["side"],
        "Index": MAIN_LINE_INDEX,
        "ContestantId": odds.get("contestantId"),
        "Line": odds["line"],
        "TeaserPoints": float(TEASER_POINTS),
        "Price": odds["price"],
        "Amount": float(stake),
        "PointsPurchased": 0.0,
        "PitcherAction": 0,
        "RoundRobinCombinations": 0,
        "UseFreePlay": False,
        "RiskOrWinType": 1,
        "IdGameType": odds["dGtId"],
        "IdSport": odds["dSpId"],
    }


def wager_body(transaction_id: str, player_id: str, wager_type_id: int, picks: list[dict]) -> list[dict]:
    """The POST body: one teaser wager risking each pick's Amount, no line changes accepted."""
    return [{
        "IdTransaction": transaction_id,
        "IdPlayer": int(player_id),
        "FillIdWager": -1,
        "WagerType": WAGER_TYPE_TEASER,
        "OpenSpots": 0,
        "Picks": picks,
        "RiskWin": 1,
        "AcceptChanges": 0,
        "IdWagerType": wager_type_id,
        "IsLive": False,
        "PhoneLine": None,
        "UserId": None,
    }]


# ---- the open bets ---------------------------------------------------------------------

def _is_teaser(wager: dict) -> bool:
    return bool(TEASER_HEADER_RE.search(str(wager.get("headerDescription") or "")))


def _detail_rotation(detail: dict) -> int | None:
    found = DETAIL_ROTATION_RE.search(str(detail.get("detailDescription") or ""))
    return int(found.group(1)) if found else None


def open_teaser_with_sides(open_wagers: list[dict], sides: set[tuple[int, int]]) -> dict | None:
    """The open teaser whose legs are exactly these (BFA game id, rotation) pairs, or None:
    the same ticket already placed."""
    for wager in open_wagers:
        if not _is_teaser(wager):
            continue
        held = {(detail.get("idGame"), _detail_rotation(detail)) for detail in wager.get("betDetails") or []}
        if held == sides:
            return wager
    return None


def new_teaser_on_games(open_wagers: list[dict], game_ids: set[int], known_ids: set) -> dict | None:
    """A teaser not open before the POST whose legs sit on exactly these BFA games, or None."""
    for wager in open_wagers:
        if wager.get("idWager") in known_ids or not _is_teaser(wager):
            continue
        if {detail.get("idGame") for detail in wager.get("betDetails") or []} == game_ids:
            return wager
    return None


def ticket_of(wager: dict) -> dict:
    return {"ticketNumber": wager.get("ticketNumber") or wager.get("idWager"), "risk": wager.get("riskAmount"),
            "toWin": wager.get("winAmount"), "placedDate": wager.get("placedDate")}


def _short(text: str) -> str:
    return text if len(text) <= MAX_ERROR_TEXT else f"{text[:MAX_ERROR_TEXT]}…"


# ---- the placer ------------------------------------------------------------------------

class BFATeaserPlacer:
    """Places one teaser at a time on the BFA source's own Keycloak session.

    `session_factory()` builds an HTTP session with .get()/.post() (tests inject a fake);
    `clock()` and `sleep()` pace the confirmation reads; `new_transaction_id()` mints the
    wager's IdTransaction. Held in memory only: the metadata (an hour) and the tickets that
    went out (RECENT_TICKET_SEC), so a double click or a stale list never bets twice.
    """

    def __init__(self, bfa_source: BFASource, session_factory: Callable[[], object] = new_session,
                 clock: Callable[[], float] = time.time, sleep: Callable[[float], None] = time.sleep,
                 new_transaction_id: Callable[[], str] = lambda: str(uuid.uuid4())):
        self._bfa = bfa_source
        self._session_factory = session_factory
        self._clock = clock
        self._sleep = sleep
        self._new_transaction_id = new_transaction_id
        self._placing = threading.Lock()
        self._metadata: dict | None = None
        self._metadata_at = 0.0
        self._sent: dict[tuple, dict] = {}

    def place(self, request: dict) -> dict:
        """One validated request (validate_place_request) -> {status, message, ticket?, openWagers?}."""
        if not self._placing.acquire(blocking=False):
            return _refused("another ticket is being placed at BFA right now")
        try:
            return self._place(request)
        finally:
            self._placing.release()

    def _place(self, request: dict) -> dict:
        key = ticket_key(request["legs"])
        recent = self._recent(key)
        if recent:
            return _refused(recent)
        headers, player_id = self._bfa.authorized()
        metadata = self._account_metadata(headers, player_id)
        teaser_type = teaser_type_of(metadata, {SPORT_CODES[leg["league"]] for leg in request["legs"]})
        if isinstance(teaser_type, str):
            return _refused(teaser_type)
        if teaser_type["maxRisk"] is not None and request["stake"] > teaser_type["maxRisk"]:
            return _refused(f"${request['stake']} is over BFA's ${teaser_type['maxRisk']:g} limit on this teaser")
        matches = self._match_legs(request["legs"], headers, player_id, metadata.get("agentId"))
        if isinstance(matches, str):
            return _refused(matches)
        open_before = self._bfa.open_wagers()
        sides = {(match["odds"].get("dGmId"), leg["rotation"]) for leg, match in zip(request["legs"], matches)}
        already = open_teaser_with_sides(open_before, sides)
        if already:
            return _refused(f"this ticket is already open at BFA (ticket {ticket_of(already)['ticketNumber']})")
        transaction_id = self._new_transaction_id()
        body = wager_body(transaction_id, player_id, teaser_type["wagerTypeId"],
                          [pick_of(match, request["stake"]) for match in matches])
        log.info("bfa teaser: posting %s, $%d, transaction %s", " / ".join(leg["label"] for leg in request["legs"]),
                 request["stake"], transaction_id)
        # Held before the POST: from here on the legs may be bet, so they never go out twice.
        self._sent[key] = {"at": self._clock(), "ticket": None}
        try:
            posted = self._post_wager(body, headers, player_id)
        except Exception as error:  # noqa: BLE001 — a lost reply is not a refusal: the wager may stand
            log.warning("bfa teaser: transaction %s POST failed: %s: %s", transaction_id, type(error).__name__, error)
            posted = {"refused": None, "note": f"the request to BFA failed ({type(error).__name__})"}
        if posted["refused"]:
            del self._sent[key]
            log.warning("bfa teaser: transaction %s refused: %s", transaction_id, posted["refused"])
            return _refused(posted["refused"])
        game_ids = {match["odds"].get("dGmId") for match in matches}
        return self._confirm(key, game_ids, {wager.get("idWager") for wager in open_before}, transaction_id, posted["note"])

    def _recent(self, key: tuple) -> str | None:
        """Why the same legs cannot go out again yet, or None."""
        now = self._clock()
        self._sent = {held: sent for held, sent in self._sent.items() if now - sent["at"] < RECENT_TICKET_SEC}
        sent = self._sent.get(key)
        if sent is None:
            return None
        if sent["ticket"]:
            return f"this ticket was placed {int(now - sent['at'])} s ago (ticket {sent['ticket']})"
        return (f"this ticket went to BFA {int(now - sent['at'])} s ago and no ticket has shown yet: "
                "check BFA's open bets before placing it again")

    def _account_metadata(self, headers: dict, player_id: str) -> dict:
        if self._metadata is not None and self._clock() - self._metadata_at < METADATA_TTL_SEC:
            return self._metadata
        response = self._session_factory().get(PLAYER_METADATA_URL, params={"playerId": player_id},
                                               headers=headers, timeout=HTTP_TIMEOUT_SEC)
        if response.status_code != 200:
            raise RuntimeError(f"BFA account metadata: HTTP {response.status_code}")
        body = response.json()
        if not isinstance(body, dict):
            raise RuntimeError(f"BFA account metadata: expected an object, got {type(body).__name__}")
        self._metadata, self._metadata_at = body, self._clock()
        return body

    def _board(self, league: str, headers: dict, player_id: str, agent_id: object) -> list[dict]:
        response = self._session_factory().get(
            BOARD_URL.format(slug=BOARD_SLUGS[league]), headers=headers, timeout=HTTP_TIMEOUT_SEC,
            params={"playerId": player_id, "agentId": agent_id or 0, "fixtureType": BOARD_FIXTURE_TYPE, "set": "Auto"})
        if response.status_code != 200:
            raise RuntimeError(f"BFA {league.upper()} board: HTTP {response.status_code}")
        games = (response.json() or {}).get("games")
        if not isinstance(games, list):
            raise RuntimeError(f"BFA {league.upper()} board carries no games list")
        return games

    def _match_legs(self, legs: list[dict], headers: dict, player_id: str, agent_id: object) -> list[dict] | str:
        boards = {league: self._board(league, headers, player_id, agent_id)
                  for league in sorted({leg["league"] for leg in legs})}
        now = datetime.fromtimestamp(self._clock(), timezone.utc)
        matches = []
        for leg in legs:
            match = match_leg(boards[leg["league"]], leg, now)
            if isinstance(match, str):
                return match
            matches.append(match)
        if len({match["fixture"]["id"] for match in matches}) != len(matches):
            return "two legs are on the same game: a teaser takes one leg per game"
        return matches

    def _post_wager(self, body: list[dict], headers: dict, player_id: str) -> dict:
        """{refused, note}: `refused` is why BFA turned the wager down — only a 4xx, which
        books nothing — else None and the open bets decide. `note` names a reply other than
        the usual {"state": 0}, for the message when no ticket shows. A network error is
        raised: the wager may have gone through."""
        response = self._session_factory().post(WAGER_URL, params={"playerId": player_id}, json=body,
                                                headers={**headers, "Content-Type": "application/json; charset=utf-8"},
                                                timeout=HTTP_TIMEOUT_SEC)
        reply = None
        try:
            reply = response.json()
        except ValueError:
            pass
        shown = _short(str(reply if reply is not None else response.text))
        if HTTP_CLIENT_ERRORS[0] <= response.status_code <= HTTP_CLIENT_ERRORS[1]:
            return {"refused": f"BFA refused the wager: HTTP {response.status_code} {shown}", "note": None}
        if response.status_code == 200 and isinstance(reply, dict) and reply.get("state") == 0:
            return {"refused": None, "note": None}
        return {"refused": None, "note": f"BFA answered HTTP {response.status_code} {shown}"}

    def _confirm(self, key: tuple, game_ids: set, known_ids: set, transaction_id: str, note: str | None) -> dict:
        deadline = self._clock() + CONFIRM_TIMEOUT_SEC
        open_wagers: list[dict] = []
        while True:
            self._sleep(CONFIRM_POLL_SEC)
            try:
                open_wagers = self._bfa.open_wagers()
            except Exception as error:  # noqa: BLE001 — a failed read is one more try, not a failed wager
                log.warning("bfa teaser: open bets read failed while confirming %s: %s", transaction_id, error)
                open_wagers = []
            placed = new_teaser_on_games(open_wagers, game_ids, known_ids)
            if placed:
                ticket = ticket_of(placed)
                self._sent[key]["ticket"] = ticket["ticketNumber"]
                log.info("bfa teaser: transaction %s placed as ticket %s", transaction_id, ticket["ticketNumber"])
                return {"status": STATUS_PLACED, "message": f"Placed · ticket {ticket['ticketNumber']}",
                        "ticket": ticket, "openWagers": open_wagers}
            if self._clock() >= deadline:
                log.warning("bfa teaser: transaction %s sent, no ticket after %.0f s (%s)", transaction_id,
                            CONFIRM_TIMEOUT_SEC, note or "BFA took it for processing")
                reply = f" ({note})" if note else ""
                return {"status": STATUS_UNCONFIRMED,
                        "message": (f"Sent to BFA{reply}, but no ticket showed in {CONFIRM_TIMEOUT_SEC:.0f} s. "
                                    "Check BFA's open bets before placing it again."),
                        "openWagers": open_wagers}


def _refused(message: str) -> dict:
    return {"status": STATUS_REFUSED, "message": f"Not placed: {message}"}
