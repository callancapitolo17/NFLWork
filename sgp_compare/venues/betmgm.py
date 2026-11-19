"""BetMGM venue adapter: NFL games, every bet-builder leg (game lines, halves,
quarters, team totals, player props), and same-game parlay prices.

Network: BetMGM's Entain "CDS" JSON endpoints over the curl_cffi Chrome-TLS
session and accessid harvest of mlb_sgp.betmgm_client.BetMGMClient (no login):
    GET  /cds-api/bettingoffer/fixtures?...&sportIds=11       NFL schedule
    GET  /cds-api/bettingoffer/fixtures?...&fixtureIds={id}&offerMapping=All
    POST /cds-api/bettingoffer/picks   {tv1Picks: [], tv2Picks: [...]}
The picks POST is the bet-builder price request: it returns a price and never
places a bet (betting needs a logged-in betslip, which this module has no code
for). No DB I/O.

BetMGM shapes (recon 2026-11-18 on Cal's Mac, IND @ HOU, tests/fixtures/betmgm_*):
  * NFL fixtures: sport 11, competition name "NFL", fixture name "Away @ Home";
    participants carry properties.type HomeTeam / AwayTeam / Player (a player's
    properties.team is its team's participantId).
  * ONE optionMarket PER LINE: "Spread" -12.5 and "Spread" -10 are separate
    markets; isBalancedLine marks the main one (isMain only on Moneyline).
  * every market's `parameters` say what it is: MarketType (2way, 2wayHandicap,
    Over/Under, Scorer, ...), MarketSubType (YesNo milestone, PlayerTD),
    Happening (Point, ReceivingYards, Touchdown, ...), Period (FullTime,
    FirstHalf, FirstQuarter), DecimalValue / DecimalHandicap (the number) and
    FixtureParticipant (the team or player).
  * milestones ("Tank Dell to record 50+ receiving yards") and player TDs
    ("Deion Burks to score 1+ TDs") have a single "Yes" option: no other side.
  * option price: {odds (decimal), americanOdds, numerator, denominator}.
  * NFL parlays price through tv2Picks (the MLB tv1Picks body returns nothing).
    The answer's betBuilderPricingGroups[group] has odds + suspensionState; an
    impossible combo comes back with an `error` (ZeroProbabilityError) and a
    one-leg request is refused (SingleSelectionRequestError).
"""
from __future__ import annotations

import re
import threading
import time
import uuid
from datetime import datetime, timezone
from typing import Any, Callable

from sgp_compare.legs import CanonicalLeg, Game, OfferedLeg, ParlayQuote, ParlayTier
from sgp_compare.odds import american_to_decimal
from sgp_compare.venues.base import VenueTransportError

MGM_FOOTBALL_SPORT_ID = 11
NFL_COMPETITION_NAME = "NFL"
PRE_MATCH_STAGE = "PreMatch"
PRICED_STATE = "MarketUnsuspended"
VISIBLE = "Visible"
REQUEST_TIMEOUT_SEC = 20
GAMES_CACHE_SEC = 300          # the schedule changes daily
MARKETS_CACHE_SEC = 20         # the leg picker; pricing always re-reads
SCHEDULE_PAGE_SIZE = 200       # football fixtures per schedule page (college fills most of them)
SCHEDULE_MAX_PAGES = 5         # 1,000 fixtures: a full college Saturday is well under that
# The only market types whose two options are each other's other side.
TWO_WAY_MARKET_TYPES = ("2way", "2wayHandicap", "Over/Under")
# Error codes in a priced group that mean "BetMGM won't price this combination"
# (recon 2026-11-18). Any other error is reported as an error, not a decline.
DECLINE_ERROR_CODES = ("ZeroProbabilityError", "SingleSelectionRequestError")
NOT_REDUNDANT = "NotRedundantLeg"
# BetMGM's 400 body when the harvested accessid went stale (the same test as
# mlb_sgp.betmgm_client._is_accessid_gate, kept here so a missing curl_cffi
# stays a VenueTransportError instead of an ImportError).
ACCESSID_GATE_TEXT = "access id"

PERIOD_OF = {
    "FullTime": "FG", "FirstHalf": "1H", "SecondHalf": "2H",
    "FirstQuarter": "1Q", "SecondQuarter": "2Q", "ThirdQuarter": "3Q", "FourthQuarter": "4Q",
}

# BetMGM Happening -> CanonicalLeg.stat for player Over/Under lines and milestones.
# Names follow the ProphetX adapter where both venues offer the stat.
PLAYER_STAT_OF = {
    "PassingYards": "passing_yards",
    "PassAttempt": "passing_attempts",
    "PassCompleted": "pass_completions",
    "TouchdownPass": "passing_touchdowns",
    "PassingRushingYards": "passing_rushing_yards",
    "RushingYards": "rushing_yards",
    "RushingAttempt": "rushing_attempts",
    "RushingReceivingYards": "rushing_receiving_yards",
    "ReceivingYards": "receiving_yards",
    "Reception": "receptions",
    "LongestReception": "longest_reception",
    "KickingPoint": "kicking_points",
    "FieldGoal": "field_goals_made",
    "ExtraPointMade": "extra_points_made",
    "Touchdown": "touchdowns",
}

# A milestone says "N+" ("to record 50+ receiving yards"): over N - 0.5.
MILESTONE_RE = re.compile(r"\b\d+\+")
# "Houston Texans (-12.5)" / "Indianapolis Colts (+3)"
SPREAD_OPTION_RE = re.compile(r"^(.+?)\s*\(([+-]?\d+(?:\.\d+)?)\)\s*$")


# ---------------------------------------------------------------- parsing (pure)

def parse_nfl_games(fixtures_payload: dict, now: datetime) -> list[Game]:
    """Not-yet-started NFL games from the schedule response, soonest first."""
    games = []
    for fixture in fixtures_payload.get("fixtures") or []:
        if _text((fixture.get("competition") or {}).get("name")) != NFL_COMPETITION_NAME:
            continue
        if fixture.get("stage") not in (None, PRE_MATCH_STAGE):
            continue
        game = _game_from_fixture(fixture)
        if game is not None and game.game_start_time > now:
            games.append(game)
    return sorted(games, key=lambda g: g.game_start_time)


def _game_from_fixture(fixture: dict) -> Game | None:
    teams = {(p.get("properties") or {}).get("type"): _text(p.get("name"))
             for p in fixture.get("participants") or []}
    home, away = teams.get("HomeTeam"), teams.get("AwayTeam")
    if not (home and away):
        return None
    start = _parse_time(fixture.get("startDate"))
    if start is None:
        return None
    return Game(venue_game_id=str(fixture["id"]), home_team=home.strip(), away_team=away.strip(),
                game_start_time=start)


def parse_legs(game: Game, markets_payload: dict) -> list[OfferedLeg]:
    """Every option of the game's visible bet-builder markets.

    Two-option markets give two legs that point at each other; one-option
    ("Yes" milestones, player TDs) and many-option markets (first TD scorer,
    correct score) give legs with no other side. A suspended option is listed
    with no price; a suspended market is left out. A fixture closed for
    betting, or one whose tree BetMGM has not built yet (it builds the full
    tree the morning before), has no legs.
    """
    fixtures = markets_payload.get("fixtures") or []
    if not fixtures:
        return []
    fixture = fixtures[0]
    if str(fixture["id"]) != game.venue_game_id:
        raise VenueTransportError(f"asked for fixture {game.venue_game_id}, BetMGM answered {fixture['id']}")
    if fixture.get("isOpenForBetting") is False or fixture.get("stage") not in (None, PRE_MATCH_STAGE):
        return []
    participants = {p.get("id"): p for p in fixture.get("participants") or []}
    markets = [m for m in fixture.get("optionMarkets") or []
               if m.get("isBetBuilder") and m.get("status") == VISIBLE and m.get("options")]
    legs: list[OfferedLeg] = []
    for market in markets:
        params = _params(market)
        options = market["options"]
        ids = [_leg_id(market, o) for o in options]
        pair = len(options) == 2 and params.get("MarketType") in TWO_WAY_MARKET_TYPES
        for index, option in enumerate(options):
            other = options[1 - index] if pair else None
            legs.append(OfferedLeg(
                leg_id=ids[index],
                market_label=_text(market.get("name")),
                category=_category(params, participants),
                side_label=_text(option.get("name")),
                canonical=_canonical(game, market, params, option, participants),
                american_odds=_american(option),
                opposite_american_odds=_american(other) if other else None,
                available_stake=None,
                opposite_leg_id=ids[1 - index] if pair else None,
                market_key=_market_key(params, market),
                is_main_line=bool(market.get("isBalancedLine") or market.get("isMain")),
                venue_ref={"fixtureId": game.venue_game_id, "optionMarketId": market.get("id"),
                           "optionId": option.get("id")},
            ))
    return legs


def _leg_id(market: dict, option: dict) -> str:
    return f"{market.get('id')}:{option.get('id')}"


def _market_key(params: dict, market: dict) -> str:
    """One key per BetMGM market FAMILY, because BetMGM posts each line as its
    own market: every "Spread" line of the game shares a key (so the page
    folds alt lines under the main one and a parlay takes one line), as do a
    player's 50+/75+/100+ yard milestones."""
    parts = [params.get(k, "") for k in ("MarketType", "MarketSubType", "Happening", "Period", "FixtureParticipant")]
    if not any(parts):
        return f"market:{market.get('id')}"
    return "|".join(parts)


def _category(params: dict, participants: dict) -> str:
    participant = participants.get(_int(params.get("FixtureParticipant")))
    kind = ((participant or {}).get("properties") or {}).get("type")
    if kind == "Player":
        return "Player Props"
    if kind in ("HomeTeam", "AwayTeam"):
        return "Team Totals"
    return "Game Lines"


def _canonical(game: Game, market: dict, params: dict, option: dict, participants: dict) -> CanonicalLeg | None:
    period = PERIOD_OF.get(params.get("Period", ""))
    if period is None:
        return None
    market_type = params.get("MarketType")
    happening = params.get("Happening")
    option_name = _text(option.get("name"))
    option_types = (option.get("parameters") or {}).get("optionTypes") or []
    participant = participants.get(_int(params.get("FixtureParticipant")))
    participant_kind = ((participant or {}).get("properties") or {}).get("type")
    subject = _text((participant or {}).get("name")) or None

    if happening == "Point" and market_type == "2way":
        side = _team_side(game, option_name)
        return CanonicalLeg(period, "moneyline", None, side, None) if side else None
    if happening == "Point" and market_type == "2wayHandicap":
        match = SPREAD_OPTION_RE.match(option_name)
        side = _team_side(game, match.group(1)) if match else None
        return CanonicalLeg(period, "spread", None, side, float(match.group(2))) if side else None
    if market_type == "Scorer" and happening == "Touchdown" and params.get("Position") == "First":
        # Each option is a player; an option with no player ("No Touchdown") stays unmatched.
        if (option.get("parameters") or {}).get("fixtureParticipant") is None:
            return None
        return CanonicalLeg(period, "first_td", option_name, "yes", None)
    if market_type != "Over/Under":
        return None

    value = _float(params.get("DecimalValue"))
    if value is None:
        return None
    sub_type = params.get("MarketSubType")
    side = "over" if "Over" in option_types else "under" if "Under" in option_types else None
    if side is None:
        return None
    if happening == "Point" and sub_type is None and "FixtureParticipant" not in params:
        return CanonicalLeg(period, "total", None, side, value)
    if happening == "Point" and sub_type is None and participant_kind in ("HomeTeam", "AwayTeam"):
        return CanonicalLeg(period, "team_total", subject, side, value)
    if participant_kind != "Player" or subject is None:
        return None
    if sub_type == "PlayerTD":
        # DecimalValue is already the half-point line: 1+ TDs = 0.5.
        if abs(value - 0.5) < 1e-9:
            return CanonicalLeg(period, "anytime_td", subject, "yes", None)
        return CanonicalLeg(period, "touchdowns", subject, "over", value)
    stat = PLAYER_STAT_OF.get(happening or "")
    if stat is None:
        return None
    if sub_type == "YesNo":
        # "N+" milestone = over N - 0.5. "Over 75 Rushing and Reception Yards"
        # does not say whether exactly 75 wins, so it stays unmatched.
        if not MILESTONE_RE.search(_text(market.get("name"))):
            return None
        return CanonicalLeg(period, stat, subject, "over", value - 0.5)
    return CanonicalLeg(period, stat, subject, side, value)


def _team_side(game: Game, team_name: str) -> str | None:
    name = team_name.strip()
    if name == game.home_team:
        return "home"
    if name == game.away_team:
        return "away"
    return None


def _american(option: dict | None) -> float | None:
    """The option's American price. BetMGM's decimal / fractional fields can
    differ from it by up to ~1 cent (-112 next to 1.88): they are rounded to a
    fraction ladder, while US BetMGM shows and pays the American price."""
    if not option or option.get("status") != VISIBLE:
        return None
    american = (option.get("price") or {}).get("americanOdds")
    if isinstance(american, (int, float)) and abs(american) >= 100:
        return float(american)
    return None


def parlay_payload(legs: list[OfferedLeg], pick_group_id: str) -> dict:
    return {"tv1Picks": [], "tv2Picks": [
        {"fixtureId": leg.venue_ref["fixtureId"], "optionMarketId": leg.venue_ref["optionMarketId"],
         "optionId": leg.venue_ref["optionId"], "useLiveFallback": False, "pickGroupId": pick_group_id}
        for leg in legs]}


def parse_parlay_quote(response: dict, legs: list[OfferedLeg], pick_group_id: str) -> ParlayQuote | None:
    """BetMGM's one price for `legs` (no stake limit), or None when it declined:
    a known decline code (impossible combo, one selection), the group or any
    leg suspended (a suspended group still carries odds that are not on
    offer), or a leg BetMGM marks redundant (it would price fewer legs than
    asked). An unknown error, or an answer that does not price exactly the
    legs sent, raises VenueTransportError rather than reading as a decline."""
    group = (response.get("betBuilderPricingGroups") or {}).get(pick_group_id)
    if group is None:
        raise VenueTransportError(f"BetMGM picks: no pricing group {pick_group_id} in the answer")
    if group.get("suspensionState") != PRICED_STATE:
        return None                  # suspended: a decline, whatever error it carries
    error = group.get("error")
    if error:
        if any(code in error for code in DECLINE_ERROR_CODES):
            return None
        raise VenueTransportError(f"BetMGM picks: {error[:200]}")
    leg_states = group["legInformation"]
    sent = sorted(_client_leg_id(leg) for leg in legs)
    if sorted(info["clientLegId"] for info in leg_states) != sent:
        raise VenueTransportError(f"BetMGM picks: priced legs {[i['clientLegId'] for i in leg_states]}, sent {sent}")
    if any(info.get("suspensionState") != PRICED_STATE or info.get("redundantState") != NOT_REDUNDANT
           for info in leg_states):
        return None
    odds = group.get("odds") or {}
    american = odds.get("americanOdds")
    if isinstance(american, (int, float)) and abs(american) >= 100:
        decimal = american_to_decimal(float(american))
    elif isinstance(odds.get("odds"), (int, float)) and odds["odds"] > 1.0:
        decimal = float(odds["odds"])
    else:
        return None
    stored = {k: v for k, v in response.items() if k != "fixturePage"}   # fixturePage repeats the whole tree
    return ParlayQuote((ParlayTier(decimal, None),), {}, stored)


def _client_leg_id(leg: OfferedLeg) -> str:
    """BetMGM echoes each pick back as "<fixtureId>:<optionMarketId>:<optionId>"."""
    ref = leg.venue_ref
    return f"{ref['fixtureId']}:{ref['optionMarketId']}:{ref['optionId']}"


def _params(market: dict) -> dict[str, str]:
    return {p.get("key"): str(p.get("value")) for p in market.get("parameters") or [] if p.get("key")}


def _text(value: Any) -> str:
    if isinstance(value, dict):
        value = value.get("value")
    return value if isinstance(value, str) else ""


def _int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _parse_time(value: Any) -> datetime | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return (parsed if parsed.tzinfo else parsed.replace(tzinfo=timezone.utc)).astimezone(timezone.utc)


# ---------------------------------------------------------------- the venue

class BetMgmVenue:
    """Live BetMGM reads. Thread-safe; caches only the schedule and the leg
    picker's market tree, never a price used to compute a parlay."""

    name = "betmgm"
    display_name = "BetMGM"

    def __init__(self, client_factory: Callable[[], Any] | None = None,
                 clock: Callable[[], float] = time.monotonic,
                 wall_clock: Callable[[], datetime] = lambda: datetime.now(timezone.utc)) -> None:
        self._client_factory = client_factory or _betmgm_client
        self._client = None
        self._clock = clock
        self._wall_clock = wall_clock
        self._lock = threading.Lock()
        self._games: tuple[float, list[Game]] | None = None
        self._markets: dict[str, tuple[float, dict]] = {}

    def list_games(self) -> list[Game]:
        """The cached schedule, re-checked against the clock on every call so a
        game drops out at kickoff, not when the cache next refills."""
        now = self._wall_clock()
        return [game for game in self._read_games() if game.game_start_time > now]

    def list_legs(self, game: Game, fresh: bool = False) -> list[OfferedLeg]:
        payload = self._read_markets(game.venue_game_id, fresh)
        return _parsed("markets", lambda: parse_legs(game, payload))

    def price_parlay(self, game: Game, legs: list[OfferedLeg]) -> ParlayQuote | None:
        pick_group_id = str(uuid.uuid4())
        response = self._request("POST", "/bettingoffer/picks", json_body=parlay_payload(legs, pick_group_id))
        return _parsed("picks", lambda: parse_parlay_quote(response, legs, pick_group_id))

    # -- reads
    def _read_games(self) -> list[Game]:
        with self._lock:
            cached = self._games
        if cached and self._clock() - cached[0] < GAMES_CACHE_SEC:
            return cached[1]
        payload = self._read_schedule()
        games = _parsed("fixtures", lambda: parse_nfl_games(payload, self._wall_clock()))
        with self._lock:
            self._games = (self._clock(), games)
        return games

    def _read_schedule(self) -> dict:
        """Every football fixture, page by page: the query cannot filter to
        the NFL, and a college Saturday alone can fill a page."""
        fixtures: dict[str, dict] = {}     # by id: a live list can shift a fixture across pages
        for page in range(SCHEDULE_MAX_PAGES):
            payload = self._request("GET", "/bettingoffer/fixtures", query=(
                "&fixtureTypes=Standard&state=Latest&offerMapping=Filtered"
                "&offerCategories=Gridable&fixtureCategories=Gridable"
                f"&sportIds={MGM_FOOTBALL_SPORT_ID}&skip={page * SCHEDULE_PAGE_SIZE}"
                f"&take={SCHEDULE_PAGE_SIZE}&sortBy=Tags"))
            page_fixtures = payload.get("fixtures") or []
            if not page_fixtures:            # stop on an empty page, in case BetMGM caps `take`
                return {"fixtures": list(fixtures.values())}
            for fixture in page_fixtures:
                fixtures.setdefault(str(fixture.get("id")), fixture)
        raise VenueTransportError(f"BetMGM schedule: more than {SCHEDULE_MAX_PAGES * SCHEDULE_PAGE_SIZE} "
                                  "football fixtures; NFL games could be missing")

    def _read_markets(self, fixture_id: str, fresh: bool) -> dict:
        with self._lock:
            cached = self._markets.get(fixture_id)
        if cached and not fresh and self._clock() - cached[0] < MARKETS_CACHE_SEC:
            return cached[1]
        payload = self._request("GET", "/bettingoffer/fixtures",
                                query=f"&fixtureIds={fixture_id}&offerMapping=All&state=Latest")
        with self._lock:
            self._markets[fixture_id] = (self._clock(), payload)
        return payload

    def _request(self, method: str, path: str, query: str = "", json_body: dict | None = None) -> dict:
        """One CDS call. A 400 "Access id ..." means the harvested id went
        stale: harvest a new one and send the call once more."""
        response = self._send(method, path, query, json_body, refresh_accessid=False)
        if response.status_code == 400 and ACCESSID_GATE_TEXT in (response.text or "").lower():
            response = self._send(method, path, query, json_body, refresh_accessid=True)
        if response.status_code != 200:
            raise VenueTransportError(f"BetMGM {method} {path}: expected HTTP 200, got "
                                      f"{response.status_code}: {(response.text or '')[:200]}")
        try:
            body = response.json()
        except ValueError as exc:
            raise VenueTransportError(f"BetMGM {method} {path}: response is not JSON") from exc
        if not isinstance(body, dict):
            raise VenueTransportError(f"BetMGM {method} {path}: expected a JSON object, got {type(body).__name__}")
        return body

    def _send(self, method: str, path: str, query: str, json_body: dict | None, refresh_accessid: bool):
        try:
            if self._client is None:
                self._client = self._client_factory()
            client = self._client
            client.accessid(refresh=refresh_accessid)
            url = f"{client._base}{path}?{client._common_q()}{query}"
            if method == "GET":
                return client.session.get(url, timeout=REQUEST_TIMEOUT_SEC)
            return client.session.post(url, json=json_body, timeout=REQUEST_TIMEOUT_SEC)
        except Exception as exc:   # curl_cffi, the client import and the accessid harvest raise their own types
            self._client = None
            raise VenueTransportError(f"BetMGM {method} {path}: {type(exc).__name__}: {exc}") from exc


def _parsed(what: str, parse: Callable[[], Any]) -> Any:
    """Run a parser over a BetMGM answer; a shape it does not expect is the
    venue answering garbage, reported as such rather than as a crash."""
    try:
        return parse()
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        raise VenueTransportError(f"BetMGM {what} answer has an unexpected shape: "
                                  f"{type(exc).__name__}: {exc}") from exc


def _betmgm_client():
    """mlb_sgp's BetMGM session: Chrome-TLS curl_cffi plus the per-state
    accessid harvest every CDS read needs (state "pa" by default)."""
    from mlb_sgp.betmgm_client import BetMGMClient
    return BetMGMClient()
