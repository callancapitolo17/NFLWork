"""Novig venue adapter: NFL games, every parlay-eligible leg (game lines and
player props), and parlay prices.

Network: anonymous endpoints on api.novig.us over a curl_cffi Chrome-TLS
session (no login; the 2026-11-19 recon on Cal's Mac read all three):
    GET  /nbx/v1/trading/NFL/page              the NFL screen's game cards
    POST /v1/graphql                           EventMarkets_Query, one game's market tree
    POST /rfq/v2/quote/unauthenticated         {outcomeIds: [...]} -> one parlay price
The quote POST is a price request: it returns a price that expires in 15-30 s
and never places a bet (betting needs a logged-in account, which this module
has no code for). No DB I/O.

The GraphQL endpoint only runs queries the novig.com app ships, so the query
text is the committed `mlb_sgp/novig_event_markets_query.json` (one copy for
the whole repo; `python3 mlb_sgp/refresh_novig_query.py` refreshes it when
Novig answers "query is not allowed").

Novig shapes (recon 2026-11-19, IND at HOU, fixtures in tests/fixtures/):
  * a market is ONE line: {type, strike, is_consensus, status, market_locks,
    player, competitor, outcomes: [two]}. Alt lines are separate markets of the
    same type; `is_consensus` marks the main one.
  * outcome: {id, description, available, competitor}. `available` is the
    price to take that side now as a chance (0.495 = pay 49.5c to win $1), or
    null when nobody offers it. The two sides of a main line sum to ~1.005.
  * SPREAD / MONEY (and their _1H forms): the team is on outcome.competitor.id;
    a spread outcome's description ends with its own signed line ("IND +12.5").
  * TOTAL / TEAM_TOTAL / player props: outcomes are "Over N" / "Under N";
    TEAM_TOTAL's team is market.competitor. FIRST_TOUCHDOWN_SCORER is Yes / No.
  * a player prop has market.player.full_name. TOUCHDOWNS Over 0.5 is the
    anytime touchdown. `<TYPE>_PROTECTED` is a twin market with its own ids and
    a higher price (believed void if the player sits out; not confirmed), kept
    apart from the plain one.
  * quote answer: HTTP 200 {price, expiresAt, serverNow}; `price` is the chance
    the parlay costs, vig included (decimal odds = 1 / price). A refused
    combination is HTTP 400 with the text "Invalid parlay combination.".
    Novig copies its parlay prices from a sportsbook; the old endpoint named
    that book per leg, this one does not, so the source book is unknown.
"""
from __future__ import annotations

import json
import logging
import re
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from sgp_compare.legs import CanonicalLeg, Game, OfferedLeg, ParlayQuote, ParlayTier
from sgp_compare.venues.base import VenueTransportError

log = logging.getLogger(__name__)

NOVIG_NFL_PAGE_URL = "https://api.novig.us/nbx/v1/trading/NFL/page"
NOVIG_GRAPHQL_URL = "https://api.novig.us/v1/graphql"
NOVIG_PARLAY_QUOTE_URL = "https://api.novig.us/rfq/v2/quote/unauthenticated"
NOVIG_WARM_UP_URL = "https://www.novig.us"
EVENT_MARKETS_QUERY_PATH = Path(__file__).resolve().parents[2] / "mlb_sgp" / "novig_event_markets_query.json"

REQUEST_TIMEOUT_SEC = 15
GAMES_CACHE_SEC = 300          # the schedule changes daily
MARKETS_CACHE_SEC = 20         # the leg picker; pricing always re-reads
GAME_CARD_TYPE = "game_event_card"
PREGAME_STATUS = "OPEN_PREGAME"
OPEN_MARKET_STATUS = "OPEN"
DECLINED_PARLAY_TEXT = "Invalid parlay combination"
FIRST_HALF_SUFFIX = "_1H"
PROTECTED_SUFFIX = "_PROTECTED"
GAME_LINES_CATEGORY = "Game Lines"
LINE_TOLERANCE = 1e-6

# Canonical names shared with the other venues (ProphetX's adapter set them).
TOUCHDOWN_STAT_OF_STRIKE = {
    0.5: "anytime_td",
    1.5: "to_score_at_least_2_touchdowns",
    2.5: "to_score_at_least_3_touchdowns",
}
TOUCHDOWN_LABEL_OF_STRIKE = {0.5: "Anytime Touchdown", 1.5: "2+ Touchdowns", 2.5: "3+ Touchdowns"}
FIRST_TOUCHDOWN_STAT = "to_score_first_touchdown"
SIGNED_LINE_AT_END = re.compile(r"([+-]\d+(?:\.\d+)?)\s*$")


# ---------------------------------------------------------------- parsing (pure)

def parse_nfl_games(page_payload: dict, now: datetime) -> list[Game]:
    """Not-yet-started NFL games from the NFL page, soonest first. An
    off-season page ({"sections": []}) gives no games."""
    games = []
    for section in page_payload.get("sections") or []:
        content = section.get("content") or {}
        for card in content.get("components") or []:
            if card.get("type") != GAME_CARD_TYPE or card.get("eventStatus") != PREGAME_STATUS:
                continue
            game = _game_from_card(card)
            if game is not None and game.game_start_time > now:
                games.append(game)
    return sorted(games, key=lambda g: g.game_start_time)


def _game_from_card(card: dict) -> Game | None:
    home, away = card.get("homeTeam") or {}, card.get("awayTeam") or {}
    if not (card.get("eventId") and card.get("scheduledStart") and home.get("name") and away.get("name")):
        return None
    start = datetime.fromisoformat(card["scheduledStart"].replace("Z", "+00:00"))
    if start.tzinfo is None:   # Novig sends UTC; a bare time must not be read as the Mac's local time
        start = start.replace(tzinfo=timezone.utc)
    return Game(venue_game_id=str(card["eventId"]), home_team=home["name"], away_team=away["name"],
                game_start_time=start.astimezone(timezone.utc))


def parse_legs(game: Game, markets_payload: dict) -> list[OfferedLeg]:
    """Every priced leg of the game's open, unlocked markets.

    A line where neither side has a price is left out (deep alt ladders have
    empty tails). An unknown game id answers {"data": {"event": []}}: no legs.
    """
    events = (markets_payload.get("data") or {}).get("event") or []
    if not events:
        return []
    event = events[0]
    teams = (event.get("game") or {})
    home_team_id = (teams.get("homeTeam") or {}).get("id")
    away_team_id = (teams.get("awayTeam") or {}).get("id")
    if not home_team_id or not away_team_id or home_team_id == away_team_id:
        raise ValueError(f"expected two different team ids on the game, got {home_team_id!r} and {away_team_id!r}")
    side_of_team_id = {home_team_id: "home", away_team_id: "away"}
    legs: list[OfferedLeg] = []
    for market in event.get("markets") or []:
        if market.get("status") != OPEN_MARKET_STATUS or market.get("market_locks"):
            continue
        legs.extend(_legs_of_market(game, market, side_of_team_id))
    return legs


def _legs_of_market(game: Game, market: dict, side_of_team_id: dict[str, str]) -> list[OfferedLeg]:
    outcomes = market.get("outcomes") or []
    if len(outcomes) != 2 or not all(o.get("id") for o in outcomes):
        return []
    prices = [american_odds_of_chance(o.get("available")) for o in outcomes]
    if prices[0] is None and prices[1] is None:
        return []
    described = [_describe(game, market, outcome, side_of_team_id) for outcome in outcomes]
    if described[0] is None or described[1] is None:
        return []
    legs = []
    for index, other_index in ((0, 1), (1, 0)):
        outcome, other = outcomes[index], outcomes[other_index]
        market_label, category, side_label, canonical, market_key = described[index]
        legs.append(OfferedLeg(
            leg_id=outcome["id"],
            market_label=market_label,
            category=category,
            side_label=side_label,
            canonical=canonical,
            american_odds=prices[index],
            opposite_american_odds=prices[other_index],
            available_stake=None,          # the market tree carries a price, not the dollars behind it
            opposite_leg_id=other["id"],
            market_key=market_key,
            is_main_line=market.get("is_consensus") is True,
            venue_ref={"outcomeId": outcome["id"], "marketId": market.get("id"), "marketType": market.get("type")},
        ))
    return legs


def american_odds_of_chance(chance: Any) -> float | None:
    """Novig's `available` (0.495) as American odds (+102.02). None for a
    missing price or one outside (0, 1), which has no payout."""
    if isinstance(chance, bool) or not isinstance(chance, (int, float)) or not 0.0 < chance < 1.0:
        return None
    if chance <= 0.5:
        return round((1.0 - chance) / chance * 100.0, 2)
    return round(-chance / (1.0 - chance) * 100.0, 2)


def _describe(game: Game, market: dict, outcome: dict,
              side_of_team_id: dict[str, str]) -> tuple[str, str, str, CanonicalLeg | None, str] | None:
    """(market_label, category, side_label, canonical, market_key) of one
    outcome, or None when the outcome cannot be read safely (the whole line is
    then left out rather than shown with a guessed side or number)."""
    market_type = market.get("type") or ""
    if market.get("player"):
        return _describe_player_prop(market, outcome)
    period, base_type = ("1H", market_type[: -len(FIRST_HALF_SUFFIX)]) if market_type.endswith(FIRST_HALF_SUFFIX) \
        else ("FG", market_type)
    period_words = "1st Half " if period == "1H" else ""
    strike = _number(market.get("strike"))
    if base_type in ("MONEY", "SPREAD"):
        side = side_of_team_id.get((outcome.get("competitor") or {}).get("id"))
        if side is None:
            return None
        team_name = game.home_team if side == "home" else game.away_team
        if base_type == "MONEY":
            return (f"{period_words}Moneyline", GAME_LINES_CATEGORY, team_name,
                    CanonicalLeg(period, "moneyline", None, side, None), market_type)
        line = _spread_line_of(outcome, strike)
        if line is None:
            return None
        return (f"{period_words}Spread", GAME_LINES_CATEGORY, f"{team_name} {line:+g}",
                CanonicalLeg(period, "spread", None, side, line), market_type)
    side = _over_under_side(outcome)
    if side is None or strike is None:
        return None
    side_label = f"{side.capitalize()} {strike:g}"
    if base_type == "TOTAL":
        return (f"{period_words}Total", GAME_LINES_CATEGORY, side_label,
                CanonicalLeg(period, "total", None, side, strike), market_type)
    if base_type == "TEAM_TOTAL":
        team = market.get("competitor") or {}
        if not team.get("name"):
            return None
        return (f"{period_words}{team['name']} Team Total", GAME_LINES_CATEGORY, side_label,
                CanonicalLeg(period, "team_total", team["name"], side, strike), f"{market_type}:{team.get('id')}")
    # An over/under market type this module has no name for: listed, not matched across venues.
    return (market.get("description") or market_type, _title(market_type), side_label, None, market_type)


def _describe_player_prop(market: dict, outcome: dict) -> tuple[str, str, str, CanonicalLeg | None, str] | None:
    market_type = market["type"]
    player = market["player"]
    player_name = player.get("full_name")
    if not player_name:
        return None
    protected = market_type.endswith(PROTECTED_SUFFIX)
    base_type = market_type[: -len(PROTECTED_SUFFIX)] if protected else market_type
    protected_words = " (Protected)" if protected else ""
    stat_suffix = "_protected" if protected else ""
    strike = _number(market.get("strike"))
    market_key = f"{market_type}:{player.get('id')}"
    description = (outcome.get("description") or "").strip().lower()

    if base_type == "FIRST_TOUCHDOWN_SCORER":
        if description not in ("yes", "no"):
            return None
        return (f"{player_name} First Touchdown{protected_words}", f"First Touchdown{protected_words}",
                description.upper(), CanonicalLeg("FG", FIRST_TOUCHDOWN_STAT + stat_suffix, player_name, description, None),
                market_key)

    side = _over_under_side(outcome)
    if side is None or strike is None:
        return None
    if base_type == "TOUCHDOWNS":
        # Each touchdown count is its own yes/no bet (anytime, 2+, 3+), as the other venues list it.
        stat = TOUCHDOWN_STAT_OF_STRIKE.get(strike)
        if stat is None:
            return None
        yes_no = "yes" if side == "over" else "no"
        label = TOUCHDOWN_LABEL_OF_STRIKE[strike]
        return (f"{player_name} {label}{protected_words}", f"{label}{protected_words}", yes_no.upper(),
                CanonicalLeg("FG", stat + stat_suffix, player_name, yes_no, None), f"{market_key}:{strike:g}")
    return (f"{player_name} {_title(base_type)}{protected_words}", f"{_title(base_type)}{protected_words}",
            f"{side.capitalize()} {strike:g}",
            CanonicalLeg("FG", base_type.lower() + stat_suffix, player_name, side, strike), market_key)


def _spread_line_of(outcome: dict, strike: float | None) -> float | None:
    """The outcome's own signed line, read off its description ("IND +12.5").
    It must agree in size with the market's strike; a line that cannot be read
    or disagrees is refused, because a wrong sign is a different bet."""
    match = SIGNED_LINE_AT_END.search(outcome.get("description") or "")
    if match is None or strike is None:
        return None
    line = float(match.group(1))
    if abs(abs(line) - abs(strike)) > LINE_TOLERANCE:
        log.warning("Novig spread outcome %r disagrees with its market strike %r; left out",
                    outcome.get("description"), strike)
        return None
    return line


def _over_under_side(outcome: dict) -> str | None:
    description = (outcome.get("description") or "").strip().lower()
    if description.startswith("over "):
        return "over"
    if description.startswith("under "):
        return "under"
    return None


def _number(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _title(market_type: str) -> str:
    """'RECEIVING_YARDS' -> 'Receiving Yards'."""
    return " ".join(word.capitalize() for word in market_type.split("_"))


def event_markets_request_body(query_text: str, event_id: str) -> str:
    """The committed EventMarkets_Query with its eventId swapped for this game."""
    payload = json.loads(query_text)
    payload["variables"]["eventId"] = event_id
    return json.dumps(payload)


def parlay_payload(legs: list[OfferedLeg]) -> dict:
    return {"outcomeIds": [leg.venue_ref["outcomeId"] for leg in legs]}


def parse_parlay_quote(response: dict) -> ParlayQuote:
    """One price, no size and no per-leg prices: Novig's quote carries neither."""
    price = float(response["price"])
    if not 0.0 < price < 1.0:
        raise ValueError(f"expected a parlay price between 0 and 1, got {response['price']!r}")
    return ParlayQuote((ParlayTier(1.0 / price, None),), {}, response)


# ---------------------------------------------------------------- the venue

class NovigVenue:
    """Live Novig reads. Thread-safe; caches only the schedule and the leg
    picker's market tree, never a price used to compute a parlay."""

    name = "novig"
    display_name = "Novig"

    def __init__(self, session_factory: Callable[[], Any] | None = None,
                 clock: Callable[[], float] = time.monotonic,
                 wall_clock: Callable[[], datetime] = lambda: datetime.now(timezone.utc),
                 event_markets_query_path: Path = EVENT_MARKETS_QUERY_PATH) -> None:
        self._session_factory = session_factory or _chrome_session
        self._session = None
        self._clock = clock
        self._wall_clock = wall_clock
        self._event_markets_query_path = event_markets_query_path
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
        response = self._send("POST", NOVIG_PARLAY_QUOTE_URL, json=parlay_payload(legs))
        if response.status_code == 400 and DECLINED_PARLAY_TEXT in response.text:
            return None
        body = _json_object("parlay quote", response)
        return _parsed("parlay quote", lambda: parse_parlay_quote(body))

    # -- reads
    def _read_games(self) -> list[Game]:
        with self._lock:
            cached = self._games
        if cached and self._clock() - cached[0] < GAMES_CACHE_SEC:
            return cached[1]
        payload = _json_object("NFL page", self._send("GET", NOVIG_NFL_PAGE_URL))
        games = _parsed("NFL page", lambda: parse_nfl_games(payload, self._wall_clock()))
        with self._lock:
            self._games = (self._clock(), games)
        return games

    def _read_markets(self, event_id: str, fresh: bool) -> dict:
        with self._lock:
            cached = self._markets.get(event_id)
        if cached and not fresh and self._clock() - cached[0] < MARKETS_CACHE_SEC:
            return cached[1]
        try:   # read per call: refresh_novig_query.py rewrites the file while the page is running
            query_text = self._event_markets_query_path.read_text()
        except OSError as exc:
            raise VenueTransportError(f"Novig markets: cannot read the committed query at "
                                      f"{self._event_markets_query_path}: {exc}") from exc
        body = _parsed("committed query", lambda: event_markets_request_body(query_text, event_id))
        payload = _json_object("markets", self._send("POST", NOVIG_GRAPHQL_URL, data=body,
                                                     headers={"Content-Type": "application/json"}))
        if payload.get("errors"):
            raise VenueTransportError(
                f"Novig markets: GraphQL refused the query: {str(payload['errors'])[:200]} "
                "(if it says 'query is not allowed', run python3 mlb_sgp/refresh_novig_query.py)")
        with self._lock:
            self._markets[event_id] = (self._clock(), payload)
        return payload

    def _send(self, method: str, url: str, **kwargs) -> Any:
        """One request; the raw response, whatever its status."""
        if self._session is None:
            self._session = self._session_factory()
        try:
            return self._session.request(method, url, timeout=REQUEST_TIMEOUT_SEC, **kwargs)
        except Exception as exc:   # curl_cffi raises its own error types
            self._session = None
            raise VenueTransportError(f"Novig {method} {url}: {type(exc).__name__}: {exc}") from exc


def _json_object(what: str, response: Any) -> dict:
    if response.status_code != 200:
        raise VenueTransportError(f"Novig {what}: expected HTTP 200, got {response.status_code}: {response.text[:200]}")
    try:
        body = response.json()
    except ValueError as exc:
        raise VenueTransportError(f"Novig {what}: response is not JSON") from exc
    if not isinstance(body, dict):
        raise VenueTransportError(f"Novig {what}: expected a JSON object, got {type(body).__name__}")
    return body


def _parsed(what: str, parse: Callable[[], Any]) -> Any:
    """Run a parser over a Novig answer; a shape it does not expect is the
    venue answering garbage, reported as such rather than as a crash."""
    try:
        return parse()
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        raise VenueTransportError(f"Novig {what} answer has an unexpected shape: "
                                  f"{type(exc).__name__}: {exc}") from exc


def _chrome_session():
    from curl_cffi import requests as curl_requests
    session = curl_requests.Session(impersonate="chrome")
    try:   # seeds Cloudflare's cookies, as mlb_sgp's Novig session does; the API usually answers without them
        session.get(NOVIG_WARM_UP_URL, timeout=10)
    except Exception as exc:
        log.info("Novig warm-up request failed (%s: %s); continuing without it", type(exc).__name__, exc)
    return session
