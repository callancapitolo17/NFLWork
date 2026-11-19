"""Caesars venue adapter: NFL games, every same-game-parlay-eligible leg (game
lines with their alt ladders, team totals, halves/quarters, player props, TD
scorers), and Caesars' correlated SGP price.

Network: api.americanwagering.com (Liberty / ZeroFlucs engine), logged out,
through mlb_sgp.caesars_client.CaesarsClient, used unchanged: it mints the
AWS-WAF token browser-free under `node` (needs node on PATH) and carries the
headers Caesars requires. Reads (state NJ):
    GET  /sb/v4/sports/americanfootball/tabs/{tab}    schedule (two tabs, unioned)
    GET  /sb/v4/events/{event_id}                     every market of one game
    POST /sb/v2/bets/details                          price request for a leg set
`bets/details` is the bet slip's preview call: it returns prices and never
places a bet (placing needs a login and a separate call this module does not
have). No DB I/O.

Caesars shapes (recon 2026-11-19 on Cal's Mac, IND at HOU; captures and notes
in sgp_compare/recon/caesars/, trimmed fixtures in tests/fixtures/):
  * There is no "Games" tab, and an unknown tab id silently returns Featured.
    Featured listed 12 NFL games and the TD Scorer tab 13, so both are read
    and unioned by event id. The NFL competition is named "NFL"; NCAAF shares
    the feed.
  * Names are pipe-wrapped ("|Spread|", "|Houston Texans|").
  * A spread market's `line` is the HOME line (-12.5 = Texans -12.5); the
    away leg is sent with the negated line. The main spread, total and team
    totals appear twice, once with `movingLines.linePrices[]`: the alt ladder,
    each rung {line (home line for spreads), selections:[{selectionType, price}]}.
    An alt rung reuses the main market's selection ids; the requested `line`
    picks the rung.
  * Player yardage / reception props are OVER-ONLY ladders: one selection per
    threshold named "|60+|", typed "over", no market line. 60+ = over 59.5.
  * TD scorer markets hold one selection per player, typed with his side;
    there is no "no" side.
  * metadata.period: MATCH_INC_OT, FIRST_HALF, FIRST_QUARTER (no 2H, Q2-Q4).
  * betLimits.maxBet is in cents ($50k on the main spread).
  * The parlay is response.parlays[0]: price {decimal, american}; a refused
    combination comes back HTTP 200 with parlays[0].errors and no price.
    bets/details reports no stake limit, so a tier's max_stake is None.
"""
from __future__ import annotations

import logging
import threading
import time
from datetime import datetime, timezone
from typing import Any, Callable, Protocol

from sgp_compare.legs import CanonicalLeg, Game, OfferedLeg, ParlayQuote, ParlayTier
from sgp_compare.venues.base import VenueTransportError

log = logging.getLogger(__name__)

CAESARS_STATE = "nj"
NFL_COMPETITION_NAME = "NFL"
SCHEDULE_TAB_PATHS = (
    "v4/sports/americanfootball/tabs/HIGHLIGHTS%7CFeatured",
    "v4/sports/americanfootball/tabs/SCHEDULE%7CTD%20Scorer%20%F0%9F%99%8C",
)
GAMES_CACHE_SEC = 300          # the schedule changes daily
EVENT_CACHE_SEC = 20           # the leg picker; pricing always re-reads
CENTS_PER_DOLLAR = 100

PERIOD_OF_CAESARS_PERIOD = {"MATCH_INC_OT": "FG", "FIRST_HALF": "1H", "FIRST_QUARTER": "1Q"}

# Named the way the ProphetX adapter names the same bets, so legs match across venues.
SCORER_STAT_OF_MARKET_TYPE = {
    "touchdown-scorer-anytime": "anytime_td",
    "touchdown-scorer-first": "to_score_first_touchdown",
    "touchdown-scorer-2-plus": "to_score_at_least_2_touchdowns",
    "touchdown-scorer-3-plus": "to_score_at_least_3_touchdowns",
}
# A different player can be added from the same scorer market, except first scorer.
SCORER_MARKETS_ALLOWING_SEVERAL_PLAYERS = {"touchdown-scorer-anytime", "touchdown-scorer-2-plus",
                                           "touchdown-scorer-3-plus"}
PLAYER_LADDER_PREFIX = "player-alt-"
PLAYER_STAT_OF_LADDER = {"receiving-rushing-yards": "rushing_receiving_yards"}
# "60+" means 60 or more, the same bet as over 59.5.
LADDER_THRESHOLD_TO_LINE = -0.5


# ---------------------------------------------------------------- parsing (pure)

def clean_name(name: str | None) -> str:
    """'|Houston Texans| |First Score|' -> 'Houston Texans First Score'."""
    return " ".join(part.strip() for part in (name or "").split("|") if part.strip())


def parse_nfl_games(tabs_payloads: list[dict], now: datetime) -> list[Game]:
    """Not-yet-started NFL games across the schedule tabs, each once, soonest first."""
    games: dict[str, Game] = {}
    for payload in tabs_payloads:
        for competition in payload.get("competitions") or []:
            if (competition.get("name") or "").strip() != NFL_COMPETITION_NAME:
                continue
            for event in competition.get("events") or []:
                game = _game_from_event(event)
                if game is not None and not event.get("started") and game.game_start_time > now:
                    games[game.venue_game_id] = game
    return sorted(games.values(), key=lambda g: g.game_start_time)


def _game_from_event(event: dict) -> Game | None:
    metadata = event.get("metadata") or {}
    home, away, start = metadata.get("homeTeamName"), metadata.get("awayTeamName"), event.get("startTime")
    if not (event.get("id") and home and away and start):
        return None
    start_time = datetime.fromisoformat(start.replace("Z", "+00:00")).astimezone(timezone.utc)
    return Game(venue_game_id=str(event["id"]), home_team=home, away_team=away, game_start_time=start_time)


def parse_legs(game: Game, event: dict) -> list[OfferedLeg]:
    """Every leg of the game's active, SGP-eligible markets.

    A market listed twice (once plain, once with its alt ladder) is read once,
    from the copy with the ladder, under the group it first appears in.
    Boosts and specials (byoEligible false) are left out: Caesars will not
    put them in a same-game parlay.
    """
    competition_id = event.get("competitionId")
    legs: list[OfferedLeg] = []
    for market, category in _unique_markets(event):
        if not (market.get("active") and market.get("display", True) and market.get("byoEligible")):
            continue
        context = _MarketContext(game, market, category, str(event.get("id")), competition_id)
        legs.extend(_legs_of_market(context))
    return legs


def _unique_markets(event: dict) -> list[tuple[dict, str]]:
    by_id: dict[str, tuple[dict, str]] = {}
    for group in event.get("keyMarketGroups") or []:
        category = group.get("marketDisplayGroupDisplayName") or "Other"
        for market in group.get("markets") or []:
            market_id = market.get("id")
            if market_id not in by_id:
                by_id[market_id] = (market, category)
            elif market.get("movingLines") and not by_id[market_id][0].get("movingLines"):
                by_id[market_id] = (market, by_id[market_id][1])
    return list(by_id.values())


class _MarketContext:
    def __init__(self, game: Game, market: dict, category: str, event_id: str, competition_id: str | None):
        self.game = game
        self.market = market
        self.category = category
        self.event_id = event_id
        self.competition_id = competition_id
        self.metadata = market.get("metadata") or {}
        self.type_code = self.metadata.get("marketTypeCode") or ""
        self.period = PERIOD_OF_CAESARS_PERIOD.get(self.metadata.get("period"), "FG")
        self.label = clean_name(market.get("name"))
        self.limit_dollars = _max_bet_dollars(market)

    def team_of_side(self, side: str | None) -> str | None:
        return {"home": self.game.home_team, "away": self.game.away_team}.get(side)


def _max_bet_dollars(market: dict) -> float | None:
    max_bet_cents = (market.get("betLimits") or {}).get("maxBet")
    return max_bet_cents / CENTS_PER_DOLLAR if isinstance(max_bet_cents, (int, float)) else None


def _legs_of_market(ctx: _MarketContext) -> list[OfferedLeg]:
    market_type = ctx.market.get("type")
    if market_type in ("two-way-handicap", "over-under"):
        return _legs_of_line_market(ctx)
    if ctx.type_code == "head-to-head":
        return _two_way_legs(ctx, _canonical_moneyline)
    if ctx.type_code in SCORER_STAT_OF_MARKET_TYPE:
        return _scorer_legs(ctx)
    if ctx.type_code.startswith(PLAYER_LADDER_PREFIX):
        return _player_ladder_legs(ctx)
    return _each_selection_legs(ctx)


# -- spreads and totals: a main line plus, when present, its alt ladder

def _legs_of_line_market(ctx: _MarketContext) -> list[OfferedLeg]:
    """Both sides at every number. Every rung is one market_key, so a parlay
    takes one number per market."""
    selection_of_type = {s.get("type"): s for s in ctx.market.get("selections") or [] if s.get("active")}
    main_line = _as_float(ctx.market.get("line"))
    rungs = _rungs(ctx.market, selection_of_type, main_line)
    legs: list[OfferedLeg] = []
    for market_line, price_of_type in rungs:
        sides = [t for t in ("home", "away", "over", "under") if t in price_of_type and t in selection_of_type]
        if len(sides) != 2:
            continue
        for side, other in ((sides[0], sides[1]), (sides[1], sides[0])):
            legs.append(_line_leg(ctx, selection_of_type, market_line, side, other, price_of_type,
                                  is_main=market_line == main_line))
    return legs


def _rungs(market: dict, selection_of_type: dict, main_line: float | None) -> list[tuple[float, dict]]:
    """[(market line, {selection type: price dict})], the main line included once."""
    rungs: dict[float, dict] = {}
    if main_line is not None:
        rungs[main_line] = {t: s.get("price") for t, s in selection_of_type.items() if _american(s.get("price"))}
    for rung in (market.get("movingLines") or {}).get("linePrices") or []:
        line = _as_float(rung.get("line"))
        if line is None or not rung.get("active") or line == main_line:
            continue
        rungs[line] = {s.get("selectionType"): s.get("price") for s in rung.get("selections") or []
                       if _american(s.get("price"))}
    return sorted(rungs.items())


def _line_leg(ctx: _MarketContext, selection_of_type: dict, market_line: float, side: str, other: str,
              price_of_type: dict, is_main: bool) -> OfferedLeg:
    selection = selection_of_type[side]
    side_line = -market_line if side == "away" else market_line   # spread line is the home line
    canonical = _canonical_line(ctx, side, side_line)
    return OfferedLeg(
        leg_id=_line_leg_id(selection["id"], market_line),
        market_label=ctx.label,
        category=ctx.category,
        side_label=_line_side_label(ctx, side, side_line),
        canonical=canonical,
        american_odds=_american(price_of_type[side]),
        opposite_american_odds=_american(price_of_type[other]),
        available_stake=ctx.limit_dollars,
        opposite_leg_id=_line_leg_id(selection_of_type[other]["id"], market_line),
        market_key=str(ctx.market.get("id")),
        is_main_line=is_main,
        venue_ref=_bet_details_leg(ctx, selection, price_of_type[side], side_line),
    )


def _line_leg_id(selection_id: str, market_line: float) -> str:
    return f"{selection_id}@{market_line:g}"


def _line_side_label(ctx: _MarketContext, side: str, side_line: float) -> str:
    if side in ("over", "under"):
        return f"{side.capitalize()} {side_line:g}"
    return f"{ctx.team_of_side(side)} {side_line:+g}"


def _canonical_line(ctx: _MarketContext, side: str, side_line: float) -> CanonicalLeg | None:
    if ctx.type_code == "line":
        return CanonicalLeg(ctx.period, "spread", None, side, side_line)
    if ctx.type_code == "total":
        return CanonicalLeg(ctx.period, "total", None, side, side_line)
    if ctx.type_code in ("total-home", "total-away"):
        team = ctx.team_of_side(ctx.type_code.removeprefix("total-"))
        return CanonicalLeg(ctx.period, "team_total", team, side, side_line)
    return None


# -- two-way markets without a number (moneyline, overtime yes/no)

def _two_way_legs(ctx: _MarketContext, canonical_of: Callable[[_MarketContext, dict], CanonicalLeg | None]
                  ) -> list[OfferedLeg]:
    selections = [s for s in ctx.market.get("selections") or [] if s.get("active") and _american(s.get("price"))]
    if len(selections) != 2:
        return _each_selection_legs(ctx)
    legs = []
    for selection, other in ((selections[0], selections[1]), (selections[1], selections[0])):
        legs.append(_selection_leg(ctx, selection, canonical_of(ctx, selection), other=other))
    return legs


def _canonical_moneyline(ctx: _MarketContext, selection: dict) -> CanonicalLeg | None:
    side = selection.get("type")
    return CanonicalLeg(ctx.period, "moneyline", None, side, None) if side in ("home", "away") else None


# -- one-sided markets

def _scorer_legs(ctx: _MarketContext) -> list[OfferedLeg]:
    stat = SCORER_STAT_OF_MARKET_TYPE[ctx.type_code]
    several = ctx.type_code in SCORER_MARKETS_ALLOWING_SEVERAL_PLAYERS
    legs = []
    for selection in _priced_selections(ctx):
        player = clean_name(selection.get("name"))
        legs.append(_selection_leg(ctx, selection, CanonicalLeg(ctx.period, stat, player, "yes", None),
                                   market_key_suffix=selection["id"] if several else None))
    return legs


def _player_ladder_legs(ctx: _MarketContext) -> list[OfferedLeg]:
    """'Nico Collins - Receiving Yards' 60+ -> receiving_yards over 59.5, one
    threshold per parlay (they share the market)."""
    ladder = ctx.type_code.removeprefix(PLAYER_LADDER_PREFIX)
    stat = PLAYER_STAT_OF_LADDER.get(ladder, ladder.replace("-", "_"))
    player = ctx.metadata.get("player") or ctx.label
    legs = []
    for selection in _priced_selections(ctx):
        threshold = _threshold(clean_name(selection.get("name")))
        canonical = None
        if threshold is not None and selection.get("type") == "over":
            canonical = CanonicalLeg(ctx.period, stat, player, "over", threshold + LADDER_THRESHOLD_TO_LINE)
        legs.append(_selection_leg(ctx, selection, canonical))
    return sorted(legs, key=lambda leg: leg.canonical.line if leg.canonical else float("inf"))


def _threshold(name: str) -> float | None:
    """'60+' -> 60.0; anything else -> None."""
    return _as_float(name[:-1]) if name.endswith("+") else None


def _each_selection_legs(ctx: _MarketContext) -> list[OfferedLeg]:
    """Markets read as a list of outcomes (3-way, winning margin, correct
    score): each selection is a yes on that outcome."""
    legs = []
    for selection in _priced_selections(ctx):
        legs.append(_selection_leg(ctx, selection, _canonical_outcome(ctx, selection)))
    return legs


def _canonical_outcome(ctx: _MarketContext, selection: dict) -> CanonicalLeg | None:
    if ctx.type_code == "win-draw-win":
        subject = ctx.team_of_side(selection.get("type")) or "Draw"
        return CanonicalLeg(ctx.period, "win_3way", subject, "yes", None)
    return None


def _priced_selections(ctx: _MarketContext) -> list[dict]:
    return [s for s in ctx.market.get("selections") or [] if s.get("active") and _american(s.get("price"))]


def _selection_leg(ctx: _MarketContext, selection: dict, canonical: CanonicalLeg | None,
                   other: dict | None = None, market_key_suffix: str | None = None) -> OfferedLeg:
    market_key = str(ctx.market.get("id"))
    return OfferedLeg(
        leg_id=str(selection["id"]),
        market_label=ctx.label,
        category=ctx.category,
        side_label=clean_name(selection.get("name")),
        canonical=canonical,
        american_odds=_american(selection.get("price")),
        opposite_american_odds=_american(other.get("price")) if other else None,
        available_stake=ctx.limit_dollars,
        opposite_leg_id=str(other["id"]) if other else None,
        market_key=f"{market_key}:{market_key_suffix}" if market_key_suffix else market_key,
        is_main_line=other is not None,
        venue_ref=_bet_details_leg(ctx, selection, selection.get("price"), None),
    )


def _bet_details_leg(ctx: _MarketContext, selection: dict, price: dict, side_line: float | None) -> dict:
    """One leg of the bets/details body, in the bet slip's own shape (the
    same shape mlb_sgp/caesars.py sends)."""
    return {
        "selectionId": selection["id"],
        "eventId": ctx.event_id,
        "marketId": ctx.market.get("id"),
        "competitionId": ctx.competition_id,
        "priceType": "fp",
        "related": True,
        "eachWay": False,
        "stakePerLine": 0,
        "line": side_line,
        "selectionType": selection.get("type"),
        "name": selection.get("name"),
        "price": price,
        "origPrice": price,
    }


def _american(price: dict | None) -> float | None:
    american = (price or {}).get("a")
    return float(american) if isinstance(american, (int, float)) and abs(american) >= 100 else None


def _as_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# -- parlay request and answer

def bet_details_body(legs: list[OfferedLeg]) -> dict:
    """`combinationSelections` must stay empty: Caesars infers the SGP from
    same-event legs (mlb_sgp recon)."""
    return {"channel": "desktop", "channelDetail": "cordova-desktop", "includeMultiLine": False,
            "combinationSelections": [], "legs": [dict(leg.venue_ref) for leg in legs]}


def parse_parlay_quote(response: dict, legs: list[OfferedLeg]) -> ParlayQuote | None:
    """None when Caesars refused the combination (no parlay, errors, or a leg
    left out). Raises VenueTransportError when the answer prices something
    other than what was asked, such as a different line."""
    parlays = response.get("parlays") or []
    if not parlays:
        return None
    parlay = parlays[0]
    decimal_odds = (parlay.get("price") or {}).get("decimal")
    if parlay.get("errors") or not isinstance(decimal_odds, (int, float)) or decimal_odds <= 1.0:
        return None
    priced_legs = {leg.get("selectionId"): leg for leg in parlay.get("legs") or []}
    for leg in legs:
        priced = priced_legs.get(leg.venue_ref["selectionId"])
        if priced is None or not priced.get("included"):
            return None
        asked_line = leg.venue_ref.get("line")
        if asked_line is not None and _as_float(priced.get("requestedLine")) != asked_line:
            raise VenueTransportError(f"Caesars priced {leg.market_label} {leg.side_label} at line "
                                      f"{priced.get('requestedLine')}, expected {asked_line:g}")
    return ParlayQuote((ParlayTier(float(decimal_odds), None),), {}, response)


# ---------------------------------------------------------------- the venue

class CaesarsTransport(Protocol):
    def get_json(self, path: str) -> dict | None:
        """GET /sb/{path}; None when Caesars answers 404. Raises VenueTransportError."""

    def post_bet_details(self, body: dict) -> dict:
        """POST /sb/v2/bets/details. Raises VenueTransportError."""


class CaesarsVenue:
    """Live Caesars reads. Thread-safe (one lock around the shared session);
    caches only the schedule and the leg picker's event read, never a price
    used to compute a parlay."""

    name = "caesars"
    display_name = "Caesars"

    def __init__(self, transport_factory: Callable[[], CaesarsTransport] | None = None,
                 clock: Callable[[], float] = time.monotonic) -> None:
        self._transport_factory = transport_factory or _LiveCaesarsTransport
        self._transport: CaesarsTransport | None = None
        self._clock = clock
        self._lock = threading.Lock()
        self._games: tuple[float, list[Game]] | None = None
        self._events: dict[str, tuple[float, dict | None]] = {}

    def list_games(self) -> list[Game]:
        cached = self._games
        if cached and self._clock() - cached[0] < GAMES_CACHE_SEC:
            return cached[1]
        payloads = [self._get(path) or {} for path in SCHEDULE_TAB_PATHS]
        games = parse_nfl_games(payloads, datetime.now(timezone.utc))
        self._games = (self._clock(), games)
        return games

    def list_legs(self, game: Game, fresh: bool = False) -> list[OfferedLeg]:
        event = self._read_event(game.venue_game_id, fresh)
        return parse_legs(game, event) if event else []

    def price_parlay(self, game: Game, legs: list[OfferedLeg]) -> ParlayQuote | None:
        with self._lock:
            response = self._transport_now().post_bet_details(bet_details_body(legs))
        return parse_parlay_quote(response, legs)

    def _read_event(self, event_id: str, fresh: bool) -> dict | None:
        cached = self._events.get(event_id)
        if cached and not fresh and self._clock() - cached[0] < EVENT_CACHE_SEC:
            return cached[1]
        payload = self._get(f"v4/events/{event_id}")
        event = (payload.get("event") or payload) if payload else None
        self._events[event_id] = (self._clock(), event)
        return event

    def _get(self, path: str) -> dict | None:
        with self._lock:
            return self._transport_now().get_json(path)

    def _transport_now(self) -> CaesarsTransport:
        if self._transport is None:
            self._transport = self._transport_factory()
        return self._transport


class _LiveCaesarsTransport:
    """mlb_sgp's CaesarsClient (token mint, headers, retry), imported only
    when the venue first goes live so tests need neither curl_cffi nor node."""

    def __init__(self) -> None:
        from mlb_sgp.caesars_client import CaesarsClient
        self._client = CaesarsClient(state=CAESARS_STATE)

    def get_json(self, path: str) -> dict | None:
        from mlb_sgp._shared import BookTransportError
        self._ensure_token()
        try:
            # The "structure" stage is the one where CaesarsClient returns None on a 404.
            return self._client._get_json(f"{self._client._sb}/{path}", "structure")
        except BookTransportError as exc:
            raise VenueTransportError(f"Caesars GET {path}: {exc}") from exc

    def post_bet_details(self, body: dict) -> dict:
        from mlb_sgp._shared import RETRY_LIVE, BookTransportError, request_with_retry
        from mlb_sgp.caesars_client import _looks_like_json
        self._ensure_token()
        client = self._client
        try:
            response = request_with_retry(
                lambda: client.session.post(f"{client._sb}/v2/bets/details", headers=client._headers(),
                                            cookies=client._cookies, json=body, impersonate="chrome", timeout=30),
                profile=RETRY_LIVE, book="caesars", stage="price", body_ok=_looks_like_json)
        except BookTransportError as exc:
            raise VenueTransportError(f"Caesars POST bets/details: {exc}") from exc
        if response.status_code != 200 or not _looks_like_json(response):
            raise VenueTransportError(f"Caesars POST bets/details: expected HTTP 200 JSON, got "
                                      f"{response.status_code}: {(response.text or '')[:200]}")
        return response.json()

    def _ensure_token(self) -> None:
        if not self._client.ensure_token():
            raise VenueTransportError("Caesars: AWS-WAF token could not be minted or validated "
                                      "(is node on PATH?)")
