"""FanDuel venue adapter: NFL games, every SGP-eligible leg (game lines, halves,
quarters, team totals, player props, touchdown scorers), and SGP prices.

Network: the same three FanDuel endpoints mlb_sgp and nfl_specials use, with
the headers / app key / PerimeterX context from mlb_sgp/scraper_fanduel_sgp.py
(one place to fix a rotated token), over a curl_cffi Chrome-TLS session:
    POST scan .../navigation/facet/v1.0/search      NFL event list
    GET  api.sportsbook.fanduel.com/sbapi/event-page?eventId=..&tab=..
    POST sib .../fixedodds/transactional/v1/implyBets?pricePolicy=SUGGESTED
implyBets is a price request: it returns what a bet slip would pay and never
places a bet (placing needs a logged-in account, which this module has no code
for). No DB I/O.

FanDuel shapes (recon 2026-11-19 on Cal's Mac, IND @ HOU, fixtures in tests/fixtures/):
  * event name "Away @ Home"; NFL competitionId 12282733.
  * an event page tab returns only that tab's markets. "same-game-parlay-"
    holds ~150 SGP markets; full-game team alt totals and 1H 2-way winner sit
    on "popular", 1H / 2H alt totals on "1st-half" / "2nd-half", quarter
    2-way winners on the quarter tabs, quarter player props on "quarter-props".
    All are read and merged by marketId.
  * only markets with sgmMarket true go in an SGP; the rest are left out.
  * runners: {selectionId, runnerName, handicap, runnerStatus,
    result.type HOME/AWAY (team runners), winRunnerOdds.trueOdds.decimalOdds}.
    selectionIds repeat across markets ("Over" is one id everywhere), so a leg
    is (marketId, selectionId).
  * main lines carry the number in `handicap` (spread signed per team; totals
    and player O/U "Name Over"/"Name Under"); alternate ladders carry it in
    the runner name ("Houston Texans (-7.5)", "Over (44.5)") with handicap 0;
    player ladders are "Name 25+ Yards" with no other side.
  * implyBets answer: betCombinations, one SINGLE per leg plus the parlay;
    the parlay entry has isSGM true, winAvgOdds.trueOdds.decimalOdds and
    betMaxStake. No isSGM entry = FanDuel will not build that SGP (a decline),
    e.g. a moneyline and spread on the same team.
"""
from __future__ import annotations

import json
import logging
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from sgp_compare.legs import CanonicalLeg, Game, OfferedLeg, ParlayQuote, ParlayTier
from sgp_compare.venues.base import VenueTransportError

log = logging.getLogger(__name__)

NFL_COMPETITION_ID = 12282733
REQUEST_TIMEOUT_SEC = 20
GAMES_CACHE_SEC = 300          # the schedule changes daily
MARKETS_CACHE_SEC = 20         # the leg picker; pricing always re-reads
PRIMARY_TAB = "same-game-parlay-"
# Tabs holding SGP markets the SGP tab leaves out. Best effort: one failing
# drops its markets (logged), the SGP tab failing fails the read.
EXTRA_TABS = ("popular", "1st-half", "2nd-half", "1st-quarter", "2nd-quarter",
              "3rd-quarter", "4th-quarter", "quarter-props")

ALT_LINE_RE = re.compile(r"^(.+?) \(([+-]?\d+(?:\.\d+)?)\)$")     # "Over (44.5)", "Houston Texans (-7.5)"
PLAYER_OU_RE = re.compile(r"^(.+) (Over|Under)$")                  # "Nico Collins Over"
PLAYER_LADDER_RE = re.compile(r"^(.+?) (\d+)\+ (.+)$")             # "Nico Collins 5+ Receptions"

# marketType prefix -> period. Checked in order; no prefix = full game.
PERIOD_PREFIXES = (
    ("1ST_HALF_", "1H"), ("FIRST_HALF_", "1H"), ("2ND_HALF_", "2H"), ("SECOND_HALF_", "2H"),
    ("1ST_QUARTER_", "1Q"), ("2ND_QUARTER_", "2Q"), ("3RD_QUARTER_", "3Q"), ("4TH_QUARTER_", "4Q"),
)
# "AWAY_-_FIRST_HALF_ALTERNATE_TOTAL_POINTS": the period sits after the team.
TEAM_PERIOD_RE = re.compile(r"^(AWAY|HOME)_-_(FIRST|SECOND)_HALF_ALTERNATE_TOTAL_POINTS$")

PLAYER_STATS = {   # FanDuel stat words (after PLAYER_X_[ALT_]) -> canonical stat
    "PASSING_YARDS": "passing_yards", "PASSING_TOUCHDOWNS": "passing_touchdowns",
    "PASSING_+_RUSHING_YARDS": "passing_rushing_yards", "RECEIVING_YARDS": "receiving_yards",
    "RECEPTIONS": "receptions", "RUSHING_YARDS": "rushing_yards",
    "RUSHING_+_RECEIVING_YARDS": "rushing_receiving_yards",
}
PLAYER_RE = re.compile(r"^PLAYER_X_(?:(1ST|2ND|3RD|4TH)_QUARTER_)?(ALT_)?(.+?)(?:_(?:HIGH|MEDIUM|LOW))?$")
QUARTER_OF = {"1ST": "1Q", "2ND": "2Q", "3RD": "3Q", "4TH": "4Q"}

# Yes-only player markets, one runner per player. Stat names match ProphetX's
# for the same bet so legs line up across venues.
TD_SCORER_STATS = {
    "ANY_TIME_TOUCHDOWN_SCORER": ("FG", "anytime_td"),
    "ANYTIME_1ST_HALF_TD_SCORER": ("1H", "anytime_td"),
    "ANYTIME_2ND_HALF_TD_SCORER": ("2H", "anytime_td"),
    "X_TOUCHDOWN_SCORER_QB": ("FG", "to_score_first_touchdown"),
    "LAST_TOUCHDOWN_SCORER": ("FG", "to_score_last_touchdown"),
    "TO_SCORE_2+_TOUCHDOWNS": ("FG", "to_score_at_least_2_touchdowns"),
    "TO_SCORE_3+_TOUCHDOWNS": ("FG", "to_score_at_least_3_touchdowns"),
}


# ---------------------------------------------------------------- parsing (pure)

def parse_nfl_games(scan_payload: dict, now: datetime) -> list[Game]:
    """Not-yet-started NFL games ("Away @ Home") from the scan response, soonest first."""
    events = (scan_payload.get("attachments") or {}).get("events") or {}
    games = []
    for event in events.values():
        if event.get("competitionId") not in (None, NFL_COMPETITION_ID) or event.get("inPlay"):
            continue
        away, sep, home = (event.get("name") or "").partition(" @ ")
        if not sep or not event.get("openDate"):
            continue
        start = datetime.fromisoformat(event["openDate"].replace("Z", "+00:00")).astimezone(timezone.utc)
        if start > now:
            games.append(Game(str(event["eventId"]), home.strip(), away.strip(), start))
    return sorted(games, key=lambda g: g.game_start_time)


def markets_of_pages(pages: list[dict]) -> list[dict]:
    """Open SGP-eligible markets across event-page tabs, first seen per marketId."""
    seen: dict[str, dict] = {}
    for page in pages:
        markets = (page.get("attachments") or {}).get("markets") or {}
        for market in markets.values():
            if market.get("sgmMarket") and market.get("marketStatus") == "OPEN" and not market.get("inPlay"):
                seen.setdefault(str(market["marketId"]), market)
    return list(seen.values())


def parse_legs(game: Game, markets: list[dict]) -> list[OfferedLeg]:
    legs: list[OfferedLeg] = []
    for market in markets:
        legs.extend(_legs_of_market(game, market))
    return legs


def leg_id(market_id: str, selection_id: int, handicap: float) -> str:
    """The handicap is in the id so a moved main line reads as a new leg
    (pricing then says "line moved") instead of silently pricing a new number."""
    return f"{market_id}:{selection_id}:{handicap:g}"


def _legs_of_market(game: Game, market: dict) -> list[OfferedLeg]:
    runners = [r for r in market.get("runners") or [] if r.get("runnerStatus") == "ACTIVE"]
    priced = [(r, _american(r)) for r in runners]
    priced = [(r, odds) for r, odds in priced if odds is not None]
    if not priced:
        return []
    market_id = str(market["marketId"])
    ids = {r["selectionId"]: leg_id(market_id, r["selectionId"], _handicap(r)) for r, _ in priced}
    odds_of = {r["selectionId"]: odds for r, odds in priced}
    canon = {r["selectionId"]: _canonical(game, market, r) for r, _ in priced}
    opposite = _opposites([r for r, _ in priced], canon)
    alt_market = "ALT" in (market.get("marketType") or "").upper()
    legs = []
    for runner, odds in priced:
        sel = runner["selectionId"]
        other = opposite.get(sel)
        legs.append(OfferedLeg(
            leg_id=ids[sel],
            market_label=market.get("marketName") or "",
            category=_category(market),
            side_label=_side_label(runner),
            canonical=canon[sel],
            american_odds=odds,
            opposite_american_odds=odds_of.get(other) if other is not None else None,
            available_stake=None,           # FanDuel states a max stake per bet slip, not per single
            opposite_leg_id=ids.get(other) if other is not None else None,
            market_key=market_id,
            is_main_line=not (alt_market or PLAYER_LADDER_RE.match(runner.get("runnerName") or "")),
            venue_ref={"marketId": market_id, "selectionId": sel},
        ))
    return legs


def _american(runner: dict) -> float | None:
    """Unrounded American odds from FanDuel's true decimal price."""
    odds = runner.get("winRunnerOdds") or {}
    decimal = ((odds.get("trueOdds") or {}).get("decimalOdds") or {}).get("decimalOdds")
    if isinstance(decimal, (int, float)) and decimal > 1.0:
        return round(_decimal_to_american(float(decimal)), 2)
    return None


def _decimal_to_american(decimal_odds: float) -> float:
    """2.5 -> +150, 1.5 -> -200 (unrounded)."""
    if decimal_odds >= 2.0:
        return (decimal_odds - 1.0) * 100.0
    return -100.0 / (decimal_odds - 1.0)


def _handicap(runner: dict) -> float:
    value = runner.get("handicap")
    return float(value) if isinstance(value, (int, float)) else 0.0


def _side_label(runner: dict) -> str:
    name = runner.get("runnerName") or ""
    hc = _handicap(runner)
    if hc and not ALT_LINE_RE.match(name):
        if name in ("Over", "Under") or name.endswith((" Over", " Under")):
            return f"{name} {hc:g}"
        return f"{name} {hc:+g}"
    return name


def _opposites(runners: list[dict], canon: dict) -> dict[int, int]:
    """selectionId -> the other side's selectionId, when the market has one.

    Pairs are found by meaning, never by runner count: Over x / Under x of the
    same stat and subject, home / away moneylines, home -x / away +x spreads,
    and a plain Yes / No market. A yes-only ladder ("25+ Yards", TD scorers,
    a drive ladder that happens to list two rungs) or a 3-way market has none."""
    flip = {"over": "under", "under": "over", "home": "away", "away": "home"}
    by_key: dict[tuple, int] = {}
    for runner in runners:
        leg = canon.get(runner["selectionId"])
        if leg is not None and leg.side in flip:
            by_key[(leg.period, leg.stat, leg.subject, leg.side, leg.line)] = runner["selectionId"]
    pairs = {}
    for (period, stat, subject, side, line), sel in by_key.items():
        other_line = -line if stat == "spread" and line is not None else line
        other = by_key.get((period, stat, subject, flip[side], other_line))
        if other is not None:
            pairs[sel] = other
    names = {(r.get("runnerName") or "").lower(): r["selectionId"] for r in runners}
    if len(runners) == 2 and set(names) == {"yes", "no"}:
        pairs[names["yes"]], pairs[names["no"]] = names["no"], names["yes"]
    return pairs


def _split_period(market_type: str) -> tuple[str, str]:
    for prefix, period in PERIOD_PREFIXES:
        if market_type.startswith(prefix):
            return period, market_type[len(prefix):]
    return "FG", market_type


def _team_side(runner: dict, game: Game) -> str | None:
    kind = ((runner.get("result") or {}).get("type") or "").upper()
    if kind in ("HOME", "AWAY"):
        return kind.lower()
    name = runner.get("runnerName") or ""
    match = ALT_LINE_RE.match(name)
    team = match.group(1) if match else name
    return {game.home_team: "home", game.away_team: "away"}.get(team)


def _canonical(game: Game, market: dict, runner: dict) -> CanonicalLeg | None:
    market_type = (market.get("marketType") or "").upper()
    name = runner.get("runnerName") or ""
    team_of = {"home": game.home_team, "away": game.away_team}

    player = PLAYER_RE.match(market_type)
    if player and player.group(3) in PLAYER_STATS:
        period = QUARTER_OF.get(player.group(1) or "", "FG")
        stat = PLAYER_STATS[player.group(3)]
        subject = (market.get("marketName") or "").split(" - ")[0].strip()
        ou = PLAYER_OU_RE.match(name)
        if ou and not player.group(2):
            return CanonicalLeg(period, stat, subject, ou.group(2).lower(), _handicap(runner))
        ladder = PLAYER_LADDER_RE.match(name)
        if ladder:   # "5+ Receptions" = over 4.5
            return CanonicalLeg(period, stat, ladder.group(1), "over", int(ladder.group(2)) - 0.5)
        return None

    if market_type in TD_SCORER_STATS:
        period, stat = TD_SCORER_STATS[market_type]
        return CanonicalLeg(period, stat, name, "yes", None)

    team_alt = TEAM_PERIOD_RE.match(market_type)
    if team_alt:
        period = "1H" if team_alt.group(2) == "FIRST" else "2H"
        return _over_under(period, "team_total", team_of[team_alt.group(1).lower()], name, runner)

    period, rest = _split_period(market_type)
    if market_type == "MONEY_LINE" or rest == "WINNER":
        side = _team_side(runner, game)
        return CanonicalLeg(period, "moneyline", None, side, None) if side else None
    if rest == "WINNER_THREE_WAY":
        side = _team_side(runner, game)
        if side:
            return CanonicalLeg(period, "win_3way", team_of[side], "yes", None)
        return CanonicalLeg(period, "win_3way", "Draw", "yes", None) if name == "Tie" else None
    if rest in ("MATCH_HANDICAP_(2-WAY)", "HANDICAP", "ALTERNATE_HANDICAP"):
        side = _team_side(runner, game)
        match = ALT_LINE_RE.match(name)
        line = float(match.group(2)) if match else _handicap(runner)
        return CanonicalLeg(period, "spread", None, side, line) if side else None
    if rest in ("TOTAL_POINTS_(OVER/UNDER)", "TOTAL", "ALTERNATE_TOTAL", "ALTERNATE_TOTAL_POINTS"):
        return _over_under(period, "total", None, name, runner)
    team_total = re.match(r"^(AWAY|HOME)_TEAM_(TOTAL_POINTS|ALTERNATE_TOTAL)$", rest)
    if team_total:
        return _over_under(period, "team_total", team_of[team_total.group(1).lower()], name, runner)
    return None


def _over_under(period: str, stat: str, subject: str | None, name: str, runner: dict) -> CanonicalLeg | None:
    match = ALT_LINE_RE.match(name)
    side, line = (match.group(1), float(match.group(2))) if match else (name, _handicap(runner))
    if side not in ("Over", "Under"):
        return None
    return CanonicalLeg(period, stat, subject, side.lower(), line)


def _category(market: dict) -> str:
    market_type = (market.get("marketType") or "").upper()
    if market_type.startswith("PLAYER_X_") or market_type.startswith(("MOST_", "PLAYERS_WITH_")):
        return "Player Props"
    if "TOUCHDOWN" in market_type or "TD_SCORER" in market_type:
        return "Touchdowns"
    if "DRIVE" in market_type:
        return "Drives"
    period, _ = _split_period(market_type)
    if TEAM_PERIOD_RE.match(market_type):
        period = "1H" if "FIRST" in market_type else "2H"
    if period in ("1H", "2H"):
        return "Halves"
    if period != "FG":
        return "Quarters"
    return "Game Lines"


def implybets_payload(legs: list[OfferedLeg]) -> dict:
    return {"betLegs": [{"legType": "SIMPLE_SELECTION",
                         "betRunners": [{"runner": {"marketId": leg.venue_ref["marketId"],
                                                    "selectionId": leg.venue_ref["selectionId"]}}]}
                        for leg in legs]}


def parse_parlay_quote(response: dict, legs: list[OfferedLeg]) -> ParlayQuote | None:
    """The SGP entry covering every leg, or None when FanDuel will not build it."""
    for combination in response.get("betCombinations") or []:
        if not combination.get("isSGM") or len(combination.get("legCombinations") or []) != len(legs):
            continue
        decimal = (((combination.get("winAvgOdds") or {}).get("trueOdds") or {})
                   .get("decimalOdds") or {}).get("decimalOdds")
        if not isinstance(decimal, (int, float)) or decimal <= 1.0:
            return None
        max_stake = combination.get("betMaxStake")
        tier = ParlayTier(float(decimal), float(max_stake) if isinstance(max_stake, (int, float)) else None)
        return ParlayQuote((tier,), {}, response)
    return None


# ---------------------------------------------------------------- the venue

class FanDuelVenue:
    """Live FanDuel reads. Thread-safe; caches only the schedule and the leg
    picker's markets, never a price used to compute a parlay."""

    name = "fanduel"
    display_name = "FanDuel"

    def __init__(self, session_factory: Callable[[], Any] | None = None,
                 api: dict[str, Any] | None = None,
                 clock: Callable[[], float] = time.monotonic) -> None:
        self._session_factory = session_factory or _chrome_session
        self._api = api
        self._local = threading.local()     # one session per thread: tabs load in parallel
        self._clock = clock
        self._lock = threading.Lock()
        self._games: tuple[float, list[Game]] | None = None
        self._markets: dict[str, tuple[float, list[dict]]] = {}

    def list_games(self) -> list[Game]:
        with self._lock:
            cached = self._games
        if cached and self._clock() - cached[0] < GAMES_CACHE_SEC:
            return cached[1]
        api = self._fd_api()
        body = {"filter": {"competitionIds": [NFL_COMPETITION_ID],
                           "contentGroup": {"language": "en", "regionCode": "NAMERICA"},
                           "marketLevels": ["AVB_EVENT"], "maxResults": 100,
                           "productTypes": ["SPORTSBOOK"], "selectBy": "FIRST_TO_START"},
                "facets": [{"type": "EVENT"}], "currencyCode": "USD"}
        payload = self._request("POST", api["scan_url"], data=json.dumps(body), json_body=True)
        games = parse_nfl_games(payload, datetime.now(timezone.utc))
        with self._lock:
            self._games = (self._clock(), games)
        return games

    def list_legs(self, game: Game, fresh: bool = False) -> list[OfferedLeg]:
        return parse_legs(game, self._read_markets(game.venue_game_id, fresh))

    def price_parlay(self, game: Game, legs: list[OfferedLeg]) -> ParlayQuote | None:
        response = self._request("POST", self._fd_api()["imply_url"],
                                 data=json.dumps(implybets_payload(legs)), json_body=True)
        return parse_parlay_quote(response, legs)

    # -- reads
    def _read_markets(self, event_id: str, fresh: bool) -> list[dict]:
        with self._lock:
            cached = self._markets.get(event_id)
        if cached and not fresh and self._clock() - cached[0] < MARKETS_CACHE_SEC:
            return cached[1]
        tabs = (PRIMARY_TAB, *EXTRA_TABS)
        with ThreadPoolExecutor(max_workers=len(tabs)) as pool:
            results = list(pool.map(lambda tab: self._read_tab(event_id, tab), tabs))
        primary = results[0]
        if isinstance(primary, VenueTransportError):
            raise primary
        pages = [primary]
        for tab, result in zip(EXTRA_TABS, results[1:]):
            if isinstance(result, VenueTransportError):
                log.warning("FanDuel tab %s of event %s left out: %s", tab, event_id, result)
            else:
                pages.append(result)
        markets = markets_of_pages(pages)
        with self._lock:
            self._markets[event_id] = (self._clock(), markets)
        return markets

    def _read_tab(self, event_id: str, tab: str) -> dict | VenueTransportError:
        api = self._fd_api()
        params = {"_ak": api["app_key"], "eventId": event_id, "tab": tab,
                  "useCombinedTouchdownsVirtualMarket": "true", "useQuickBets": "true"}
        try:
            return self._request("GET", api["event_page_url"], params=params)
        except VenueTransportError as exc:
            return exc

    def _request(self, method: str, url: str, json_body: bool = False, **kwargs) -> dict:
        api = self._fd_api()
        headers = {**api["headers"], **({"Content-Type": "application/json"} if json_body else {})}
        session = getattr(self._local, "session", None)
        if session is None:
            session = self._local.session = self._session_factory()
        try:
            response = session.request(method, url, headers=headers, timeout=REQUEST_TIMEOUT_SEC, **kwargs)
        except Exception as exc:   # curl_cffi raises its own error types
            self._local.session = None
            raise VenueTransportError(f"FanDuel {method} {url.split('?')[0]}: {type(exc).__name__}: {exc}") from exc
        if response.status_code != 200:
            raise VenueTransportError(f"FanDuel {method} {url.split('?')[0]}: expected HTTP 200, got "
                                      f"{response.status_code}: {response.text[:200]}")
        try:
            return response.json()
        except ValueError as exc:
            raise VenueTransportError(f"FanDuel {method} {url.split('?')[0]}: response is not JSON") from exc

    def _fd_api(self) -> dict[str, Any]:
        if self._api is None:
            self._api = _load_fd_api()
        return self._api


def _load_fd_api() -> dict[str, Any]:
    """Headers and URLs from mlb_sgp/scraper_fanduel_sgp.py, imported on first
    use: that module pulls in the MLB stack (DuckDB paths, Answer Keys), which
    the page should not load unless FanDuel is actually called."""
    mlb_sgp = Path(__file__).resolve().parents[2] / "mlb_sgp"
    if str(mlb_sgp) not in sys.path:
        sys.path.insert(0, str(mlb_sgp))
    import scraper_fanduel_sgp as fd   # noqa: E402  (bare name, as nfl_specials/fd_book.py does)
    return {"headers": dict(fd.FD_HEADERS), "app_key": fd.FD_AK, "scan_url": fd.FD_SCAN_URL,
            "event_page_url": fd.FD_EVENT_PAGE_URL, "imply_url": fd.FD_IMPLY_BETS_URL}


def _chrome_session():
    from curl_cffi import requests as curl_requests
    return curl_requests.Session(impersonate="chrome")
