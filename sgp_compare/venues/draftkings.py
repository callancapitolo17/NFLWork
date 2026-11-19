"""DraftKings venue adapter: NFL games, every same-game-parlay leg (game lines,
halves, quarters, team totals, player props) and SGP prices.

Network: the same three plain-HTTP endpoints nfl_specials/dk_book.py uses
(curl_cffi Chrome TLS, no login), reused from there:
    league listing   dk_book.DK_LEAGUE_URL       (NFL league 88808)
    SGP leg list     dk_book.DK_SGP_EVENT_URL    /parlays/v1/sgp/events/{event_id}
    SGP price        dk_book.DK_SGP_PRICE_URL    GET, the call DK's SGP builder widget makes
The price call is a read: it returns prices and never places a bet (this module
has no betslip or login code). DraftKings refuses data-center IPs, so this
runs from the Mac's home IP, or through the Mac's nfl_specials dk_relay when
SGP_COMPARE_DK_PROXY is set (it tunnels to dk_book.DK_HOSTS, which covers every
host here). No DB I/O.

DraftKings shapes (recon 2026-11-19 on Cal's Mac, IND at HOU, fixtures in tests/fixtures/):
  * The SGP leg list holds every SGP market of the event (468: game lines,
    periods, team totals, 282 player-prop markets). Each market has
    `marketType.name` ("Receiving Yards O/U - 1st Half"), `category` (the
    shape: HomeAway, HomeAwayByPoints, OverUnderSingleLine, ScorerAnytime ...),
    `groups` (ids of the event's marketGroups, e.g. "Receiving Props").
  * Selections carry `outcomeType` (Home/Away/Tie/Over/Under/Yes/No/ToScoreAnyTime),
    `points` (signed spread from that side, or the total), `milestoneValue`
    ("70+" ladders), `players` / `participants` with venueRole, `displayOdds.american`
    (with a U+2212 minus) and `trueOdds` (unrounded decimal).
  * The price call takes base selections + one candidate market and answers the
    SGP price of base + EACH selection of that market ("trueOdds"). DK does not
    refuse a base it cannot combine: it drops legs (`selectionsNotMapped`) and
    prices the rest, so a base not taken whole prices nothing
    (dk_book.parse_added_leg_odds). A candidate market it will not combine is
    missing from `compatibleMarkets`.
  * DK prices a single leg a bit worse inside the builder than on its listing
    (HOU ML -750 vs -700), so singles come from the listing and parlays only
    from the builder.
Any payload this module cannot parse raises VenueTransportError, so one odd DK
answer shows as a DraftKings error instead of failing every venue's page.
"""
from __future__ import annotations

import os
import re
import threading
import time
from datetime import datetime, timezone
from typing import Any, Callable

from nfl_specials.books import parse_start
from nfl_specials.config import parse_dk_proxy
from nfl_specials.dk_book import (DK_CONNECT_TIMEOUT_SECONDS, DK_FOOTBALL_SPORT_ID, DK_LEAGUE_IDS, DK_LEAGUE_URL,
                                  DK_SGP_EVENT_URL, DK_SGP_PRICE_URL, DK_WARMUP_URL,
                                  MIN_SECONDS_BETWEEN_PRICE_CALLS, added_leg_params, dk_team_to_abbr,
                                  league_listing_params, parse_added_leg_odds)
from nfl_specials.teams import TEAMS
from sgp_compare.legs import CanonicalLeg, Game, OfferedLeg, ParlayQuote, ParlayTier
from sgp_compare.venues.base import VenueTransportError

PROXY_ENV = "SGP_COMPARE_DK_PROXY"
REQUEST_TIMEOUT_SEC = 30          # dk_book's read timeouts: 30 s, and 60 s for the ~4.6 MB SGP leg list
LEG_LIST_TIMEOUT_SEC = 60
HTTP_FORBIDDEN = 403
GAMES_CACHE_SEC = 300
MARKETS_CACHE_SEC = 20            # the leg picker; pricing always re-reads
# One builder answer prices base + every side of a market, so the partition
# cells that differ only in their last leg share it. A fresh leg read (the start
# of every Price click) drops the game's answers, so they never outlive one
# click; this is only the backstop.
ADDED_LEG_CACHE_SEC = 10

# Period words in DK's marketType names: "Spread 1st Half", "Total Alternate - 1st Half",
# "1st Half (3 Way)".
PERIOD_OF_WORDS = {"1st Half": "1H", "2nd Half": "2H", "1st Quarter": "1Q",
                   "2nd Quarter": "2Q", "3rd Quarter": "3Q", "4th Quarter": "4Q"}
_PERIOD_SUFFIX_RE = re.compile(r"\s*(?:-\s*)?(" + "|".join(PERIOD_OF_WORDS) + r")\s*$")

# Scorer categories -> stat, named as ProphetX names them so legs can be matched later.
SCORER_STATS = {
    "ScorerAnytime": ("FG", "anytime_td"),
    "ScorerAnytime1stHalf": ("1H", "anytime_td"),
    "ScorerAnytime2ndHalf": ("2H", "anytime_td"),
    "Scorer2OrMore": ("FG", "to_score_at_least_2_touchdowns"),
    "Scorer3OrMore": ("FG", "to_score_at_least_3_touchdowns"),
    "ScorerFirst": ("FG", "to_score_first_touchdown"),
    "ScorerFirst2ndHalf": ("2H", "to_score_first_touchdown"),
    "ScorerLast": ("FG", "to_score_last_touchdown"),
}
# DK words one stat two ways ("Passing Touchdown Milestones" / "Passing Touchdowns O/U").
STAT_ALIASES = {"passing_touchdown": "passing_touchdowns",
                "rushing_and_receiving_yards": "rushing_receiving_yards"}
GAME_LINES_GROUP = "Game Lines"
POPULAR_GROUP = "Popular"
OTHER_GROUP = "Other"
SIDE_OF_OUTCOME = {"Home": "home", "Away": "away", "Over": "over", "Under": "under", "Yes": "yes", "No": "no"}
OPPOSITE_OUTCOME = {"Home": "Away", "Away": "Home", "Over": "Under", "Under": "Over", "Yes": "No", "No": "Yes"}


# ---------------------------------------------------------------- parsing (pure)

def dk_proxy_from_env() -> str | None:
    """SGP_COMPARE_DK_PROXY, the Mac's dk_relay ("http://<tailnet ip>:8097"); None = straight to DK."""
    return parse_dk_proxy(os.environ.get(PROXY_ENV), env_name=PROXY_ENV)


def full_team_name(dk_name: str) -> str:
    """'HOU Texans' -> 'Houston Texans' (ProphetX's naming); DK's name when unknown."""
    abbr = dk_team_to_abbr(dk_name) if dk_name else None
    return TEAMS[abbr][0] if abbr in TEAMS else dk_name


def american_of(selection: dict) -> float | None:
    """DK's listed American price ('−118' with a U+2212 minus), None when unusable."""
    text = ((selection.get("displayOdds") or {}).get("american") or "").replace("−", "-").replace("+", "")
    try:
        odds = float(text)
    except ValueError:
        return None
    return odds if abs(odds) >= 100 else None


def parse_nfl_games(league_payload: dict, now: datetime) -> list[Game]:
    """Not-yet-started NFL games from the league listing, soonest first.
    Raises KeyError/ValueError/TypeError on a shape it does not know."""
    games = []
    for event in league_payload.get("events") or []:
        roles = {p.get("venueRole"): p.get("name", "") for p in event.get("participants") or []}
        if not roles.get("Home") or not roles.get("Away") or not event.get("startEventDate"):
            continue
        if event.get("status") not in (None, "NOT_STARTED"):
            continue
        start = parse_start(event["startEventDate"]).astimezone(timezone.utc)
        if start <= now:
            continue
        games.append(Game(venue_game_id=str(event["id"]), home_team=full_team_name(roles["Home"]),
                          away_team=full_team_name(roles["Away"]), game_start_time=start))
    return sorted(games, key=lambda g: g.game_start_time)


def parse_legs(sgp_payload: dict) -> list[OfferedLeg]:
    """One leg per live selection of every SGP market of the event."""
    data = sgp_payload.get("data") or {}
    events = data.get("events") or []
    group_names = {g.get("id"): g.get("name", "") for g in (events[0].get("marketGroups") if events else None) or []}
    legs: list[OfferedLeg] = []
    for market in data.get("markets") or []:
        if market.get("isSuspended") or "SGP" not in (market.get("tags") or []):
            continue
        live = [s for s in market.get("selections") or []
                if not s.get("isDisabled") and s.get("status") in (None, "Active") and american_of(s) is not None]
        category = _category(market, group_names)
        for selection in live:
            opposite = _opposite(market, selection, live)
            legs.append(OfferedLeg(
                leg_id=selection["id"],
                market_label=_market_label(market),
                category=category,
                side_label=_side_label(market, selection),
                canonical=_canonical(market, selection),
                american_odds=american_of(selection),
                opposite_american_odds=american_of(opposite) if opposite else None,
                available_stake=None,
                opposite_leg_id=opposite["id"] if opposite else None,
                market_key=str(market["id"]),
                is_main_line=_is_main_line(market),
                venue_ref={"selectionId": selection["id"], "marketId": str(market["id"])},
            ))
    return legs


def _category(market: dict, group_names: dict[str, str]) -> str:
    """The picker tab: Game Lines when DK lists the market there, else its first
    group other than DK's cross-cutting "Popular"."""
    names = [group_names.get(group_id, "") for group_id in market.get("groups") or []]
    if GAME_LINES_GROUP in names:
        return GAME_LINES_GROUP
    return next((name for name in names if name and name != POPULAR_GROUP), OTHER_GROUP)


def _type_name(market: dict) -> str:
    return ((market.get("marketType") or {}).get("name") or market.get("name") or "").strip()


def _period_and_base(type_name: str) -> tuple[str, str]:
    """'Receiving Yards O/U - 1st Half' -> ('1H', 'Receiving Yards O/U');
    '1st Half (3 Way)' -> ('1H', '(3 Way)')."""
    found = _PERIOD_SUFFIX_RE.search(type_name)
    if found:
        return PERIOD_OF_WORDS[found.group(1)], type_name[:found.start()].strip()
    for words, period in PERIOD_OF_WORDS.items():
        if type_name.startswith(words + " "):
            return period, type_name[len(words):].strip()
    return "FG", type_name


def _slug(words: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", words.lower()).strip("_")
    return STAT_ALIASES.get(slug, slug)


def _team_names(market: dict, selection: dict) -> list[str]:
    participants = selection.get("participants") or market.get("participants") or []
    return [p.get("name", "") for p in participants]


def _players(selection: dict) -> list[str]:
    return [p.get("name", "") for p in selection.get("players") or []]


def _market_label(market: dict) -> str:
    name = market.get("name") or _type_name(market)
    team = re.match(r"^([A-Z]{2,3} [A-Za-z0-9 ]+?): (.+)$", name)
    return f"{full_team_name(team.group(1))}: {team.group(2)}" if team else name


def _side_label(market: dict, selection: dict) -> str:
    outcome = selection.get("outcomeType")
    points = selection.get("points")
    name = selection.get("name") or ""
    if outcome in ("Home", "Away"):
        team = full_team_name(name)
        return f"{team} {points:+g}" if points is not None else team
    if outcome in ("Over", "Under") and points is not None:
        return f"{outcome} {points:g}"
    return name


def _is_main_line(market: dict) -> bool:
    category = market.get("category")
    type_name = _type_name(market)
    return not ("Alternate" in type_name or "Milestones" in type_name or category == "OverUnderMultipleLine")


def _opposite(market: dict, selection: dict, live: list[dict]) -> dict | None:
    """The selection that wins exactly when `selection` loses, when the market
    has one. Only for legs whose canonical is understood: the pricer's push and
    tie checks read the canonical, so an unknown market gets no other side and
    no vig-removed factor rather than a partition that may be missing an outcome."""
    outcome = selection.get("outcomeType")
    wanted = OPPOSITE_OUTCOME.get(outcome)
    if wanted is None or market.get("category") == "HomeTieAway" or _canonical(market, selection) is None:
        return None
    points = selection.get("points")
    candidates = [s for s in live if s is not selection and s.get("outcomeType") == wanted
                  and _players(s) == _players(selection)]
    if outcome in ("Home", "Away") and points is not None:      # spread: the other team at the opposite number
        candidates = [s for s in candidates if s.get("points") is not None and abs(s["points"] + points) < 1e-6]
    elif points is not None:                                      # total: same number
        candidates = [s for s in candidates if s.get("points") is not None and abs(s["points"] - points) < 1e-6]
    elif len(market.get("selections") or []) != 2:                # judged on every listed outcome, not the live ones
        return None
    return candidates[0] if len(candidates) == 1 else None


def _canonical(market: dict, selection: dict) -> CanonicalLeg | None:
    category = market.get("category") or ""
    outcome = selection.get("outcomeType")
    points = selection.get("points")
    period, base = _period_and_base(_type_name(market))
    players = _players(selection)

    if category in SCORER_STATS:
        if len(players) != 1:
            return None
        scorer_period, stat = SCORER_STATS[category]
        return CanonicalLeg(period if period != "FG" else scorer_period, stat, players[0], "yes", None)

    if category == "HomeAway" and base == "Moneyline" and outcome in ("Home", "Away"):
        return CanonicalLeg(period, "moneyline", None, SIDE_OF_OUTCOME[outcome], None)
    if (category == "HomeAwayByPoints" and base in ("Spread", "Spread Alternate") and points is not None
            and outcome in ("Home", "Away")):
        return CanonicalLeg(period, "spread", None, SIDE_OF_OUTCOME[outcome], float(points))
    if category == "HomeTieAway" and base == "(3 Way)":
        subject = "Draw" if outcome == "Tie" else full_team_name(selection.get("name") or "")
        return CanonicalLeg(period, "win_3way", subject, "yes", None)

    if outcome in ("Over", "Under") and points is not None:
        side, line = SIDE_OF_OUTCOME[outcome], float(points)
        if base in ("Total", "Total Alternate"):
            return CanonicalLeg(period, "total", None, side, line)
        if base == "Team Total Points":
            teams = _team_names(market, selection)
            return CanonicalLeg(period, "team_total", full_team_name(teams[0]), side, line) if len(teams) == 1 else None
        if base.endswith(" O/U") and len(players) == 1:
            return CanonicalLeg(period, _slug(base[:-4]), players[0], side, line)
        if not players and not market.get("participants"):
            return CanonicalLeg(period, _slug(base), None, side, line)          # "Total Touchdowns"
        teams = _team_names(market, selection)
        if not players and len(teams) == 1:                                     # "Team Total Touchdowns"
            return CanonicalLeg(period, _slug(base), full_team_name(teams[0]), side, line)
        return None

    milestone = selection.get("milestoneValue")
    if base.endswith(" Milestones") and milestone is not None and len(players) == 1:
        # "70+" receiving yards = over 69.5.
        return CanonicalLeg(period, _slug(base[: -len(" Milestones")]), players[0], "over", float(milestone) - 0.5)

    if category == "YesNo" and outcome in ("Yes", "No") and "[" not in base:
        teams = _team_names(market, selection)
        subject = full_team_name(teams[0]) if len(teams) == 1 else None
        return CanonicalLeg(period, _slug(base), subject, SIDE_OF_OUTCOME[outcome], None)
    return None


def quote_from_added_leg_odds(added_leg_odds: dict[str, float], last: OfferedLeg, raw: dict) -> ParlayQuote | None:
    """DK's price of base + `last`, or None when DK will not combine them. DK
    quotes one price with no stated limit and no per-leg prices inside it."""
    decimal_odds = added_leg_odds.get(last.leg_id)
    if decimal_odds is None or decimal_odds <= 1.0:
        return None
    return ParlayQuote((ParlayTier(decimal_odds, None),), {}, raw)


# ---------------------------------------------------------------- the venue

class DraftKingsVenue:
    """Live DraftKings reads. Thread-safe. Caches the schedule (re-checked
    against the clock on every call), the leg picker's market list, and each
    builder answer until the next fresh leg read of its game."""

    name = "draftkings"
    display_name = "DraftKings"

    def __init__(self, proxy_url: str | None = None,
                 session_factory: Callable[[], Any] | None = None,
                 clock: Callable[[], float] = time.monotonic,
                 sleep: Callable[[float], None] = time.sleep,
                 now: Callable[[], datetime] = lambda: datetime.now(timezone.utc)) -> None:
        self._session_factory = session_factory or (lambda: _chrome_session(proxy_url))
        self._session = None
        self._clock = clock
        self._sleep = sleep
        self._now = now
        self._lock = threading.Lock()
        self._session_lock = threading.Lock()
        self._price_lock = threading.Lock()
        self._last_price_call_at = float("-inf")
        self._games: tuple[float, list[Game]] | None = None
        self._markets: dict[str, tuple[float, dict]] = {}
        self._added: dict[tuple[str, frozenset[str], str], tuple[float, dict[str, float], dict]] = {}

    def list_games(self) -> list[Game]:
        with self._lock:
            cached = self._games
        if not cached or self._clock() - cached[0] >= GAMES_CACHE_SEC:
            payload = self._get(DK_LEAGUE_URL, "league listing", params=league_listing_params(DK_LEAGUE_IDS["NFL"]))
            games = _parsed("league listing", parse_nfl_games, payload, self._now())
            cached = (self._clock(), games)
            with self._lock:
                self._games = cached
        now = self._now()
        return [game for game in cached[1] if game.game_start_time > now]

    def list_legs(self, game: Game, fresh: bool = False) -> list[OfferedLeg]:
        event_id = game.venue_game_id
        with self._lock:
            cached = self._markets.get(event_id)
        if cached and not fresh and self._clock() - cached[0] < MARKETS_CACHE_SEC:
            return _parsed("SGP leg list", parse_legs, cached[1])
        payload = self._get(DK_SGP_EVENT_URL.format(event_id=event_id), "SGP leg list",
                            read_timeout=LEG_LIST_TIMEOUT_SEC)
        legs = _parsed("SGP leg list", parse_legs, payload)
        with self._lock:
            now = self._clock()
            # One ~4.6 MB payload per game: keep only the ones still young enough to serve.
            self._markets = {k: v for k, v in self._markets.items() if now - v[0] < MARKETS_CACHE_SEC}
            self._markets[event_id] = (now, payload)
            if fresh:
                # A fresh read starts a Price click: its builder answers start fresh too.
                self._added = {k: v for k, v in self._added.items() if k[0] != event_id}
        return legs

    def price_parlay(self, game: Game, legs: list[OfferedLeg]) -> ParlayQuote | None:
        """The SGP price of all `legs`: the builder prices the last leg's market
        as an addition to the others."""
        if len(legs) < 2:
            raise ValueError(f"an SGP needs at least 2 legs, got {len(legs)}")
        base = frozenset(leg.leg_id for leg in legs[:-1])
        last = legs[-1]
        key = (game.venue_game_id, base, last.market_key)
        with self._lock:
            cached = self._added.get(key)
        if cached is None or self._clock() - cached[0] >= ADDED_LEG_CACHE_SEC:
            raw = self._price_call(*key)
            cached = (self._clock(), _parsed("SGP price", parse_added_leg_odds, raw, base), raw)
            with self._lock:
                now = self._clock()
                self._added = {k: v for k, v in self._added.items() if now - v[0] < ADDED_LEG_CACHE_SEC}
                self._added[key] = cached
        return quote_from_added_leg_odds(cached[1], last, cached[2])

    # -- network
    def _price_call(self, event_id: str, base: frozenset[str], candidate_market_id: str) -> dict:
        """One builder GET, never faster than DK's own widget (one per 500 ms)."""
        with self._price_lock:
            wait = MIN_SECONDS_BETWEEN_PRICE_CALLS - (self._clock() - self._last_price_call_at)
            if wait > 0:
                self._sleep(wait)
            try:
                return self._get(DK_SGP_PRICE_URL, "SGP price", headers={"X-SportId": DK_FOOTBALL_SPORT_ID},
                                 params=added_leg_params(event_id, base, candidate_market_id))
            finally:
                self._last_price_call_at = self._clock()

    def _warm_session(self):
        """The shared session, made and warmed once (DK sets its cookies on the page)."""
        with self._session_lock:
            if self._session is None:
                session = self._session_factory()
                try:
                    session.get(DK_WARMUP_URL, timeout=(DK_CONNECT_TIMEOUT_SECONDS, REQUEST_TIMEOUT_SEC))
                except Exception as exc:   # curl_cffi raises its own error types
                    raise VenueTransportError(f"DraftKings warmup: {type(exc).__name__}: {exc}") from exc
                self._session = session
            return self._session

    def _drop_session(self, session) -> None:
        """Forget and close `session` (not one another thread already replaced it with)."""
        with self._session_lock:
            if self._session is session:
                self._session = None
        try:
            session.close()
        except Exception:   # closing a broken curl handle must not mask the error being raised
            pass

    def _get(self, url: str, what: str, read_timeout: float = REQUEST_TIMEOUT_SEC, **kwargs) -> dict:
        session = self._warm_session()
        try:
            response = session.get(url, timeout=(DK_CONNECT_TIMEOUT_SECONDS, read_timeout), **kwargs)
        except Exception as exc:
            self._drop_session(session)
            raise VenueTransportError(f"DraftKings {what}: {type(exc).__name__}: {exc}") from exc
        if response.status_code != 200:
            if response.status_code == HTTP_FORBIDDEN:
                # Akamai may have flagged this cookie jar: start the next call on a new one.
                self._drop_session(session)
            raise VenueTransportError(f"DraftKings {what}: expected HTTP 200, got "
                                      f"{response.status_code}: {response.text[:200]}")
        try:
            body = response.json()
        except ValueError as exc:
            raise VenueTransportError(f"DraftKings {what}: response is not JSON") from exc
        if not isinstance(body, dict):
            raise VenueTransportError(f"DraftKings {what}: expected a JSON object, got {type(body).__name__}")
        return body


def _parsed(what: str, parse: Callable[..., Any], *args) -> Any:
    """parse(*args), with a shape DK changed reported as this venue's error."""
    try:
        return parse(*args)
    except (KeyError, IndexError, TypeError, ValueError, AttributeError) as exc:
        raise VenueTransportError(f"DraftKings {what}: unexpected shape: {type(exc).__name__}: {exc}") from exc


def _chrome_session(proxy_url: str | None):
    from curl_cffi import requests as curl_requests
    return curl_requests.Session(impersonate="chrome", proxy=proxy_url)
