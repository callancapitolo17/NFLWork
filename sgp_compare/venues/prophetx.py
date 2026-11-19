"""ProphetX venue adapter: NFL games, every parlay-eligible leg (game lines and
player props), and parlay prices.

Network: anonymous public web endpoints on www.prophetx.co over a curl_cffi
Chrome-TLS session (the same session mlb_sgp uses; the 2026-11-18 recon on Cal's
Mac read all three without a login):
    GET  /trade/public/api/v1/tournaments?expand=events&type=highlight&limit=150
    GET  /trade/public/api/v2/events/{event_id}/markets
    POST /parlay/public/api/v1/user/request    {marketLines: [...], stake: 1}
The parlay POST is a price request: it returns offers and never places a bet
(betting needs a separate accept on a logged-in account, which this module has
no code for). No DB I/O.

ProphetX shapes (recon 2026-11-18, IND at HOU, fixtures in tests/fixtures/):
  * NFL tournament id 31, name "NFL"; "NFL Win Totals" is a different tournament.
  * competitor seq 0 is HOME.
  * a spread/total market keeps its numbers in `marketLines`, each with two
    `outcomes` [{id, name, line, lineID, isMain}] and `selections` = one
    best-first price ladder per side ([] or null when nobody offers it).
  * a moneyline or yes/no market has top-level `outcomes` and `selections`.
  * outcome ids name the side: spread home 1714 / away 1715, total over 12 /
    under 13, moneyline home 4 / away 5, yes 1 / no 2. `lineID` is unique per
    side per number, so it is the leg id.
  * ladder entries: `odds` (American, the price to take that side), `stake`
    (dollars offered at that price).
  * parlay answer: {success, data: {offers: [{odds (American, fractional),
    stake (max dollars at that price), estimatedPrices: [{lineId, odds}]}]}};
    ProphetX folds its correlation charge into each leg's estimatedPrices.
"""
from __future__ import annotations

import logging
import re
import threading
import time
from datetime import datetime, timezone
from typing import Any, Callable

from sgp_compare.legs import CanonicalLeg, Game, OfferedLeg, ParlayQuote, ParlayTier
from sgp_compare.odds import american_to_decimal
from sgp_compare.venues.base import VenueTransportError

log = logging.getLogger(__name__)

PROPHETX_BASE_URL = "https://www.prophetx.co"
NFL_TOURNAMENT_NAME = "NFL"
REQUEST_TIMEOUT_SEC = 15
GAMES_CACHE_SEC = 300          # the schedule changes daily
MARKETS_CACHE_SEC = 20         # the leg picker; pricing always re-reads
PARLAY_PRICE_REQUEST_STAKE = 1

SIDE_OF_OUTCOME_ID = {
    1714: "home", 1715: "away",      # spread
    4: "home", 5: "away",            # moneyline
    12: "over", 13: "under",         # total
    1: "yes", 2: "no",               # yes/no props
}

# subType prefix -> period. Checked in order; no prefix = full game.
PERIOD_PREFIXES = (
    ("first_half_", "1H"), ("second_half_", "2H"),
    ("1st_quarter_", "1Q"), ("2nd_quarter_", "2Q"), ("3rd_quarter_", "3Q"), ("4th_quarter_", "4Q"),
)

# Player prop market names are "<Player> <words>"; the words depend on the
# subType. Unknown subTypes keep the whole market name as the subject.
PLAYER_NAME_SUFFIXES = {
    "to_score_a_touchdown": " To Score a Touchdown",
    "player_to_score_first_touchdown": " To Score First Touchdown",
    "player_to_score_last_touchdown": " To Score Last Touchdown",
    "player_to_score_at_least_2_touchdowns": " To Score 2+ Touchdowns",
    "player_to_score_at_least_3_touchdowns": " To Score 3+ Touchdowns",
}


# ---------------------------------------------------------------- parsing (pure)

def parse_nfl_games(tournaments_payload: dict, now: datetime) -> list[Game]:
    """Not-yet-started NFL games from the tournaments response, soonest first."""
    tournaments = (tournaments_payload.get("data") or {}).get("tournaments") or []
    games = []
    for tournament in tournaments:
        if tournament.get("name") != NFL_TOURNAMENT_NAME:
            continue
        for event in tournament.get("sportEvents") or []:
            game = _game_from_event(event)
            if game is not None and event.get("status") == "not_started" and game.game_start_time > now:
                games.append(game)
    return sorted(games, key=lambda g: g.game_start_time)


def _game_from_event(event: dict) -> Game | None:
    competitors = {c.get("seq"): c for c in event.get("competitors") or []}
    home, away = competitors.get(0), competitors.get(1)
    if home is None or away is None or not event.get("scheduled"):
        return None
    start = datetime.fromisoformat(event["scheduled"].replace("Z", "+00:00"))
    return Game(venue_game_id=str(event["id"]), home_team=home["name"], away_team=away["name"],
                game_start_time=start.astimezone(timezone.utc))


def parse_legs(game: Game, markets_payload: dict, home_competitor_id: int | None) -> list[OfferedLeg]:
    """Every priced leg of the game's active markets.

    A market line where neither side has a price is left out (ProphetX lists
    deep alt ladders whose tails are empty). `home_competitor_id` decides the
    side of moneyline outcomes that carry only a competitor id.
    """
    markets = (markets_payload.get("data") or {}).get("markets") or []
    legs: list[OfferedLeg] = []
    for market in markets:
        if market.get("status") != "active":
            continue
        if market.get("marketLines"):
            for market_line in market["marketLines"]:
                legs.extend(_legs_of_line(game, market, market_line, home_competitor_id))
        elif market.get("outcomes"):
            legs.extend(_legs_of_line(game, market, market, home_competitor_id))
    return legs


def _legs_of_line(game: Game, market: dict, holder: dict, home_competitor_id: int | None) -> list[OfferedLeg]:
    """Both sides of one two-outcome line. `holder` is the marketLine, or the
    market itself for moneyline / yes-no markets."""
    outcomes = holder.get("outcomes") or []
    if len(outcomes) != 2:
        return []
    best = _best_price_by_line_id(holder.get("selections") or [])
    if not any(o.get("lineID") in best for o in outcomes):
        return []
    legs = []
    for outcome, other in ((outcomes[0], outcomes[1]), (outcomes[1], outcomes[0])):
        price, stake = best.get(outcome.get("lineID"), (None, None))
        other_price, _ = best.get(other.get("lineID"), (None, None))
        legs.append(OfferedLeg(
            leg_id=outcome["lineID"],
            game_id=game.venue_game_id,
            market_label=market.get("name") or "",
            category=market.get("categoryName") or "",
            side_label=_side_label(game, outcome, home_competitor_id),
            canonical=_canonical(game, market, outcome, home_competitor_id),
            american_odds=price,
            opposite_american_odds=other_price,
            available_stake=stake,
            opposite_leg_id=other.get("lineID"),
            market_key=str(market.get("id")),
            is_main_line=bool(outcome.get("isMain")),
            venue_ref={"sportEventId": int(game.venue_game_id), "marketId": market.get("id"),
                       "outcomeId": outcome.get("id"), "lineId": outcome.get("lineID"),
                       "line": outcome.get("line", 0)},
        ))
    return legs


def _best_price_by_line_id(selections: list) -> dict[str, tuple[float, float | None]]:
    """lineID -> (best American price, dollars at it). A ladder is best first."""
    best: dict[str, tuple[float, float | None]] = {}
    for ladder in selections:
        if not ladder:
            continue
        top = ladder[0]
        odds = top.get("odds")
        if top.get("lineID") and isinstance(odds, (int, float)) and abs(odds) >= 100:
            best[top["lineID"]] = (float(odds), top.get("stake"))
    return best


def _side_of(outcome: dict, home_competitor_id: int | None) -> str | None:
    side = SIDE_OF_OUTCOME_ID.get(outcome.get("id"))
    if side is None and outcome.get("competitorId") is not None and home_competitor_id is not None:
        side = "home" if outcome["competitorId"] == home_competitor_id else "away"
    return side


def _side_label(game: Game, outcome: dict, home_competitor_id: int | None) -> str:
    name = outcome.get("name")
    if name:
        return name[:1].upper() + name[1:]
    side = _side_of(outcome, home_competitor_id)
    return {"home": game.home_team, "away": game.away_team}.get(side, side or "?")


def _period_and_stat(sub_type: str) -> tuple[str, str]:
    for prefix, period in PERIOD_PREFIXES:
        if sub_type.startswith(prefix):
            return period, sub_type[len(prefix):]
    return "FG", sub_type


def _canonical(game: Game, market: dict, outcome: dict, home_competitor_id: int | None) -> CanonicalLeg | None:
    sub_type = market.get("subType") or ""
    side = _side_of(outcome, home_competitor_id)
    if not sub_type or side is None:
        return None
    period, stat = _period_and_stat(sub_type)
    name = market.get("name") or ""
    line = outcome.get("line")
    line = float(line) if isinstance(line, (int, float)) and stat not in ("moneyline",) else None
    if stat in ("moneyline", "spread", "total"):
        return CanonicalLeg(period, stat, None, side, line)
    if stat.endswith("team_total_points"):
        return CanonicalLeg(period, "team_total", name.split(":")[0].strip(), side, line)
    if stat.startswith("moneyline_3_way"):
        # "Houston Texans To Win (2nd Quarter)" / "Draw (2nd Quarter)", each a yes/no market.
        subject = re.sub(r"\s*(To Win)?\s*\(.*\)$", "", name).strip()
        return CanonicalLeg(period, "win_3way", subject, side, None)
    if sub_type in PLAYER_NAME_SUFFIXES or sub_type.startswith("player_"):
        suffix = PLAYER_NAME_SUFFIXES.get(sub_type)
        subject = name[: -len(suffix)] if suffix and name.endswith(suffix) else _strip_category(name, market)
        player_stat = re.sub(r"^(player_)?(total_)?", "", sub_type)
        return CanonicalLeg(period, player_stat, subject, side, line if side in ("over", "under") else None)
    return CanonicalLeg(period, stat, name, side, line if side in ("over", "under") else None)


def _strip_category(name: str, market: dict) -> str:
    """'Nico Collins Receiving Yards' -> 'Nico Collins' when the stat words
    match the subType ('player_total_receiving_yards')."""
    words = re.sub(r"^(player_)?(total_)?", "", market.get("subType") or "").split("_")
    tail = " " + " ".join(w.capitalize() for w in words)
    return name[: -len(tail)] if name.lower().endswith(tail.lower()) else name


def parlay_payload(legs: list[OfferedLeg]) -> dict:
    return {"marketLines": [dict(leg.venue_ref) for leg in legs], "stake": PARLAY_PRICE_REQUEST_STAKE}


def parse_parlay_quote(response: dict, legs: list[OfferedLeg]) -> ParlayQuote | None:
    """None when ProphetX declined (success false, or no offers)."""
    data = response.get("data") or {}
    offers = data.get("offers") or []
    if not response.get("success") or not offers:
        return None
    tiers = []
    for offer in offers:
        odds = offer.get("odds")
        if isinstance(odds, (int, float)) and abs(odds) >= 100:
            tiers.append(ParlayTier(american_to_decimal(float(odds)), offer.get("stake")))
    if not tiers:
        return None
    tiers.sort(key=lambda tier: tier.decimal_odds, reverse=True)
    leg_id_of_line = {leg.venue_ref["lineId"]: leg.leg_id for leg in legs}
    best_offer = max(offers, key=lambda o: o.get("odds") or float("-inf"))
    leg_prices = {leg_id_of_line[p["lineId"]]: float(p["odds"])
                  for p in best_offer.get("estimatedPrices") or [] if p.get("lineId") in leg_id_of_line}
    return ParlayQuote(tuple(tiers), leg_prices, response)


# ---------------------------------------------------------------- the venue

class ProphetXVenue:
    """Live ProphetX reads. Thread-safe; caches only the schedule and the leg
    picker's market list, never a price used to compute a parlay."""

    name = "prophetx"
    display_name = "ProphetX"

    def __init__(self, session_factory: Callable[[], Any] | None = None,
                 clock: Callable[[], float] = time.monotonic) -> None:
        self._session_factory = session_factory or _chrome_session
        self._session = None
        self._clock = clock
        self._lock = threading.Lock()
        self._games: tuple[float, list[Game], dict[str, int]] | None = None
        self._markets: dict[str, tuple[float, dict]] = {}

    def list_games(self) -> list[Game]:
        return self._read_games()[0]

    def list_legs(self, game: Game, fresh: bool = False) -> list[OfferedLeg]:
        payload = self._read_markets(game.venue_game_id, fresh)
        home_id = self._read_games()[1].get(game.venue_game_id)
        return parse_legs(game, payload, home_id)

    def price_parlay(self, game: Game, legs: list[OfferedLeg]) -> ParlayQuote | None:
        response = self._request("POST", "/parlay/public/api/v1/user/request", json=parlay_payload(legs))
        return parse_parlay_quote(response, legs)

    # -- reads
    def _read_games(self) -> tuple[list[Game], dict[str, int]]:
        with self._lock:
            cached = self._games
        if cached and self._clock() - cached[0] < GAMES_CACHE_SEC:
            return cached[1], cached[2]
        payload = self._request("GET", "/trade/public/api/v1/tournaments",
                                params={"expand": "events", "type": "highlight", "limit": 150})
        games = parse_nfl_games(payload, datetime.now(timezone.utc))
        home_ids = _home_competitor_ids(payload)
        with self._lock:
            self._games = (self._clock(), games, home_ids)
        return games, home_ids

    def _read_markets(self, event_id: str, fresh: bool) -> dict:
        with self._lock:
            cached = self._markets.get(event_id)
        if cached and not fresh and self._clock() - cached[0] < MARKETS_CACHE_SEC:
            return cached[1]
        payload = self._request("GET", f"/trade/public/api/v2/events/{event_id}/markets")
        with self._lock:
            self._markets[event_id] = (self._clock(), payload)
        return payload

    def _request(self, method: str, path: str, **kwargs) -> dict:
        if self._session is None:
            self._session = self._session_factory()
        url = PROPHETX_BASE_URL + path
        try:
            response = self._session.request(method, url, timeout=REQUEST_TIMEOUT_SEC, **kwargs)
        except Exception as exc:   # curl_cffi raises its own error types
            self._session = None
            raise VenueTransportError(f"ProphetX {method} {path}: {type(exc).__name__}: {exc}") from exc
        if response.status_code != 200:
            raise VenueTransportError(f"ProphetX {method} {path}: expected HTTP 200, got "
                                      f"{response.status_code}: {response.text[:200]}")
        try:
            return response.json()
        except ValueError as exc:
            raise VenueTransportError(f"ProphetX {method} {path}: response is not JSON") from exc


def _home_competitor_ids(tournaments_payload: dict) -> dict[str, int]:
    home_ids = {}
    for tournament in (tournaments_payload.get("data") or {}).get("tournaments") or []:
        for event in tournament.get("sportEvents") or []:
            for competitor in event.get("competitors") or []:
                if competitor.get("seq") == 0:
                    home_ids[str(event.get("id"))] = competitor.get("id")
    return home_ids


def _chrome_session():
    from curl_cffi import requests as curl_requests
    return curl_requests.Session(impersonate="chrome")
