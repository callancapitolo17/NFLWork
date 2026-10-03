"""FanDuel adapter (books.SgpBook) for NFL fecta legs. Plain HTTP, no browser.

Reads:  scan (event list), event-page tab "same-game-parlay-" (every market a
        trifecta needs, one GET per game), implyBets (SGP price).
No DB writes. Headers / URLs / the PerimeterX context come from the MLB
scraper so a rotated token is fixed in one place.

Leg -> FD market (only sgmMarket=True markets are used):
    win GM         "Moneyline"                                  (home, away)
    win Q1 / H1    "1st Quarter Winner (3-Way)" / "1st Half Winner (3-Way)"
                   (home, Tie, away). The 2-way "Winner" markets push a tie
                   and FD refuses them in an SGP.
    spread         main period spread at that line, else the alternate ladder
    scores_first   not offered: "Team to Score First" is not SGP-eligible,
                   so FD cannot price a superfecta.
"""
from __future__ import annotations

import json
import re
import sys
from typing import Hashable

from curl_cffi import requests as cffi_requests

from nfl_specials.books import BookGame, LegMarket, Outcome, next_game, is_half_point, widen_with_tie
from nfl_specials.config import REPO_ROOT
from nfl_specials.special_parser import Leg
from nfl_specials.teams import ODDS_API_NAME_TO_ABBR

# scraper_fanduel_sgp imports its siblings by bare name.
sys.path.insert(0, str(REPO_ROOT / "mlb_sgp"))
from scraper_fanduel_sgp import (FD_AK, FD_EVENT_PAGE_URL, FD_HEADERS,  # noqa: E402
                                 FD_IMPLY_BETS_URL, FD_SCAN_URL)

FD_NFL_COMPETITION_ID = 12282733
FD_SGP_TAB = "same-game-parlay-"
SPREAD_MAIN = {"GM": "Spread", "Q1": "1st Quarter Spread", "H1": "1st Half Spread"}
SPREAD_ALT = {"GM": "Alternate Spread", "H1": "1st Half Alternate Spread"}
THREE_WAY_WIN = {"Q1": "1st Quarter Winner (3-Way)", "H1": "1st Half Winner (3-Way)"}
ALT_RUNNER_RE = re.compile(r"^(.+?) \(([+-]\d+(?:\.\d+)?)\)$")


class FanDuelBook:
    name = "FanDuel"

    def __init__(self) -> None:
        self.session = cffi_requests.Session(impersonate="chrome")
        self._games = self._list_games()
        self._markets: dict[str, dict[str, dict]] = {}

    # --- discovery ---------------------------------------------------------
    def _list_games(self) -> list[BookGame]:
        body = {"filter": {"competitionIds": [FD_NFL_COMPETITION_ID],
                           "contentGroup": {"language": "en", "regionCode": "NAMERICA"},
                           "marketLevels": ["AVB_EVENT"], "maxResults": 100,
                           "productTypes": ["SPORTSBOOK"], "selectBy": "FIRST_TO_START"},
                "facets": [{"type": "EVENT"}], "currencyCode": "USD"}
        resp = self.session.post(FD_SCAN_URL, headers={**FD_HEADERS, "Content-Type": "application/json"},
                                 data=json.dumps(body), timeout=30)
        if resp.status_code != 200:
            raise RuntimeError(f"FanDuel event scan returned HTTP {resp.status_code}")
        games = []
        for event in _walk_dicts(resp.json(), lambda d: "eventId" in d and "openDate" in d and "name" in d):
            away_name, sep, home_name = event["name"].partition(" @ ")
            home, away = ODDS_API_NAME_TO_ABBR.get(home_name), ODDS_API_NAME_TO_ABBR.get(away_name)
            if not sep or home is None or away is None:
                continue
            games.append(BookGame(book_event_id=str(event["eventId"]), home=home, away=away,
                                  game_start_time=event["openDate"]))
        return games

    def find_game(self, team: str) -> BookGame | None:
        return next_game(self._games, team)

    def _event_markets(self, game: BookGame) -> dict[str, dict]:
        if game.book_event_id not in self._markets:
            url = (f"{FD_EVENT_PAGE_URL}?_ak={FD_AK}&eventId={game.book_event_id}"
                   f"&tab={FD_SGP_TAB}&useCombinedTouchdownsVirtualMarket=true&useQuickBets=true")
            resp = self.session.get(url, headers=FD_HEADERS, timeout=30)
            if resp.status_code != 200:
                raise RuntimeError(f"FanDuel event page {game.book_event_id} returned HTTP {resp.status_code}")
            markets = {}
            for market in _walk_dicts(resp.json(), lambda d: "marketId" in d and "runners" in d and "marketName" in d):
                if market.get("sgmMarket") and market.get("marketStatus") == "OPEN":
                    markets.setdefault(market["marketName"], market)
            self._markets[game.book_event_id] = markets
        return self._markets[game.book_event_id]

    # --- legs --------------------------------------------------------------
    def leg_market(self, game: BookGame, role: str, leg: Leg) -> LegMarket | None:
        markets = self._event_markets(game)
        team = game.home if role == "home" else game.away
        if leg.kind == "scores_first":
            return None
        if leg.kind == "win" and leg.period == "GM":
            return _team_outcome_market(markets.get("Moneyline"), team)
        if leg.kind == "win":
            return _team_outcome_market(markets.get(THREE_WAY_WIN[leg.period]), team)
        if leg.kind == "spread":
            found = _main_spread(markets.get(SPREAD_MAIN[leg.period]), team, leg.line)
            if found is None and leg.period in SPREAD_ALT:
                found = _alt_spread(markets.get(SPREAD_ALT[leg.period]), team, leg.line)
            if found is None and leg.period in THREE_WAY_WIN and is_half_point(leg.line):
                # FD has no 1Q alternate spread; +/-0.5 is the 3-way winner.
                three_way = _team_outcome_market(markets.get(THREE_WAY_WIN[leg.period]), team)
                found = three_way if leg.line < 0 else widen_with_tie(three_way)
            return found
        raise ValueError(f"unknown leg kind {leg.kind!r}")

    # --- pricing -----------------------------------------------------------
    def price(self, game: BookGame, refs: tuple[Hashable, ...]) -> float | None:
        body = {"betLegs": [{"legType": "SIMPLE_SELECTION",
                             "betRunners": [{"runner": {"marketId": market_id, "selectionId": selection_id}}]}
                            for market_id, selection_id in refs]}
        resp = self.session.post(FD_IMPLY_BETS_URL, headers={**FD_HEADERS, "Content-Type": "application/json"},
                                 data=json.dumps(body), timeout=20)
        if resp.status_code != 200:
            raise RuntimeError(f"FanDuel implyBets returned HTTP {resp.status_code}: {resp.text[:200]}")
        data = resp.json()
        for combination in data.get("betCombinations", []):
            if combination.get("isSGM"):
                return float(combination["winAvgOdds"]["trueOdds"]["decimalOdds"]["decimalOdds"])
        # INVALID_SGM_COMBINATION = FD will not build this parlay (a decline);
        # INVALID_COMBINATION rides on every same-game request and means nothing.
        return None


def _walk_dicts(node, keep) -> list[dict]:
    found, stack = [], [node]
    while stack:
        item = stack.pop()
        if isinstance(item, dict):
            if keep(item):
                found.append(item)
            stack.extend(item.values())
        elif isinstance(item, list):
            stack.extend(item)
    return found


def _active(market: dict | None) -> list[dict]:
    if market is None:
        return []
    return [r for r in market.get("runners", []) if r.get("runnerStatus") == "ACTIVE"]


def _team_outcome_market(market: dict | None, team: str) -> LegMarket | None:
    runners = _active(market)
    group = tuple(Outcome(ref=(market["marketId"], r["selectionId"]), label=r["runnerName"]) for r in runners)
    chosen = next((o for o, r in zip(group, runners)
                   if ODDS_API_NAME_TO_ABBR.get(r["runnerName"]) == team), None)
    if chosen is None or len(group) < 2:
        return None
    return LegMarket.single(chosen, group)


def _main_spread(market: dict | None, team: str, line: float) -> LegMarket | None:
    runners = _active(market)
    if len(runners) != 2:
        return None
    mine = next((r for r in runners if ODDS_API_NAME_TO_ABBR.get(r["runnerName"]) == team), None)
    other = next((r for r in runners if r is not mine), None)
    if mine is None or abs(float(mine.get("handicap", 0)) - line) > 1e-6:
        return None
    pick = Outcome(ref=(market["marketId"], mine["selectionId"]), label=f"{mine['runnerName']} {line:+g}")
    rest = Outcome(ref=(market["marketId"], other["selectionId"]), label=f"{other['runnerName']} {-line:+g}")
    return LegMarket.single(pick, (pick, rest))


def _alt_spread(market: dict | None, team: str, line: float) -> LegMarket | None:
    """Alternate ladders name the line in the runner: 'Seattle Seahawks (-7.5)'."""
    pick = other = None
    for runner in _active(market):
        match = ALT_RUNNER_RE.match(runner["runnerName"])
        if not match:
            continue
        runner_team, runner_line = ODDS_API_NAME_TO_ABBR.get(match.group(1)), float(match.group(2))
        outcome = Outcome(ref=(market["marketId"], runner["selectionId"]), label=runner["runnerName"])
        if runner_team == team and abs(runner_line - line) < 1e-6:
            pick = outcome
        elif runner_team not in (team, None) and abs(runner_line + line) < 1e-6:
            other = outcome
    if pick is None or other is None:
        return None
    return LegMarket.single(pick, (pick, other))
