"""FanDuel adapter (books.SgpBook) for NFL fecta legs and NFL / college 1H
spread + total legs. Plain HTTP, no browser.

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
    total          "1st Half Total" at that line (Over, Under); FD posts no
                   1H alternate totals, so another line is not offered
    scores_first   not offered: "Team to Score First" is not SGP-eligible,
                   so FD cannot price a superfecta.
"""
from __future__ import annotations

import json
import re
import sys
from typing import Hashable

from curl_cffi import requests as cffi_requests

from nfl_specials.books import (NFL, BookGame, LegMarket, Outcome, RequestPacer, TeamKey, is_half_point,
                                next_game, widen_with_tie)
from nfl_specials.config import REPO_ROOT
from nfl_specials.special_parser import Leg
from nfl_specials.teams import ODDS_API_NAME_TO_ABBR

# scraper_fanduel_sgp imports its siblings by bare name.
sys.path.insert(0, str(REPO_ROOT / "mlb_sgp"))
from scraper_fanduel_sgp import (FD_AK, FD_EVENT_PAGE_URL, FD_HEADERS,  # noqa: E402
                                 FD_IMPLY_BETS_URL, FD_SCAN_URL)

# College: FBS and FCS are separate competitions.
FD_COMPETITION_IDS = {"NFL": (12282733,), "CFB": (12529073, 12623176)}
FD_SGP_TAB = "same-game-parlay-"
SPREAD_MAIN = {"GM": "Spread", "Q1": "1st Quarter Spread", "H1": "1st Half Spread"}
SPREAD_ALT = {"GM": "Alternate Spread", "H1": "1st Half Alternate Spread"}
THREE_WAY_WIN = {"Q1": "1st Quarter Winner (3-Way)", "H1": "1st Half Winner (3-Way)"}
TOTAL_MAIN = {"H1": "1st Half Total"}
ALT_RUNNER_RE = re.compile(r"^(.+?) \(([+-]\d+(?:\.\d+)?)\)$")
# FanDuel's WAF blocks an address for ~5 min after 100 implyBets in a
# trailing 5 min (VM probe 2026-11-13: 163 calls at 1.5/s, the 164th 403,
# unblocked as the window fell under 100; a burst of 234 in 40 s got through
# first only because the WAF counts with a lag). So FanDuel prices to a
# window budget, with margin; its scan and event pages are other hosts and do
# not count. The 1 s gap keeps a full budget under the probe's 1.5/s and
# costs nothing (90 s of sends in a 300 s window).
MIN_SECONDS_BETWEEN_PRICE_CALLS = 1.0
FD_PRICE_CALLS_PER_WINDOW = 90
FD_PRICE_WINDOW_SECONDS = 300.0


class FanDuelBook:
    name = "FanDuel"

    def __init__(self, league: str = NFL, pacer: RequestPacer | None = None) -> None:
        """`pacer`: shared by every FanDuelBook of one refresh; None = its own."""
        self.league = league
        self.pacer = pacer if pacer is not None else RequestPacer(self.name, MIN_SECONDS_BETWEEN_PRICE_CALLS)
        self.session = cffi_requests.Session(impersonate="chrome")
        self._games = self._list_games()
        self._markets: dict[str, dict[str, dict]] = {}

    # --- discovery ---------------------------------------------------------
    def _list_games(self) -> list[BookGame]:
        body = {"filter": {"competitionIds": list(FD_COMPETITION_IDS[self.league]),
                           "contentGroup": {"language": "en", "regionCode": "NAMERICA"},
                           "marketLevels": ["AVB_EVENT"], "maxResults": 100,
                           "productTypes": ["SPORTSBOOK"], "selectBy": "FIRST_TO_START"},
                "facets": [{"type": "EVENT"}], "currencyCode": "USD"}
        resp = self.pacer.guard(lambda: self.session.post(
            FD_SCAN_URL, headers={**FD_HEADERS, "Content-Type": "application/json"},
            data=json.dumps(body), timeout=30), "event scan")
        if resp.status_code != 200:
            raise RuntimeError(f"FanDuel event scan returned HTTP {resp.status_code}")
        games = []
        for event in _walk_dicts(resp.json(), lambda d: "eventId" in d and "openDate" in d and "name" in d):
            away_name, sep, home_name = event["name"].partition(" @ ")
            home, away = self._team_key(home_name), self._team_key(away_name)
            if not sep or home is None or away is None:
                continue
            games.append(BookGame(book_event_id=str(event["eventId"]), home=home, away=away,
                                  game_start_time=event["openDate"]))
        return games

    def _team_key(self, fd_name: str) -> str | None:
        """NFL: 'Seattle Seahawks' -> 'SEA'. College: FD's name as is."""
        return ODDS_API_NAME_TO_ABBR.get(fd_name) if self.league == NFL else fd_name

    def games(self) -> list[BookGame]:
        return list(self._games)

    def find_game(self, team: str) -> BookGame | None:
        return next_game(self._games, team)

    def _event_markets(self, game: BookGame) -> dict[str, dict]:
        if game.book_event_id not in self._markets:
            url = (f"{FD_EVENT_PAGE_URL}?_ak={FD_AK}&eventId={game.book_event_id}"
                   f"&tab={FD_SGP_TAB}&useCombinedTouchdownsVirtualMarket=true&useQuickBets=true")
            resp = self.pacer.guard(lambda: self.session.get(url, headers=FD_HEADERS, timeout=30), "event page")
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
        key = self._team_key
        if leg.kind == "scores_first":
            return None
        if leg.kind == "total":
            return _total(markets.get(TOTAL_MAIN.get(leg.period, "")), leg.line)
        if leg.kind == "win" and leg.period == "GM":
            return _team_outcome_market(markets.get("Moneyline"), team, key)
        if leg.kind == "win":
            return _team_outcome_market(markets.get(THREE_WAY_WIN[leg.period]), team, key)
        if leg.kind == "spread":
            found = _main_spread(markets.get(SPREAD_MAIN[leg.period]), team, leg.line, key)
            if found is None and leg.period in SPREAD_ALT:
                found = _alt_spread(markets.get(SPREAD_ALT[leg.period]), team, leg.line, key)
            if found is None and leg.period in THREE_WAY_WIN and is_half_point(leg.line):
                # FD has no 1Q alternate spread; +/-0.5 is the 3-way winner.
                three_way = _team_outcome_market(markets.get(THREE_WAY_WIN[leg.period]), team, key)
                found = three_way if leg.line < 0 else widen_with_tie(three_way)
            return found
        raise ValueError(f"unknown leg kind {leg.kind!r}")

    # --- pricing -----------------------------------------------------------
    def price(self, game: BookGame, refs: tuple[Hashable, ...]) -> float | None:
        body = {"betLegs": [{"legType": "SIMPLE_SELECTION",
                             "betRunners": [{"runner": {"marketId": market_id, "selectionId": selection_id}}]}
                            for market_id, selection_id in refs]}
        resp = self.pacer.call(lambda: self.session.post(
            FD_IMPLY_BETS_URL, headers={**FD_HEADERS, "Content-Type": "application/json"},
            data=json.dumps(body), timeout=20), "implyBets")
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


def _team_outcome_market(market: dict | None, team: str, team_key: TeamKey) -> LegMarket | None:
    runners = _active(market)
    group = tuple(Outcome(ref=(market["marketId"], r["selectionId"]), label=r["runnerName"]) for r in runners)
    chosen = next((o for o, r in zip(group, runners) if team_key(r["runnerName"]) == team), None)
    if chosen is None or len(group) < 2:
        return None
    return LegMarket.single(chosen, group)


def _main_spread(market: dict | None, team: str, line: float, team_key: TeamKey) -> LegMarket | None:
    runners = _active(market)
    if len(runners) != 2:
        return None
    mine = next((r for r in runners if team_key(r["runnerName"]) == team), None)
    other = next((r for r in runners if r is not mine), None)
    if mine is None or abs(float(mine.get("handicap", 0)) - line) > 1e-6:
        return None
    pick = Outcome(ref=(market["marketId"], mine["selectionId"]), label=f"{mine['runnerName']} {line:+g}")
    rest = Outcome(ref=(market["marketId"], other["selectionId"]), label=f"{other['runnerName']} {-line:+g}")
    return LegMarket.single(pick, (pick, rest))


def _alt_spread(market: dict | None, team: str, line: float, team_key: TeamKey) -> LegMarket | None:
    """Alternate ladders name the line in the runner: 'Seattle Seahawks (-7.5)'."""
    pick = other = None
    for runner in _active(market):
        match = ALT_RUNNER_RE.match(runner["runnerName"])
        if not match:
            continue
        runner_team, runner_line = team_key(match.group(1)), float(match.group(2))
        outcome = Outcome(ref=(market["marketId"], runner["selectionId"]), label=runner["runnerName"])
        if runner_team == team and abs(runner_line - line) < 1e-6:
            pick = outcome
        elif runner_team not in (team, None) and abs(runner_line + line) < 1e-6:
            other = outcome
    if pick is None or other is None:
        return None
    return LegMarket.single(pick, (pick, other))


def _total(market: dict | None, line: float) -> LegMarket | None:
    """A main total market at exactly `line`: the leg is the Over."""
    runners = _active(market)
    over = next((r for r in runners if r["runnerName"] == "Over"), None)
    under = next((r for r in runners if r["runnerName"] == "Under"), None)
    if over is None or under is None or len(runners) != 2:
        return None
    if abs(float(over.get("handicap", 0)) - line) > 1e-6 or abs(float(under.get("handicap", 0)) - line) > 1e-6:
        return None
    pick = Outcome(ref=(market["marketId"], over["selectionId"]), label=f"Over {line:g}")
    other = Outcome(ref=(market["marketId"], under["selectionId"]), label=f"Under {line:g}")
    return LegMarket.single(pick, (pick, other))
