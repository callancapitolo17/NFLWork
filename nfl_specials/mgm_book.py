"""BetMGM adapter (books.SgpBook) for NFL fecta legs. Plain HTTP, no browser.

Reuses mlb_sgp.betmgm_client.BetMGMClient for the session and the accessid
harvest (Entain CDS). NFL differs from MLB in three ways: sport id 11 with
competition "NFL", fixture names "<Away> @ <Home>", and V2 fixtures that
price through `tv2Picks` (the MLB `tv1Picks` body returns an empty group).
No DB writes.

Leg -> BetMGM market (only isBetBuilder markets):
    win GM         "Moneyline"                                  (home, away)
    win Q1 / H1    "1st quarter spread" / "1st half spread" at -0.5 for the
                   team (+0.5 opponent). The period moneylines are 2-way with a
                   tie pushed, and no 3-way period market can go in an SGP —
                   with integer scores, -0.5 IS "wins the period, a tie loses".
    spread         the period's spread market at that exact line
    scores_first   not offered, so BetMGM cannot price a superfecta.
"""
from __future__ import annotations

import re
import uuid
from typing import Hashable

from mlb_sgp.betmgm_client import BetMGMClient
from nfl_specials.books import BookGame, LegMarket, Outcome, next_game
from nfl_specials.special_parser import Leg
from nfl_specials.teams import ODDS_API_NAME_TO_ABBR

MGM_FOOTBALL_SPORT_ID = 11
MGM_NFL_COMPETITION_NAME = "NFL"
SPREAD_MARKET = {"GM": "Spread", "Q1": "1st quarter spread", "H1": "1st half spread"}
# "Seattle Seahawks -7.5" / "Los Angeles Chargers (+6)"
OPTION_RE = re.compile(r"^(.+?)\s*\(?([+-]\d+(?:\.\d+)?)\)?\s*$")
PRICED_STATE = "MarketUnsuspended"


class BetMgmBook:
    name = "BetMGM"

    def __init__(self) -> None:
        self.client = BetMGMClient()
        self._games = self._list_games()
        self._markets: dict[str, list[dict]] = {}

    # --- discovery ---------------------------------------------------------
    def _list_games(self) -> list[BookGame]:
        url = (f"{self.client._base}/bettingoffer/fixtures?{self.client._common_q()}"
               f"&fixtureTypes=Standard&state=Latest&offerMapping=Filtered"
               f"&offerCategories=Gridable&fixtureCategories=Gridable"
               f"&sportIds={MGM_FOOTBALL_SPORT_ID}&skip=0&take=200&sortBy=Tags")
        resp = self.client.session.get(url, timeout=60)
        if resp.status_code != 200:
            raise RuntimeError(f"BetMGM fixture listing returned HTTP {resp.status_code}")
        games = []
        for fixture in resp.json().get("fixtures", []):
            competition = ((fixture.get("competition") or {}).get("name") or {}).get("value")
            if competition != MGM_NFL_COMPETITION_NAME:
                continue
            away_name, sep, home_name = ((fixture.get("name") or {}).get("value") or "").partition(" @ ")
            home, away = ODDS_API_NAME_TO_ABBR.get(home_name), ODDS_API_NAME_TO_ABBR.get(away_name)
            if not sep or home is None or away is None:
                continue
            games.append(BookGame(book_event_id=str(fixture["id"]), home=home, away=away,
                                  game_start_time=fixture.get("startDate", "")))
        return games

    def find_game(self, team: str) -> BookGame | None:
        return next_game(self._games, team)

    def _option_markets(self, game: BookGame) -> list[dict]:
        if game.book_event_id not in self._markets:
            url = (f"{self.client._base}/bettingoffer/fixtures?{self.client._common_q()}"
                   f"&fixtureIds={game.book_event_id}&offerMapping=All&state=Latest")
            resp = self.client.session.get(url, timeout=60)
            if resp.status_code != 200:
                raise RuntimeError(f"BetMGM market tree {game.book_event_id} returned HTTP {resp.status_code}")
            fixtures = resp.json().get("fixtures") or [{}]
            self._markets[game.book_event_id] = [
                m for m in fixtures[0].get("optionMarkets", [])
                if m.get("isBetBuilder") and m.get("status") == "Visible"]
        return self._markets[game.book_event_id]

    # --- legs --------------------------------------------------------------
    def leg_market(self, game: BookGame, role: str, leg: Leg) -> LegMarket | None:
        markets = self._option_markets(game)
        team = game.home if role == "home" else game.away
        if leg.kind == "scores_first":
            return None
        if leg.kind == "win" and leg.period == "GM":
            return _moneyline(markets, team)
        if leg.kind == "win":
            return _spread(markets, SPREAD_MARKET[leg.period], team, -0.5)
        if leg.kind == "spread":
            return _spread(markets, SPREAD_MARKET[leg.period], team, leg.line)
        raise ValueError(f"unknown leg kind {leg.kind!r}")

    # --- pricing -----------------------------------------------------------
    def price(self, game: BookGame, refs: tuple[Hashable, ...]) -> float | None:
        pick_group = str(uuid.uuid4())
        body = {"tv1Picks": [], "tv2Picks": [
            {"fixtureId": game.book_event_id, "optionMarketId": market_id, "optionId": option_id,
             "useLiveFallback": False, "pickGroupId": pick_group}
            for market_id, option_id in refs]}
        url = f"{self.client._base}/bettingoffer/picks?{self.client._common_q()}"
        resp = self.client.session.post(url, json=body, timeout=30)
        if resp.status_code != 200:
            raise RuntimeError(f"BetMGM picks returned HTTP {resp.status_code}: {resp.text[:200]}")
        groups = resp.json().get("betBuilderPricingGroups") or {}
        if not groups:
            return None
        group = next(iter(groups.values()))
        # A suspended group still carries odds (a +500 that is not on offer);
        # only an unsuspended, error-free group is a price.
        if group.get("error") or group.get("suspensionState") != PRICED_STATE:
            return None
        decimal = (group.get("odds") or {}).get("odds")
        return float(decimal) if decimal and decimal > 1.0 else None


def _visible_options(market: dict) -> list[dict]:
    return [o for o in market.get("options", []) if o.get("status") == "Visible"]


def _moneyline(markets: list[dict], team: str) -> LegMarket | None:
    for market in markets:
        if market["name"]["value"] != "Moneyline":
            continue
        options = _visible_options(market)
        group = tuple(Outcome(ref=(market["id"], o["id"]), label=o["name"]["value"]) for o in options)
        chosen = next((out for out, o in zip(group, options)
                       if ODDS_API_NAME_TO_ABBR.get(o["name"]["value"]) == team), None)
        if chosen is not None and len(group) == 2:
            return LegMarket.single(chosen, group)
    return None


def _spread(markets: list[dict], market_name: str, team: str, line: float) -> LegMarket | None:
    """Each BetMGM spread market is ONE line pair, e.g. SEA -0.5 / LAC +0.5."""
    for market in markets:
        if market["name"]["value"] != market_name:
            continue
        options = _visible_options(market)
        if len(options) != 2:
            continue
        parsed = [OPTION_RE.match(o["name"]["value"]) for o in options]
        if not all(parsed):
            continue
        for option, match, other_option in ((options[0], parsed[0], options[1]),
                                            (options[1], parsed[1], options[0])):
            if ODDS_API_NAME_TO_ABBR.get(match.group(1)) == team and abs(float(match.group(2)) - line) < 1e-6:
                pick = Outcome(ref=(market["id"], option["id"]), label=option["name"]["value"])
                other = Outcome(ref=(market["id"], other_option["id"]), label=other_option["name"]["value"])
                return LegMarket.single(pick, (pick, other))
    return None
