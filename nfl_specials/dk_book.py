"""DraftKings adapter (books.SgpBook) for NFL fecta legs.

Plain HTTP throughout (curl_cffi Chrome TLS; none of these is Akamai-gated):
    league listing   sportsbook-nash .../leagueSubcategory/v1/markets
    SGP payload      sportsbook-nash .../parlays/v1/sgp/events/{event_id}
    SGP price        sportsbook-nash .../sgp/dkuswv/sportsdata/v2/sgp   (GET)
The price endpoint is the one DK's own SGP builder widget calls as you click
legs (dk-same-game-parlay 9.0.1). Given base legs and a candidate market it
returns the SGP price of base + EACH outcome of that market, so one request
prices every cell that differs only in its last leg. Verified 2026-10-02:
10/10 superfecta cells identical to the betslip's calculateBets (the call
issue #102 locked behind a real Chrome). No DB writes.

Leg -> DK market:
    scores_first   "1st to Score"                     (home, away)
    win GM         "Moneyline"                        (home, away)
    win Q1 / H1    "1st Quarter (3 Way)" / "1st Half (3 Way)"  (home, Tie, away)
    spread         main or alternate spread for the period at that exact line,
                   paired with the opponent at the opposite line
"""
from __future__ import annotations

import logging
import time
from typing import Hashable

from curl_cffi import requests as cffi_requests

from nfl_specials.books import BookGame, LegMarket, Outcome, next_game, is_half_point, widen_with_tie
from nfl_specials.special_parser import Leg
from nfl_specials.teams import NICKNAME_TO_ABBR

log = logging.getLogger("nfl_specials.dk_book")

DK_NFL_LEAGUE_ID = "88808"
DK_GAME_LINES_SUBCATEGORY = "4518"
DK_LEAGUE_URL = ("https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/"
                 "controldata/league/leagueSubcategory/v1/markets")
DK_SGP_EVENT_URL = ("https://sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/"
                    "parlays/v1/sgp/events/{event_id}")
DK_SGP_PRICE_URL = ("https://sportsbook-nash.draftkings.com/sites/US-WV-SB/api/sportscontent/"
                    "sgp/dkuswv/sportsdata/v2/sgp")
DK_WARMUP_URL = "https://sportsbook.draftkings.com/leagues/football/nfl"
# DK's sport id for football. The price endpoint answers 404 without the header.
DK_FOOTBALL_SPORT_ID = "3"
# DK's own SGP widget throttles its price requests to one per 500 ms; never
# go faster than a person clicking legs.
MIN_SECONDS_BETWEEN_PRICE_CALLS = 0.5

SPREAD_MAIN = {"GM": "Spread", "Q1": "Spread 1st Quarter", "H1": "Spread 1st Half"}
SPREAD_ALT = {"GM": "Spread Alternate", "Q1": "Spread Alternate - 1st Quarter",
              "H1": "Spread Alternate - 1st Half"}
THREE_WAY_WIN = {"Q1": "1st Quarter (3 Way)", "H1": "1st Half (3 Way)"}
ROLE_TO_DK = {"home": "Home", "away": "Away"}


def dk_team_to_abbr(dk_name: str) -> str | None:
    """'SEA Seahawks' -> 'SEA'. The nickname decides: DK's city prefix is
    ambiguous for the two LA and two NY teams."""
    return NICKNAME_TO_ABBR.get(dk_name.split()[-1].upper())


class DraftKingsBook:
    name = "DraftKings"

    def __init__(self) -> None:
        self.session = cffi_requests.Session(impersonate="chrome")
        self.session.get(DK_WARMUP_URL, timeout=30)
        self._games = self._list_games()
        self._markets: dict[str, dict[str, dict]] = {}
        self._market_of_selection: dict[str, str] = {}
        # (event id, base selection ids, candidate market id) -> decimal odds of
        # base + each candidate selection DK can combine with it
        self._added_leg_odds: dict[tuple[str, frozenset[str], str], dict[str, float]] = {}
        self._last_price_call_at = 0.0

    # --- discovery ---------------------------------------------------------
    def _list_games(self) -> list[BookGame]:
        params = {
            "isBatchable": "false",
            "templateVars": DK_NFL_LEAGUE_ID,
            "eventsQuery": (f"$filter=leagueId eq '{DK_NFL_LEAGUE_ID}' AND clientMetadata/"
                            f"Subcategories/any(s: s/Id eq '{DK_GAME_LINES_SUBCATEGORY}')"),
            "marketsQuery": (f"$filter=clientMetadata/subCategoryId eq '{DK_GAME_LINES_SUBCATEGORY}' "
                             "AND tags/all(t: t ne 'SportcastBetBuilder')"),
            "include": "Events",
            "entity": "events",
        }
        resp = self.session.get(DK_LEAGUE_URL, params=params, timeout=30)
        if resp.status_code != 200:
            raise RuntimeError(f"DK event listing returned HTTP {resp.status_code}")
        games = []
        for event in resp.json().get("events", []):
            roles = {p.get("venueRole"): p.get("name", "") for p in event.get("participants", [])}
            home, away = dk_team_to_abbr(roles.get("Home", "")), dk_team_to_abbr(roles.get("Away", ""))
            if home is None or away is None:
                log.warning("DK event %s: unmapped teams %s", event.get("id"), roles)
                continue
            games.append(BookGame(book_event_id=str(event["id"]), home=home, away=away,
                                  game_start_time=event.get("startEventDate", "")))
        return games

    def find_game(self, team: str) -> BookGame | None:
        return next_game(self._games, team)

    def _event_markets(self, game: BookGame) -> dict[str, dict]:
        if game.book_event_id not in self._markets:
            url = DK_SGP_EVENT_URL.format(event_id=game.book_event_id)
            resp = self.session.get(url, timeout=60)
            if resp.status_code != 200:
                raise RuntimeError(f"DK SGP payload for event {game.book_event_id} "
                                   f"returned HTTP {resp.status_code}")
            markets = (resp.json().get("data") or {}).get("markets") or []
            self._markets[game.book_event_id] = {
                m["name"]: m for m in markets
                if "SGP" in (m.get("tags") or []) and m.get("selections")}
            for market in self._markets[game.book_event_id].values():
                for selection in market["selections"]:
                    self._market_of_selection[selection["id"]] = market["id"]
        return self._markets[game.book_event_id]

    # --- legs --------------------------------------------------------------
    def leg_market(self, game: BookGame, role: str, leg: Leg) -> LegMarket | None:
        markets = self._event_markets(game)
        dk_role = ROLE_TO_DK[role]
        if leg.kind == "scores_first":
            return _outcome_market(markets.get("1st to Score"), dk_role)
        if leg.kind == "win" and leg.period == "GM":
            return _outcome_market(markets.get("Moneyline"), dk_role)
        if leg.kind == "win":
            # Wins the period (a tie loses) = -0.5. A 2-way -0.5 spread, when DK
            # posts one, halves the cells this leg adds to a partition.
            for name in (SPREAD_MAIN[leg.period], SPREAD_ALT[leg.period]):
                found = _spread_market(markets.get(name), dk_role, -0.5)
                if found:
                    return found
            return _outcome_market(markets.get(THREE_WAY_WIN[leg.period]), dk_role)
        if leg.kind == "spread":
            for name in (SPREAD_MAIN[leg.period], SPREAD_ALT[leg.period]):
                found = _spread_market(markets.get(name), dk_role, leg.line)
                if found:
                    return found
            if leg.period in THREE_WAY_WIN and is_half_point(leg.line):
                three_way = _outcome_market(markets.get(THREE_WAY_WIN[leg.period]), dk_role)
                return three_way if leg.line < 0 else widen_with_tie(three_way)
            return None
        raise ValueError(f"unknown leg kind {leg.kind!r}")

    # --- pricing -----------------------------------------------------------
    def price(self, game: BookGame, refs: tuple[Hashable, ...]) -> float | None:
        """DK's SGP price for `refs`, or None when DK will not combine them.

        The last ref is priced as an addition to the others, so cells that
        share every leg but the last (consecutive cells of a partition) cost
        one request between them.
        """
        if len(refs) < 2:
            raise ValueError(f"an SGP needs at least 2 legs, got {len(refs)}")
        base, last = frozenset(refs[:-1]), refs[-1]
        self._event_markets(game)                    # fills _market_of_selection
        key = (game.book_event_id, base, self._market_of_selection[last])
        if key not in self._added_leg_odds:
            self._added_leg_odds[key] = self._fetch_added_leg_odds(*key)
        return self._added_leg_odds[key].get(last)

    def _fetch_added_leg_odds(self, event_id: str, base: frozenset[str],
                              candidate_market_id: str) -> dict[str, float]:
        """One GET: decimal odds of base + each selection of the candidate
        market that DK can combine with the base. Raises on HTTP failure."""
        wait = MIN_SECONDS_BETWEEN_PRICE_CALLS - (time.monotonic() - self._last_price_call_at)
        if wait > 0:
            time.sleep(wait)
        try:
            resp = self.session.get(DK_SGP_PRICE_URL, timeout=30,
                                    headers={"X-SportId": DK_FOOTBALL_SPORT_ID},
                                    params={"eventId": event_id, "selections": ",".join(sorted(base)),
                                            "marketCandidates": candidate_market_id,
                                            "oddsStyle": "american"})
        finally:
            self._last_price_call_at = time.monotonic()
        if resp.status_code != 200:
            raise RuntimeError(f"DK SGP price for event {event_id} returned HTTP "
                               f"{resp.status_code}: {resp.text[:200]}")
        return parse_added_leg_odds(resp.json(), base)


def parse_added_leg_odds(body: dict, base: frozenset[str]) -> dict[str, float]:
    """Map DK's SGP-builder answer to {candidate selection id: decimal odds}.

    DK does not refuse a base it cannot combine: it drops the offending legs
    (`selectionsNotMapped`) and prices the rest. A price on a smaller base
    would be a silent wrong number, so a base DK did not take whole prices
    nothing. A candidate absent from `compatibleMarkets` is one DK declines.
    """
    data = body.get("data") or {}
    if data.get("selectionsNotMapped") or set(data.get("selectionsMapped") or []) != base:
        return {}
    return {selection["id"]: float(selection["trueOdds"])
            for market in data.get("compatibleMarkets") or []
            for selection in market.get("selections") or []
            if selection.get("trueOdds")}


def _live(market: dict | None) -> list[dict]:
    if market is None:
        return []
    return [s for s in market.get("selections", [])
            if not s.get("isDisabled") and (s.get("displayOdds") or {}).get("american")]


def _role_of(selection: dict) -> str | None:
    participants = selection.get("participants") or []
    return participants[0].get("venueRole") if participants else None


def _outcome_market(market: dict | None, dk_role: str) -> LegMarket | None:
    """2- or 3-way team market: every live selection is one outcome."""
    selections = _live(market)
    group = tuple(Outcome(ref=s["id"], label=s.get("name", "")) for s in selections)
    chosen = next((Outcome(ref=s["id"], label=s.get("name", ""))
                   for s in selections if s.get("outcomeType") == dk_role), None)
    if chosen is None or len(group) < 2:
        return None
    return LegMarket.single(chosen, group)


def _spread_market(market: dict | None, dk_role: str, line: float) -> LegMarket | None:
    selections = _live(market)
    chosen = next((s for s in selections if _role_of(s) == dk_role
                   and s.get("points") is not None and abs(s["points"] - line) < 1e-6), None)
    opponent = next((s for s in selections if _role_of(s) not in (dk_role, None)
                     and s.get("points") is not None and abs(s["points"] + line) < 1e-6), None)
    if chosen is None or opponent is None:
        return None
    pick = Outcome(ref=chosen["id"], label=f"{chosen.get('name', '')} {line:+g}")
    other = Outcome(ref=opponent["id"], label=f"{opponent.get('name', '')} {-line:+g}")
    return LegMarket.single(pick, (pick, other))
