"""Wagerzon side of the 1H parlays: read every NFL and college 1st-half
spread + total, price each same-game spread + total parlay, and place one.

Reads:  ActiveLeaguesHelper (resolve WZ_HALF_LEAGUES by Description), then
        NewScheduleHelper for those leagues (JSON; needs the XHR header).
        Each 1H line is its own Wagerzon game (idgmtyp 15) with its own idgm.
Prices: ConfirmWagerHelper with RiskWin=2 (a preview: no balance check, never
        placed), through wagerzon_odds/parlay_pricer.get_combined_parlay_price.
Places: wagerzon_odds/parlay_placer.place_parlays with an explicit account
        (ConfirmWagerHelper preflight, $0.01 win-drift check, then
        PostWagerMultipleHelper).
No DB writes here. Wagerzon refuses quarter lines in any parlay
(WAGERTYPENOALLOWGAME, checked live 2026-11-12), so only halves are read.
"""
from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime

from nfl_specials import wz
from nfl_specials.books import CFB, NFL
from nfl_specials.teams import resolve_wz_team

import parlay_placer  # noqa: E402  (wz put wagerzon_odds on sys.path)
from parlay_pricer import get_combined_parlay_price  # noqa: E402
from scraper_v2 import _wz_game_start_time, fetch_active_leagues  # noqa: E402

log = logging.getLogger("nfl_specials.wz_halves")

# Wagerzon catalog Description -> league (resolved to an id on every read).
WZ_HALF_LEAGUES = {"NFL - 1ST HALVES": NFL, "NCAA FOOTBALL - 1ST HALVES": CFB}
SCHEDULE_URL = wz.config.WZ_BASE_URL + "/wager/NewScheduleHelper.aspx?WT=0&lg={league_ids}"
WZ_FIRST_HALF_GAME_TYPE = 15
# Play codes ConfirmWagerHelper takes (parlay_pricer).
PLAY_AWAY_SPREAD, PLAY_HOME_SPREAD, PLAY_OVER, PLAY_UNDER = 0, 1, 2, 3
# What parlay_pricer sends; Wagerzon priced football parlays with it (2026-11-12).
WZ_PITCHER_FIELD = 3
# Preview stakes, largest first: $500 prices to +/-0.001 in decimal; $100 is
# the fallback when a parlay's win would pass Wagerzon's max (MLB measured).
PREVIEW_STAKES = (500, 100)
PREVIEW_WORKERS = 8
# The preview at the stake may differ from the board price by Wagerzon's
# whole-dollar rounding; past this many dollars of to-win the price moved.
PRICE_MOVE_TOLERANCE_DOLLARS = 1.0

SPREAD_SIDES = ("away", "home")
TOTAL_SIDES = ("over", "under")
# corner = "<spread side>_<total side>", e.g. "away_over"
CORNERS = tuple(f"{side}_{total}" for side in SPREAD_SIDES for total in TOTAL_SIDES)


@dataclass(frozen=True)
class WzHalfGame:
    league: str                  # 'NFL' | 'CFB'
    wz_game_id: int              # idgm of the 1H line, what the parlay is placed against
    away: str                    # NFL: nflverse abbreviation; CFB: Wagerzon's name, e.g. 'MEMPHIS'
    home: str
    game_start_time: datetime    # UTC
    away_spread: float           # away team's 1H line; the home line is its negative
    away_spread_american: int
    home_spread_american: int
    total: float
    over_american: int
    under_american: int

    @property
    def home_spread(self) -> float:
        return -self.away_spread


@dataclass(frozen=True)
class WzParlayPrice:
    decimal: float | None        # Wagerzon's parlay payout as decimal odds
    preview_stake: int | None
    error: str | None = None     # Wagerzon's ErrorMsgKey or why there is no price


def whole_number_skip_reason(game: WzHalfGame) -> str | None:
    """The one place whole-number lines are kept off the board (v1, Cal
    2026-11-12): a pushed leg makes Wagerzon grade the parlay as the other leg
    alone, and the books' prices do not carry the push chance."""
    if game.away_spread == int(game.away_spread) or game.total == int(game.total):
        return "whole-number line, push not priced"
    return None


def corner_legs(game: WzHalfGame, corner: str) -> list[dict]:
    """ConfirmWagerHelper legs for one corner: the spread side at its line and
    price, then the total (Wagerzon wants -total for the Over, +total for the Under)."""
    spread_side, total_side = corner.split("_")
    spread_leg = ({"play": PLAY_AWAY_SPREAD, "points": game.away_spread, "odds": game.away_spread_american}
                  if spread_side == "away" else
                  {"play": PLAY_HOME_SPREAD, "points": game.home_spread, "odds": game.home_spread_american})
    total_leg = ({"play": PLAY_OVER, "points": -game.total, "odds": game.over_american}
                 if total_side == "over" else
                 {"play": PLAY_UNDER, "points": game.total, "odds": game.under_american})
    return [{"idgm": game.wz_game_id, **spread_leg}, {"idgm": game.wz_game_id, **total_leg}]


def fetch_half_games() -> list[WzHalfGame]:
    """Every NFL and college 1H spread + total Wagerzon posts now."""
    return wz._retry_reads(_fetch_half_games_once)


def _fetch_half_games_once() -> list[WzHalfGame]:
    session = wz._session()
    league_of_id = {int(row["IdLeague"]): WZ_HALF_LEAGUES[(row.get("Description") or "").strip()]
                    for row in fetch_active_leagues(session) or []
                    if (row.get("Description") or "").strip() in WZ_HALF_LEAGUES and row.get("IdLeague")}
    if not league_of_id:
        raise RuntimeError(f"none of {sorted(WZ_HALF_LEAGUES)} in Wagerzon's catalog")
    resp = session.get(SCHEDULE_URL.format(league_ids=",".join(str(i) for i in sorted(league_of_id))),
                       timeout=30, headers={"Accept": "application/json", "X-Requested-With": "XMLHttpRequest"})
    resp.raise_for_status()
    games: dict[int, WzHalfGame] = {}
    for league in ((resp.json().get("result") or {}).get("listLeagues") or [[]])[0]:
        league_name = WZ_HALF_LEAGUES.get((league.get("Description") or "").strip())
        if league_name is None:
            continue
        for raw in league.get("Games", []):
            game = parse_half_game(raw, league_name)
            if game is not None:
                games[game.wz_game_id] = game
    return list(games.values())


def parse_half_game(raw: dict, league: str) -> WzHalfGame | None:
    """One NewScheduleHelper 1H game, or None when it is not a 1H line with
    both a spread and a total posted, or its NFL teams do not resolve."""
    if raw.get("idgmtyp") != WZ_FIRST_HALF_GAME_TYPE:
        return None
    line = (raw.get("GameLines") or [{}])[0] or {}
    numbers = [_number(line.get(key)) for key in ("vsprdt", "unt")]
    prices = [_american(line.get(key)) for key in ("vsprdoddst", "hsprdoddst", "ovoddst", "unoddst")]
    start = _wz_game_start_time(raw.get("gmdt", ""), raw.get("gmtm", ""))
    if None in numbers or None in prices or start is None:
        return None
    away, home = _team(raw.get("vtm"), league), _team(raw.get("htm"), league)
    if away is None or home is None:
        log.warning("Wagerzon 1H game %s: unmapped teams %r @ %r", raw.get("idgm"), raw.get("vtm"), raw.get("htm"))
        return None
    return WzHalfGame(league=league, wz_game_id=int(raw["idgm"]), away=away, home=home, game_start_time=start,
                      away_spread=numbers[0], away_spread_american=prices[0], home_spread_american=prices[1],
                      total=numbers[1], over_american=prices[2], under_american=prices[3])


def _team(name: str | None, league: str) -> str | None:
    """'1H WAS COMMANDERS' -> 'WAS' (NFL, by nickname) or '1H MEMPHIS' -> 'MEMPHIS' (CFB)."""
    text = (name or "").strip().upper()
    if text.startswith("1H "):
        text = text[3:].strip()
    if not text:
        return None
    return resolve_wz_team(text.split()[-1]) if league == NFL else text


def _number(value) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except ValueError:
        return None


def _american(value) -> int | None:
    number = _number(value)
    return int(number) if number is not None and number != 0 else None


def price_corners(games: list[WzHalfGame]) -> dict[tuple[int, str], WzParlayPrice]:
    """Wagerzon's preview price of every corner of every game, PREVIEW_WORKERS at a time."""
    session = wz._session()
    work = [(game, corner) for game in games for corner in CORNERS]
    with ThreadPoolExecutor(max_workers=PREVIEW_WORKERS) as pool:
        prices = list(pool.map(lambda item: _price_corner(session, *item), work))
    return {(game.wz_game_id, corner): price for (game, corner), price in zip(work, prices)}


def _price_corner(session, game: WzHalfGame, corner: str) -> WzParlayPrice:
    legs = corner_legs(game, corner)
    for stake in PREVIEW_STAKES:
        quote = get_combined_parlay_price(session, legs, amount=stake)
        if quote is not None:
            return WzParlayPrice(decimal=quote["decimal"], preview_stake=stake)
    return WzParlayPrice(None, None, error="Wagerzon gave no preview price (see the log for its error)")


def place_half_parlay(account_label: str, game: WzHalfGame, corner: str, risk: float,
                      board_decimal: float) -> parlay_placer.PlacementResult:
    """Place `risk` dollars on one corner at the board's price.

    Previews the parlay at `risk` first (RiskWin=2, nothing placed) and
    refuses with status 'price_moved' when the to-win differs from
    risk x (board_decimal - 1) by more than PRICE_MOVE_TOLERANCE_DOLLARS;
    otherwise places it with that previewed to-win, which place_parlays'
    own preflight must match to the cent. Never retried here: a resent
    submission can place the bet twice.
    """
    account = next((a for a in wz.list_accounts() if a.label == account_label), None)
    if account is None:
        raise ValueError(f"unknown Wagerzon account {account_label!r}")
    legs = corner_legs(game, corner)
    parlay_hash = f"wz-1h-{game.wz_game_id}-{corner}"
    quote = get_combined_parlay_price(wz.wagerzon_auth.get_session(account), legs, amount=risk)
    if quote is None:
        return parlay_placer.PlacementResult(parlay_hash, "rejected", error_msg="Wagerzon refused the preview")
    expected_win = risk * (board_decimal - 1.0)
    if abs(quote["win"] - expected_win) > PRICE_MOVE_TOLERANCE_DOLLARS:
        return parlay_placer.PlacementResult(
            parlay_hash, "price_moved", actual_win=quote["win"],
            error_msg=f"board to-win ${expected_win:.2f}, Wagerzon now ${quote['win']:.2f}")
    spec = parlay_placer.ParlaySpec(
        parlay_hash=parlay_hash,
        legs=[parlay_placer.Leg(idgm=leg["idgm"], play=leg["play"], points=leg["points"], odds=leg["odds"],
                                pitcher=WZ_PITCHER_FIELD) for leg in legs],
        amount=risk, expected_win=float(quote["win"]), expected_risk=risk)
    return parlay_placer.place_parlays([spec], account)[0]
