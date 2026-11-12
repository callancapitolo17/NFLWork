"""Wagerzon side of the 1H parlays: read every NFL and college 1st-half
spread + total, price each same-game spread + total parlay, and place one.

Reads:  ActiveLeaguesHelper (resolve WZ_HALF_LEAGUES by Description), then
        NewScheduleHelper for those leagues (JSON; needs the XHR header).
        Each 1H line is its own Wagerzon game (idgmtyp 15) with its own idgm.
Prices: ConfirmWagerHelper with RiskWin=2 (a preview: no balance check, never
        placed), through wagerzon_odds/parlay_pricer.get_combined_parlay_price.
Places: re-reads the line, then parlay_placer's ConfirmWagerHelper preflight
        (RiskWin=0) and ONE PostWagerMultipleHelper submit on an explicit
        account. Never parlay_placer.place_parlays: it re-submits after a
        logged-out reply to the submit, which can place the bet twice.
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
import requests  # noqa: E402
from parlay_pricer import (PLAY_AWAY_SPREAD, PLAY_HOME_SPREAD, PLAY_OVER, PLAY_UNDER,  # noqa: E402
                           get_combined_parlay_price)
from scraper_v2 import _wz_game_start_time, fetch_active_leagues  # noqa: E402

log = logging.getLogger("nfl_specials.wz_halves")

# Wagerzon catalog Description -> league (resolved to an id on every read).
WZ_HALF_LEAGUES = {"NFL - 1ST HALVES": NFL, "NCAA FOOTBALL - 1ST HALVES": CFB}
WZ_FIRST_HALF_GAME_TYPE = 15
# What parlay_pricer sends; Wagerzon priced football parlays with it (2026-11-12).
WZ_PITCHER_FIELD = 3
# Preview stakes, largest first: $500 prices to +/-0.001 in decimal; $100 is
# the fallback when a parlay's win would pass Wagerzon's max (MLB measured).
PREVIEW_STAKES = (500, 100)
PREVIEW_WORKERS = 8
# Wagerzon pays whole dollars (wz.wagerzon_to_win), so the to-win at a stake
# can sit up to $0.50 under the board's price, plus the board decimal's own
# rounding (a whole-dollar to-win on a $100 preview, the fallback: +/-0.005).
WIN_ROUNDING_DOLLARS = 0.5
BOARD_DECIMAL_PRECISION = 0.005

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
    error: str | None = None     # Wagerzon's ErrorMsgKey or why there is no price


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
    games: dict[int, WzHalfGame] = {}
    # One request per league: a combined lg=430,540 answers only the first sport (2026-11-12).
    for league_id, league_name in sorted(league_of_id.items()):
        resp = session.get(wz.SCHEDULE_URL.format(league_id=league_id), timeout=30,
                           headers={"Accept": "application/json", "X-Requested-With": "XMLHttpRequest"})
        resp.raise_for_status()
        for league in ((resp.json().get("result") or {}).get("listLeagues") or [[]])[0]:
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
            return WzParlayPrice(decimal=quote["decimal"])
    return WzParlayPrice(None, error="Wagerzon gave no preview price (see the log for its error)")


@dataclass(frozen=True)
class HalfPlacement:
    status: str               # 'placed' | 'unknown' | 'network_error' (both: may have placed) |
                              # 'line_moved' | 'price_moved' | 'rejected' | 'auth_error' (not placed)
    win: float | None         # the to-win Wagerzon quoted at the stake
    ticket_number: str | None = None
    error: str | None = None


# Placement statuses after which the bet may be on Wagerzon.
MAYBE_PLACED = ("placed", "unknown", "network_error")
CHECK_OPEN_BETS = "check Wagerzon's open bets before placing again"


def place_half_parlay(account_label: str, game: WzHalfGame, corner: str, risk: float,
                      board_decimal: float) -> HalfPlacement:
    """Place `risk` dollars on one corner at the board's price, or refuse.

    1. Re-reads Wagerzon's 1H lines: a game whose numbers or prices differ
       from `game` is 'line_moved' (the fair was priced at those numbers).
    2. ConfirmWagerHelper preflight at `risk` (RiskWin=0, balance-checked;
       one re-login on a logged-out reply, nothing is placed by it). A to-win
       under risk x (board_decimal - 1) by more than the two roundings
       (WIN_ROUNDING_DOLLARS, BOARD_DECIMAL_PRECISION) is 'price_moved'.
       A failure in steps 1-2 is 'rejected': nothing was sent.
    3. ONE PostWagerMultipleHelper submit, never re-sent. A logged-out,
       failed or unreadable reply is 'unknown' / 'network_error': the bet
       may be on Wagerzon.
    """
    account = next((a for a in wz.list_accounts() if a.label == account_label), None)
    if account is None:
        raise ValueError(f"unknown Wagerzon account {account_label!r}")
    try:
        current = {g.wz_game_id: g for g in fetch_half_games()}.get(game.wz_game_id)
    except Exception as exc:  # nothing sent yet: a refusal, never a maybe-placed bet
        return HalfPlacement("rejected", None, error=f"couldn't re-read Wagerzon's 1H lines ({type(exc).__name__})")
    if current != game:
        return HalfPlacement("line_moved", None, error="Wagerzon's 1H line changed since the board priced it; refresh")
    legs = [parlay_placer.Leg(idgm=leg["idgm"], play=leg["play"], points=leg["points"], odds=leg["odds"],
                              pitcher=WZ_PITCHER_FIELD) for leg in corner_legs(game, corner)]
    try:
        win, _risk = _preflight(legs, risk, account)
    except parlay_placer.PreflightRejected as exc:
        return HalfPlacement("rejected", None, error=exc.error_msg_key or exc.error_msg or "rejected")
    except parlay_placer.AuthExpired:
        return HalfPlacement("auth_error", None, error="Wagerzon login failed at the preflight")
    except Exception as exc:  # nothing sent yet: a refusal, never a maybe-placed bet
        return HalfPlacement("rejected", None, error=f"preflight failed: {type(exc).__name__}")
    if win < risk * (board_decimal - 1.0) - WIN_ROUNDING_DOLLARS - risk * BOARD_DECIMAL_PRECISION:
        return HalfPlacement("price_moved", win,
                             error=f"board to-win ${risk * (board_decimal - 1.0):.2f}, Wagerzon now ${win:.2f}")
    spec = parlay_placer.ParlaySpec(parlay_hash=f"wz-1h-{game.wz_game_id}-{corner}", legs=legs,
                                    amount=risk, expected_win=win, expected_risk=risk)
    try:
        result = parlay_placer._post_wagers([spec], account)[0]
    except parlay_placer.AuthExpired:
        return HalfPlacement("unknown", win, error=f"Wagerzon answered the submit logged out; {CHECK_OPEN_BETS}")
    except requests.RequestException as exc:
        return HalfPlacement("network_error", win, error=f"submit failed ({type(exc).__name__}); {CHECK_OPEN_BETS}")
    except (KeyError, TypeError, ValueError, IndexError) as exc:
        return HalfPlacement("unknown", win, error=f"unreadable submit reply ({type(exc).__name__}); {CHECK_OPEN_BETS}")
    if result.status == "placed":
        return HalfPlacement("placed", result.actual_win or win, ticket_number=result.ticket_number)
    return HalfPlacement(result.status, win, error=result.error_msg or result.error_msg_key or None)


def _preflight(legs: list, risk: float, account) -> tuple[float, float]:
    """parlay_placer's ConfirmWagerHelper preflight, logged in again once when
    Wagerzon answers logged out (a preflight places nothing)."""
    try:
        return parlay_placer._confirm_preflight(legs, risk, account)
    except parlay_placer.AuthExpired:
        wz.wagerzon_auth.clear_session_cache(account.label)
        return parlay_placer._confirm_preflight(legs, risk, account)
