"""Kalshi venue adapter: NFL games, every leg Kalshi lets into a combo, and a
parlay price from a request for quotes (RFQ).

Kalshi posts no parlay price. A parlay is priced by asking market makers:
    POST   /multivariate_event_collections/KXMVECROSSCATEGORY-R   mint (or look up) the combo
    POST   /communications/rfqs        {market_ticker, target_cost_dollars: 100}
    GET    /communications/quotes?rfq_id=...   polled for up to QUOTE_WINDOW_MAX_SEC
    DELETE /communications/rfqs/{id}   in a finally block; after an unclear create
                                       reply, any open RFQ of ours on the combo
This module has no code that accepts a quote, so it can never bet. Market
makers do see each RFQ (Cal OK'd that, 2026-11-18). Reads:
    GET /events?series_ticker=KXNFLGAME&status=open&with_nested_markets=true   the games
    GET /multivariate_event_collections/KXMVECROSSCATEGORY-R   which events can go in a combo
    GET /events/{event_ticker}?with_nested_markets=true         each such event of the game
All requests are signed (kalshi_common.auth_client) with the key in the
environment or kalshi_mlb_rfq/.env. No DB I/O.

Kalshi shapes (recon 2026-11-19 on Cal's Mac, IND at HOU; trimmed copies in
tests/fixtures/kalshi_*.json):
  * game key "26NOV19INDHOU" = date + AWAY code + HOME code; the KXNFLGAME
    markets end "-IND" / "-HOU" and their yes_sub_title is the team name.
  * a game market's close_time is the scheduled start + 48 h (the rules'
    postponement window); that is the only start time the event list gives.
  * every market is a binary YES/NO contract; `yes_ask_dollars` /
    `no_ask_dollars` are the prices to buy each side, `yes_ask_size_fp` the
    contracts at the YES ask and `yes_bid_size_fp` those at the NO ask.
  * spreads, totals, team totals and player props are one market per number
    ("Houston wins by over 4.5 points", "Nico Collins: 150+"), strike
    `floor_strike` (4.5, 149.5); YES = over / the named team covers.
  * a full-game moneyline tie pays 0.50 a side, so NO on one team pays exactly
    what YES on the other does. 1H / 1Q winner markets have a third "Tie"
    market and a tie makes both team markets NO.
  * the collection's associated event with `size_max` 1 takes one market per
    combo; the player-prop and team-total events (`size_max` absent) take
    several, one per player / team here.
  * an RFQ quote's `no_bid_dollars` is what the maker pays for NO, so buying
    the combo YES costs 1 - no_bid; `no_contracts_fp` is the YES contracts it
    gives for the RFQ's dollars. Makers answered within ~1 s (20 quotes).
  * taker fee 0.07 x p x (1 - p) per contract (series fee_type
    quadratic_with_maker_fees, multiplier 1). Every price here includes it.
"""
from __future__ import annotations

import logging
import os
import threading
import time
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, replace
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

from sgp_compare.legs import CanonicalLeg, Game, OfferedLeg, ParlayQuote, ParlayTier
from sgp_compare.venues.base import VenueTransportError

log = logging.getLogger(__name__)

COLLECTION_TICKER = "KXMVECROSSCATEGORY-R"
GAME_SERIES = "KXNFLGAME"
LEAGUE_PREFIX = "KXNFL"
KALSHI_BASE_URL = "https://api.elections.kalshi.com/trade-api/v2"
TAKER_FEE_RATE = 0.07
GAME_CLOSE_AFTER_START = timedelta(hours=48)
RFQ_TARGET_COST_DOLLARS = 100.0       # makers skip tiny RFQs; never accepted, so never spent
QUOTE_POLL_SEC = 0.5
QUOTE_WINDOW_MIN_SEC = 1.5            # wait at least this long for quotes
QUOTE_WINDOW_MAX_SEC = 4.0            # and stop early once the count holds steady
GAMES_CACHE_SEC = 300
COLLECTION_CACHE_SEC = 300
EVENTS_CACHE_SEC = 20                 # the leg picker; pricing always re-reads
EVENT_READ_WORKERS = 4                # ~33 events per game; the key's read bucket is shared with the MLB taker
MINT_DECLINE_STATUSES = (400, 422)    # Kalshi refuses legs that can't go together
GAMES_MAX_PAGES = 5
MIN_QUOTE_DECIMAL_ODDS = 1.01         # a 0.01c no_bid "quote" returns nothing; not an offer
DEFAULT_ENV_FILE = Path(__file__).resolve().parents[2] / "kalshi_mlb_rfq" / ".env"

# Signed Kalshi call: (method, path, body) -> (status, body, headers).
KalshiApi = Callable[..., tuple[int, Any, dict]]


@dataclass(frozen=True)
class SeriesKind:
    """How one Kalshi series reads as legs."""
    period: str
    kind: str          # "moneyline" | "win_3way" | "spread" | "total" | "team_total"
                       # | "player_over" | "anytime_td" | "player_yes_no" | "team_yes_no" | "other"
    category: str
    label: str         # market label; player props prefix the player's name
    stat: str | None = None


SERIES_KINDS: dict[str, SeriesKind] = {
    "KXNFLGAME": SeriesKind("FG", "moneyline", "Game Lines", "Moneyline", "moneyline"),
    "KXNFLSPREAD": SeriesKind("FG", "spread", "Game Lines", "Spread", "spread"),
    "KXNFLTOTAL": SeriesKind("FG", "total", "Game Lines", "Total", "total"),
    "KXNFLTEAMTOTAL": SeriesKind("FG", "team_total", "Team Totals", "Team Total", "team_total"),
    "KXNFLWINMARGIN": SeriesKind("FG", "other", "Game Lines", "Winning Margin"),
    "KXNFL1H": SeriesKind("1H", "win_3way", "1st Half", "1H Winner (tie loses)", "win_3way"),
    "KXNFL1HSPREAD": SeriesKind("1H", "spread", "1st Half", "1H Spread", "spread"),
    "KXNFL1HTOTAL": SeriesKind("1H", "total", "1st Half", "1H Total", "total"),
    "KXNFL1HTEAMTOTAL": SeriesKind("1H", "team_total", "1st Half", "1H Team Total", "team_total"),
    "KXNFL1Q": SeriesKind("1Q", "win_3way", "1st Quarter", "1Q Winner (tie loses)", "win_3way"),
    "KXNFL1QSPREAD": SeriesKind("1Q", "spread", "1st Quarter", "1Q Spread", "spread"),
    "KXNFL1QTOTAL": SeriesKind("1Q", "total", "1st Quarter", "1Q Total", "total"),
    "KXNFLTD": SeriesKind("FG", "anytime_td", "TD Scorers", "Touchdowns", "anytime_td"),
    "KXNFLFIRSTTD": SeriesKind("FG", "player_yes_no", "TD Scorers", "First TD", "first_td"),
    "KXNFLFIRSTTDTEAM": SeriesKind("FG", "team_yes_no", "TD Scorers", "First TD Team", "first_td_team"),
    "KXNFLTEAMFIRSTTD": SeriesKind("FG", "other", "TD Scorers", "Team First TD"),
    "KXNFLPASSYDS": SeriesKind("FG", "player_over", "Passing", "Passing Yards", "passing_yards"),
    "KXNFLPASSTDS": SeriesKind("FG", "player_over", "Passing", "Passing TDs", "passing_touchdowns"),
    "KXNFLPASSATT": SeriesKind("FG", "player_over", "Passing", "Pass Attempts", "passing_attempts"),
    "KXNFLPASSCOMP": SeriesKind("FG", "player_over", "Passing", "Completions", "passing_completions"),
    "KXNFLPASSINT": SeriesKind("FG", "player_over", "Passing", "Interceptions", "interceptions"),
    "KXNFLRSHYDS": SeriesKind("FG", "player_over", "Rushing", "Rushing Yards", "rushing_yards"),
    "KXNFLRSHATT": SeriesKind("FG", "player_over", "Rushing", "Rush Attempts", "rushing_attempts"),
    "KXNFLRECYDS": SeriesKind("FG", "player_over", "Receiving", "Receiving Yards", "receiving_yards"),
    "KXNFLREC": SeriesKind("FG", "player_over", "Receiving", "Receptions", "receptions"),
    "KXNFLRRYDS": SeriesKind("FG", "player_over", "Rushing", "Rush + Rec Yards", "rushing_and_receiving_yards"),
}
# "...PROT" series: the same prop, "protected" (a player who leaves at halftime
# settles at the pre-game fair price). A different bet, so its own stat name.
PROTECTED_SUFFIX = "PROT"


# ---------------------------------------------------------------- parsing (pure)

def all_in_price(price_dollars: float) -> float:
    """What one contract costs a taker, fee included."""
    return price_dollars + TAKER_FEE_RATE * price_dollars * (1.0 - price_dollars)


def cost_to_american(cost: float) -> float:
    """Cost per $1 payout (0 < cost < 1) -> American odds: 0.40 -> +150, 0.60 -> -150."""
    if not 0.0 < cost < 1.0:
        raise ValueError(f"expected a cost in (0, 1), got {cost!r}")
    return (1.0 - cost) / cost * 100.0 if cost <= 0.5 else -cost / (1.0 - cost) * 100.0


def _american_all_in(price_dollars: float | None) -> float | None:
    """American odds to buy at `price_dollars` with the fee, or None when
    nobody is offering (ask of 0 or 1) or the fee makes it a sure loss."""
    if price_dollars is None or not 0.0 < price_dollars < 1.0:
        return None
    cost = all_in_price(price_dollars)
    return cost_to_american(cost) if cost < 1.0 else None


def _dollars(value: Any) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def parse_nfl_games(events_payload: dict, now: datetime) -> list[Game]:
    """Not-yet-started NFL games from the KXNFLGAME event list, soonest first."""
    return [game for game, _ in parse_nfl_games_with_teams(events_payload, now)]


def parse_nfl_games_with_teams(events_payload: dict, now: datetime) -> list[tuple[Game, dict[str, str]]]:
    """Each not-yet-started game with its Kalshi team codes -> names
    ({"IND": "Indianapolis", "HOU": "Houston"}), soonest first."""
    games = []
    for event in events_payload.get("events") or []:
        read = _game_from_event(event)
        if read is not None and read[0].game_start_time > now:
            games.append(read)
    return sorted(games, key=lambda read: read[0].game_start_time)


def _game_from_event(event: dict) -> tuple[Game, dict[str, str]] | None:
    """None (logged) when the event does not read as away-code + home-code."""
    event_ticker = event.get("event_ticker") or ""
    game_key = event_ticker.split("-", 1)[1] if "-" in event_ticker else ""
    teams_part = game_key[7:]                     # "26NOV19" + "INDHOU"
    names = {}
    closes = []
    for market in event.get("markets") or []:
        code = (market.get("ticker") or "").rsplit("-", 1)[-1]
        names[code] = market.get("yes_sub_title") or code
        if market.get("close_time"):
            closes.append(datetime.fromisoformat(market["close_time"].replace("Z", "+00:00")))
    # Every pair, not the first prefix: "LACLA" is LAC at LA, and LA is also a prefix.
    pairs = [(away, home) for away in names for home in names if away != home and away + home == teams_part]
    if len(pairs) != 1 or not closes:
        log.warning("kalshi: %s does not read as away+home teams (%s); skipped", event_ticker, sorted(names))
        return None
    away, home = pairs[0]
    start = (min(closes) - GAME_CLOSE_AFTER_START).astimezone(timezone.utc)
    game = Game(venue_game_id=game_key, home_team=names[home], away_team=names[away], game_start_time=start)
    return game, {away: names[away], home: names[home]}


def combo_event_tickers(collection_payload: dict, game: Game) -> dict[str, int | None]:
    """The collection's NFL events for `game` -> their size_max (None = no
    per-event limit). Events missing here cannot go in a combo."""
    contract = collection_payload.get("multivariate_contract", collection_payload)
    events = {}
    for associated in contract.get("associated_events") or []:
        ticker = associated.get("ticker") or ""
        if ticker.startswith(LEAGUE_PREFIX) and f"-{game.venue_game_id}" in ticker:
            events[ticker] = associated.get("size_max")
    return events


def parse_legs(game: Game, team_names: dict[str, str], event_payloads: list[dict],
               size_max_of_event: dict[str, int | None]) -> list[OfferedLeg]:
    """Both sides of every active, priced market of the game's combo events.
    `team_names` is the game's Kalshi team codes -> names (from the game list)."""
    legs: list[OfferedLeg] = []
    for payload in event_payloads:
        event = payload.get("event") or {}
        event_ticker = event.get("event_ticker") or ""
        series = event.get("series_ticker") or event_ticker.split("-", 1)[0]
        markets = payload.get("markets") or event.get("markets") or []
        event_legs = []
        for market in markets:
            if market.get("status") == "active":
                event_legs.extend(_legs_of_market(game, series, event_ticker, market, team_names,
                                                  size_max_of_event.get(event_ticker)))
        legs.extend(_mark_main_lines(event_legs))
    return legs


def _side_of_team(game: Game, team_name: str) -> str | None:
    if team_name == game.home_team:
        return "home"
    if team_name == game.away_team:
        return "away"
    return None


def _legs_of_market(game: Game, series: str, event_ticker: str, market: dict,
                    team_names: dict[str, str], size_max: int | None) -> list[OfferedLeg]:
    ticker = market["ticker"]
    yes_ask = _dollars(market.get("yes_ask_dollars"))
    no_ask = _dollars(market.get("no_ask_dollars"))
    yes_odds, no_odds = _american_all_in(yes_ask), _american_all_in(no_ask)
    if yes_odds is None and no_odds is None:
        return []
    protected = series.endswith(PROTECTED_SUFFIX)
    kind = SERIES_KINDS.get(series[: -len(PROTECTED_SUFFIX)] if protected else series)
    reading = _read_market(game, kind, market, team_names, protected)
    strike = market.get("custom_strike") or {}
    group_id = strike.get("football_player") or strike.get("football_team") or ticker
    market_key = event_ticker if size_max == 1 else f"{event_ticker}|{group_id}"
    yes_size, no_size = _dollars(market.get("yes_ask_size_fp")), _dollars(market.get("yes_bid_size_fp"))
    sides = (
        ("yes", reading.yes_label, reading.yes_canonical, yes_odds, no_odds, yes_ask, yes_size),
        ("no", reading.no_label, reading.no_canonical, no_odds, yes_odds, no_ask, no_size),
    )
    legs = []
    for side, side_label, canonical, odds, other_odds, ask, size in sides:
        legs.append(OfferedLeg(
            leg_id=f"{ticker}|{side}",
            market_label=reading.market_label,
            category=reading.category,
            side_label=side_label,
            canonical=canonical,
            american_odds=odds,
            opposite_american_odds=other_odds,
            available_stake=size * all_in_price(ask) if size is not None and odds is not None else None,
            opposite_leg_id=f"{ticker}|{'no' if side == 'yes' else 'yes'}",
            market_key=market_key,
            is_main_line=False,
            venue_ref={"market_ticker": ticker, "event_ticker": event_ticker, "side": side,
                       "yes_ask_dollars": yes_ask, "group": f"{event_ticker}|{group_id}|{reading.market_label}"},
        ))
    return legs


@dataclass(frozen=True)
class MarketReading:
    market_label: str
    category: str
    yes_label: str
    no_label: str
    yes_canonical: CanonicalLeg | None
    no_canonical: CanonicalLeg | None


def _read_market(game: Game, kind: SeriesKind | None, market: dict, team_names: dict[str, str],
                 protected: bool) -> MarketReading:
    """Labels and canonical legs of one market's YES and NO."""
    title = market.get("title") or market["ticker"]
    yes_sub = market.get("yes_sub_title") or title
    line = market.get("floor_strike")
    line = float(line) if isinstance(line, (int, float)) else None
    # "-HOU", "-HOU14" (spread at 13.5), "-IND17" (team total at 16.5) -> the team code.
    team_code = market["ticker"].rsplit("-", 1)[-1].rstrip("0123456789")
    team = team_names.get(team_code)
    unread = MarketReading(kind.label if kind else title, kind.category if kind else "Other",
                           yes_sub, f"Not {yes_sub}", None, None)
    if kind is None or kind.kind == "other":
        return unread

    period, category, stat = kind.period, kind.category, kind.stat
    side = _side_of_team(game, team) if team else None
    other = {"home": "away", "away": "home"}.get(side)
    other_name = game.home_team if other == "home" else game.away_team
    if kind.kind == "moneyline":
        return MarketReading(kind.label, category, yes_sub, f"No on {yes_sub} (= {other_name})",
                             CanonicalLeg(period, "moneyline", None, side, None) if side else None,
                             CanonicalLeg(period, "moneyline", None, other, None) if other else None)
    if kind.kind == "win_3way":
        subject = "Tie" if team_code == "TIE" else team
        return MarketReading(kind.label, category, yes_sub, f"Not {yes_sub}",
                             CanonicalLeg(period, "win_3way", subject, "yes", None) if subject else None,
                             CanonicalLeg(period, "win_3way", subject, "no", None) if subject else None)
    if kind.kind == "spread" and line is not None:
        if side is None:
            return unread
        return MarketReading(kind.label, category, f"{team} {-line:+g}", f"{other_name} {line:+g}",
                             CanonicalLeg(period, "spread", None, side, -line),
                             CanonicalLeg(period, "spread", None, other, line))
    if kind.kind == "total" and line is not None:
        return MarketReading(kind.label, category, f"Over {line:g}", f"Under {line:g}",
                             CanonicalLeg(period, "total", None, "over", line),
                             CanonicalLeg(period, "total", None, "under", line))
    if kind.kind == "team_total" and line is not None:
        team = team or yes_sub.split(" over ")[0]
        return MarketReading(f"{team} {kind.label}", category, f"Over {line:g}", f"Under {line:g}",
                             CanonicalLeg(period, "team_total", team, "over", line),
                             CanonicalLeg(period, "team_total", team, "under", line))
    player = yes_sub.split(":")[0].strip()
    suffix = " (Protected)" if protected else ""
    stat_name = f"{stat}_protected" if protected else stat
    if kind.kind == "player_over" and line is not None:
        return MarketReading(f"{player} {kind.label}{suffix}", category + suffix, f"Over {line:g}",
                             f"Under {line:g}",
                             CanonicalLeg(period, stat_name, player, "over", line),
                             CanonicalLeg(period, stat_name, player, "under", line))
    if kind.kind == "anytime_td" and line is not None:
        if line == 0.5:
            stat_name = "anytime_td_protected" if protected else "anytime_td"
            return MarketReading(f"{player} Anytime TD{suffix}", category + suffix, "YES", "NO",
                                 CanonicalLeg(period, stat_name, player, "yes", None),
                                 CanonicalLeg(period, stat_name, player, "no", None))
        stat_name = "touchdowns_protected" if protected else "touchdowns"
        return MarketReading(f"{player} Touchdowns{suffix}", category + suffix, f"Over {line:g}",
                             f"Under {line:g}",
                             CanonicalLeg(period, stat_name, player, "over", line),
                             CanonicalLeg(period, stat_name, player, "under", line))
    if kind.kind in ("player_yes_no", "team_yes_no"):
        subject = player if kind.kind == "player_yes_no" else (team or yes_sub)
        return MarketReading(f"{subject} {kind.label}{suffix}", category + suffix, "YES", "NO",
                             CanonicalLeg(period, stat_name, subject, "yes", None),
                             CanonicalLeg(period, stat_name, subject, "no", None))
    return unread


def _mark_main_lines(legs: list[OfferedLeg]) -> list[OfferedLeg]:
    """Per ladder of numbers (one team's spreads, one player's yards...), the
    number whose YES ask is nearest 50c is the main line, both sides marked.
    A market with no number (moneyline, winner, yes/no prop) is always main."""
    best: dict[str, tuple[float, str]] = {}
    for leg in legs:
        ask = leg.venue_ref["yes_ask_dollars"]
        if leg.venue_ref["side"] != "yes" or ask is None or not _is_ladder(leg):
            continue
        ladder = _ladder_of(leg)
        distance = abs(ask - 0.5)
        if ladder not in best or distance < best[ladder][0]:
            best[ladder] = (distance, leg.venue_ref["market_ticker"])
    main_tickers = {ticker for _, ticker in best.values()}
    return [_with_main(leg, not _is_ladder(leg) or leg.venue_ref["market_ticker"] in main_tickers)
            for leg in legs]


def _is_ladder(leg: OfferedLeg) -> bool:
    return leg.canonical is not None and leg.canonical.line is not None


def _ladder_of(leg: OfferedLeg) -> str:
    """Spreads and totals ladder per event (both teams' spread numbers are one
    ladder); everything else per event + player/team + label."""
    if leg.canonical.stat in ("spread", "total"):
        return leg.venue_ref["event_ticker"]
    return leg.venue_ref["group"]


def _with_main(leg: OfferedLeg, is_main: bool) -> OfferedLeg:
    return leg if leg.is_main_line == is_main else replace(leg, is_main_line=is_main)


def selected_markets(legs: list[OfferedLeg]) -> list[dict]:
    """The combo's legs as Kalshi's mint endpoint takes them."""
    return [{"market_ticker": leg.venue_ref["market_ticker"], "event_ticker": leg.venue_ref["event_ticker"],
             "side": leg.venue_ref["side"]} for leg in legs]


def parse_quotes(quotes: list[dict], rfq_id: str, combo_ticker: str) -> tuple[ParlayTier, ...]:
    """Every open quote on our RFQ as a price to buy the combo YES, best first;
    quotes at the same price merge. A quote paying under MIN_QUOTE_DECIMAL_ODDS
    (YES cost 1 - no_bid, plus fee, near $1) is dropped. A tier's max_stake is
    what one quote fills: about the RFQ's $100, since makers size to it; one
    RFQ accepts one quote, so same-price quotes take the max, not the sum."""
    stake_at: dict[float, float] = {}
    for quote in quotes:
        if quote.get("status") != "open" or quote.get("rfq_id") != rfq_id or quote.get("market_ticker") != combo_ticker:
            continue
        no_bid = _dollars(quote.get("no_bid_dollars"))
        contracts = _dollars(quote.get("no_contracts_fp"))
        if no_bid is None or not 0.0 < no_bid < 1.0:
            continue
        cost = all_in_price(1.0 - no_bid)
        if cost >= 1.0 or 1.0 / cost < MIN_QUOTE_DECIMAL_ODDS:
            continue
        decimal_odds = round(1.0 / cost, 6)
        stake_at[decimal_odds] = max(stake_at.get(decimal_odds, 0.0), (contracts or 0.0) * cost)
    return tuple(ParlayTier(odds, round(stake, 2) or None) for odds, stake in sorted(stake_at.items(), reverse=True))


def _trim_quote(quote: dict) -> dict:
    """The quote fields worth storing; never the maker's or our account ids."""
    keys = ("id", "no_bid_dollars", "yes_bid_dollars", "no_contracts_fp", "yes_contracts_fp", "status", "created_ts")
    return {key: quote.get(key) for key in keys}


# ---------------------------------------------------------------- credentials

def _signed_api_from_env(env_file: Path = DEFAULT_ENV_FILE) -> tuple[KalshiApi, str]:
    """Kalshi's signed client on the key in the environment, else in
    kalshi_mlb_rfq/.env (the MLB taker's key, the one the 2026-11-19 recon
    used). Raises VenueTransportError naming what is missing."""
    repo_root = Path(__file__).resolve().parents[2]
    file_env = _read_env_file(env_file)
    key_id = os.environ.get("KALSHI_API_KEY_ID") or file_env.get("KALSHI_API_KEY_ID")
    key_path = os.environ.get("KALSHI_PRIVATE_KEY_PATH") or file_env.get("KALSHI_PRIVATE_KEY_PATH")
    user_id = os.environ.get("KALSHI_USER_ID") or file_env.get("KALSHI_USER_ID")
    missing = [name for name, value in (("KALSHI_API_KEY_ID", key_id), ("KALSHI_PRIVATE_KEY_PATH", key_path),
                                        ("KALSHI_USER_ID", user_id)) if not value]
    if missing:
        raise VenueTransportError(f"Kalshi credentials missing: set {', '.join(missing)} "
                                  "in the environment or kalshi_mlb_rfq/.env")
    from kalshi_common import auth_client
    auth_client.configure(key_id, key_path, os.environ.get("KALSHI_BASE_URL", KALSHI_BASE_URL), repo_root)
    return auth_client.api, user_id


def _read_env_file(path: Path) -> dict[str, str]:
    env: dict[str, str] = {}
    if not path.exists():
        return env
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line and not line.startswith("#") and "=" in line:
            key, value = line.split("=", 1)
            env[key.strip()] = value.strip().strip('"').strip("'")
    return env


# ---------------------------------------------------------------- the venue

class KalshiVenue:
    """Live Kalshi reads and RFQ price requests. Caches only the schedule,
    the collection and the leg picker's events, never a price used to compute
    a parlay."""

    name = "kalshi"
    display_name = "Kalshi"

    def __init__(self, api_factory: Callable[[], tuple[KalshiApi, str]] = _signed_api_from_env,
                 clock: Callable[[], float] = time.monotonic,
                 sleep: Callable[[float], None] = time.sleep) -> None:
        self._api_factory = api_factory
        self._api: KalshiApi | None = None
        self._user_id = ""
        self._clock = clock
        self._sleep = sleep
        self._lock = threading.Lock()
        self._games: tuple[float, list[tuple[Game, dict[str, str]]]] | None = None
        self._collection: tuple[float, dict] | None = None
        self._events: dict[str, tuple[float, list[dict]]] = {}

    def list_games(self) -> list[Game]:
        return [game for game, _ in self._read_games()]

    def list_legs(self, game: Game, fresh: bool = False) -> list[OfferedLeg]:
        team_names = next((teams for listed, teams in self._read_games()
                           if listed.venue_game_id == game.venue_game_id), None)
        if team_names is None:
            raise VenueTransportError(f"Kalshi: {game.label} is not an upcoming KXNFLGAME event")
        size_max_of_event = combo_event_tickers(self._read_collection(), game)
        payloads = self._read_events(game, sorted(size_max_of_event), fresh)
        return parse_legs(game, team_names, payloads, size_max_of_event)

    def price_parlay(self, game: Game, legs: list[OfferedLeg]) -> ParlayQuote | None:
        """Mint the combo, ask for quotes, read them, cancel. Never accepts."""
        api = self._signed_api()
        status, minted, _ = api("POST", f"/multivariate_event_collections/{COLLECTION_TICKER}",
                                body={"selected_markets": selected_markets(legs)})
        if status in MINT_DECLINE_STATUSES:
            log.info("kalshi: combo refused (%s): %s", status, str(minted)[:200])
            return None
        if status != 200 or not isinstance(minted, dict) or not minted.get("market_ticker"):
            raise VenueTransportError(f"Kalshi combo mint: expected 200 with market_ticker, got {status}: {str(minted)[:200]}")
        combo_ticker = minted["market_ticker"]

        try:
            status, created, _ = api("POST", "/communications/rfqs", body={
                "market_ticker": combo_ticker, "rest_remainder": False, "replace_existing": True,
                "target_cost_dollars": f"{RFQ_TARGET_COST_DOLLARS:.2f}"})
        except VenueTransportError:
            # Timeout or reset: the POST may have landed. Cancel what we hold on this combo.
            self._cancel_open_rfqs_on(api, combo_ticker)
            raise
        if status not in (200, 201) or not isinstance(created, dict) or not created.get("id"):
            # The combo minted, so this is not "these legs can't combine". A 5xx
            # or an answer without an id may still have created the RFQ; a 4xx
            # did not, and cleaning up then could cancel another click's live RFQ.
            if status >= 500 or status in (200, 201):
                self._cancel_open_rfqs_on(api, combo_ticker)
            raise VenueTransportError(f"Kalshi create RFQ: expected 201 with id, got {status}: {str(created)[:200]}")
        rfq_id = created["id"]

        try:
            quotes, polls, waited = self._collect_quotes(api, rfq_id)
        finally:
            delete_status = self._delete_rfq(api, rfq_id)
        tiers = parse_quotes(quotes, rfq_id, combo_ticker)
        if not tiers:
            # Raised, not None, so the page says why instead of "would not price".
            if quotes:
                raise VenueTransportError(f"Kalshi: {len(quotes)} quote(s) in {waited:.1f} s, none a usable price")
            raise VenueTransportError(f"Kalshi: no market maker quoted within {waited:.1f} s")
        raw = {"combo_market_ticker": combo_ticker, "rfq_id": rfq_id,
               "rfq_target_cost_dollars": RFQ_TARGET_COST_DOLLARS, "quote_polls": polls,
               "seconds_waited": round(waited, 2), "delete_status": delete_status,
               "selected_markets": selected_markets(legs), "quotes": [_trim_quote(q) for q in quotes]}
        return ParlayQuote(tiers, {}, raw)

    # -- RFQ steps
    def _collect_quotes(self, api: KalshiApi, rfq_id: str) -> tuple[list[dict], int, float]:
        """Poll until QUOTE_WINDOW_MAX_SEC, or sooner once past the minimum
        window with quotes in and the count unchanged since the last poll."""
        path = f"/communications/quotes?rfq_id={rfq_id}&rfq_creator_user_id={self._user_id}"
        start = self._clock()
        quotes: list[dict] = []
        polls = 0
        while True:
            self._sleep(QUOTE_POLL_SEC)
            status, body, _ = api("GET", path)
            if status != 200 or not isinstance(body, dict):
                raise VenueTransportError(f"Kalshi quotes for RFQ {rfq_id}: expected 200, got {status}: {str(body)[:200]}")
            previous_count = len(quotes)
            quotes = [q for q in body.get("quotes") or [] if q.get("rfq_id") == rfq_id and q.get("status") == "open"]
            polls += 1
            waited = self._clock() - start
            settled = quotes and len(quotes) == previous_count and waited >= QUOTE_WINDOW_MIN_SEC
            if settled or waited >= QUOTE_WINDOW_MAX_SEC:
                return quotes, polls, waited

    def _delete_rfq(self, api: KalshiApi, rfq_id: str) -> int | None:
        """Cancel our RFQ. A failure is logged loudly, not raised: the RFQ
        cannot fill without an accept this module never sends, and Kalshi
        expires it."""
        try:
            status, body, _ = api("DELETE", f"/communications/rfqs/{rfq_id}")
        except VenueTransportError as exc:
            log.error("kalshi: could not cancel RFQ %s: %s: %s", rfq_id, type(exc).__name__, exc)
            return None
        error = body.get("error") if isinstance(body, dict) else None
        expired = status == 400 and isinstance(error, dict) and error.get("code") == "expired"
        if status not in (200, 204) and not expired:
            log.error("kalshi: cancel RFQ %s: expected 204, got %s: %s", rfq_id, status, str(body)[:200])
        return status

    def _cancel_open_rfqs_on(self, api: KalshiApi, combo_ticker: str) -> None:
        """Delete every open RFQ of ours on `combo_ticker` (after an unclear create reply)."""
        try:
            status, body, _ = api("GET", f"/communications/rfqs?status=open&creator_user_id={self._user_id}&limit=100")
        except VenueTransportError as exc:
            log.error("kalshi: could not list open RFQs to clean up %s: %s", combo_ticker, exc)
            return
        if status != 200 or not isinstance(body, dict):
            log.error("kalshi: list open RFQs to clean up %s: expected 200, got %s", combo_ticker, status)
            return
        for rfq in body.get("rfqs") or []:
            if rfq.get("market_ticker") == combo_ticker and rfq.get("id"):
                self._delete_rfq(api, rfq["id"])

    # -- reads
    def _read_games(self) -> list[tuple[Game, dict[str, str]]]:
        """Upcoming games with their team codes. The cache is re-filtered on
        every read, so a game drops out at its start, not up to 5 min later."""
        with self._lock:
            cached = self._games
        if not cached or self._clock() - cached[0] >= GAMES_CACHE_SEC:
            events: list[dict] = []
            cursor = ""
            for _ in range(GAMES_MAX_PAGES):
                page = self._get(f"/events?series_ticker={GAME_SERIES}&status=open&with_nested_markets=true"
                                 f"&limit=200" + (f"&cursor={cursor}" if cursor else ""))
                events.extend(page.get("events") or [])
                cursor = page.get("cursor") or ""
                if not cursor:
                    break
            else:
                log.warning("kalshi: game list cut at %d pages (%d events); more remain", GAMES_MAX_PAGES, len(events))
            games = parse_nfl_games_with_teams({"events": events}, datetime.now(timezone.utc))
            cached = (self._clock(), games)
            with self._lock:
                self._games = cached
        now = datetime.now(timezone.utc)
        return [(game, teams) for game, teams in cached[1] if game.game_start_time > now]

    def _read_collection(self) -> dict:
        with self._lock:
            cached = self._collection
        if cached and self._clock() - cached[0] < COLLECTION_CACHE_SEC:
            return cached[1]
        payload = self._get(f"/multivariate_event_collections/{COLLECTION_TICKER}")
        with self._lock:
            self._collection = (self._clock(), payload)
        return payload

    def _read_events(self, game: Game, event_tickers: list[str], fresh: bool) -> list[dict]:
        with self._lock:
            cached = self._events.get(game.venue_game_id)
        if cached and not fresh and self._clock() - cached[0] < EVENTS_CACHE_SEC:
            return cached[1]
        with ThreadPoolExecutor(max_workers=EVENT_READ_WORKERS) as pool:
            results = list(pool.map(self._read_event, event_tickers))
        payloads = [payload for payload in results if payload is not None]
        with self._lock:
            self._events[game.venue_game_id] = (self._clock(), payloads)
        return payloads

    def _read_event(self, event_ticker: str) -> dict | None:
        """None when Kalshi has no such event (a collection entry with no
        market for this game, e.g. KXNFLTEAMFIRSTTD on 2026-11-19)."""
        status, body, _ = self._signed_api()("GET", f"/events/{event_ticker}?with_nested_markets=true")
        if status == 404:
            return None
        if status != 200 or not isinstance(body, dict):
            raise VenueTransportError(f"Kalshi GET event {event_ticker}: expected 200, got {status}: {str(body)[:200]}")
        return body

    def _get(self, path: str) -> dict:
        status, body, _ = self._signed_api()("GET", path)
        if status != 200 or not isinstance(body, dict):
            raise VenueTransportError(f"Kalshi GET {path.split('?')[0]}: expected 200, got {status}: {str(body)[:200]}")
        return body

    def _signed_api(self) -> KalshiApi:
        with self._lock:
            if self._api is None:
                self._api, self._user_id = self._api_factory()
            api = self._api

        def guarded(method: str, path: str, body: dict | None = None) -> tuple[int, Any, dict]:
            try:
                return api(method, path, body=body)
            except VenueTransportError:
                raise
            except Exception as exc:   # urllib / socket / signing errors
                raise VenueTransportError(f"Kalshi {method} {path.split('?')[0]}: {type(exc).__name__}: {exc}") from exc
        return guarded
