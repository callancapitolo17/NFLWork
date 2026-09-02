"""Discover open KXMLBRFI markets and turn them into quotable games.

Inputs: Kalshi REST (one GET per refresh). Outputs: list[RfiGame].
Side effects: none (no DB writes; auth_client must already be configured).

KXMLBRFI grammar (live-verified 2026-08-13, issue #87): one binary market
per game, market ticker == event ticker == KXMLBRFI-{suffix} where suffix is
YYMMMDDHHMM{AwayCode}{HomeCode} and HHMM is the first pitch in US/Eastern
(kalshi_close_time_not_start: close_time is first pitch + 72h — never use it
as the game start).
"""
import logging
from dataclasses import dataclass
from datetime import datetime, timezone

from kalshi_common import auth_client
from kalshi_common.leg_types import (_ET, _MLB_CODE_TO_TEAM,
                                     _parse_event_suffix,
                                     parse_suffix_start_utc)

log = logging.getLogger(__name__)

# parse_suffix_start_utc moved to kalshi_common.leg_types (the maker's leg
# surface keys games on the same suffix). Re-exported so this module's public
# surface is unchanged.
__all__ = ["RfiGame", "parse_suffix_start_utc", "parse_market",
           "drop_doubleheaders", "fetch_open_rfi_games"]


@dataclass(frozen=True)
class RfiGame:
    ticker: str            # market ticker, e.g. KXMLBRFI-26AUG262105MINATH
    suffix: str            # 26AUG262105MINATH — the family-independent game code
    home_team: str         # Odds-API-style full name ("Minnesota Twins")
    away_team: str
    commence_utc: datetime  # naive UTC (GameRef convention across the bots)
    yes_bid_cents: int | None
    yes_ask_cents: int | None
    status: str
    exchange_index: int | None = None   # Kalshi shard; None = let it auto-route


def _cents(market: dict, key: str) -> int | None:
    """Price field -> integer cents. The v2 API serves string decimal-dollar
    fields (`yes_bid_dollars: "0.4700"`, live-verified 2026-08-25); older
    payloads carried bare integer cents (`yes_bid: 47`) — accept both."""
    v = market.get(f"{key}_dollars")
    if v is not None:
        try:
            return round(float(v) * 100)
        except (TypeError, ValueError):
            return None
    v = market.get(key)
    try:
        return int(v) if v is not None else None
    except (TypeError, ValueError):
        return None


def _exchange_index(market: dict) -> int | None:
    """Kalshi's exchange shard for this market (sharding announced
    2026-08-24: baseball is 3, NFL/NBA are still 0). Read it, never assume
    it — order cancels must target the right shard or they 404 (bug 2)."""
    v = market.get("exchange_index")
    try:
        return int(v) if v is not None else None
    except (TypeError, ValueError):
        return None


def parse_market(market: dict) -> RfiGame | None:
    """One /markets row -> RfiGame, or None if unparseable."""
    ticker = str(market.get("ticker", ""))
    if not ticker.startswith("KXMLBRFI-"):
        return None
    suffix = ticker[len("KXMLBRFI-"):]
    away_code, home_code = _parse_event_suffix(suffix)
    if away_code is None or home_code is None:
        return None
    home_team = _MLB_CODE_TO_TEAM.get(home_code)
    away_team = _MLB_CODE_TO_TEAM.get(away_code)
    commence = parse_suffix_start_utc(suffix)
    if not home_team or not away_team or commence is None:
        return None
    return RfiGame(ticker=ticker, suffix=suffix,
                   home_team=home_team, away_team=away_team,
                   commence_utc=commence,
                   yes_bid_cents=_cents(market, "yes_bid"),
                   yes_ask_cents=_cents(market, "yes_ask"),
                   status=str(market.get("status", "")),
                   exchange_index=_exchange_index(market))


def drop_doubleheaders(games: list[RfiGame]) -> list[RfiGame]:
    """Fail-closed doubleheader exclusion (unabated_edge MLB precedent).

    Two open RFI markets for the same team pair on the same ET date mean a
    doubleheader; the books' event matchers key on team names and can pick
    the wrong game of the pair, so neither game is quoted.
    """
    def day_key(g: RfiGame):
        et_date = (g.commence_utc.replace(tzinfo=timezone.utc)
                   .astimezone(_ET).date())
        return (g.home_team, g.away_team, et_date)

    counts: dict = {}
    for g in games:
        counts[day_key(g)] = counts.get(day_key(g), 0) + 1
    kept = [g for g in games if counts[day_key(g)] == 1]
    dropped = [g.ticker for g in games if counts[day_key(g)] > 1]
    if dropped:
        log.info("discovery: doubleheader fail-closed, dropped %s", dropped)
    return kept


def fetch_open_rfi_games() -> list[RfiGame] | None:
    """All open KXMLBRFI markets as RfiGames (doubleheaders excluded).

    Returns None on an API failure — the caller must skip the cycle, NOT
    treat it as an empty board (an empty list would read as "every market
    vanished" and mass-cancel resting quotes as market_gone).

    limit=1000 is the API max and far above a 3-day listing window
    (~45 markets), so no cursor pagination is needed.
    """
    status, body, _ = auth_client.api(
        "GET", "/markets?series_ticker=KXMLBRFI&status=open&limit=1000")
    if status != 200 or not isinstance(body, dict):
        log.warning("discovery: markets fetch failed status=%s", status)
        return None
    games = []
    for m in body.get("markets", []):
        g = parse_market(m)
        if g is not None:
            games.append(g)
    return drop_doubleheaders(games)
