"""Kalshi slate + leg-ladder discovery for the leg surface (issue #96).

Inputs: Kalshi REST (auth_client must already be configured).
Outputs: list[SurfaceGame] — one per in-window game, each carrying every
CanonicalLeg the surface will price.
Side effects: NONE (no DB writes, no BOOK requests — Kalshi API only).

Why the surface enumerates its own slate instead of reading
``mlb_target_lines``: that table's game_id is the ODDS API id, resolved by
``WHERE home_team=? AND away_team=? LIMIT 1``. Team names alone collapse a
doubleheader onto one row and pick whichever — a wrong number, not a decline
(#95 hit the identical trap on FanDuel's two PHI @ SEA rows). The Kalshi
event-ticker suffix encodes date, ET first pitch and both team codes, is
unique per doubleheader game, and is exactly what ``CanonicalLeg.game_id``
carries, so #98's lookup needs no translation layer.

Legs are built by handing synthetic ``{event_ticker, market_ticker, side}``
dicts to ``legset.parse_leg`` — the SAME function the quote path parses real
RFQ legs with. The surface therefore cannot disagree with the router about
what a leg means: F5-winner -> +-0.5 F5 spread (#86), RFI -> I1 total at 0.5
(#87), spread signs home-perspective (#70).
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import datetime, timezone

from kalshi_common import auth_client, legset
from kalshi_common.leg_types import (_MLB_CODE_TO_TEAM, _parse_event_suffix,
                                     parse_suffix_start_utc)

log = logging.getLogger(__name__)

# Every Kalshi MLB series a combo leg can come from. KXMLBF5's TIE market is
# fetched and parsed but drops out below: it types as (ml, F5), which no book
# leg can express, and classify_subcombo already keeps such combos
# unpriceable.
LEG_SERIES = ("KXMLBGAME", "KXMLBSPREAD", "KXMLBTOTAL",
              "KXMLBF5", "KXMLBF5SPREAD", "KXMLBF5TOTAL", "KXMLBRFI")


@dataclass(frozen=True)
class SurfaceGame:
    """One in-window game and every leg the surface prices for it."""
    game_id: str              # Kalshi event-ticker suffix, e.g. 26AUG252138CLELAA
    home_team: str            # Odds-API canonical name
    away_team: str
    start_utc: datetime       # naive UTC first pitch (GameRef convention)
    legs: tuple               # tuple[CanonicalLeg, ...]

    def game_ref(self):
        """The GameRef the per-book event matchers consume. commence_time is
        the SUFFIX-derived first pitch, which is what makes doubleheaders
        resolve correctly — the books' match_events helpers already bucket on
        UTC date+hour, they were just being fed a team-resolved time."""
        from mlb_sgp._shared import GameRef
        return GameRef(game_id=self.game_id, home_team=self.home_team,
                       away_team=self.away_team, commence_time=self.start_utc)


def _fetch_open_game_events() -> list[str]:
    """Open KXMLBGAME event tickers, or [] on any API failure."""
    status, body, _ = auth_client.api(
        "GET", "/events?series_ticker=KXMLBGAME&status=open&limit=50")
    if status != 200 or not isinstance(body, dict):
        log.warning("surface slate: events fetch failed status=%s", status)
        return []
    return [str(e.get("event_ticker", "")) for e in body.get("events", [])
            if str(e.get("event_ticker", "")).startswith("KXMLBGAME-")]


def _fetch_series_tickers(series: str, suffix: str) -> list[str]:
    """Market tickers listed under one series for one game, or []."""
    event_ticker = f"{series}-{suffix}"
    status, body, _ = auth_client.api(
        "GET", f"/markets?event_ticker={event_ticker}&limit=100")
    if status != 200 or not isinstance(body, dict):
        return []
    return [str(m.get("ticker", "")) for m in body.get("markets", [])
            if m.get("ticker")]


def legs_for_market(series: str, suffix: str, market_ticker: str) -> list:
    """Both sides of one Kalshi market as CanonicalLegs (0, 1 or 2 of them).

    Unparseable sides and the F5 TIE market drop out here rather than
    downstream: a leg the surface cannot type is a leg the router cannot look
    up, so carrying it would only inflate the miss counts.
    """
    event_ticker = f"{series}-{suffix}"
    out = []
    for side in ("yes", "no"):
        leg = legset.parse_leg({"event_ticker": event_ticker,
                                "market_ticker": market_ticker,
                                "side": side})
        if leg is None:
            continue
        if leg.market_type == "ml" and leg.period == "F5":
            continue          # KXMLBF5 TIE — no book leg expresses it
        out.append(leg)
    return out


def _dedupe_legs(legs: list) -> tuple:
    """Distinct legs in a stable order. Two Kalshi markets CAN produce the
    same canonical leg (a team's -0.5 F5 market and the other team's, seen
    from opposite sides), and the surface stores one row per leg."""
    seen, out = set(), []
    for leg in legset.canonical_legs(legs):
        k = (leg.game_id, leg.period, leg.market_type, leg.line, leg.side)
        if k in seen:
            continue
        seen.add(k)
        out.append(leg)
    return tuple(out)


def in_window(start_utc: datetime, now_utc: datetime, *,
              min_minutes: float, max_hours: float) -> bool:
    """True if first pitch is far enough out to quote and near enough to be
    worth fetching. Games inside min_minutes are already being cancelled by
    the tipoff gate; games past max_hours have book ladders that are mostly
    main-line-only (#95: DK went ~14 spread lines/game -> exactly 1 at T-15h,
    and BetMGM posted no F5 or YRFI markets at all for next-day games)."""
    delta_sec = (start_utc - now_utc).total_seconds()
    return (min_minutes * 60.0) <= delta_sec <= (max_hours * 3600.0)


def discover_slate(*, min_minutes: float, max_hours: float,
                   now_utc: datetime | None = None) -> list[SurfaceGame]:
    """Every in-window game with its full leg ladder. Kalshi API only.

    Cost: 1 events call + ``len(LEG_SERIES)`` market-list calls per in-window
    game (~6 x 15 = 90 Kalshi calls per refresh on a full slate). ZERO book
    requests.

    Fail-safe: a game whose ladder cannot be read is skipped, not raised —
    the surface simply carries no rows for it and the router declines.
    """
    now = now_utc or datetime.now(timezone.utc).replace(tzinfo=None)
    games: list[SurfaceGame] = []
    skipped_window = 0
    for event_ticker in _fetch_open_game_events():
        suffix = event_ticker[len("KXMLBGAME-"):]
        away_code, home_code = _parse_event_suffix(suffix)
        home_team = _MLB_CODE_TO_TEAM.get(home_code or "")
        away_team = _MLB_CODE_TO_TEAM.get(away_code or "")
        start_utc = parse_suffix_start_utc(suffix)
        if not (home_team and away_team and start_utc):
            log.warning("surface slate: unparseable event %s", event_ticker)
            continue
        if not in_window(start_utc, now, min_minutes=min_minutes,
                         max_hours=max_hours):
            skipped_window += 1
            continue
        legs = []
        for series in LEG_SERIES:
            for market_ticker in _fetch_series_tickers(series, suffix):
                legs.extend(legs_for_market(series, suffix, market_ticker))
        if not legs:
            log.warning("surface slate: %s listed no priceable legs", suffix)
            continue
        games.append(SurfaceGame(game_id=suffix, home_team=home_team,
                                 away_team=away_team, start_utc=start_utc,
                                 legs=_dedupe_legs(legs)))
    log.info("surface slate: %d games in window (%d out), %d legs total",
             len(games), skipped_window, sum(len(g.legs) for g in games))
    return games


def rungs(legs) -> dict:
    """Group legs into two-way rungs: (period, market_type, line) -> legs.

    A rung is what gets devigged once and produces both of its legs' rows, so
    it is also the unit the exclusion counts are expressed in.
    """
    out: dict = {}
    for leg in legs:
        out.setdefault((leg.period, leg.market_type, leg.line), []).append(leg)
    return out
