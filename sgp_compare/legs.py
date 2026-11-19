"""The shapes every venue adapter speaks: a game, a leg it offers, a parlay quote.

Pure data, no network. A leg is described twice:
  * `CanonicalLeg` — venue-neutral (period, stat, subject, side, line), so the
    same bet can later be matched across venues (Novig, Kalshi, sportsbooks).
  * `OfferedLeg` — what one venue offers: its own id for the leg, the single
    price to take it, the price of the other side, and the venue's opaque
    reference needed to put it in a parlay request.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any


@dataclass(frozen=True)
class Game:
    venue_game_id: str
    home_team: str
    away_team: str
    game_start_time: datetime      # timezone-aware UTC

    @property
    def label(self) -> str:
        return f"{self.away_team} at {self.home_team}"


@dataclass(frozen=True)
class CanonicalLeg:
    period: str          # "FG" (full game, overtime included) | "1H" | "2H" | "1Q" | "2Q" | "3Q" | "4Q"
    stat: str            # "moneyline" | "spread" | "total" | "team_total" | a player stat
                         # such as "receiving_yards" | a yes/no prop such as "anytime_td"
    subject: str | None  # team name or player name; None for a game-wide market
    side: str            # "home"/"away" (moneyline, spread), "over"/"under", "yes"/"no"
    line: float | None   # signed from `side`'s view for spreads; None when there is no number


@dataclass(frozen=True)
class OfferedLeg:
    """One side of one market line at one venue.

    `american_odds` is the price to take this side alone right now (None when
    nobody is offering it); `opposite_american_odds` is the other side's price,
    needed to remove the vig. `available_stake` is the dollars on offer at that
    price. `opposite_leg_id` names the other side, so a parlay can flip a leg
    to price every combination of outcomes.
    """
    leg_id: str
    market_label: str                 # "Nico Collins Receiving Yards"
    category: str                     # the venue's grouping, e.g. "Rec Yards", "Game Lines"
    side_label: str                   # "Over 70.5", "Houston Texans -12.5", "YES"
    canonical: CanonicalLeg | None    # None when the venue's market is not understood
    american_odds: float | None
    opposite_american_odds: float | None
    available_stake: float | None
    opposite_leg_id: str | None
    market_key: str                   # legs sharing a market_key cannot go in one parlay
    is_main_line: bool
    venue_ref: dict[str, Any] = field(default_factory=dict, compare=False, hash=False)


@dataclass(frozen=True)
class ParlayTier:
    """One price a venue offers for the whole parlay, good up to `max_stake` dollars."""
    decimal_odds: float
    max_stake: float | None


@dataclass(frozen=True)
class ParlayQuote:
    """A venue's answer to "what would you pay on these legs together?".

    `tiers` is best price first. `leg_prices_american` holds any per-leg
    prices the venue reports inside the parlay (ProphetX shades each leg),
    keyed by leg_id. `raw` is the venue's response, stored for later study.
    """
    tiers: tuple[ParlayTier, ...]
    leg_prices_american: dict[str, float]
    raw: dict[str, Any] = field(compare=False, hash=False)

    @property
    def best(self) -> ParlayTier:
        return self.tiers[0]
