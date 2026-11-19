"""The interface every venue adapter (ProphetX now; Novig, Kalshi, Polymarket,
sportsbooks next) fills in. Adding a venue = one module implementing `Venue`
plus one line in `sgp_compare.venues.build_venues()`."""
from __future__ import annotations

from typing import Protocol

from sgp_compare.legs import Game, OfferedLeg, ParlayQuote


class VenueTransportError(RuntimeError):
    """The venue could not be reached or answered garbage (as opposed to a
    clean "I won't price this", which is a None quote)."""


class Venue(Protocol):
    name: str            # "prophetx"
    display_name: str    # "ProphetX"

    def list_games(self) -> list[Game]:
        """Upcoming NFL games, soonest first. Raises VenueTransportError."""

    def list_legs(self, game: Game, fresh: bool = False) -> list[OfferedLeg]:
        """Every leg of `game` the venue can put in a parlay, with its single
        prices; `fresh` skips any cache (pricing always passes it). Raises
        VenueTransportError."""

    def price_parlay(self, game: Game, legs: list[OfferedLeg]) -> ParlayQuote | None:
        """The venue's parlay price for `legs` together, or None when it
        declines the combination. Never places a bet. Raises VenueTransportError."""
