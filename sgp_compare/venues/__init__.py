"""Venue adapters. Each implements sgp_compare.venues.base.Venue."""
from __future__ import annotations

from sgp_compare.venues.base import Venue
from sgp_compare.venues.novig import NovigVenue
from sgp_compare.venues.prophetx import ProphetXVenue


def build_venues() -> dict[str, Venue]:
    """Every venue the page prices at, keyed by Venue.name."""
    venues: list[Venue] = [ProphetXVenue(), NovigVenue()]
    return {venue.name: venue for venue in venues}
