"""Price one parlay at one venue and say how correlated the venue thinks it is.

`price_at_venue` is the only function that talks to a venue: it re-reads the
game's legs (fresh single prices), asks for the parlay price, and for 2-3 legs
also prices every other combination of the legs' outcomes so the vig can be
removed. No DB I/O; the caller stores the result (store.py).
"""
from __future__ import annotations

import logging
from dataclasses import asdict, dataclass

from sgp_compare.correlation import (MAX_PARTITION_LEGS, CorrelationReport, correlation_report,
                                     leg_chance, partition_flips)
from sgp_compare.legs import Game, OfferedLeg, ParlayQuote
from sgp_compare.odds import american_to_decimal
from sgp_compare.venues.base import Venue, VenueTransportError

log = logging.getLogger(__name__)

MIN_PARLAY_LEGS = 2


@dataclass(frozen=True)
class VenueResult:
    venue: str
    display_name: str
    status: str                     # "priced" | "declined" | "missing_legs" | "error"
    detail: str | None
    legs: list[OfferedLeg]
    quote: ParlayQuote | None
    partition_decimal_odds: list[float | None] | None
    report: CorrelationReport | None

    def to_json(self) -> dict:
        return {
            "venue": self.venue,
            "display_name": self.display_name,
            "status": self.status,
            "detail": self.detail,
            "legs": [_leg_json(leg, self.quote) for leg in self.legs],
            "tiers": [asdict(tier) for tier in self.quote.tiers] if self.quote else [],
            "partition_decimal_odds": self.partition_decimal_odds,
            "report": asdict(self.report) if self.report else None,
        }


def _leg_json(leg: OfferedLeg, quote: ParlayQuote | None) -> dict:
    return {
        "leg_id": leg.leg_id,
        "market_label": leg.market_label,
        "side_label": leg.side_label,
        "american_odds": leg.american_odds,
        "opposite_american_odds": leg.opposite_american_odds,
        "available_stake": leg.available_stake,
        "price_inside_parlay": quote.leg_prices_american.get(leg.leg_id) if quote else None,
    }


def validate_leg_set(legs: list[OfferedLeg]) -> None:
    """Raise ValueError for a set no venue would take as one parlay."""
    if len(legs) < MIN_PARLAY_LEGS:
        raise ValueError(f"expected >= {MIN_PARLAY_LEGS} legs, got {len(legs)}")
    market_keys = [leg.market_key for leg in legs]
    if len(set(market_keys)) != len(market_keys):
        raise ValueError("two legs come from the same market; pick one line per market")


def price_at_venue(venue: Venue, game: Game, leg_ids: list[str]) -> VenueResult:
    """Never raises for a venue problem: a transport failure is status 'error'."""
    try:
        offered = {leg.leg_id: leg for leg in venue.list_legs(game, fresh=True)}
    except VenueTransportError as exc:
        return _result(venue, "error", str(exc))
    missing = [leg_id for leg_id in leg_ids if leg_id not in offered]
    if missing:
        return _result(venue, "missing_legs", f"{len(missing)} leg(s) no longer offered (line moved or pulled)")
    legs = [offered[leg_id] for leg_id in leg_ids]
    validate_leg_set(legs)

    try:
        quote = venue.price_parlay(game, legs)
    except VenueTransportError as exc:
        return _result(venue, "error", str(exc), legs)
    if quote is None:
        return _result(venue, "declined", f"{venue.display_name} would not price this combination", legs)

    partition = _price_partition(venue, game, legs, offered, quote)
    chances = [leg_chance(leg.leg_id, _decimal(leg.american_odds), _decimal(leg.opposite_american_odds))
               for leg in legs]
    report = correlation_report(quote.best.decimal_odds, chances, partition)
    return VenueResult(venue.name, venue.display_name, "priced", None, legs, quote, partition, report)


def _price_partition(venue: Venue, game: Game, legs: list[OfferedLeg], offered: dict[str, OfferedLeg],
                     quote: ParlayQuote) -> list[float | None] | None:
    """Best decimal price of every outcome combination (index 0 = the parlay),
    or None when the partition is out of reach (too many legs, a leg with no
    other side). A combination the venue declines or fails on is a None entry."""
    if len(legs) > MAX_PARTITION_LEGS:
        return None
    if any(leg.opposite_leg_id not in offered for leg in legs):
        return None
    prices: list[float | None] = [quote.best.decimal_odds]
    for flips in partition_flips(len(legs))[1:]:
        cell = [offered[leg.opposite_leg_id] if flip else leg for leg, flip in zip(legs, flips)]
        try:
            cell_quote = venue.price_parlay(game, cell)
        except VenueTransportError as exc:
            log.warning("%s: partition cell failed: %s", venue.name, exc)
            cell_quote = None
        prices.append(cell_quote.best.decimal_odds if cell_quote else None)
    return prices


def _decimal(american_odds: float | None) -> float | None:
    return american_to_decimal(american_odds) if american_odds is not None else None


def _result(venue: Venue, status: str, detail: str, legs: list[OfferedLeg] | None = None) -> VenueResult:
    return VenueResult(venue.name, venue.display_name, status, detail, legs or [], None, None, None)
