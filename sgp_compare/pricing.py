"""Price one parlay at one venue and say how correlated the venue thinks it is.

`price_at_venue` is the only function that talks to a venue: it re-reads the
game's legs (fresh single prices), asks for the parlay price, and for 2-3 legs
also prices every other combination of the legs' outcomes so the vig can be
removed (stopping at the first combination that fails, so a blocked venue gets
one request, not seven). No DB I/O; the caller stores the result (store.py).
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
# Over/under and spread stats whose whole-number lines can push (a third
# outcome), so their two sides are not every outcome of the leg.
PUSHABLE_SIDES = ("over", "under", "home", "away")


@dataclass(frozen=True)
class VenueResult:
    venue: str
    display_name: str
    status: str                     # "priced" | "declined" | "missing_legs" | "error"
    detail: str | None
    requested_leg_ids: list[str]
    legs: list[OfferedLeg]
    quote: ParlayQuote | None
    partition_decimal_odds: list[float] | None   # partition_flips order; flipped legs = opposite_leg_id
    report: CorrelationReport | None

    def to_json(self) -> dict:
        return {
            "venue": self.venue,
            "display_name": self.display_name,
            "status": self.status,
            "detail": self.detail,
            "requested_leg_ids": self.requested_leg_ids,
            "legs": [_leg_json(leg, self.quote) for leg in self.legs],
            "tiers": [asdict(tier) for tier in self.quote.tiers] if self.quote else [],
            "partition_decimal_odds": self.partition_decimal_odds,
            "report": asdict(self.report) if self.report else None,
        }


def _leg_json(leg: OfferedLeg, quote: ParlayQuote | None) -> dict:
    return {
        "leg_id": leg.leg_id,
        "market_label": leg.market_label,
        "category": leg.category,
        "side_label": leg.side_label,
        "canonical": asdict(leg.canonical) if leg.canonical else None,
        "market_key": leg.market_key,
        "opposite_leg_id": leg.opposite_leg_id,
        "american_odds": leg.american_odds,
        "opposite_american_odds": leg.opposite_american_odds,
        "available_stake": leg.available_stake,
        "price_inside_parlay": quote.leg_prices_american.get(leg.leg_id) if quote else None,
        "venue_ref": leg.venue_ref,
    }


def validate_leg_set(leg_ids: list[str]) -> None:
    """Raise ValueError for a request no venue could take as one parlay."""
    if len(leg_ids) < MIN_PARLAY_LEGS:
        raise ValueError(f"expected >= {MIN_PARLAY_LEGS} legs, got {len(leg_ids)}")
    if len(set(leg_ids)) != len(leg_ids):
        raise ValueError("the same leg was picked twice")


def price_at_venue(venue: Venue, game: Game, leg_ids: list[str]) -> VenueResult:
    """Never raises for a venue problem: a failure is status 'error', with the reason."""
    validate_leg_set(leg_ids)
    try:
        offered = {leg.leg_id: leg for leg in venue.list_legs(game, fresh=True)}
    except VenueTransportError as exc:
        return _result(venue, leg_ids, "error", str(exc))
    missing = [leg_id for leg_id in leg_ids if leg_id not in offered]
    if missing:
        return _result(venue, leg_ids, "missing_legs",
                       f"{len(missing)} of {len(leg_ids)} leg(s) no longer offered (line moved or pulled); "
                       "reload the game and pick again")
    legs = [offered[leg_id] for leg_id in leg_ids]
    if len({leg.market_key for leg in legs}) != len(legs):
        raise ValueError("two legs come from the same market; pick one line per market")

    try:
        quote = venue.price_parlay(game, legs)
    except VenueTransportError as exc:
        return _result(venue, leg_ids, "error", str(exc), legs)
    if quote is None:
        return _result(venue, leg_ids, "declined", f"{venue.display_name} would not price this combination", legs)

    partition, partition_reason = _price_partition(venue, game, legs, offered, quote)
    chances = [leg_chance(leg.leg_id, _decimal(leg.american_odds), _decimal(leg.opposite_american_odds))
               for leg in legs]
    report = correlation_report(quote.best.decimal_odds, chances, partition, partition_reason)
    return VenueResult(venue.name, venue.display_name, "priced", None, list(leg_ids), legs, quote, partition, report)


def _price_partition(venue: Venue, game: Game, legs: list[OfferedLeg], offered: dict[str, OfferedLeg],
                     quote: ParlayQuote) -> tuple[list[float] | None, str | None]:
    """Best decimal price of every outcome combination (index 0 = the parlay),
    or (None, why not). Stops at the first combination that fails."""
    if len(legs) > MAX_PARTITION_LEGS:
        return None, f"vig removal prices every outcome combination; only done up to {MAX_PARTITION_LEGS} legs"
    one_sided = [leg for leg in legs if leg.opposite_leg_id not in offered]
    if one_sided:
        return None, f"{one_sided[0].side_label} has no other side at {venue.display_name} to flip to"
    pushable = [leg for leg in legs if _can_push(leg)]
    if pushable:
        return None, (f"{pushable[0].market_label} {pushable[0].side_label} can push or tie, "
                      "so its two sides are not every outcome")
    prices = [quote.best.decimal_odds]
    for flips in partition_flips(len(legs))[1:]:
        cell = [offered[leg.opposite_leg_id] if flip else leg for leg, flip in zip(legs, flips)]
        label = " + ".join(leg.side_label for leg in cell)
        try:
            cell_quote = venue.price_parlay(game, cell)
        except VenueTransportError as exc:
            log.warning("%s: outcome combination %s failed: %s", venue.name, label, exc)
            return None, f"pricing the combination {label} failed: {exc}"[:300]
        if cell_quote is None:
            return None, f"{venue.display_name} would not price the combination {label}"
        prices.append(cell_quote.best.decimal_odds)
    return prices, None


def _can_push(leg: OfferedLeg) -> bool:
    """A whole-number line can land exactly; a half or quarter moneyline can be
    tied (a full-game NFL tie is rare enough to ignore)."""
    canonical = leg.canonical
    if canonical is None or canonical.side not in PUSHABLE_SIDES:
        return False
    if canonical.stat == "moneyline":
        return canonical.period != "FG"
    return canonical.line is not None and float(canonical.line).is_integer()


def _decimal(american_odds: float | None) -> float | None:
    return american_to_decimal(american_odds) if american_odds is not None else None


def _result(venue: Venue, leg_ids: list[str], status: str, detail: str,
            legs: list[OfferedLeg] | None = None) -> VenueResult:
    return VenueResult(venue.name, venue.display_name, status, detail, list(leg_ids), legs or [], None, None, None)
