"""Rung -> surface rows: the two-way devig and its exclusion accounting.

Pure functions. No network, no DB, no config — the workers pass the envelope
in, which keeps this fully unit-testable and keeps the band a knob.

A RUNG is one line at one book (period, market_type, line) and carries both
its legs. It is the unit that is devigged once and the unit an exclusion is
counted in.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime

from kalshi_common.fair_value import two_way_fair
from kalshi_mlb_mm.leg_surface.store import SurfaceRow


@dataclass
class ExclusionCounts:
    """Why rungs did not make it onto the surface. Fixed, small vocabulary —
    each becomes a column in surface_refresh_log, so counting one is SQL."""
    crossed: int = 0          # raw implied sum < 1.0
    overround: int = 0        # sum outside the vig envelope
    one_sided: int = 0        # book posts only one side of the line
    unresolved: int = 0       # book does not offer the line at all
    game_unmatched: int = 0   # book does not list the game
    game_ambiguous: int = 0   # 2+ of the book's games matched; fail closed

    def add(self, other: "ExclusionCounts") -> None:
        self.crossed += other.crossed
        self.overround += other.overround
        self.one_sided += other.one_sided
        self.unresolved += other.unresolved
        self.game_unmatched += other.game_unmatched
        self.game_ambiguous += other.game_ambiguous

    def bump(self, reason: str) -> None:
        setattr(self, reason, getattr(self, reason) + 1)


@dataclass
class RungResult:
    rows: list = field(default_factory=list)
    reason: str | None = None     # None when the rung priced


def devig_rung(*, book: str, route: str, game, legs, decimals: dict,
               built_at: datetime, band_min: float,
               band_max: float) -> RungResult:
    """One rung's two legs -> up to two SurfaceRows, or an exclusion reason.

    ``legs`` are the rung's CanonicalLegs (normally both sides of one market).
    ``decimals`` maps a ``side`` to the decimal price the book posted; the
    caller fills in a side it learned from the book's OPPOSITE price rather
    than from a second resolved leg.

    With only one side priced the rung is excluded ``one_sided``. It is NEVER
    haircut by a vig fallback — that constant is calibrated for cancelling a
    COMBO's compounded margin, and applied to a lone leg it invents a number
    where declining costs nothing (#95: 4-5 books price a typical leg, so one
    dropping out never approaches quorum).
    """
    by_side = {leg.side: leg for leg in legs}
    priced_sides = [s for s in by_side if decimals.get(s) is not None]
    if not priced_sides:
        return RungResult(reason="unresolved")
    if len(priced_sides) < 2:
        return RungResult(reason="one_sided")

    side_a, side_b = priced_sides[0], priced_sides[1]
    dec_a, dec_b = decimals[side_a], decimals[side_b]
    fair_a, reason = two_way_fair(dec_a, dec_b, band_min=band_min,
                                  band_max=band_max)
    if fair_a is None:
        # 'bad_price' folds into 'unresolved': from the surface's point of
        # view a price it cannot use is a line the book did not usably post.
        return RungResult(reason=reason if reason in ("crossed", "overround")
                          else "unresolved")
    overround = 1.0 / dec_a + 1.0 / dec_b
    fairs = {side_a: fair_a, side_b: 1.0 - fair_a}
    opposite_of = {side_a: side_b, side_b: side_a}
    rows = []
    for side in (side_a, side_b):
        leg = by_side[side]
        rows.append(SurfaceRow(
            book=book, game_id=game.game_id,
            game_start_time=game.start_utc, period=leg.period,
            market_type=leg.market_type, line=leg.line, side=side,
            fair_prob=fairs[side], raw_decimal=decimals[side],
            raw_decimal_opp=decimals[opposite_of[side]],
            raw_overround=overround, route=route, built_at=built_at))
    return RungResult(rows=rows)
