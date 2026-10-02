"""Fair value, edge and stake for one Wagerzon fecta.

Fair value at one book = probit devig over the FULL partition of the legs'
outcomes. Every combination of outcomes (SEA/Tie/LAC 1Q x SEA/Tie/LAC 1H x
SEA/LAC game, ...) is priced as an SGP at the book; the implied
probabilities sum to the book's SGP overround (~1.28 for 2 legs per the MLB
maker's measurements, far above the ~1.10 compounded single-leg vig), and the
devig spreads that margin across the cells. Pricing only the special's own
SGP and dividing out single-leg vig would leave most of that margin in and
overstate every edge.

Consensus = mean of the books that priced the full partition. Stake = Kelly
fraction x full Kelly at Wagerzon's price; within one game only the best
special per team gets a stake, because a team's fectas win together.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

from kalshi_common.fair_value import _probit_devig_n as probit_devig_n
from nfl_specials.books import BookGame, SgpBook
from nfl_specials.special_parser import Fecta

# A partition's implied sum must sit in [1, 1 + this x legs]; outside it a
# cell is stale or mispriced and the devig would be garbage (same envelope as
# kalshi_common.fair_value.PARTITION_OVERROUND_PER_LEG).
PARTITION_OVERROUND_PER_LEG = 0.25


@dataclass(frozen=True)
class BookFair:
    book: str
    fair_prob: float | None
    sgp_decimal: float | None      # the book's own price for the special, when it is one cell
    overround: float | None
    n_cells: int
    reason: str | None = None      # why fair_prob is None


def american_to_decimal(american: int) -> float:
    if american > 0:
        return 1.0 + american / 100.0
    return 1.0 + 100.0 / abs(american)


def decimal_to_american(decimal: float) -> int:
    if decimal >= 2.0:
        return int(round((decimal - 1.0) * 100))
    return int(round(-100.0 / (decimal - 1.0)))


def price_fecta_at_book(book: SgpBook, game: BookGame, role: str, fecta: Fecta) -> BookFair:
    """Price every cell of the fecta's partition at `book` (one SGP call each)."""
    leg_markets = []
    for leg in fecta.legs:
        market = book.leg_market(game, role, leg)
        if market is None:
            return BookFair(book.name, None, None, None, 0,
                            reason=f"no SGP market for '{leg.describe(fecta.team)}'")
        leg_markets.append(market)

    cells = list(itertools.product(*(m.group for m in leg_markets)))
    winning_cells = [index for index, cell in enumerate(cells)
                     if all(outcome in market.winning for outcome, market in zip(cell, leg_markets))]
    implied = []
    for cell in cells:
        decimal = book.price(game, tuple(outcome.ref for outcome in cell))
        if decimal is None:
            labels = " + ".join(outcome.label for outcome in cell)
            return BookFair(book.name, None, None, None, len(cells), reason=f"declined cell: {labels}")
        implied.append(1.0 / decimal)

    # The book's own price for the special exists only when it is one cell
    # (a "+0.5 off a 3-way market" leg makes it a union of cells).
    sgp_decimal = 1.0 / implied[winning_cells[0]] if len(winning_cells) == 1 else None
    overround = sum(implied)
    max_overround = 1.0 + PARTITION_OVERROUND_PER_LEG * len(fecta.legs)
    if not 1.0 <= overround <= max_overround:
        return BookFair(book.name, None, sgp_decimal, overround, len(cells),
                        reason=f"partition overround {overround:.3f} outside [1, {max_overround:.2f}]")
    devigged = probit_devig_n(implied)
    fair = sum(devigged[index] for index in winning_cells)
    return BookFair(book.name, float(fair), sgp_decimal, overround, len(cells))


def consensus_fair(book_fairs: list[BookFair]) -> float | None:
    priced = [b.fair_prob for b in book_fairs if b.fair_prob is not None]
    if not priced:
        return None
    return sum(priced) / len(priced)


def expected_value(fair_prob: float, wz_american: int) -> float:
    """EV per $1 risked at Wagerzon's price."""
    return fair_prob * american_to_decimal(wz_american) - 1.0


def kelly_stake(fair_prob: float, wz_american: int, bankroll: float, kelly_fraction: float) -> float:
    net_odds = american_to_decimal(wz_american) - 1.0
    full_kelly = expected_value(fair_prob, wz_american) / net_odds
    return max(0.0, full_kelly) * kelly_fraction * bankroll


def log_growth(fair_prob: float, wz_american: int, stake: float, bankroll: float) -> float:
    """Expected log bankroll growth of one bet; ranks overlapping specials."""
    if stake <= 0 or bankroll <= 0:
        return 0.0
    fraction = stake / bankroll
    net_odds = american_to_decimal(wz_american) - 1.0
    return fair_prob * math.log1p(net_odds * fraction) + (1 - fair_prob) * math.log1p(-fraction)
