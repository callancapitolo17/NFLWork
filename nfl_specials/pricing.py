"""Fair value, edge and stake for one Wagerzon fecta.

Fair value at one book = probit devig over the FULL partition of the legs'
outcomes. Every combination of outcomes (SEA/Tie/LAC 1Q x SEA/Tie/LAC 1H x
SEA/LAC game, ...) is priced as an SGP at the book; the implied
probabilities sum to the book's SGP overround (~1.28 for 2 legs per the MLB
maker's measurements, far above the ~1.10 compounded single-leg vig), and the
devig spreads that margin across the cells. Pricing only the special's own
SGP and dividing out single-leg vig would leave most of that margin in and
overstate every edge.

Superfectas add "scores first", which only DraftKings lets into an SGP, so
they price at DraftKings alone: P(trifecta part) x P(scores first | trifecta
part). The first factor is DK's own partition fair for the other three legs;
the second comes from TWO DK SGPs — (team scores first + the trifecta legs)
and (opponent scores first + the same legs) — which share every other leg,
so DK's margin cancels in their ratio. That is ~20 DK calls per superfecta
instead of a 36-cell partition, at a book that denies a page after ~6 calls.

Fair = the WORST case (lowest probability) among the books that priced the
full partition — a special is only ever backed. Stake = Kelly fraction x full
Kelly at Wagerzon's price, fitted to the Wagerzon balance available
(budgeted_stakes); within one game only the best special per team gets a
stake, because a team's fectas win together.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

from scipy.optimize import brentq

from kalshi_common.fair_value import _probit_devig_n as probit_devig_n
from nfl_specials.books import BookGame, SgpBook
from nfl_specials.special_parser import Fecta

# A partition's implied sum must sit in [1, 1 + this x legs]; outside it a
# cell is stale or mispriced and the devig would be garbage (same envelope as
# kalshi_common.fair_value.PARTITION_OVERROUND_PER_LEG).
PARTITION_OVERROUND_PER_LEG = 0.25


@dataclass(frozen=True)
class ScoresFirstShare:
    book: str
    share: float | None            # P(team scores first | the special's other legs)
    n_calls: int
    reason: str | None = None      # why share is None
    sgp_decimal: float | None = None  # the book's own (vigged) price for the whole
                                      # special, when it is a single SGP


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


def worst_case_fair(book_fairs: list[BookFair]) -> float | None:
    """The LOWEST fair among books that priced the full partition.

    A Wagerzon special is only ever backed, never laid, so the book giving it
    the smallest chance is the conservative price (user decision 2026-10-02,
    replacing the mean)."""
    priced = [b.fair_prob for b in book_fairs if b.fair_prob is not None]
    if not priced:
        return None
    return min(priced)


def expected_value(fair_prob: float, wz_american: int) -> float:
    """EV per $1 risked at Wagerzon's price."""
    return fair_prob * american_to_decimal(wz_american) - 1.0


def kelly_stake(fair_prob: float, wz_american: int, bankroll: float, kelly_fraction: float) -> float:
    net_odds = american_to_decimal(wz_american) - 1.0
    full_kelly = expected_value(fair_prob, wz_american) / net_odds
    return max(0.0, full_kelly) * kelly_fraction * bankroll


def _kelly_fraction_at_hurdle(fair_prob: float, net_odds: float, hurdle: float) -> float:
    """Full-Kelly fraction f at which the bet's marginal log growth
    p*b/(1+b*f) - q/(1-f) has fallen to `hurdle` (hurdle 0 = plain Kelly)."""
    q = 1.0 - fair_prob
    marginal = lambda f: fair_prob * net_odds / (1 + net_odds * f) - q / (1 - f) - hurdle
    full_kelly = (fair_prob * net_odds - q) / net_odds
    if full_kelly <= 0 or marginal(0.0) <= 0:
        return 0.0
    if hurdle <= 0:
        return full_kelly
    return brentq(marginal, 0.0, full_kelly)


def budgeted_stakes(bets: list[tuple[float, int]], bankroll: float, kelly_fraction: float,
                    budget: float | None, min_stake: float) -> list[float]:
    """Whole-dollar stakes for (fair_prob, wz_american) bets that fit `budget`.

    Uncapped, each stake is kelly_fraction x full Kelly x bankroll. When those
    add up past the budget, every bet's marginal log growth must clear one
    common hurdle, raised until the stakes fit: weaker edges shrink first and
    drop to zero, and the budget lands on the strongest. A stake under
    `min_stake` cannot be placed, so the smallest such bet is dropped and the
    rest re-fit. budget=None means no cap.
    """
    active = [i for i, (p, a) in enumerate(bets) if expected_value(p, a) > 0]
    while True:
        stakes = _fit(bets, active, bankroll, kelly_fraction, budget)
        too_small = [i for i in active if stakes[i] < min_stake]
        if not too_small:
            return [float(stakes.get(i, 0)) for i in range(len(bets))]
        active.remove(min(too_small, key=lambda i: stakes[i]))


def _fit(bets, active, bankroll, kelly_fraction, budget) -> dict[int, int]:
    def stakes_at(hurdle: float) -> dict[int, float]:
        return {i: kelly_fraction * bankroll *
                   _kelly_fraction_at_hurdle(bets[i][0], american_to_decimal(bets[i][1]) - 1.0, hurdle)
                for i in active}

    uncapped = stakes_at(0.0)
    if budget is None or sum(uncapped.values()) <= budget:
        return {i: math.floor(s) for i, s in uncapped.items()}
    low, high = 0.0, max(expected_value(*bets[i]) for i in active)
    for _ in range(60):   # bisection: the total is decreasing in the hurdle
        mid = (low + high) / 2
        if sum(stakes_at(mid).values()) > budget:
            low = mid
        else:
            high = mid
    return {i: math.floor(s) for i, s in stakes_at(high).items()}


def log_growth(fair_prob: float, wz_american: int, stake: float, bankroll: float) -> float:
    """Expected log bankroll growth of one bet; ranks overlapping specials."""
    if stake <= 0 or bankroll <= 0:
        return 0.0
    fraction = stake / bankroll
    net_odds = american_to_decimal(wz_american) - 1.0
    return fair_prob * math.log1p(net_odds * fraction) + (1 - fair_prob) * math.log1p(-fraction)


def trifecta_part(fecta: Fecta) -> Fecta:
    """A superfecta without its leading 'scores first' leg."""
    if fecta.prop_type != "SUPERFECTA" or fecta.legs[0].kind != "scores_first":
        raise ValueError(f"expected a superfecta led by scores_first, got {fecta}")
    return Fecta(team=fecta.team, prop_type="TRIFECTA", legs=fecta.legs[1:])


def scores_first_share(book: SgpBook, game: BookGame, role: str, fecta: Fecta) -> ScoresFirstShare:
    """P(team scores first | the superfecta's other legs) from `book`'s SGPs.

    Prices (team scores first + rest) and (opponent scores first + rest);
    a rest leg that wins on two outcomes (a +0.5 read off a 3-way market)
    sums its cells. Assumes the book's margin is the same multiple on both
    SGPs — they differ only in the scores-first side — so it cancels.
    A 0-0 game (no first score) is ignored (~0.1% of NFL games).
    """
    first_score = book.leg_market(game, role, fecta.legs[0])
    if first_score is None or len(first_score.group) != 2:
        return ScoresFirstShare(book.name, None, 0, reason="no SGP 'scores first' market")
    team_scores_first = next(iter(first_score.winning))
    opponent_scores_first = next(o for o in first_score.group if o != team_scores_first)

    rest_markets = []
    for leg in fecta.legs[1:]:
        market = book.leg_market(game, role, leg)
        if market is None:
            return ScoresFirstShare(book.name, None, 0, reason=f"no SGP market for '{leg.describe(fecta.team)}'")
        rest_markets.append(market)
    rest_cells = list(itertools.product(*(sorted(m.winning, key=lambda o: o.label) for m in rest_markets)))

    implied = {team_scores_first: 0.0, opponent_scores_first: 0.0}
    calls = 0
    for first in (team_scores_first, opponent_scores_first):
        for cell in rest_cells:
            decimal = book.price(game, (first.ref,) + tuple(outcome.ref for outcome in cell))
            calls += 1
            if decimal is None:
                labels = " + ".join([first.label] + [outcome.label for outcome in cell])
                return ScoresFirstShare(book.name, None, calls, reason=f"declined: {labels}")
            implied[first] += 1.0 / decimal
    share = implied[team_scores_first] / (implied[team_scores_first] + implied[opponent_scores_first])
    own_price = 1.0 / implied[team_scores_first] if len(rest_cells) == 1 else None
    return ScoresFirstShare(book.name, share, calls, sgp_decimal=own_price)
