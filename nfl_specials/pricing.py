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
so DK's margin cancels in their ratio.

Fair = the WORST case (lowest probability) among the books that priced the
full partition — a special is only ever backed. Stake = Kelly fraction x full
Kelly at Wagerzon's price, fitted to the Wagerzon balance available
(budgeted_stakes). A team's fectas in one game are sized jointly with any
bet already held on that team (ExposureGroup): they win together, so each
stake is the Kelly stake GIVEN the others, not a blanket one-per-team rule.
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


@dataclass(frozen=True)
class Candidate:
    """A special that may take a new stake."""
    fair_prob: float
    wz_american: int
    max_stake: float               # the book's limit less what is already on this special


@dataclass(frozen=True)
class HeldBet:
    """A bet already placed: fixed exposure the new stakes are sized around."""
    fair_prob: float
    wz_american: int
    risk: float


@dataclass(frozen=True)
class ExposureGroup:
    """One team's fectas in one game, sized jointly.

    A team's fectas are nested: the stricter one winning means the looser one
    won too (GM -7.5 wins => GM wins; a superfecta wins => its trifecta part
    did). So the group is modelled as ONE chain — sorted by fair, the k-th
    most likely bet wins only if every likelier one did. That is exact for
    nested specials and the maximum positive correlation the fairs allow
    otherwise, the conservative case for bets that are only ever backed.
    """
    candidates: tuple[Candidate, ...]
    held: tuple[HeldBet, ...] = ()


def _chain_states(fair_probs: list[float]) -> list[tuple[float, int]]:
    """(probability, number of winners) per outcome of a nested chain whose
    bets are sorted by fair, likeliest first: state k = the first k win."""
    states = []
    upper = 1.0
    for winners, fair_prob in enumerate(fair_probs):
        states.append((max(0.0, upper - fair_prob), winners))
        upper = fair_prob
    states.append((upper, len(fair_probs)))
    return states


def _marginal_growth(index: int, fractions: list[float], net_odds: list[float],
                     states: list[tuple[float, int]]) -> float:
    """d E[log wealth] / d fraction[index] for a chain sorted likeliest first."""
    marginal = 0.0
    for probability, winners in states:
        if probability == 0.0:
            continue
        wealth = 1.0 + sum(net_odds[i] * f if i < winners else -f for i, f in enumerate(fractions))
        payoff = net_odds[index] if index < winners else -1.0
        marginal += probability * payoff / wealth
    return marginal


def _group_fractions_at_hurdle(fair_probs: list[float], net_odds: list[float], fixed: list[float | None],
                               caps: list[float], hurdle: float) -> list[float]:
    """Full-Kelly fractions of a chain (sorted likeliest first) at which each
    free bet's marginal log growth has fallen to `hurdle` or its cap/zero
    bound binds. `fixed[i]` is a held bet's fraction (None = free). Coordinate
    ascent on a concave objective: each step solves one bet exactly given the
    rest. A single free bet with nothing held is plain hurdle Kelly."""
    states = _chain_states(fair_probs)
    fractions = [f if f is not None else 0.0 for f in fixed]
    free = [i for i, f in enumerate(fixed) if f is None]
    for _sweep in range(500):
        largest_move = 0.0
        for i in free:
            others_at_risk = sum(f for j, f in enumerate(fractions) if j != i)
            # Wealth when every bet loses must stay positive (all-lose has
            # probability 1 - top fair > 0, so the optimum is interior).
            upper = min(caps[i], 1.0 - others_at_risk - 1e-9)

            def excess(f: float) -> float:
                trial = fractions[:i] + [f] + fractions[i + 1:]
                return _marginal_growth(i, trial, net_odds, states) - hurdle

            if upper <= 0 or excess(0.0) <= 0:
                best = 0.0
            elif excess(upper) >= 0:
                best = upper
            else:
                best = brentq(excess, 0.0, upper, xtol=1e-12)
            largest_move = max(largest_move, abs(best - fractions[i]))
            fractions[i] = best
        if largest_move < 1e-10:
            break
    return fractions


def budgeted_stakes(groups: list[ExposureGroup], bankroll: float, kelly_fraction: float,
                    budget: float | None, min_stake: float) -> list[list[float]]:
    """Whole-dollar stakes, per group and candidate, that fit `budget`.

    Uncapped, a group's stakes are joint fractional Kelly given its held
    bets: full Kelly on a bankroll of kelly_fraction x bankroll, the same
    scaling as a lone bet's kelly_fraction x full Kelly. A lone candidate
    with nothing held is exactly that. When the stakes add up past the
    budget, every stake's marginal log growth must clear one common hurdle, raised until they fit: weaker edges shrink
    first and drop to zero. A stake under `min_stake` cannot be placed, so
    the smallest such stake is dropped and the rest re-fit. Each stake is
    capped at its candidate's max_stake; the hurdle then hands what a capped
    bet can't take to the next-best edge. budget=None means no cap.
    """
    if budget is not None:
        budget = max(0.0, budget)   # Wagerzon can report a negative 'available'
    active = {(g, c) for g, group in enumerate(groups) for c, cand in enumerate(group.candidates)
              if expected_value(cand.fair_prob, cand.wz_american) > 0}
    while True:
        stakes = _fit(groups, active, bankroll * kelly_fraction, budget)
        too_small = [key for key in active if stakes[key] < min_stake]
        if not too_small:
            return [[float(stakes.get((g, c), 0)) for c in range(len(group.candidates))]
                    for g, group in enumerate(groups)]
        active.remove(min(too_small, key=lambda key: stakes[key]))


def _group_stakes_at_hurdle(group: ExposureGroup, active_candidates: list[int], kelly_bankroll: float,
                            hurdle: float) -> dict[int, float]:
    """Stakes (dollars) of a group's active candidates at `hurdle`."""
    if not active_candidates or kelly_bankroll <= 0:
        return {}
    chain = ([(h.fair_prob, h.wz_american, h.risk / kelly_bankroll, math.inf, None) for h in group.held] +
             [(group.candidates[c].fair_prob, group.candidates[c].wz_american, None,
               group.candidates[c].max_stake / kelly_bankroll, c) for c in active_candidates])
    chain.sort(key=lambda bet: -bet[0])
    fractions = _group_fractions_at_hurdle(
        [bet[0] for bet in chain], [american_to_decimal(bet[1]) - 1.0 for bet in chain],
        [bet[2] for bet in chain], [bet[3] for bet in chain], hurdle)
    return {bet[4]: fraction * kelly_bankroll for bet, fraction in zip(chain, fractions) if bet[4] is not None}


def _fit(groups: list[ExposureGroup], active: set[tuple[int, int]], kelly_bankroll: float,
         budget: float | None) -> dict[tuple[int, int], int]:
    def stakes_at(hurdle: float) -> dict[tuple[int, int], float]:
        stakes = {}
        for g, group in enumerate(groups):
            members = sorted(c for gg, c in active if gg == g)
            for c, stake in _group_stakes_at_hurdle(group, members, kelly_bankroll, hurdle).items():
                stakes[(g, c)] = stake
        return stakes

    uncapped = stakes_at(0.0)
    if budget is None or sum(uncapped.values()) <= budget:
        return {key: math.floor(s) for key, s in uncapped.items()}
    low = 0.0
    high = max(expected_value(groups[g].candidates[c].fair_prob, groups[g].candidates[c].wz_american)
               for g, c in active)
    for _ in range(60):   # bisection: the total is decreasing in the hurdle
        mid = (low + high) / 2
        if sum(stakes_at(mid).values()) > budget:
            low = mid
        else:
            high = mid
    return {key: math.floor(s) for key, s in stakes_at(high).items()}


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
            # Scores first goes LAST: DraftKings prices every outcome of the
            # last leg in one request, so both sides cost one call.
            decimal = book.price(game, tuple(outcome.ref for outcome in cell) + (first.ref,))
            calls += 1
            if decimal is None:
                labels = " + ".join([first.label] + [outcome.label for outcome in cell])
                return ScoresFirstShare(book.name, None, calls, reason=f"declined: {labels}")
            implied[first] += 1.0 / decimal
    share = implied[team_scores_first] / (implied[team_scores_first] + implied[opponent_scores_first])
    return ScoresFirstShare(book.name, share, calls)
