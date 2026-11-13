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

import numpy as np
from scipy.optimize import brentq

from kalshi_common.fair_value import _probit_devig_n as probit_devig_n
from nfl_specials.books import BookGame, LegMarket, SgpBook
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
    reason: str | None = None      # why share is None, or why a fixed share stands in


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


@dataclass(frozen=True)
class PartitionFair:
    """Every cell of a leg partition at one book, devigged."""
    cells: tuple[tuple, ...]               # one Outcome per leg, in itertools.product order
    cell_fairs: tuple[float, ...] | None   # devigged probability per cell; None when unpriced
    implied: tuple[float, ...]             # 1 / the book's decimal per cell priced
    overround: float | None
    reason: str | None = None              # why cell_fairs is None


def price_partition(book: SgpBook, game: BookGame, leg_markets: list[LegMarket]) -> PartitionFair:
    """Price every cell of the legs' partition at `book` (one SGP call each)
    and probit-devig across them. A declined cell or an implied sum outside
    [1, 1 + PARTITION_OVERROUND_PER_LEG x legs] prices nothing."""
    cells = tuple(itertools.product(*(m.group for m in leg_markets)))
    implied = []
    for cell in cells:
        decimal = book.price(game, tuple(outcome.ref for outcome in cell))
        if decimal is None:
            labels = " + ".join(outcome.label for outcome in cell)
            return PartitionFair(cells, None, tuple(implied), None, reason=f"declined cell: {labels}")
        implied.append(1.0 / decimal)
    overround = sum(implied)
    max_overround = 1.0 + PARTITION_OVERROUND_PER_LEG * len(leg_markets)
    if not 1.0 <= overround <= max_overround:
        return PartitionFair(cells, None, tuple(implied), overround,
                             reason=f"partition overround {overround:.3f} outside [1, {max_overround:.2f}]")
    devigged = probit_devig_n(implied)
    return PartitionFair(cells, tuple(float(p) for p in devigged), tuple(implied), overround)


def price_fecta_at_book(book: SgpBook, game: BookGame, role: str, fecta: Fecta) -> BookFair:
    """The fecta's fair at `book`: the devigged cells on which every leg wins."""
    leg_markets = []
    for leg in fecta.legs:
        market = book.leg_market(game, role, leg)
        if market is None:
            return BookFair(book.name, None, None, None, 0,
                            reason=f"no SGP market for '{leg.describe(fecta.team)}'")
        leg_markets.append(market)

    partition = price_partition(book, game, leg_markets)
    winning_cells = [index for index, cell in enumerate(partition.cells)
                     if all(outcome in market.winning for outcome, market in zip(cell, leg_markets))]
    # The book's own price for the special exists only when it is one cell
    # (a "+0.5 off a 3-way market" leg makes it a union of cells).
    sgp_decimal = None
    if len(winning_cells) == 1 and winning_cells[0] < len(partition.implied):
        sgp_decimal = 1.0 / partition.implied[winning_cells[0]]
    if partition.cell_fairs is None:
        if partition.overround is None:   # a declined cell: no own price either
            sgp_decimal = None
        return BookFair(book.name, None, sgp_decimal, partition.overround, len(partition.cells),
                        reason=partition.reason)
    fair = sum(partition.cell_fairs[index] for index in winning_cells)
    return BookFair(book.name, float(fair), sgp_decimal, partition.overround, len(partition.cells))


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


def free_play_value(fair_prob: float, wz_american: int) -> float:
    """What $1 of free play is worth at Wagerzon's price: a win pays only the
    to-win and a loss costs no cash, so fair x net odds."""
    return fair_prob * (american_to_decimal(wz_american) - 1.0)


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

    def as_outcome_group(self) -> "OutcomeGroup":
        """The same chain written as explicit outcome states."""
        bets = ([(h.fair_prob, h.wz_american) for h in self.held] +
                [(c.fair_prob, c.wz_american) for c in self.candidates])
        order = sorted(range(len(bets)), key=lambda i: -bets[i][0])   # likeliest first
        states = chain_states([bets[i][0] for i in order])
        rank = {bet: position for position, bet in enumerate(order)}
        payoffs = []
        for i, (_fair, american) in enumerate(bets):
            net_odds = american_to_decimal(american) - 1.0
            payoffs.append(tuple(net_odds if rank[i] < winners else -1.0 for _p, winners in states))
        held_count = len(self.held)
        return OutcomeGroup(
            state_probs=tuple(p for p, _winners in states),
            candidates=tuple(OutcomeBet(payoffs[held_count + c], max_stake=cand.max_stake)
                             for c, cand in enumerate(self.candidates)),
            held=tuple(OutcomeBet(payoffs[h], risk=bet.risk) for h, bet in enumerate(self.held)))


@dataclass(frozen=True)
class OutcomeBet:
    """One bet inside an OutcomeGroup: what it returns in each outcome state."""
    payoffs: tuple[float, ...]     # $ won (+) or lost (-) per $1 staked, one per state; 0 = stake back
    max_stake: float = math.inf    # candidate: the book's limit less what is already on it
    risk: float = 0.0              # held bet: the dollars already staked


@dataclass(frozen=True)
class OutcomeGroup:
    """Bets whose results all follow from one set of exclusive outcome states,
    sized jointly. The general form of ExposureGroup: any mix of bets that
    win, lose, push or pay partly in each state (a game's four 1H corners and
    both teams' fectas, say)."""
    state_probs: tuple[float, ...]
    candidates: tuple[OutcomeBet, ...]
    held: tuple[OutcomeBet, ...] = ()

    def __post_init__(self):
        if any(p < 0 for p in self.state_probs) or abs(sum(self.state_probs) - 1.0) > 1e-9:
            raise ValueError(f"state probabilities must be >= 0 and sum to 1, got sum {sum(self.state_probs)}")
        for bet in self.candidates + self.held:
            if len(bet.payoffs) != len(self.state_probs):
                raise ValueError(f"expected {len(self.state_probs)} payoffs per bet, got {len(bet.payoffs)}")

    def expected_value(self, candidate: int) -> float:
        """EV per $1 of a candidate on its own."""
        return sum(p * r for p, r in zip(self.state_probs, self.candidates[candidate].payoffs))


def chain_states(fair_probs: list[float]) -> list[tuple[float, int]]:
    """(probability, number of winners) per outcome of a nested chain whose
    bets are sorted by fair, likeliest first: state k = the first k win."""
    states = []
    upper = 1.0
    for winners, fair_prob in enumerate(fair_probs):
        states.append((max(0.0, upper - fair_prob), winners))
        upper = fair_prob
    states.append((upper, len(fair_probs)))
    return states


def _group_fractions_at_hurdle(state_probs: tuple[float, ...], payoffs: list[tuple[float, ...]],
                               fixed: list[float | None], caps: list[float], hurdle: float) -> list[float]:
    """Full-Kelly fractions at which each free bet's marginal log growth has
    fallen to `hurdle` or its cap/zero bound binds. `fixed[i]` is a held
    bet's fraction (None = free). Coordinate ascent on a concave objective:
    each step solves one bet exactly given the rest. A single free bet with
    nothing held is plain hurdle Kelly."""
    possible = np.array(state_probs) > 0.0           # impossible states never constrain or count
    probs = np.array(state_probs)[possible]
    returns = np.array(payoffs, dtype=float).reshape(len(payoffs), -1)[:, possible]
    fractions = np.array([f if f is not None else 0.0 for f in fixed])
    wealth = 1.0 + fractions @ returns               # per state
    free = [i for i, f in enumerate(fixed) if f is None]
    for _sweep in range(500):
        largest_move = 0.0
        for i in free:
            bet_returns = returns[i]
            base = wealth - fractions[i] * bet_returns   # wealth from every other bet
            losing = bet_returns < 0
            if not losing.any():
                raise ValueError("a candidate that never loses cannot be Kelly-sized")
            # Wealth must stay positive in every possible state; a bet that
            # can lose has an interior optimum below that bound.
            upper = min(caps[i], float(np.min(base[losing] / -bet_returns[losing])) - 1e-9)

            def excess(f: float) -> float:
                return float(np.dot(probs, bet_returns / (base + f * bet_returns))) - hurdle

            if upper <= 0 or excess(0.0) <= 0:
                best = 0.0
            elif excess(upper) >= 0:
                best = upper
            else:
                best = brentq(excess, 0.0, upper, xtol=1e-12)
            largest_move = max(largest_move, abs(best - fractions[i]))
            fractions[i] = best
            wealth = base + best * bet_returns
        if largest_move < 1e-10:
            break
    return [float(f) for f in fractions]


def budgeted_stakes(groups: list[ExposureGroup | OutcomeGroup], bankroll: float, kelly_fraction: float,
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
    outcome_groups = [g.as_outcome_group() if isinstance(g, ExposureGroup) else g for g in groups]
    active = {(g, c) for g, group in enumerate(outcome_groups) for c in range(len(group.candidates))
              if group.expected_value(c) > 0}
    solved: dict[tuple, dict[int, float]] = {}   # (group, its active candidates, hurdle) -> stakes
    while True:
        stakes = _fit(outcome_groups, active, bankroll * kelly_fraction, budget, solved)
        too_small = [key for key in active if stakes[key] < min_stake]
        if not too_small:
            return [[float(stakes.get((g, c), 0)) for c in range(len(group.candidates))]
                    for g, group in enumerate(outcome_groups)]
        active.remove(min(too_small, key=lambda key: stakes[key]))


def _group_stakes_at_hurdle(group: OutcomeGroup, active_candidates: list[int], kelly_bankroll: float,
                            hurdle: float) -> dict[int, float]:
    """Stakes (dollars) of a group's active candidates at `hurdle`; the
    inactive ones stake nothing."""
    if not active_candidates or kelly_bankroll <= 0:
        return {}
    bets = ([(h.payoffs, h.risk / kelly_bankroll, math.inf, None) for h in group.held] +
            [(group.candidates[c].payoffs, None, group.candidates[c].max_stake / kelly_bankroll, c)
             for c in active_candidates])
    fractions = _group_fractions_at_hurdle(
        group.state_probs, [bet[0] for bet in bets], [bet[1] for bet in bets], [bet[2] for bet in bets], hurdle)
    return {bet[3]: fraction * kelly_bankroll for bet, fraction in zip(bets, fractions) if bet[3] is not None}


# Raising the hurdle this many times over without fitting the budget means
# the stakes do not shrink with it: a bug, not a board to size.
MAX_HURDLE_DOUBLINGS = 60


def _fit(groups: list[OutcomeGroup], active: set[tuple[int, int]], kelly_bankroll: float,
         budget: float | None, solved: dict[tuple, dict[int, float]]) -> dict[tuple[int, int], int]:
    """`solved` caches each group's stakes across calls: dropping one stake
    under the minimum re-fits every group, but only that group changed."""
    def stakes_at(hurdle: float) -> dict[tuple[int, int], float]:
        stakes = {}
        for g, group in enumerate(groups):
            members = tuple(sorted(c for gg, c in active if gg == g))
            key = (g, members, hurdle)
            if key not in solved:
                solved[key] = _group_stakes_at_hurdle(group, list(members), kelly_bankroll, hurdle)
            for c, stake in solved[key].items():
                stakes[(g, c)] = stake
        return stakes

    uncapped = stakes_at(0.0)
    if budget is None or sum(uncapped.values()) <= budget:
        return {key: math.floor(s) for key, s in uncapped.items()}
    low = 0.0
    high = max(groups[g].expected_value(c) for g, c in active)
    # A bet that hedges a held one can grow faster at the margin than its own
    # EV, so the first guess may still be too low a hurdle.
    for _ in range(MAX_HURDLE_DOUBLINGS):
        if sum(stakes_at(high).values()) <= budget:
            break
        low, high = high, high * 2
    else:
        raise RuntimeError(f"no hurdle up to {high:g} fits the stakes into the ${budget:g} budget")
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
