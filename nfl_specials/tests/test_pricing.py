import itertools

import pytest
from scipy.stats import norm

from nfl_specials.books import BookGame, LegMarket, Outcome, next_game, parse_start, widen_with_tie
from nfl_specials.pricing import (american_to_decimal, decimal_to_american, expected_value,
                                  kelly_stake, price_fecta_at_book)
from nfl_specials.special_parser import Fecta, Leg

GAME = BookGame("evt", home="SEA", away="LAC", game_start_time="2099-01-01T00:00:00Z")


class FakeBook:
    """Two-outcome independent legs with known true probabilities. Each SGP
    cell's implied probability is the true one shifted up by `vig_z` in
    probit space — the margin shape probit devig removes exactly."""

    name = "Fake"

    def __init__(self, leg_probs: list[float], vig_z: float = 0.1, declined: set = frozenset()):
        self.leg_probs = leg_probs
        self.vig_z = vig_z
        self.declined = declined

    def leg_market(self, game, role, leg):
        index = {"Q1": 0, "H1": 1, "GM": 2}[leg.period]
        win, lose = Outcome((index, "win"), f"win {index}"), Outcome((index, "lose"), f"lose {index}")
        return LegMarket.single(win, (win, lose))

    def price(self, game, refs):
        if refs in self.declined:
            return None
        prob = 1.0
        for index, side in refs:
            p = self.leg_probs[index]
            prob *= p if side == "win" else 1 - p
        return 1.0 / norm.cdf(norm.ppf(prob) + self.vig_z)


TRIFECTA = Fecta("SEA", "TRIFECTA", (Leg("win", "Q1"), Leg("win", "H1"), Leg("win", "GM")))


def test_partition_devig_recovers_true_joint_probability():
    fair = price_fecta_at_book(FakeBook([0.6, 0.65, 0.75]), GAME, "home", TRIFECTA)
    assert fair.n_cells == 8
    assert 1.1 < fair.overround < 1.5
    assert fair.fair_prob == pytest.approx(0.6 * 0.65 * 0.75, abs=1e-6)
    assert fair.fair_prob < 1.0 / fair.sgp_decimal


def test_declined_cell_means_no_fair():
    declined = {((0, "lose"), (1, "lose"), (2, "lose"))}
    fair = price_fecta_at_book(FakeBook([0.6, 0.65, 0.75], declined=declined), GAME, "home", TRIFECTA)
    assert fair.fair_prob is None
    assert "declined cell" in fair.reason


def test_overround_outside_envelope_is_rejected():
    fair = price_fecta_at_book(FakeBook([0.6, 0.65, 0.75], vig_z=0.8), GAME, "home", TRIFECTA)
    assert fair.fair_prob is None
    assert "overround" in fair.reason


def test_widen_with_tie_adds_the_tie():
    """A '+0.5' leg read off a 3-way market wins on the team OR the tie."""
    team, tie, opp = Outcome("t", "SEA"), Outcome("x", "Tie"), Outcome("o", "LAC")
    three_way = LegMarket.single(team, (team, tie, opp))
    widened = widen_with_tie(three_way)
    assert widened.winning == frozenset({team, tie})
    assert widen_with_tie(LegMarket.single(team, (team, opp))) is None


def test_ev_and_kelly_at_wagerzon_price():
    assert american_to_decimal(150) == pytest.approx(2.5)
    assert american_to_decimal(-200) == pytest.approx(1.5)
    assert decimal_to_american(2.5) == 150
    assert decimal_to_american(1.5) == -200
    assert expected_value(0.5, 150) == pytest.approx(0.25)
    # full Kelly = EV / net odds = 0.25 / 1.5; quarter Kelly on $1000
    assert kelly_stake(0.5, 150, 1000, 0.25) == pytest.approx(1000 * 0.25 * 0.25 / 1.5)
    assert kelly_stake(0.3, 150, 1000, 0.25) == 0.0


def test_parse_start_handles_seven_fractional_digits():
    parsed = parse_start("2026-10-04T20:25:00.0000000Z")
    assert parsed.utcoffset().total_seconds() == 0
    assert (parsed.hour, parsed.minute) == (20, 25)


def test_next_game_skips_started_and_later_weeks():
    games = [BookGame("past", "SEA", "LAC", "2000-01-01T00:00:00Z"),
             BookGame("later", "SEA", "SF", "2099-02-01T00:00:00Z"),
             BookGame("next", "ARI", "SEA", "2099-01-01T00:00:00Z")]
    assert next_game(games, "SEA").book_event_id == "next"
    assert next_game(games, "KC") is None


class ScoresFirstBook:
    """Independent legs; the SGP price is the true joint probability times a
    flat margin, so the scores-first ratio must come back exact."""

    name = "Fake"
    TRUE = {"sf": 0.55, "Q1": 0.6, "H1": 0.65, "GM": 0.75}

    def leg_market(self, game, role, leg):
        key = "sf" if leg.kind == "scores_first" else leg.period
        win, lose = Outcome((key, True), f"{key} yes"), Outcome((key, False), f"{key} no")
        return LegMarket.single(win, (win, lose))

    def __init__(self):
        self.priced = []

    def price(self, game, refs):
        self.priced.append(refs)
        prob = 1.0
        for key, side in refs:
            prob *= self.TRUE[key] if side else 1 - self.TRUE[key]
        return 1.0 / (prob * 1.3)


def test_scores_first_share_cancels_the_margin():
    from nfl_specials.pricing import scores_first_share
    superfecta = Fecta("SEA", "SUPERFECTA", (Leg("scores_first"),) + TRIFECTA.legs)
    book = ScoresFirstBook()
    result = scores_first_share(book, GAME, "home", superfecta)
    assert result.n_calls == 2
    assert result.share == pytest.approx(0.55)
    # Scores first is the last leg, so DraftKings prices both sides in one request.
    assert [refs[-1][0] for refs in book.priced] == ["sf", "sf"]


def test_trifecta_part_drops_scores_first():
    from nfl_specials.pricing import trifecta_part
    superfecta = Fecta("SEA", "SUPERFECTA", (Leg("scores_first"),) + TRIFECTA.legs)
    assert trifecta_part(superfecta) == TRIFECTA
    with pytest.raises(ValueError):
        trifecta_part(TRIFECTA)


def _lone_stakes(bets, bankroll, kelly_fraction, budget, min_stake, max_stake):
    """budgeted_stakes with every bet in its own group, flattened."""
    from nfl_specials.pricing import Candidate, ExposureGroup, budgeted_stakes
    groups = [ExposureGroup((Candidate(p, a, max_stake),)) for p, a in bets]
    return [stakes[0] for stakes in budgeted_stakes(groups, bankroll, kelly_fraction, budget, min_stake)]


def test_budgeted_stakes_are_plain_fractional_kelly_when_they_fit():
    bets = [(0.12, 1730), (0.17, 855)]
    uncapped = _lone_stakes(bets, 35000, 0.25, None, 20, 1e9)
    assert uncapped == [float(int(kelly_stake(p, a, 35000, 0.25))) for p, a in bets]
    assert _lone_stakes(bets, 35000, 0.25, 10_000, 20, 1e9) == uncapped


def test_budgeted_stakes_trim_the_weakest_first_and_fit():
    strong, weak = (0.115, 1730), (0.081, 1600)        # EV +111% vs +38%
    stakes = _lone_stakes([strong, weak], 35000, 0.25, 300, 20, 1e9)
    assert sum(stakes) <= 300
    assert stakes[0] > 0 and stakes[1] == 0.0


def test_stakes_under_the_minimum_are_dropped_and_refit():
    assert _lone_stakes([(0.12, 1730), (0.081, 1600)], 1000, 0.25, None, 20, 250) == [0.0, 0.0]
    assert _lone_stakes([(0.12, 1730)], 35000, 0.25, 15, 20, 250) == [0.0]


def test_negative_budget_means_no_stakes_not_a_crash():
    assert _lone_stakes([(0.12, 1730), (0.17, 855)], 35000, 0.25, -50, 20, 250) == [0.0, 0.0]


def test_stakes_cap_at_the_max_and_the_rest_flows_to_the_next_edge():
    stakes = _lone_stakes([(1 / 7.92, 1730), (1 / 6.04, 855)], 35000, 0.25, 500, 20, 250)
    assert stakes[0] == 250.0
    assert 240 <= stakes[1] <= 250


def test_exclusive_outcomes_size_bigger_than_independent_bets():
    """31% at +250 and 28.35% at +275 that cannot both win hedge each other."""
    from nfl_specials.pricing import Candidate, ExposureGroup, OutcomeBet, OutcomeGroup, budgeted_stakes
    exclusive = OutcomeGroup((0.31, 0.2835, 0.4065), (OutcomeBet((2.5, -1, -1)), OutcomeBet((-1, 2.75, -1))))
    assert budgeted_stakes([exclusive], 35000, 0.25, None, 0) == [[442.0, 361.0]]
    independent = [ExposureGroup((Candidate(0.31, 250, 1e9),)), ExposureGroup((Candidate(0.2835, 275, 1e9),))]
    assert budgeted_stakes(independent, 35000, 0.25, None, 0) == [[297.0], [200.0]]


def test_outcome_group_kelly_matches_brute_force_with_held_bets():
    import numpy as np
    from scipy.optimize import minimize
    from nfl_specials.pricing import OutcomeBet, OutcomeGroup, budgeted_stakes
    probs = (0.12, 0.08, 0.2, 0.15, 0.25, 0.2)
    candidates = (OutcomeBet((4.5, 4.5, -1, -1, -1, -1)), OutcomeBet((-1, -1, 2.8, 0.9, -1, -1)),
                  OutcomeBet((-1, -1, -1, 0.0, 2.8, -1)))
    held = (OutcomeBet((9.0, -1, -1, -1, -1, -1), risk=150.0),)
    kelly_bankroll = 10000 * 0.25
    group = OutcomeGroup(probs, candidates, held)
    assert all(group.expected_value(c) > 0 for c in range(3))   # only +EV bets are staked
    stakes = budgeted_stakes([group], 10000, 0.25, None, 0)[0]

    payoffs = np.array([c.payoffs for c in candidates]).T
    held_wealth = np.array(held[0].payoffs) * held[0].risk / kelly_bankroll

    def neg_growth(f):
        return -float(np.dot(probs, np.log(1 + held_wealth + payoffs @ f)))

    best = minimize(neg_growth, np.full(3, 0.01), bounds=[(0, 0.5)] * 3, method="L-BFGS-B",
                    options={"ftol": 1e-14, "gtol": 1e-12})
    assert stakes == pytest.approx(list(np.floor(best.x * kelly_bankroll)), abs=1.0)


def test_outcome_group_refuses_probabilities_that_do_not_sum_to_one():
    from nfl_specials.pricing import OutcomeBet, OutcomeGroup
    with pytest.raises(ValueError, match="sum to 1"):
        OutcomeGroup((0.5, 0.4), (OutcomeBet((1.0, -1.0)),))
    with pytest.raises(ValueError, match="payoffs per bet"):
        OutcomeGroup((0.5, 0.5), (OutcomeBet((1.0,)),))


def test_a_stake_a_hair_under_a_dollar_is_that_dollar():
    from nfl_specials.pricing import _whole_dollars
    assert _whole_dollars(49.9999999) == 50 and _whole_dollars(49.99) == 49 and _whole_dollars(50.0) == 50


def test_the_leg_with_the_most_outcomes_goes_last_in_every_request():
    """DraftKings prices every outcome of the LAST leg in one request, so a
    3 x 3 x 2 partition with a 3-way leg last needs 6 distinct bases, not 9."""
    from nfl_specials.pricing import price_partition

    def market(name, size):
        outcomes = tuple(Outcome((name, i), f"{name}{i}") for i in range(size))
        return LegMarket.single(outcomes[0], outcomes)

    class RecordingBook:
        name = "Recorder"

        def __init__(self):
            self.requests = []

        def price(self, game, refs):
            self.requests.append(refs)
            return 1.0 / (1.05 / 18)

    book = RecordingBook()
    partition = price_partition(book, GAME, [market("q1", 3), market("h1", 3), market("ml", 2)])
    assert len(partition.cells) == 18 and partition.cell_fairs is not None
    assert all(refs[-1][0] == "h1" for refs in book.requests)
    assert len({frozenset(refs[:-1]) for refs in book.requests}) == 6
    # The cells keep the legs' own order.
    assert [outcome.ref[0] for outcome in partition.cells[0]] == ["q1", "h1", "ml"]


# --- PriceMemory: reuse across refreshes -------------------------------------------
class SinglesBook(FakeBook):
    """FakeBook that counts SGP calls and shows each outcome's single price."""

    def __init__(self, leg_probs, singles: dict | None = None, **kwargs):
        super().__init__(leg_probs, **kwargs)
        self.singles = singles if singles is not None else {}
        self.calls = 0

    def price(self, game, refs):
        self.calls += 1
        return super().price(game, refs)

    def single_odds(self, game, refs):
        odds = tuple(self.singles.get(ref, "-110") for ref in refs)
        return None if None in odds else odds


class FakeClock:
    def __init__(self):
        self.now = 0.0

    def __call__(self):
        return self.now


def test_an_unchanged_partition_is_reused_with_its_priced_at():
    from nfl_specials.pricing import PriceMemory
    memory = PriceMemory()
    first_book, second_book = SinglesBook([0.6, 0.55, 0.7]), SinglesBook([0.6, 0.55, 0.7])
    first = price_fecta_at_book(first_book, GAME, "home", TRIFECTA, memory)
    second = price_fecta_at_book(second_book, GAME, "home", TRIFECTA, memory)
    assert first_book.calls == 8 and second_book.calls == 0
    assert second == first and second.priced_at is not None


def test_a_moved_single_price_prices_the_partition_again():
    from nfl_specials.pricing import PriceMemory
    memory = PriceMemory()
    price_fecta_at_book(SinglesBook([0.6, 0.55, 0.7]), GAME, "home", TRIFECTA, memory)
    moved = SinglesBook([0.6, 0.55, 0.72], singles={(2, "win"): "-125"})
    fair = price_fecta_at_book(moved, GAME, "home", TRIFECTA, memory)
    assert moved.calls == 8
    assert fair.fair_prob == pytest.approx(0.6 * 0.55 * 0.72, abs=1e-6)


def test_a_partition_older_than_the_age_limit_prices_again():
    from nfl_specials.pricing import MAX_REUSE_AGE_SECONDS, PriceMemory
    clock = FakeClock()
    memory = PriceMemory(clock=clock)
    price_fecta_at_book(SinglesBook([0.6, 0.55, 0.7]), GAME, "home", TRIFECTA, memory)
    clock.now = MAX_REUSE_AGE_SECONDS - 1
    fresh = SinglesBook([0.6, 0.55, 0.7])
    price_fecta_at_book(fresh, GAME, "home", TRIFECTA, memory)
    clock.now = MAX_REUSE_AGE_SECONDS + 1    # age counts from the first pricing, not the reuse
    stale = SinglesBook([0.6, 0.55, 0.7])
    price_fecta_at_book(stale, GAME, "home", TRIFECTA, memory)
    assert (fresh.calls, stale.calls) == (0, 8)


def test_no_single_prices_or_a_declined_cell_is_never_reused():
    from nfl_specials.pricing import PriceMemory
    memory = PriceMemory()
    declined = SinglesBook([0.6, 0.55, 0.7], declined={((0, "win"), (1, "win"), (2, "win"))})
    assert price_fecta_at_book(declined, GAME, "home", TRIFECTA, memory).fair_prob is None
    again = SinglesBook([0.6, 0.55, 0.7])
    assert price_fecta_at_book(again, GAME, "home", TRIFECTA, memory).fair_prob is not None
    assert again.calls == 8
    unknown = SinglesBook([0.6, 0.55, 0.7], singles={(0, "win"): None})
    price_fecta_at_book(unknown, GAME, "home", TRIFECTA, memory)
    assert unknown.calls == 8


def test_an_unchanged_scores_first_share_is_reused():
    from nfl_specials.pricing import PriceMemory, scores_first_share

    class SinglesScoresFirstBook(ScoresFirstBook):
        def single_odds(self, game, refs):
            return tuple("-110" for _ in refs)

    memory = PriceMemory()
    superfecta = Fecta("SEA", "SUPERFECTA", (Leg("scores_first"),) + TRIFECTA.legs)
    first_book, second_book = SinglesScoresFirstBook(), SinglesScoresFirstBook()
    first = scores_first_share(first_book, GAME, "home", superfecta, memory)
    second = scores_first_share(second_book, GAME, "home", superfecta, memory)
    assert (len(first_book.priced), len(second_book.priced)) == (2, 0)
    assert second == first


def test_request_pacer_keeps_the_gap_and_stops_after_halt():
    import time
    from nfl_specials.books import BookBlocked, RequestPacer
    pacer = RequestPacer("FanDuel", 0.05)
    stamps = []
    for _ in range(3):
        pacer.call(lambda: stamps.append(time.monotonic()), "implyBets")
    assert min(b - a for a, b in zip(stamps, stamps[1:])) >= 0.045

    class Answer:
        status_code = 403
    pacer.guard(lambda: Answer(), "event page")          # the blocked answer itself comes back
    with pytest.raises(BookBlocked, match="FanDuel blocked this refresh: event page HTTP 403"):
        pacer.call(lambda: stamps.append("sent"), "implyBets")
    with pytest.raises(BookBlocked):
        pacer.guard(lambda: stamps.append("sent"), "event page")
    assert "sent" not in stamps


def test_a_line_moved_under_the_same_ids_and_prices_is_priced_again():
    """FanDuel re-centres a main spread under the same market and selection
    ids at the same price: only the label carries the line, so the label is
    part of what is remembered."""
    from nfl_specials.pricing import PriceMemory, price_partition

    def markets(spread_line):
        away = Outcome(("m1", "away"), f"away {spread_line:+g}")
        home = Outcome(("m1", "home"), f"home {-spread_line:+g}")
        over, under = Outcome(("m2", "over"), "Over 23.5"), Outcome(("m2", "under"), "Under 23.5")
        return [LegMarket.single(away, (away, home)), LegMarket.single(over, (over, under))]

    class SameIdsBook:
        name = "FanDuel"
        calls = 0

        def price(self, game, refs):
            self.calls += 1
            return 3.6

        def single_odds(self, game, refs):
            return tuple(1.909 for _ in refs)

    memory = PriceMemory()
    book = SameIdsBook()
    first = price_partition(book, GAME, markets(-2.5), memory)
    moved = price_partition(book, GAME, markets(-3.5), memory)
    assert book.calls == 8
    assert moved.cells == tuple(itertools.product(*(m.group for m in markets(-3.5))))
    assert first.cells != moved.cells
