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
