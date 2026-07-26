# tests/test_maker_ledger.py
from unabated_edge.maker import ledger


def test_empty_book_is_flat():
    assert ledger.worst_case([]) == 0.0
    assert ledger.pnl_grid([]) == [0.0]


def test_single_yes_fill_worst_case():
    # 40x YES Over-2.5 @ 0.42: loses 16.8 when g<=2, wins 23.2 when g>=3
    fills = [(2.5, "yes", 40, 0.42)]
    grid = ledger.pnl_grid(fills)
    assert round(grid[0], 2) == -16.8 and round(grid[1], 2) == 23.2
    assert round(ledger.worst_case(fills), 2) == -16.8


def test_offset_under25_over15_worst_and_middle():
    # From the design conversation: 20x NO Over-2.5 @0.55 + 20x YES Over-1.5 @0.70
    fills = [(2.5, "no", 20, 0.55), (1.5, "yes", 20, 0.70)]
    grid = ledger.pnl_grid(fills)
    assert round(grid[0], 2) == -5.0          # low-goal worlds
    assert round(grid[1], 2) == 15.0          # the middle: both win at exactly 2
    assert round(grid[2], 2) == -5.0          # high-goal worlds
    assert round(ledger.worst_case(fills), 2) == -5.0


def test_grid_constant_above_top_rung():
    # No Kalshi rung above 5.5 -> pnl must be identical above 5.5
    fills = [(5.5, "yes", 10, 0.05), (2.5, "no", 7, 0.6)]
    grid = ledger.pnl_grid(fills)
    # With interval ledger, above the top rung (5.5) there's only one interval
    # so all outcomes > 5.5 map to the same grid value
    assert len(grid) == 3  # t<2.5, 2.5<t<=5.5, t>5.5
    assert round(grid[2], 9) == round(ledger.pnl(fills, 6), 9)


def test_max_contracts_shrinks_to_budget():
    # existing 40x YES Over-2.5 @0.42 (worst -16.8); candidate YES Over-1.5 @0.70,
    # budget 30 -> allowed n = floor((30-16.8)/0.70) = 18
    fills = [(2.5, "yes", 40, 0.42)]
    n = ledger.max_contracts(fills, 1.5, "yes", 0.70, 30.0)
    assert n == 18
    assert ledger.worst_case(fills + [(1.5, "yes", 18, 0.70)]) >= -30.0
    assert ledger.worst_case(fills + [(1.5, "yes", 19, 0.70)]) < -30.0


def test_max_contracts_opposite_direction_gets_room():
    # candidate NO Over-2.5 offsets the existing YES: plenty of room
    fills = [(2.5, "yes", 40, 0.42)]
    n = ledger.max_contracts(fills, 2.5, "no", 0.58, 30.0)
    assert n > 40


def test_max_contracts_empty_book():
    assert ledger.max_contracts([], 2.5, "yes", 0.42, 400.0) == int(400 / 0.42)


def test_max_contracts_zero_when_already_breached_or_no_budget():
    fills = [(2.5, "yes", 100, 0.42)]          # worst -42
    assert ledger.max_contracts(fills, 1.5, "yes", 0.70, 30.0) == 0
    assert ledger.max_contracts([], 2.5, "yes", 0.42, 0.0) == 0


def test_mark_to_fair_yes_profit_when_fair_rises():
    val = ledger.mark_to_fair([(2.5, "yes", 10, 0.40)], {2.5: 0.50})
    assert round(val, 9) == round(10 * (0.50 - 0.40), 9) == 1.0


def test_mark_to_fair_no_side():
    val = ledger.mark_to_fair([(2.5, "no", 10, 0.55)], {2.5: 0.50})
    assert round(val, 9) == round(10 * ((1 - 0.50) - 0.55), 9) == -0.5


def test_mark_to_fair_unmarked_line_contributes_zero():
    assert ledger.mark_to_fair([(2.5, "yes", 10, 0.40)], {}) == 0.0


# NEW TESTS FOR INTERVAL LEDGER
import itertools
import random

from unabated_edge.maker import ledger


def _worst_case_bruteforce(fills, t_max=40):
    """Reference: exhaustive integer totals 0..t_max (superset of any old grid)."""
    return min(ledger.pnl(fills, t) for t in range(t_max + 1))


def test_outcome_points_cover_every_interval():
    pts = ledger.outcome_points([2.5, 8.5])
    # intervals: t<2.5, 2.5<t<8.5, t>8.5 -> one representative each
    assert pts == [2.0, 3.0, 9.0]


def test_outcome_points_extra_line_included():
    assert ledger.outcome_points([2.5], extra_line=8.5) == [2.0, 3.0, 9.0]
    assert ledger.outcome_points([], extra_line=None) == [0.0]


def test_high_line_no_longer_truncated():
    # A NO fill at 12.5 loses only when t >= 13 — beyond the old g=0..10 grid,
    # which valued this book at riskless +0.05/contract. Interval ledger sees it.
    fills = [(12.5, "no", 10, 0.95)]
    assert ledger.worst_case(fills) == -9.5  # t>=13: lose 0.95 x 10


def test_worst_case_matches_bruteforce_on_random_books():
    rng = random.Random(20260725)
    lines = [0.5, 1.5, 2.5, 6.5, 7.5, 8.5, 9.5, 12.5]
    for _ in range(200):
        fills = [(rng.choice(lines), rng.choice(["yes", "no"]),
                  rng.randint(1, 20), round(rng.uniform(0.05, 0.95), 2))
                 for _ in range(rng.randint(1, 6))]
        assert abs(ledger.worst_case(fills) - _worst_case_bruteforce(fills)) < 1e-9


def test_max_contracts_binds_beyond_old_grid():
    # Candidate NO at 12.5, price 0.95: worst case is -0.95/contract (t>=13).
    # budget 9.5 -> exactly 10 contracts. The old grid said unlimited (bound None).
    assert ledger.max_contracts([], 12.5, "no", 0.95, 9.5) == 10


def test_max_contracts_agrees_with_bruteforce_greedy():
    rng = random.Random(99)
    lines = [0.5, 2.5, 8.5, 12.5]
    for _ in range(100):
        fills = [(rng.choice(lines), rng.choice(["yes", "no"]),
                  rng.randint(1, 5), round(rng.uniform(0.1, 0.9), 2))
                 for _ in range(rng.randint(0, 3))]
        line, side = rng.choice(lines), rng.choice(["yes", "no"])
        price, budget = round(rng.uniform(0.1, 0.9), 2), rng.uniform(1, 20)
        n = ledger.max_contracts(fills, line, side, price, budget)

        # Skip Property checks if existing fills already exceed budget (fail-closed case)
        existing_worst = _worst_case_bruteforce(fills) if fills else 0.0
        if existing_worst >= -budget - 1e-9:
            # Property 1: the returned size respects the budget.
            assert _worst_case_bruteforce(fills + [(line, side, n, price)]) >= -budget - 1e-6
            # Property 2 (maximality): with price in (0,1) every candidate has a
            # losing region, so one more contract must break the budget.
            assert _worst_case_bruteforce(fills + [(line, side, n + 1, price)]) < -budget + 1e-6
