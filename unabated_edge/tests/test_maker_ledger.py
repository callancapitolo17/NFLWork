# unabated_edge/tests/test_maker_ledger.py
from unabated_edge.maker import ledger


def test_empty_book_is_flat():
    assert ledger.worst_case([]) == 0.0
    assert ledger.pnl_grid([]) == [0.0] * (ledger.G_MAX + 1)


def test_single_yes_fill_worst_case():
    # 40x YES Over-2.5 @ 0.42: loses 16.8 when g<=2, wins 23.2 when g>=3
    fills = [(2.5, "yes", 40, 0.42)]
    grid = ledger.pnl_grid(fills)
    assert round(grid[0], 2) == -16.8 and round(grid[2], 2) == -16.8
    assert round(grid[3], 2) == 23.2
    assert round(ledger.worst_case(fills), 2) == -16.8


def test_offset_under25_over15_worst_and_middle():
    # From the design conversation: 20x NO Over-2.5 @0.55 + 20x YES Over-1.5 @0.70
    fills = [(2.5, "no", 20, 0.55), (1.5, "yes", 20, 0.70)]
    grid = ledger.pnl_grid(fills)
    assert round(grid[0], 2) == -5.0          # low-goal worlds
    assert round(grid[2], 2) == 15.0          # the middle: both win at exactly 2
    assert round(grid[3], 2) == -5.0          # high-goal worlds
    assert round(ledger.worst_case(fills), 2) == -5.0


def test_grid_constant_above_top_rung():
    # No Kalshi rung above 5.5 -> pnl must be identical for g=6..10
    fills = [(5.5, "yes", 10, 0.05), (2.5, "no", 7, 0.6)]
    grid = ledger.pnl_grid(fills)
    assert len(set(round(v, 9) for v in grid[6:])) == 1


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
