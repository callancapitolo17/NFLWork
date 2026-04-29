import numpy as np

from kalshi_mlb_rfq import kelly


def test_single_kelly_no_existing_positions():
    contracts = kelly.kelly_size_combo(
        outcome_vec=np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0]),  # 30% hit rate
        existing_positions=[],
        effective_price=0.26,
        bankroll=1000.0,
        kelly_fraction=0.25,
    )
    assert 40 <= contracts <= 65   # ballpark check


def test_kelly_zero_for_negative_edge():
    # Fair 0.20 vs price 0.30 → -EV → 0 contracts
    contracts = kelly.kelly_size_combo(
        outcome_vec=np.array([1, 1, 0, 0, 0, 0, 0, 0, 0, 0]),  # 20% hit rate
        existing_positions=[],
        effective_price=0.30,
        bankroll=1000.0,
        kelly_fraction=0.25,
    )
    assert contracts == 0


def test_conditional_kelly_shrinks_for_correlated_position():
    n = 1000
    rng = np.random.default_rng(42)
    outcome = (rng.random(n) < 0.30).astype(int)

    single = kelly.kelly_size_combo(
        outcome_vec=outcome, existing_positions=[],
        effective_price=0.25, bankroll=1000.0, kelly_fraction=0.25,
    )

    placed = [{"outcome_vec": outcome.copy(), "contracts": single, "effective_price": 0.25}]
    conditional = kelly.kelly_size_combo(
        outcome_vec=outcome, existing_positions=placed,
        effective_price=0.25, bankroll=1000.0, kelly_fraction=0.25,
    )
    # Quarter-Kelly with one perfectly-correlated existing position: conditional should
    # shrink (math says by ~kelly_fraction = 25% — but with quarter-scaled f_placed the
    # actual reduction per cycle is ~7-15%). Just verify it's strictly smaller.
    assert conditional < single
    assert conditional > 0
