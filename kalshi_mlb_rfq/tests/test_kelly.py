from kalshi_mlb_rfq import kelly


BANK = 1000.0
KF = 0.25


def test_no_positions_single_bet_kelly():
    # +EV: fair 0.30 vs price 0.20
    n = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=[],
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    assert n > 0


def test_negative_ev_returns_zero():
    n = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.20, existing_positions=[],
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    assert n == 0


def test_degenerate_price_returns_zero():
    assert kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.5, existing_positions=[],
        effective_price=0.0, bankroll=BANK, kelly_fraction=KF) == 0
    assert kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.5, existing_positions=[],
        effective_price=1.0, bankroll=BANK, kelly_fraction=KF) == 0


def test_perfectly_correlated_position_downsizes():
    """A held, perfectly correlated position must shrink the new size vs the
    no-position case — the core anti-overbet guarantee."""
    base = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=[],
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    # comonotonic: cov_return positive and large
    pos = [{"cov_return": (0.30 - 0.30 * 0.30) / (0.20 * 0.20),
            "contracts": 50, "effective_price": 0.20}]
    corr = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=pos,
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    assert corr < base


def test_independent_position_matches_base():
    base = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=[],
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    pos = [{"cov_return": 0.0, "contracts": 50, "effective_price": 0.20}]
    same = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=pos,
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    assert same == base


def test_negative_correlation_upsizes():
    base = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=[],
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    pos = [{"cov_return": -0.5, "contracts": 50, "effective_price": 0.20}]
    hedge = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=pos,
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    assert hedge >= base
