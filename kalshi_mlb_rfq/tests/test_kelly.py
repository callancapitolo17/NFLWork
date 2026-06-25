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


# ---------------------------------------------------------------------------
# Legacy model-path tests (USE_MODEL=True shape: outcome_vec, no cov_return)
# ---------------------------------------------------------------------------

import numpy as np


def _make_outcome_vec(p: float, n: int = 10_000, seed: int = 42) -> np.ndarray:
    """Binary (0/1) sample path with hit-rate ~p."""
    rng = np.random.default_rng(seed)
    return rng.binomial(1, p, size=n).astype(float)


def test_legacy_position_no_keyerror():
    """Legacy-shape position (outcome_vec, no cov_return) must not crash."""
    new_vec = _make_outcome_vec(0.30)
    pos_vec = _make_outcome_vec(0.30, seed=99)
    legacy_pos = [{"outcome_vec": pos_vec, "contracts": 50, "effective_price": 0.20}]
    n = kelly.kelly_size_combo(
        outcome_vec=new_vec, blended_fair=0.30, existing_positions=legacy_pos,
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    assert isinstance(n, int) and n >= 0


def test_legacy_perfectly_correlated_downsizes():
    """A perfectly-correlated legacy position must shrink size vs no-position base."""
    # Use the same vector for both new bet and existing position → max correlation.
    shared_vec = _make_outcome_vec(0.30)
    base = kelly.kelly_size_combo(
        outcome_vec=shared_vec, blended_fair=0.30, existing_positions=[],
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    legacy_pos = [{"outcome_vec": shared_vec, "contracts": 50, "effective_price": 0.20}]
    corr = kelly.kelly_size_combo(
        outcome_vec=shared_vec, blended_fair=0.30, existing_positions=legacy_pos,
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    assert corr < base


def test_legacy_position_outcome_vec_none_no_crash():
    """Legacy-shape position present but new outcome_vec=None → cov=0 (independent), no crash."""
    pos_vec = _make_outcome_vec(0.30)
    legacy_pos = [{"outcome_vec": pos_vec, "contracts": 50, "effective_price": 0.20}]
    # With no new outcome_vec, falls back to cov=0 → result equals independent base.
    base = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=[],
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    n = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=legacy_pos,
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    assert n == base  # cov=0 → identical to no-position case
