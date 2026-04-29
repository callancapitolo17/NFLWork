"""Conditional Kelly sizing on combo outcome vectors.

Adapted from kalshi_mm/kelly.py. Uses return-space (small-f Taylor) approximation
which underbets vs exact binary Kelly by ~10% — conservative, fine for v1.

For a new bet with mean outcome μ_new (in returns) and existing positions on the
same game, conditional Kelly:
    f_new* = Σ_nn⁻¹ × (μ_new − Σ_np × f_placed)
"""

import math

import numpy as np


def kelly_size_combo(
    outcome_vec: np.ndarray,
    existing_positions: list[dict],
    effective_price: float,
    bankroll: float,
    kelly_fraction: float,
) -> int:
    """Return contract count to take.

    Args:
        outcome_vec: 0/1 array, len = sample count.
        existing_positions: list of {outcome_vec, contracts, effective_price}.
        effective_price: post-fee price per contract for this new bet.
        bankroll: dollars.
        kelly_fraction: e.g., 0.25 for quarter-Kelly.

    Returns:
        Floored, non-negative contract count.
    """
    if effective_price <= 0 or effective_price >= 1:
        return 0
    p = float(outcome_vec.mean())
    if p <= effective_price:
        return 0  # -EV

    # Return-space formulation: r_i = X_i / price - 1 (multi-asset Kelly's natural units).
    new_returns = outcome_vec / effective_price - 1.0
    mu_new = float(new_returns.mean())            # = (p - price) / price
    var_new = float(np.var(new_returns))           # = p(1-p) / price²
    if var_new <= 1e-12:
        return 0

    base_frac = mu_new / var_new   # single-bet Kelly fraction (return-space approx)

    if not existing_positions:
        contracts = math.floor(kelly_fraction * base_frac * bankroll / effective_price)
        return max(0, contracts)

    f_placed = []
    np_cov_terms = []
    new_centered = new_returns - mu_new
    for pos in existing_positions:
        pos_vec = pos["outcome_vec"]
        if len(pos_vec) != len(outcome_vec):
            continue
        pos_price = float(pos["effective_price"])
        if pos_price <= 0 or pos_price >= 1:
            continue
        pos_returns = pos_vec / pos_price - 1.0
        pos_mean = float(pos_returns.mean())
        cov = float(((pos_returns - pos_mean) * new_centered).mean())
        np_cov_terms.append(cov)
        # Stake fraction the existing position represents (kelly_fraction-scaled bankroll).
        f_placed.append(
            kelly_fraction * pos["contracts"] * pos_price / bankroll
        )

    if not f_placed:
        contracts = math.floor(kelly_fraction * base_frac * bankroll / effective_price)
        return max(0, contracts)

    f_placed_arr = np.array(f_placed)
    np_cov = np.array(np_cov_terms)

    # Full-Kelly target after accounting for existing positions:
    f_full_new = max(0.0, (mu_new - float(np.dot(np_cov, f_placed_arr))) / var_new)

    # Apply Kelly fraction to the residual.
    contracts = math.floor(kelly_fraction * f_full_new * bankroll / effective_price)
    return max(0, contracts)
