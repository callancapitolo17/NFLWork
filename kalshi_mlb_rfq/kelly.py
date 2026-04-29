"""Conditional Kelly sizing on combo outcome vectors.

Adapted from kalshi_mm/kelly.py. The new bet's outcome is a boolean vector
(1 if combo hits in that simulation path, else 0). Existing positions on the
same game contribute their own outcome vectors; we compute joint covariance
and apply the conditional Kelly formula:

    f_new* = Σ_nn⁻¹ × (μ_new − Σ_np × f_placed)

where Σ_nn is the new bet's variance, Σ_np is the covariance vector with
existing positions, and f_placed is each existing position's stake fraction.
"""

import math

import numpy as np


def _single_bet_kelly_fraction(p: float, effective_price: float) -> float:
    """Standard Kelly: (bp - q) / b where b = (1 - price) / price."""
    if effective_price <= 0 or effective_price >= 1:
        return 0.0
    b = (1 - effective_price) / effective_price
    q = 1 - p
    f = (b * p - q) / b
    return max(0.0, f)


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
    p = float(outcome_vec.mean())
    base_frac = _single_bet_kelly_fraction(p, effective_price)
    if base_frac <= 0:
        return 0

    if not existing_positions:
        contracts = math.floor(kelly_fraction * base_frac * bankroll / effective_price)
        return max(0, contracts)

    # Conditional Kelly. Build the new outcome's centered vector and Σ_np.
    n_samples = len(outcome_vec)
    new_centered = outcome_vec - p
    var_new = float(np.var(outcome_vec))
    if var_new <= 1e-12:
        return 0

    # f_placed in stake-fraction units (contracts * price / bankroll).
    f_placed = []
    np_cov_terms = []
    for pos in existing_positions:
        pos_vec = pos["outcome_vec"]
        if len(pos_vec) != n_samples:
            continue
        pos_mean = float(pos_vec.mean())
        cov = float(((pos_vec - pos_mean) * new_centered).mean())
        f_placed.append(pos["contracts"] * pos["effective_price"] / bankroll)
        np_cov_terms.append(cov)

    if not f_placed:
        contracts = math.floor(kelly_fraction * base_frac * bankroll / effective_price)
        return max(0, contracts)

    f_placed = np.array(f_placed)
    np_cov = np.array(np_cov_terms)

    mu_new = p - effective_price  # excess return over price
    f_new_star = (mu_new - np.dot(np_cov, f_placed)) / var_new
    f_new_star = max(0.0, f_new_star)

    contracts = math.floor(kelly_fraction * f_new_star * bankroll / effective_price)
    return max(0, contracts)
