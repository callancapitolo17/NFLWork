"""Conditional Kelly sizing for combo bets.

Sizing uses the *blended fair* (median of model + sportsbooks) for the new
bet's mean and variance — keeping EV gate and Kelly gate consistent.
Covariance with existing positions is supplied by the caller as a precomputed
`cov_return` float on each position dict (book-implied, from the correlation
engine). Model sample paths (`outcome_vec`) are no longer used.

Return-space (small-f Taylor) approximation; underbets exact binary Kelly
by ~10% — conservative, fine for v1.
"""

import math

import numpy as np


def kelly_size_combo(
    outcome_vec,
    blended_fair: float,
    existing_positions: list[dict],
    effective_price: float,
    bankroll: float,
    kelly_fraction: float,
) -> int:
    """Return contract count to take.

    Args:
        outcome_vec: ignored (kept for call-site compatibility; pass None).
        blended_fair: median(model, *book fairs) — drives mu_new and var_new.
        existing_positions: list of {cov_return, contracts, effective_price}.
            `cov_return` is the precomputed return-space covariance between
            the existing position and the new bet (computed by the correlation
            engine in the caller).
        effective_price: post-fee price per contract for this new bet.
        bankroll: dollars.
        kelly_fraction: e.g., 0.25 for quarter-Kelly.

    Returns:
        Floored, non-negative contract count.
    """
    if effective_price <= 0 or effective_price >= 1:
        return 0
    p = float(blended_fair)
    if not 0.0 < p < 1.0 or p <= effective_price:
        return 0  # -EV or degenerate fair

    # Binary outcome at price: both moments derived from blended fair.
    mu_new = (p - effective_price) / effective_price
    var_new = p * (1.0 - p) / (effective_price ** 2)
    if var_new <= 1e-12:
        return 0

    base_frac = mu_new / var_new   # single-bet Kelly fraction (return-space approx)

    if not existing_positions:
        contracts = math.floor(kelly_fraction * base_frac * bankroll / effective_price)
        return max(0, contracts)

    # Book-implied covariance: each position carries a precomputed return-space
    # covariance with the new bet (correlation.cov_returns). Full-Kelly residual
    # after the existing book, then apply the Kelly fraction.
    np_cov = []
    f_placed = []
    for pos in existing_positions:
        cov = float(pos["cov_return"])
        pos_price = float(pos["effective_price"])
        if pos_price <= 0 or pos_price >= 1:
            continue
        np_cov.append(cov)
        f_placed.append(kelly_fraction * pos["contracts"] * pos_price / bankroll)

    if not f_placed:
        contracts = math.floor(kelly_fraction * base_frac * bankroll / effective_price)
        return max(0, contracts)

    f_full_new = max(0.0, (mu_new - float(np.dot(np.array(np_cov), np.array(f_placed)))) / var_new)
    contracts = math.floor(kelly_fraction * f_full_new * bankroll / effective_price)
    return max(0, contracts)
