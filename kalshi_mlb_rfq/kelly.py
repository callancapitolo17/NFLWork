"""Conditional Kelly sizing for combo bets.

Sizing uses the *blended fair* (median of model + sportsbooks) for the new
bet's mean and variance — keeping EV gate and Kelly gate consistent.

Two position dict shapes are supported:
- Book path (USE_MODEL=False): {cov_return, contracts, effective_price}
  `cov_return` is the precomputed return-space covariance from the correlation
  engine.
- Legacy model path (USE_MODEL=True): {outcome_vec, contracts, effective_price}
  Covariance is computed on-the-fly from the numpy sample arrays.

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
        outcome_vec: numpy array of outcome samples for the new bet (model
            path), or None (book path — ignored).
        blended_fair: median(model, *book fairs) — drives mu_new and var_new.
        existing_positions: list of position dicts. Two shapes accepted:
            - {cov_return, contracts, effective_price}  (book path)
            - {outcome_vec, contracts, effective_price} (legacy model path)
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

    # Precompute new-bet return series once (only needed for legacy model path).
    new_centered = None
    if outcome_vec is not None:
        new_returns = np.asarray(outcome_vec, dtype=float) / effective_price - 1.0
        new_centered = new_returns - new_returns.mean()

    np_cov = []
    f_placed = []
    for pos in existing_positions:
        pos_price = float(pos["effective_price"])
        if pos_price <= 0 or pos_price >= 1:
            continue

        if "cov_return" in pos:
            # Book path: caller precomputed return-space covariance.
            cov = float(pos["cov_return"])
        elif "outcome_vec" in pos and new_centered is not None:
            # Legacy model path: compute from sample arrays.
            pos_vec = np.asarray(pos["outcome_vec"], dtype=float)
            if len(pos_vec) != len(new_centered):
                continue  # mismatched sample count — skip
            pos_returns = pos_vec / pos_price - 1.0
            cov = float(((pos_returns - pos_returns.mean()) * new_centered).mean())
        else:
            # No covariance info available — treat as independent.
            cov = 0.0

        np_cov.append(cov)
        f_placed.append(kelly_fraction * pos["contracts"] * pos_price / bankroll)

    if not f_placed:
        contracts = math.floor(kelly_fraction * base_frac * bankroll / effective_price)
        return max(0, contracts)

    f_full_new = max(0.0, (mu_new - float(np.dot(np.array(np_cov), np.array(f_placed)))) / var_new)
    contracts = math.floor(kelly_fraction * f_full_new * bankroll / effective_price)
    return max(0, contracts)
