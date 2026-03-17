"""
Conditional Kelly criterion sizing for the Kalshi MM bot.

Ports the multivariate Kelly logic from Tools.R to Python/numpy.
Uses game simulation samples (exported from CBB.R to cbb.duckdb)
to compute covariances between correlated positions on the same game.
"""

import math
import logging
import numpy as np
import duckdb
from config import (
    CBB_DB_PATH, BANKROLL, KELLY_FRACTION,
    MIN_KELLY_CONTRACTS, MAX_KELLY_CONTRACTS, CONTRACT_SIZE
)

import db as _db

log = logging.getLogger(__name__)

# In-memory cache: game_id -> DataFrame with (home_margin_h1, total_h1)
_sample_cache = {}


def clear_sample_cache():
    """Clear cached samples (call on pipeline refresh)."""
    _sample_cache.clear()
    log.info("Kelly sample cache cleared")


def load_game_samples(game_id):
    """Load simulation samples for a game from cbb.duckdb.

    Returns numpy array of shape (n_sims, 2) with columns
    [home_margin_h1, total_h1], or None if unavailable.
    """
    if game_id in _sample_cache:
        return _sample_cache[game_id]

    try:
        con = duckdb.connect(str(CBB_DB_PATH), read_only=True)
        try:
            result = con.execute(
                "SELECT home_margin_h1, total_h1 FROM cbb_game_samples "
                "WHERE game_id = ? ORDER BY sim_idx",
                [game_id]
            ).fetchnumpy()
        finally:
            con.close()

        if result is None or len(result["home_margin_h1"]) == 0:
            _sample_cache[game_id] = None
            return None

        samples = np.column_stack([
            result["home_margin_h1"].astype(np.float64),
            result["total_h1"].astype(np.float64),
        ])
        _sample_cache[game_id] = samples
        return samples

    except Exception as e:
        log.warning("Failed to load samples for %s: %s", game_id, e)
        _sample_cache[game_id] = None
        return None


def evaluate_outcomes(samples, market_type, side, line=None):
    """Evaluate binary outcomes for each simulation row.

    Args:
        samples: (n_sims, 2) array with [home_margin_h1, total_h1].
        market_type: "spreads", "totals", or "moneyline".
        side: "home"/"away" for spreads/ML, "over"/"under" for totals,
              "tie" for ML tie.
        line: spread or total line (float). None for ML.

    Returns:
        (n_sims,) float array with 1.0 (hit), 0.0 (miss), NaN (push).
    """
    margin = samples[:, 0]
    total = samples[:, 1]
    n = len(margin)
    out = np.empty(n, dtype=np.float64)

    if market_type == "spreads":
        # Home -3.5 → home covers if margin > 3.5
        # line is stored as the home spread (negative = favorite)
        threshold = -line
        if side == "home":
            out[:] = np.where(margin > threshold, 1.0,
                              np.where(margin < threshold, 0.0, np.nan))
        else:  # away
            out[:] = np.where(margin < threshold, 1.0,
                              np.where(margin > threshold, 0.0, np.nan))

    elif market_type == "totals":
        if side == "over":
            out[:] = np.where(total > line, 1.0,
                              np.where(total < line, 0.0, np.nan))
        else:  # under
            out[:] = np.where(total < line, 1.0,
                              np.where(total > line, 0.0, np.nan))

    elif market_type == "moneyline":
        if side == "home":
            out[:] = np.where(margin > 0, 1.0,
                              np.where(margin < 0, 0.0, np.nan))
        elif side == "away":
            out[:] = np.where(margin < 0, 1.0,
                              np.where(margin > 0, 0.0, np.nan))
        else:  # tie
            out[:] = np.where(margin == 0, 1.0, 0.0)
    else:
        out[:] = np.nan

    return out


def build_outcome_matrix(samples, positions):
    """Build outcome matrix for multiple positions on the same game.

    Args:
        samples: (n_sims, 2) array from load_game_samples().
        positions: list of dicts with market_type, side, line.

    Returns:
        (n_valid, n_positions) float array with push rows dropped,
        or None if <30 valid rows.
    """
    n_pos = len(positions)
    cols = []
    for pos in positions:
        col = evaluate_outcomes(
            samples, pos["market_type"], pos["side"], pos.get("line")
        )
        cols.append(col)

    matrix = np.column_stack(cols)

    # Drop rows with any NaN (pushes on any leg)
    valid_mask = ~np.isnan(matrix).any(axis=1)
    matrix = matrix[valid_mask]

    if len(matrix) < 30:
        return None
    return matrix


def _compute_covariance(outcome_matrix, fair_probs, kalshi_prices):
    """Compute return covariance matrix from binary outcome matrix.

    Args:
        outcome_matrix: (n_sims, n_positions) of 0/1.
        fair_probs: array of model probabilities.
        kalshi_prices: array of Kalshi prices in cents.

    Returns:
        Covariance matrix Σ, or None if ill-conditioned.
    """
    n = outcome_matrix.shape[1]
    if n < 2:
        # Single position — return 1x1 variance
        p = fair_probs[0]
        price = kalshi_prices[0]
        b = (100 - price) / price  # decimal odds - 1
        sigma = math.sqrt(p * (1 - p)) * (b + 1)
        return np.array([[sigma ** 2 + 0.01]])

    R = np.corrcoef(outcome_matrix, rowvar=False)
    if np.any(np.isnan(R)):
        return None

    # σ_i = sqrt(p_i * (1-p_i)) * (b_i + 1)
    sigmas = np.array([
        math.sqrt(p * (1 - p)) * ((100 - price) / price + 1)
        for p, price in zip(fair_probs, kalshi_prices)
    ])

    Sigma = R * np.outer(sigmas, sigmas)

    # Ridge regularization
    Sigma += np.eye(n) * 0.01

    if np.linalg.cond(Sigma) > 100:
        return None

    return Sigma


def kelly_size_single(fair_prob, kalshi_price_cents, bankroll, kelly_mult):
    """Single-bet Kelly for an uncorrelated position.

    Returns integer contract count, clamped to [MIN, MAX].
    """
    price = kalshi_price_cents / 100.0
    if fair_prob <= price or price <= 0 or price >= 1:
        return 0

    b = (1 - price) / price  # decimal odds - 1
    f = (b * fair_prob - (1 - fair_prob)) / b
    if f <= 0:
        return 0

    # Convert fraction to contracts: each contract costs price dollars, pays $1
    cost_per_contract = price
    contracts = f * kelly_mult * bankroll / cost_per_contract
    contracts = max(MIN_KELLY_CONTRACTS, min(MAX_KELLY_CONTRACTS, round(contracts)))
    return int(contracts)


def conditional_kelly_sizes(new_positions, placed_positions, samples,
                            bankroll, kelly_mult):
    """Compute Kelly sizes for new positions, accounting for placed exposure.

    Args:
        new_positions: list of position dicts (market_type, side, line,
                       fair_prob, kalshi_price).
        placed_positions: list of position dicts with additional 'size' field.
        samples: (n_sims, 2) array from load_game_samples().
        bankroll: total bankroll in dollars.
        kelly_mult: fractional Kelly multiplier.

    Returns:
        list of int contract sizes, one per new_position.
        Falls back to single-bet Kelly with rho scaling if covariance fails.
    """
    n_new = len(new_positions)
    n_placed = len(placed_positions)
    all_positions = placed_positions + new_positions

    # Build outcome matrix for all positions
    outcome_matrix = build_outcome_matrix(samples, all_positions)

    if outcome_matrix is not None:
        fair_probs = np.array([p["fair_prob"] for p in all_positions])
        kalshi_prices = np.array([p["kalshi_price"] for p in all_positions])
        Sigma = _compute_covariance(outcome_matrix, fair_probs, kalshi_prices)

        if Sigma is not None:
            # EV vector: edge as fraction of cost
            mu = np.array([
                (p["fair_prob"] - p["kalshi_price"] / 100.0)
                for p in all_positions
            ])

            if n_placed > 0:
                # Conditional Kelly: f_new* = Σ_nn⁻¹ × (μ_new − Σ_np × f_placed)
                p_idx = slice(0, n_placed)
                n_idx = slice(n_placed, n_placed + n_new)

                Sigma_nn = Sigma[n_idx, n_idx]
                Sigma_np = Sigma[n_idx, p_idx]

                # f_placed = contracts * cost / (kelly_mult * bankroll)
                f_placed = np.array([
                    p["size"] * (p["kalshi_price"] / 100.0) / (kelly_mult * bankroll)
                    for p in placed_positions
                ])

                mu_new = mu[n_idx]
                adjusted_mu = mu_new - Sigma_np @ f_placed

                try:
                    f_new = np.linalg.solve(Sigma_nn, adjusted_mu)
                except np.linalg.LinAlgError:
                    f_new = None
            else:
                # Standard multivariate Kelly: f* = Σ⁻¹ × μ
                try:
                    f_star = np.linalg.solve(Sigma, mu)
                    f_new = f_star[n_placed:]  # just the new positions
                except np.linalg.LinAlgError:
                    f_new = None

            if f_new is not None:
                f_new = np.maximum(f_new, 0)
                sizes = []
                for i, pos in enumerate(new_positions):
                    cost = pos["kalshi_price"] / 100.0
                    contracts = f_new[i] * kelly_mult * bankroll / cost
                    contracts = max(MIN_KELLY_CONTRACTS,
                                    min(MAX_KELLY_CONTRACTS, round(contracts)))
                    sizes.append(int(contracts))
                return sizes

    # --- Fallback: per-bet average ρ scaling ---
    return _fallback_rho_scaling(new_positions, placed_positions, samples,
                                 bankroll, kelly_mult)


def _fallback_rho_scaling(new_positions, placed_positions, samples,
                          bankroll, kelly_mult):
    """Fallback sizing using average correlation scaling.

    scale_i = 1 / sqrt(1 + (n-1) * avg_ρ_i)
    If placed bet has ρ > 0.90 with new bet → size = 0.
    """
    n_new = len(new_positions)
    n_placed = len(placed_positions)
    all_positions = placed_positions + new_positions

    outcome_matrix = build_outcome_matrix(samples, all_positions)

    if outcome_matrix is not None and outcome_matrix.shape[1] >= 2:
        R = np.corrcoef(outcome_matrix, rowvar=False)
        if not np.any(np.isnan(R)):
            sizes = []
            n_all = len(all_positions)
            for j in range(n_new):
                i = n_placed + j  # index in full group
                pos = new_positions[j]

                # Check if any placed bet is highly correlated
                if n_placed > 0:
                    max_placed_rho = np.max(R[i, :n_placed])
                    if max_placed_rho > 0.90:
                        sizes.append(0)
                        continue

                # Average correlation with all other positions
                other_indices = [k for k in range(n_all) if k != i]
                avg_rho = np.mean(R[i, other_indices])

                denom = 1 + (n_all - 1) * avg_rho
                if denom <= 0:
                    scale = 1.5
                else:
                    scale = min(1.5, 1.0 / math.sqrt(denom))
                scale = max(scale, 0)

                base = kelly_size_single(
                    pos["fair_prob"], pos["kalshi_price"],
                    bankroll, kelly_mult
                )
                adjusted = max(MIN_KELLY_CONTRACTS,
                               min(MAX_KELLY_CONTRACTS, round(base * scale)))
                sizes.append(int(adjusted))
            return sizes

    # Final fallback: independent single-bet Kelly
    return [
        kelly_size_single(p["fair_prob"], p["kalshi_price"], bankroll, kelly_mult)
        for p in new_positions
    ]


def kelly_size_for_quote(market, side, placed_positions, bankroll, kelly_mult):
    """Compute Kelly size for one side of a quote on a market.

    Args:
        market: quotable market dict (must have game_id, market_type, fair_prob, etc.)
        side: "bid" (buying YES) or "ask" (selling YES = buying NO)
        placed_positions: list of placed position dicts for this game
        bankroll: total bankroll
        kelly_mult: fractional Kelly

    Returns:
        int contract count
    """
    game_id = market.get("game_id")
    if not game_id:
        return CONTRACT_SIZE  # no game_id → can't load samples

    samples = load_game_samples(game_id)
    if samples is None:
        return CONTRACT_SIZE  # no samples → fall back to fixed size

    # Build the new position dict
    if side == "bid":
        # Buying YES: we profit when the event happens
        pos = _market_to_position(market, "yes")
        pos["kalshi_price"] = market.get("book_bid", int(market["fair_prob"] * 100))
    else:
        # Selling YES (buying NO): we profit when the event doesn't happen
        pos = _market_to_position(market, "no")
        pos["kalshi_price"] = 100 - market.get("book_ask",
                                                int((1 - market["fair_prob"]) * 100))

    if pos["kalshi_price"] <= 0 or pos["kalshi_price"] >= 100:
        return CONTRACT_SIZE

    if not placed_positions:
        return kelly_size_single(pos["fair_prob"], pos["kalshi_price"],
                                 bankroll, kelly_mult)

    sizes = conditional_kelly_sizes([pos], placed_positions, samples,
                                    bankroll, kelly_mult)
    return sizes[0] if sizes else CONTRACT_SIZE


def _market_to_position(market, yes_or_no):
    """Convert a quotable market dict to a Kelly position dict."""
    mt = market["market_type"]

    if yes_or_no == "yes":
        fair_prob = market["fair_prob"]
        if mt == "spreads":
            side = "home" if market.get("is_home_contract") else "away"
        elif mt == "totals":
            side = "over"
        elif mt == "moneyline":
            # Determine from contract_team or yes_sub_title
            ct = market.get("contract_team", "")
            if ct == market.get("home_team"):
                side = "home"
            elif ct == market.get("away_team"):
                side = "away"
            else:
                side = "tie"
        else:
            side = "home"
    else:
        # NO side: flip the probability and side
        fair_prob = 1 - market["fair_prob"]
        if mt == "spreads":
            side = "away" if market.get("is_home_contract") else "home"
        elif mt == "totals":
            side = "under"
        elif mt == "moneyline":
            ct = market.get("contract_team", "")
            if ct == market.get("home_team"):
                side = "away"
            elif ct == market.get("away_team"):
                side = "home"
            else:
                # Tie NO is complex (home or away wins) — use simple Kelly
                return {
                    "market_type": mt, "side": "home",
                    "line": None, "fair_prob": fair_prob,
                    "kalshi_price": 50,
                }
        else:
            side = "away"

    return {
        "market_type": mt,
        "side": side,
        "line": market.get("strike"),
        "fair_prob": fair_prob,
        "kalshi_price": 0,  # caller fills this in
    }


def get_placed_positions_for_game(game_key, game_id):
    """Build placed position dicts for Kelly from DB positions on this game.

    Returns list of position dicts with market_type, side, line, fair_prob,
    kalshi_price, and size fields.
    """
    positions = _db.get_all_positions()
    placed = []
    for pos in positions:
        pos_home = pos.get("home_team", "")
        pos_away = pos.get("away_team", "")
        if (pos_home, pos_away) != game_key:
            continue
        net = pos.get("net_yes", 0)
        if net == 0:
            continue

        mt = pos.get("market_type", "spreads")
        if mt == "spreads":
            side = "home" if net > 0 else "away"
        elif mt == "totals":
            side = "over" if net > 0 else "under"
        elif mt == "moneyline":
            side = "home" if net > 0 else "away"
        else:
            continue

        placed.append({
            "market_type": mt,
            "side": side,
            "line": pos.get("line_value"),
            "fair_prob": pos.get("fair_prob", 0.5),
            "kalshi_price": int(pos.get("avg_entry_price", 50)),
            "size": abs(net),
        })
    return placed
