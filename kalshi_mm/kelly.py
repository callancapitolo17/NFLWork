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
)

import db as _db

log = logging.getLogger(__name__)

# In-memory cache: game_id -> DataFrame with (home_margin_h1, total_h1)
_sample_cache = {}

# Per-cycle cache for positions from Kalshi API (source of truth)
_positions_cache = None
_positions_cache_time = 0
_positions_markets_ref = None  # quotable_markets for ticker→metadata mapping
_POSITIONS_CACHE_TTL = 5  # seconds — prevents taker from hammering API


def _load_kalshi_positions():
    """Load positions from Kalshi API and map to internal format.

    Returns list of dicts matching db.get_all_positions() format, or None on failure.
    """
    import orders as _orders
    raw = _orders.get_kalshi_positions()
    if raw is None:
        log.warning("Kalshi positions API failed — cannot size")
        return None

    # Build ticker→market lookup from quotable_markets
    market_lookup = {}
    if _positions_markets_ref:
        for m in _positions_markets_ref:
            market_lookup[m["ticker"]] = m

    positions = []
    unmapped = []
    for kp in raw:
        ticker = kp.get("ticker", "")
        position = int(float(kp.get("position_fp", "0")))
        if position == 0:
            continue

        market = market_lookup.get(ticker)
        if not market:
            unmapped.append(ticker)
            continue

        # Compute avg_entry_price from Kalshi's exposure data
        exposure_dollars = float(kp.get("market_exposure_dollars", "0"))
        avg_price = (exposure_dollars / abs(position) * 100) if position != 0 else 0

        mt = market.get("market_type", "spreads")
        strike = market.get("strike")
        if mt == "spreads" and strike is not None:
            line_value = -strike if market.get("is_home_contract") else strike
        elif mt == "totals" and strike is not None:
            line_value = strike
        else:
            line_value = None

        positions.append({
            "ticker": ticker,
            "event_ticker": market.get("event_ticker", ""),
            "home_team": market.get("home_team", ""),
            "away_team": market.get("away_team", ""),
            "market_type": mt,
            "line_value": line_value,
            "net_yes": position,
            "avg_entry_price": avg_price,
            "fair_prob": market.get("fair_prob", 0.5),
            "contract_team": market.get("contract_team", ""),
        })

    if unmapped:
        log.info("Kalshi positions: %d mapped, %d unmapped (settled/removed games)",
                 len(positions), len(unmapped))
    return positions


def clear_positions_cache(current_markets=None):
    """Mark positions cache as stale. Next read re-fetches from Kalshi API.

    Args:
        current_markets: quotable_markets list for ticker→metadata mapping.
            Only updated if provided (not None).
    """
    global _positions_cache, _positions_cache_time, _positions_markets_ref
    _positions_cache = None
    _positions_cache_time = 0
    if current_markets is not None:
        _positions_markets_ref = current_markets


_positions_api_ok = True  # Track if last API call succeeded


def _ensure_positions_loaded():
    """Lazy-load positions from Kalshi if cache is stale (>TTL seconds old)."""
    global _positions_cache, _positions_cache_time, _positions_api_ok
    import time as _t
    now = _t.time()
    if _positions_cache is not None and (now - _positions_cache_time) < _POSITIONS_CACHE_TTL:
        return  # Cache is fresh
    loaded = _load_kalshi_positions()
    if loaded is not None:
        _positions_cache = loaded
        _positions_cache_time = now
        _positions_api_ok = True
    else:
        _positions_api_ok = False
        if _positions_cache is None:
            _positions_cache = []  # API failed and no prior cache — empty


def get_cached_positions():
    """Return the current positions cache (for use by risk.py exposure cap)."""
    _ensure_positions_loaded()
    return _positions_cache


def positions_api_healthy():
    """Check if the last Kalshi positions API call succeeded."""
    _ensure_positions_loaded()
    return _positions_api_ok


def get_net_position(ticker):
    """Get net YES position for a ticker from Kalshi positions.

    Returns 0 if no position exists.
    """
    _ensure_positions_loaded()
    for pos in _positions_cache:
        if pos.get("ticker") == ticker:
            return pos.get("net_yes", 0)
    return 0


def prewarm_sample_cache():
    """Bulk-load all game samples into cache in one query.

    Replaces clear_sample_cache() — clears stale data and reloads fresh.
    Called at startup and after prediction refresh to avoid cold-cache
    Kelly cycles (~15 min → seconds).
    """
    _sample_cache.clear()
    try:
        con = duckdb.connect(str(CBB_DB_PATH), read_only=True)
        try:
            result = con.execute(
                "SELECT game_id, home_margin_h1, total_h1 "
                "FROM cbb_game_samples ORDER BY game_id, sim_idx"
            ).fetchall()
        finally:
            con.close()

        if not result:
            log.info("Sample cache prewarm: no samples found")
            return 0

        # Group rows by game_id
        from itertools import groupby
        from operator import itemgetter
        for game_id, rows in groupby(result, key=itemgetter(0)):
            data = [(r[1], r[2]) for r in rows]
            _sample_cache[game_id] = np.array(data, dtype=np.float64)

        log.info("Sample cache prewarmed: %d games", len(_sample_cache))
        return len(_sample_cache)

    except Exception as e:
        log.warning("Sample cache prewarm failed: %s", e)
        return 0


def clear_sample_cache():
    """Clear and reload sample cache. Alias for prewarm_sample_cache()."""
    return prewarm_sample_cache()


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
        elif side == "not_home":
            # NOT(home wins) = away wins or tie
            out[:] = np.where(margin > 0, 0.0, 1.0)
        elif side == "not_away":
            # NOT(away wins) = home wins or tie
            out[:] = np.where(margin < 0, 0.0, 1.0)
        elif side == "not_tie":
            out[:] = np.where(margin != 0, 1.0, 0.0)
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
    contracts = max(0, round(contracts))
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
            # EV vector: expected return per dollar wagered
            mu = np.array([
                (p["fair_prob"] * 100 - p["kalshi_price"]) / p["kalshi_price"]
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
                    contracts = max(0, round(contracts))
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
                adjusted = max(0, round(base * scale))
                sizes.append(int(adjusted))
            return sizes

    # Final fallback: independent single-bet Kelly
    return [
        kelly_size_single(p["fair_prob"], p["kalshi_price"], bankroll, kelly_mult)
        for p in new_positions
    ]


def batch_kelly_sizes_for_game(markets, placed_positions, bankroll, kelly_mult):
    """Compute Kelly sizes for all markets on a game with cross-side awareness.

    Bids are batched into one conditional_kelly_sizes call (all YES positions
    are correlated but not perfectly, so the matrix is well-conditioned).

    Asks are computed per-ticker. Each ask includes bid results from OTHER
    tickers as "placed" positions, so Kelly sees the cross-side correlation
    (e.g., YES SLU and NO MICH are both directional bets on Saint Louis).
    Same-ticker bid results are excluded to avoid the ρ=-1 singularity.

    Args:
        markets: list of quotable market dicts for this game
        placed_positions: existing positions on this game
        bankroll, kelly_mult: Kelly params

    Returns:
        dict: ticker -> {"bid_size": int, "ask_size": int}
    """
    result = {m["ticker"]: {"bid_size": 0, "ask_size": 0} for m in markets}
    if not markets:
        return result

    game_id = markets[0].get("game_id")
    if not game_id:
        return result

    samples = load_game_samples(game_id)
    if samples is None:
        return result

    # Build bid and ask positions
    bid_positions = []
    bid_tickers = []
    ask_positions = []
    ask_tickers = []

    for market in markets:
        ticker = market["ticker"]
        fair_cents = market["fair_prob"] * 100

        # Bid (buying YES)
        bid_pos = _market_to_position(market, "yes")
        book_bid = market.get("book_bid", 0)
        if book_bid > 0:
            bid_pos["kalshi_price"] = book_bid
        else:
            bid_pos["kalshi_price"] = int(math.floor(fair_cents / 1.05))

        if 0 < bid_pos["kalshi_price"] < 100:
            bid_positions.append(bid_pos)
            bid_tickers.append(ticker)

        # Ask (selling YES = buying NO)
        ask_pos = _market_to_position(market, "no")
        book_ask = market.get("book_ask", 0)
        if book_ask > 0:
            ask_pos["kalshi_price"] = 100 - book_ask
        else:
            ask_pos["kalshi_price"] = 100 - int(math.ceil((fair_cents + 5) / 1.05))

        if 0 < ask_pos["kalshi_price"] < 100:
            ask_positions.append(ask_pos)
            ask_tickers.append(ticker)

    # Diagnostic: log placed position count per game
    if placed_positions:
        total_size = sum(p["size"] for p in placed_positions)
        log.info("Game %s: %d placed positions (%d total contracts), %d bid + %d ask markets",
                 markets[0].get("game_key", "?"), len(placed_positions), total_size,
                 len(bid_positions), len(ask_positions))

    # Pass 1: Compute bid sizes (one batch call — all YES, well-conditioned)
    bid_results = {}  # ticker -> (position_dict, size)
    if bid_positions:
        bid_sizes = conditional_kelly_sizes(bid_positions, placed_positions,
                                            samples, bankroll, kelly_mult)
        for ticker, pos, size in zip(bid_tickers, bid_positions, bid_sizes):
            result[ticker]["bid_size"] = size
            if size > 0:
                bid_results[ticker] = (pos, size)

    # Pass 2: Compute ask sizes per-ticker with cross-ticker bid awareness.
    # For each ask, include bid results from OTHER tickers as placed positions
    # so Kelly sees the correlation. Exclude same-ticker bids (ρ=-1 singularity).
    if ask_positions:
        for ask_pos, ask_ticker in zip(ask_positions, ask_tickers):
            # Build placed = existing positions + bid results from other tickers
            cross_placed = list(placed_positions)
            for bid_ticker, (bid_pos, bid_size) in bid_results.items():
                if bid_ticker != ask_ticker:
                    cross_placed.append({
                        **bid_pos,
                        "size": bid_size,
                    })

            sizes = conditional_kelly_sizes([ask_pos], cross_placed,
                                            samples, bankroll, kelly_mult)
            result[ask_ticker]["ask_size"] = sizes[0] if sizes else 0

    return result


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
        return 0  # no game_id → can't compute Kelly → don't quote

    samples = load_game_samples(game_id)
    if samples is None:
        return 0  # no samples → can't compute Kelly → don't quote

    # Build the new position dict
    fair_cents = market["fair_prob"] * 100
    if side == "bid":
        # Buying YES: we profit when the event happens
        pos = _market_to_position(market, "yes")
        book_bid = market.get("book_bid", 0)
        if book_bid > 0:
            pos["kalshi_price"] = book_bid
        else:
            # Empty book — use EV-floor price (same formula as quoter)
            pos["kalshi_price"] = int(math.floor(fair_cents / 1.05))
    else:
        # Selling YES (buying NO): we profit when the event doesn't happen
        pos = _market_to_position(market, "no")
        book_ask = market.get("book_ask", 0)
        if book_ask > 0:
            pos["kalshi_price"] = 100 - book_ask
        else:
            # Empty book — use EV-floor price for NO side
            pos["kalshi_price"] = 100 - int(math.ceil((fair_cents + 5) / 1.05))

    if pos["kalshi_price"] <= 0 or pos["kalshi_price"] >= 100:
        return 0

    if not placed_positions:
        return kelly_size_single(pos["fair_prob"], pos["kalshi_price"],
                                 bankroll, kelly_mult)

    sizes = conditional_kelly_sizes([pos], placed_positions, samples,
                                    bankroll, kelly_mult)
    return sizes[0] if sizes else 0


def _market_to_position(market, yes_or_no):
    """Convert a quotable market dict to a Kelly position dict.

    For spreads, converts the Kalshi strike to home-spread convention
    (negative = home favorite) so evaluate_outcomes computes correctly:
      threshold = -line → home covers when margin > threshold.
    """
    mt = market["market_type"]
    strike = market.get("strike")

    # Convert strike to home-spread convention for evaluate_outcomes
    if mt == "spreads" and strike is not None:
        # Home contract strike=3 ("home wins by >3") → home_spread = -3
        # Away contract strike=5 ("away wins by >5") → home_spread = +5
        line = -strike if market.get("is_home_contract") else strike
    elif mt == "totals":
        line = strike  # total line is used as-is
    else:
        line = None

    if yes_or_no == "yes":
        fair_prob = market["fair_prob"]
        if mt == "spreads":
            side = "home" if market.get("is_home_contract") else "away"
        elif mt == "totals":
            side = "over"
        elif mt == "moneyline":
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
                side = "not_home"
            elif ct == market.get("away_team"):
                side = "not_away"
            else:
                side = "not_tie"
        else:
            side = "away"

    return {
        "market_type": mt,
        "side": side,
        "line": line,
        "fair_prob": fair_prob,
        "kalshi_price": 0,  # caller fills this in
    }


def kelly_size_for_take(market, take_side, execution_price, fee_cents,
                       placed_positions, bankroll, kelly_mult):
    """Compute Kelly size for a taker order at the actual execution price.

    Unlike kelly_size_for_quote (which uses the book bid/ask for resting
    orders), this uses the confirmed execution price and accounts for
    the taker fee — producing smaller sizes that reflect true edge.

    Args:
        market: quotable market dict
        take_side: "yes" or "no"
        execution_price: price in cents we'll pay (ask for YES, 100-bid for NO)
        fee_cents: taker fee in cents per contract
        placed_positions: existing positions for conditional adjustment
        bankroll: total bankroll
        kelly_mult: fractional Kelly

    Returns:
        int contract count
    """
    game_id = market.get("game_id")
    if not game_id:
        return 0

    samples = load_game_samples(game_id)
    if samples is None:
        return 0

    # Build position at fee-adjusted price
    if take_side == "yes":
        pos = _market_to_position(market, "yes")
    else:
        pos = _market_to_position(market, "no")

    # Effective price includes the taker fee
    effective_price = int(round(execution_price + fee_cents))
    effective_price = max(1, min(99, effective_price))
    pos["kalshi_price"] = effective_price

    if pos["fair_prob"] <= effective_price / 100.0:
        return 0  # No edge after fees

    if not placed_positions:
        return kelly_size_single(pos["fair_prob"], effective_price,
                                 bankroll, kelly_mult)

    sizes = conditional_kelly_sizes([pos], placed_positions, samples,
                                    bankroll, kelly_mult)
    return sizes[0] if sizes else 0


def get_placed_positions_for_game(game_key, game_id, current_markets=None):
    """Build placed position dicts for Kelly from Kalshi positions on this game.

    Args:
        game_key: (home_team, away_team) tuple.
        game_id: game ID for sample lookup.
        current_markets: list of quotable market dicts with current fair_prob.
            Used to override stale fair_prob with live model values.

    Returns list of position dicts with market_type, side, line, fair_prob,
    kalshi_price, and size fields.
    """
    _ensure_positions_loaded()
    positions = _positions_cache

    # Build ticker → current fair_prob lookup from live markets
    live_fair = {}
    if current_markets:
        for m in current_markets:
            live_fair[m["ticker"]] = m["fair_prob"]

    placed = []
    skipped_team = []
    skipped_line = []
    for pos in positions:
        pos_home = pos.get("home_team", "")
        pos_away = pos.get("away_team", "")
        if (pos_home, pos_away) != game_key:
            # Track skipped non-zero positions for diagnostics
            if pos.get("net_yes", 0) != 0:
                skipped_team.append((pos.get("ticker", "?"), pos_home, pos_away))
            continue
        net = pos.get("net_yes", 0)
        if net == 0:
            continue

        ticker = pos.get("ticker", "")
        mt = pos.get("market_type", "spreads")
        ct = pos.get("contract_team", "")

        if mt == "spreads":
            side = "home" if net > 0 else "away"
        elif mt == "totals":
            side = "over" if net > 0 else "under"
        elif mt == "moneyline":
            # Use contract_team to determine correct side
            if ct == pos_home:
                side = "home" if net > 0 else "not_home"
            elif ct == pos_away:
                side = "away" if net > 0 else "not_away"
            elif ct:  # tie contract
                side = "tie" if net > 0 else "not_tie"
            else:
                # No contract_team stored (legacy position) — skip
                continue
        else:
            continue

        # line_value is stored in home-spread convention (see _market_to_position)
        line = pos.get("line_value")

        # Spreads/totals require a line for evaluate_outcomes (threshold = -line).
        # Legacy positions without line_value would crash — skip them.
        if mt in ("spreads", "totals") and line is None:
            skipped_line.append((ticker, mt, net))
            continue

        # Use live fair_prob if available, fall back to DB value
        fair_prob = live_fair.get(ticker, pos.get("fair_prob", 0.5))
        # For NO positions (net < 0), the DB stores the YES-side ticker.
        # The live_fair lookup uses the same ticker, giving YES fair_prob.
        # We need the side's fair_prob: if we're short YES, our position
        # is effectively long NO, so fair_prob for the position = 1 - fair_yes.
        if net < 0:
            fair_prob = 1 - fair_prob

        placed.append({
            "market_type": mt,
            "side": side,
            "line": line,
            "fair_prob": fair_prob,
            "kalshi_price": int(pos.get("avg_entry_price", 50)),
            "size": abs(net),
        })

    # Diagnostic: warn if positions exist but were filtered out
    if skipped_line:
        log.warning("Skipped %d positions with NULL line_value for %s: %s",
                     len(skipped_line), game_key,
                     [(t, mt, n) for t, mt, n in skipped_line])
    # Only log team mismatches if we have few placed positions
    # (otherwise it's just other games, not a bug)
    if not placed and skipped_team:
        # Check if any skipped position's ticker matches this game's market tickers
        game_ticker_prefix = None
        if current_markets:
            game_ticker_prefix = current_markets[0].get("ticker", "")[:30]
        matching_skips = [s for s in skipped_team
                          if game_ticker_prefix and game_ticker_prefix[:20] in s[0]]
        if matching_skips:
            log.warning("NO placed positions for %s but found %d with matching tickers "
                        "and different team names: %s",
                        game_key, len(matching_skips), matching_skips[:3])

    return placed
