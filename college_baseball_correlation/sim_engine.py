#!/usr/bin/env python3
"""
College Baseball Inning-Level Simulation Engine

Prices ML+Total parlays by simulating games inning-by-inning with walk-off rules.
The walk-off truncation mechanism naturally produces ML+Total correlation:
  - Away wins → full 9+ innings → higher totals → away_over positively correlated
  - Home wins → game can end early (walk-off) → lower totals → home_under positively correlated

Models:
  - Fit empirical half-inning run distributions from ESPN linescores
  - Derive per-team scoring rates from ML odds + total line
  - Simulate games inning-by-inning with proper baseball rules
  - Price any ML+Total parlay by counting simulated outcomes

Usage:
    python sim_engine.py --fit              # Fit inning-level distributions from linescores
    python sim_engine.py --validate         # Backtest on 486 games with odds
    python sim_engine.py --price            # Price today's games (production)
"""

import argparse
import json
import math
import os
import sys
from pathlib import Path

import duckdb
import numpy as np
from scipy import stats, optimize

DB_PATH = Path(__file__).parent / "college_baseball.duckdb"
PARAMS_PATH = Path(__file__).parent / "model_params.json"
N_SIMS = 200_000  # simulations per game


# =============================================================================
# UTILITIES
# =============================================================================

def devig_ml(home_ml: float, away_ml: float) -> tuple[float, float]:
    """Devig American moneyline odds to true probabilities."""
    def implied(odds):
        if odds > 0:
            return 100 / (odds + 100)
        else:
            return abs(odds) / (abs(odds) + 100)
    p_h = implied(home_ml)
    p_a = implied(away_ml)
    total = p_h + p_a
    return p_h / total, p_a / total


def prob_to_american(p: float) -> int:
    """Convert probability to American odds."""
    if p <= 0 or p >= 1:
        return 0
    if p >= 0.5:
        return round(-(p / (1 - p)) * 100)
    else:
        return round(((1 - p) / p) * 100)


def american_to_decimal(odds: int) -> float:
    """Convert American odds to decimal odds."""
    if odds > 0:
        return 1 + odds / 100
    else:
        return 1 + 100 / abs(odds)


# =============================================================================
# INNING-LEVEL DISTRIBUTION FITTING
# =============================================================================

def fit_distributions():
    """Fit half-inning run distributions from linescore data.

    Learns:
      - Empirical PMF of runs per half-inning (0, 1, 2, 3, ...)
      - NegBin parameters (r, p) for half-innings
      - Separate distributions by inning position (early/mid/late)
      - Home vs away scoring differences
      - Walk-off truncation frequency
    """
    con = duckdb.connect(str(DB_PATH), read_only=True)

    # Check linescores exist
    n_ls = con.execute("SELECT COUNT(DISTINCT game_id) FROM linescores").fetchone()[0]
    if n_ls == 0:
        print("No linescores in database. Run: python validate.py --pull-linescores 2025 2026")
        con.close()
        return

    # Get all half-inning data for regulation games (no extras, no mercy)
    rows = con.execute("""
        SELECT l.game_id, l.inning, l.half, l.runs,
               g.home_score, g.away_score, g.went_extras, g.mercy_rule
        FROM linescores l
        INNER JOIN games g ON l.game_id = g.game_id
        WHERE NOT g.went_extras AND NOT g.mercy_rule
        ORDER BY l.game_id, l.inning, CASE WHEN l.half = 'top' THEN 0 ELSE 1 END
    """).fetchall()
    con.close()

    if not rows:
        print("No matching linescore data for regulation games.")
        return

    # Organize by game
    games = {}
    for game_id, inning, half, runs, h_score, a_score, extras, mercy in rows:
        if game_id not in games:
            games[game_id] = {"innings": [], "home_score": h_score, "away_score": a_score}
        games[game_id]["innings"].append((inning, half, runs))

    n_games = len(games)
    print(f"=== INNING-LEVEL DISTRIBUTION FITTING ({n_games} games, {len(rows)} half-innings) ===\n")

    # Collect all half-inning runs
    all_runs = np.array([r[3] for r in rows if r[3] >= 0])  # filter out ESPN -2 glitch
    top_runs = np.array([r[3] for r in rows if r[2] == "top" and r[3] >= 0])
    bot_runs = np.array([r[3] for r in rows if r[2] == "bottom" and r[3] >= 0])

    print(f"Half-inning observations: {len(all_runs)}")
    print(f"  Top (away batting): {len(top_runs)}  mean={np.mean(top_runs):.3f}")
    print(f"  Bot (home batting): {len(bot_runs)}  mean={np.mean(bot_runs):.3f}")

    # PMF of runs per half-inning
    max_runs = 10
    pmf_counts = np.bincount(np.clip(all_runs.astype(int), 0, max_runs), minlength=max_runs + 1)
    pmf = pmf_counts / pmf_counts.sum()
    print(f"\nEmpirical half-inning PMF:")
    for k in range(max_runs + 1):
        bar = "#" * int(pmf[k] * 100)
        label = f"{k}+" if k == max_runs else str(k)
        print(f"  {label:>3}: {pmf[k]:.4f} ({pmf_counts[k]:>6}) {bar}")

    # Fit NegBin to half-inning runs
    mean_runs = np.mean(all_runs)
    var_runs = np.var(all_runs, ddof=1)
    vr = var_runs / mean_runs
    print(f"\nMean runs/half-inning: {mean_runs:.4f}")
    print(f"Variance:              {var_runs:.4f}")
    print(f"Var/Mean ratio:        {vr:.3f} (Poisson=1.0)")

    # NegBin: r = mu^2 / (var - mu)
    if var_runs > mean_runs:
        r_half = mean_runs ** 2 / (var_runs - mean_runs)
        p_half = r_half / (r_half + mean_runs)
    else:
        r_half = 999.0  # Poisson-like
        p_half = r_half / (r_half + mean_runs)

    # Log-likelihoods
    poisson_ll = np.sum(stats.poisson.logpmf(all_runs.astype(int), mean_runs))
    nb_ll = np.sum(stats.nbinom.logpmf(all_runs.astype(int), r_half, p_half))
    print(f"\nNegBin: r={r_half:.3f}, p={p_half:.4f}")
    print(f"Poisson log-lik:  {poisson_ll:.1f}")
    print(f"NegBin log-lik:   {nb_ll:.1f}")
    print(f"NegBin improvement: {nb_ll - poisson_ll:.1f} log-units")

    # Scoring by inning position
    print(f"\n--- Scoring by Inning ---")
    print(f"{'Inning':>8} {'N_top':>8} {'Top Mean':>10} {'N_bot':>8} {'Bot Mean':>10}")
    print("-" * 50)
    inning_params = {}
    for inn in range(1, 10):
        top = np.array([r[3] for r in rows if r[1] == inn and r[2] == "top" and r[3] >= 0])
        bot = np.array([r[3] for r in rows if r[1] == inn and r[2] == "bottom" and r[3] >= 0])
        if len(top) > 0 and len(bot) > 0:
            inning_params[inn] = {
                "top_mean": float(np.mean(top)),
                "bot_mean": float(np.mean(bot)),
                "top_var": float(np.var(top, ddof=1)),
                "bot_var": float(np.var(bot, ddof=1)),
                "n_top": len(top),
                "n_bot": len(bot),
            }
            print(f"  {inn:>5} {len(top):>8} {np.mean(top):>10.3f} {len(bot):>8} {np.mean(bot):>10.3f}")

    # Walk-off detection: games where home team has fewer half-innings
    walkoff_count = 0
    nine_inning_count = 0
    for gid, gdata in games.items():
        innings = gdata["innings"]
        top_innings = [i for i in innings if i[1] == "top"]
        bot_innings = [i for i in innings if i[1] == "bottom"]
        if len(top_innings) == 9 and len(bot_innings) == 8:
            # Home didn't bat in 9th (was ahead)
            walkoff_count += 1
        if len(top_innings) == 9:
            nine_inning_count += 1

    truncation_rate = walkoff_count / nine_inning_count if nine_inning_count > 0 else 0
    print(f"\n--- Walk-off Truncation ---")
    print(f"9-inning games: {nine_inning_count}")
    print(f"Home skipped bottom 9th: {walkoff_count} ({truncation_rate*100:.1f}%)")
    print(f"(This is what drives ML+Total correlation)")

    # Game-level validation: sum of half-innings vs reported scores
    score_mismatches = 0
    for gid, gdata in games.items():
        innings = gdata["innings"]
        top_total = sum(r for _, h, r in innings if h == "top" and r >= 0)
        bot_total = sum(r for _, h, r in innings if h == "bottom" and r >= 0)
        if top_total != gdata["away_score"] or bot_total != gdata["home_score"]:
            score_mismatches += 1
    print(f"Score mismatches (linescore sum ≠ final): {score_mismatches}/{n_games}")

    # --- Game-level pace factor ---
    # Each game has a latent "pace" that scales scoring for both teams.
    # Estimate pace per game as total_runs / expected_total (based on average rates).
    # The distribution of pace factors tells us how much game-level variance exists
    # beyond the half-inning level.
    print(f"\n--- Game-Level Pace Factor ---")
    mean_top = np.mean(top_runs)  # avg away batting per half-inning
    mean_bot = np.mean(bot_runs)  # avg home batting per half-inning

    pace_factors = []
    for gid, gdata in games.items():
        innings = gdata["innings"]
        n_top = sum(1 for _, h, r in innings if h == "top" and r >= 0)
        n_bot = sum(1 for _, h, r in innings if h == "bottom" and r >= 0)
        expected = n_top * mean_top + n_bot * mean_bot
        actual = gdata["home_score"] + gdata["away_score"]
        if expected > 0:
            pace_factors.append(actual / expected)

    pace_arr = np.array(pace_factors)
    pace_mean = np.mean(pace_arr)
    pace_std = np.std(pace_arr, ddof=1)
    pace_cv = pace_std / pace_mean  # coefficient of variation

    print(f"Pace factor: mean={pace_mean:.3f}, std={pace_std:.3f}, CV={pace_cv:.3f}")
    print(f"  P10={np.percentile(pace_arr, 10):.2f}  P50={np.percentile(pace_arr, 50):.2f}  P90={np.percentile(pace_arr, 90):.2f}")

    # Fit Gamma distribution to pace factors (must be > 0)
    # Gamma shape a, scale s: mean = a*s, var = a*s^2
    # So CV^2 = 1/a → a = 1/CV^2
    gamma_shape = 1 / (pace_cv ** 2)
    gamma_scale = pace_mean / gamma_shape
    print(f"Gamma fit: shape={gamma_shape:.3f}, scale={gamma_scale:.4f}")

    # Verify: simulate pace factors and check match
    rng = np.random.default_rng(42)
    sim_pace = rng.gamma(gamma_shape, gamma_scale, size=10000)
    print(f"  Simulated: mean={np.mean(sim_pace):.3f}, std={np.std(sim_pace):.3f}")

    # Save parameters
    params = {
        "n_games": n_games,
        "n_half_innings": len(all_runs),
        "half_inning": {
            "mean": float(mean_runs),
            "var": float(var_runs),
            "var_mean_ratio": float(vr),
            "negbin_r": float(r_half),
            "negbin_p": float(p_half),
        },
        "home_batting_mean": float(np.mean(bot_runs)),
        "away_batting_mean": float(np.mean(top_runs)),
        "home_advantage_per_inning": float(np.mean(bot_runs) - np.mean(top_runs)),
        "truncation_rate": float(truncation_rate),
        "pace_factor": {
            "gamma_shape": float(gamma_shape),
            "gamma_scale": float(gamma_scale),
            "mean": float(pace_mean),
            "std": float(pace_std),
            "cv": float(pace_cv),
        },
        "inning_params": {str(k): v for k, v in inning_params.items()},
        "empirical_pmf": {str(k): float(pmf[k]) for k in range(max_runs + 1)},
    }

    with open(PARAMS_PATH, "w") as f:
        json.dump(params, f, indent=2)
    print(f"\nParameters saved to {PARAMS_PATH}")
    return params


# =============================================================================
# INNING-LEVEL SIMULATION
# =============================================================================

def load_params() -> dict:
    """Load fitted model parameters."""
    if not PARAMS_PATH.exists():
        raise FileNotFoundError(f"Run --fit first. No params at {PARAMS_PATH}")
    with open(PARAMS_PATH) as f:
        return json.load(f)


def _nb_game_win_prob(mu_home_inn: float, mu_away_inn: float, r: float,
                      n_innings: int = 9, max_k: int = 25) -> float:
    """Analytically compute P(home wins) from per-inning NegBin rates.

    Game total for each team = sum of 9 NegBin(r, p) half-innings.
    Sum of 9 iid NegBin(r, p) = NegBin(9*r, p).
    This ignores walk-off truncation but gives a fast approximation.
    """
    r9 = n_innings * r
    p_h = r / (r + mu_home_inn)
    p_a = r / (r + mu_away_inn)

    # PMFs for 9-inning totals
    ks = np.arange(max_k + 1)
    pmf_h = stats.nbinom.pmf(ks, r9, p_h)
    pmf_a = stats.nbinom.pmf(ks, r9, p_a)

    # P(home > away) = sum h>a of P(H=h)*P(A=a)
    prob = 0.0
    for h in range(1, max_k + 1):
        prob += pmf_h[h] * np.sum(pmf_a[:h])
    return prob


def solve_team_rates(p_home: float, total_line: float, params: dict) -> tuple[float, float]:
    """Solve for per-inning scoring rates (mu_home, mu_away) from odds.

    Uses analytical NegBin PMF convolution for speed. The pace factor
    (mean ≈ 1.0) doesn't affect the marginal rates, only the joint distribution,
    so we solve rates analytically and apply pace in the simulation step.

    Constraints:
        9 * mu_away + 8.5 * mu_home ≈ total_line  (accounting for truncation)
        P(home wins | rates, r) ≈ p_home
    """
    r = params["half_inning"]["negbin_r"]
    trunc_rate = params.get("truncation_rate", 0.6)

    # Expected home innings ≈ 9 - truncation_rate (games where bot 9th skipped)
    # But truncation_rate depends on home win prob... use rough approximation
    home_eff_innings = 9 - 0.5  # ~8.5 on average

    def error(log_ratio):
        ratio = np.exp(log_ratio)
        # total_line ≈ 9 * mu_away + home_eff_innings * mu_home
        mu_away = total_line / (home_eff_innings * ratio + 9)
        mu_home = mu_away * ratio
        if mu_home < 0.01 or mu_away < 0.01:
            return 10.0
        wp = _nb_game_win_prob(mu_home, mu_away, r)
        return (wp - p_home) ** 2

    result = optimize.minimize_scalar(
        error,
        bounds=(-2.0, 2.0),
        method="bounded",
        options={"xatol": 0.005, "maxiter": 60}
    )

    ratio = np.exp(result.x)
    mu_away = total_line / (home_eff_innings * ratio + 9)
    mu_home = mu_away * ratio
    return mu_home, mu_away


def simulate_games(mu_home: float, mu_away: float, params: dict,
                   total_line: float, n_sims: int = N_SIMS,
                   seed: int = None) -> dict:
    """Simulate n_sims baseball games inning-by-inning with walk-off rules.

    Key features:
    1. NegBin half-inning scoring with game-specific rates
    2. Walk-off truncation: bottom 9th skipped if home leads after top 9
    3. Walk-off run capping: in bottom 9th, runs capped at deficit+1

    The simulation captures the walk-off truncation mechanism (~1.05 CF),
    which is then scaled by empirical calibration factors for production pricing.

    Returns dict with arrays: home_score, away_score, home_won, total, home_truncated
    """
    rng = np.random.default_rng(seed)
    r = params["half_inning"]["negbin_r"]

    p_h = r / (r + mu_home)
    p_a = r / (r + mu_away)

    # Simulate all 9 innings for both sides
    away_by_inning = rng.negative_binomial(r, p_a, size=(n_sims, 9))
    home_by_inning = rng.negative_binomial(r, p_h, size=(n_sims, 9))

    away_total = away_by_inning.sum(axis=1)
    home_thru_8 = home_by_inning[:, :8].sum(axis=1)

    # Walk-off truncation: bottom 9th skipped if home leads after top 9
    home_ahead_after_8 = home_thru_8 > away_total

    # Walk-off run capping for bottom 9th
    deficit = away_total - home_thru_8
    raw_bot_9 = home_by_inning[:, 8]
    walkoff_cap = np.maximum(deficit + 1, 0)
    capped_bot_9 = np.where(
        raw_bot_9 > walkoff_cap,
        walkoff_cap,
        raw_bot_9
    )

    home_total = np.where(
        home_ahead_after_8,
        home_thru_8,
        home_thru_8 + capped_bot_9
    )

    home_won = home_total > away_total
    away_won = away_total > home_total
    total = home_total + away_total

    return {
        "home_score": home_total,
        "away_score": away_total,
        "home_won": home_won,
        "away_won": away_won,
        "total": total,
        "home_truncated": home_ahead_after_8,
    }


def price_parlay(sims: dict, total_line: float) -> dict:
    """Price all 4 ML+Total parlay combos from simulation results."""
    n = len(sims["home_score"])
    home_won = sims["home_won"]
    away_won = sims["away_won"]
    over = sims["total"] > total_line
    under = sims["total"] < total_line
    push = sims["total"] == total_line

    # Marginal probs (excluding pushes for total)
    p_home = np.mean(home_won)
    p_away = np.mean(away_won)
    p_over = np.mean(over)
    p_under = np.mean(under)

    combos = {
        "home_over": (home_won & over, p_home * p_over),
        "home_under": (home_won & under, p_home * p_under),
        "away_over": (away_won & over, p_away * p_over),
        "away_under": (away_won & under, p_away * p_under),
    }

    results = {}
    for name, (joint_mask, indep_prob) in combos.items():
        joint_prob = np.mean(joint_mask)
        cf = joint_prob / indep_prob if indep_prob > 0 else 1.0
        results[name] = {
            "joint_prob": float(joint_prob),
            "independent_prob": float(indep_prob),
            "correlation_factor": float(cf),
            "fair_american": prob_to_american(joint_prob),
            "fair_decimal": float(1 / joint_prob) if joint_prob > 0 else 999,
        }

    results["marginals"] = {
        "p_home": float(p_home),
        "p_away": float(p_away),
        "p_over": float(p_over),
        "p_under": float(p_under),
        "truncation_rate": float(np.mean(sims["home_truncated"])),
    }

    return results


# =============================================================================
# EMPIRICAL CORRELATION FACTORS
# =============================================================================

# Empirical CFs from 486 matched games (ESPN scores + Odds API closing lines).
# These capture ALL correlation sources: walk-off truncation, pitcher quality,
# game dynamics, within-inning rally effects.
#
# The simulation model captures only the walk-off mechanism (~1.05 CF).
# The empirical CFs are the best available estimate for production pricing.
#
# Format: {total_bucket: {combo: CF}}
EMPIRICAL_CFS = {
    "low": {       # total_line < 12
        "away_over": 0.935, "away_under": 1.068,
        "home_over": 1.053, "home_under": 0.945,
    },
    "medium": {    # 12 <= total_line <= 13.5
        "away_over": 1.280, "away_under": 0.799,
        "home_over": 0.791, "home_under": 1.150,
    },
    "high": {      # total_line > 13.5
        "away_over": 1.097, "away_under": 0.956,
        "home_over": 0.948, "home_under": 1.023,
    },
}


def get_total_bucket(total_line: float) -> str:
    """Classify total line into bucket for CF lookup."""
    if total_line < 12:
        return "low"
    elif total_line <= 13.5:
        return "medium"
    else:
        return "high"


def price_parlay_with_cfs(p_home: float, p_over: float, total_line: float) -> dict:
    """Price ML+Total parlays using empirical correlation factors.

    This is the production pricing function. It:
    1. Takes devigged marginal probabilities from market odds
    2. Applies empirical CFs to compute fair joint probabilities
    3. Returns fair odds for all 4 combos

    Args:
        p_home: True probability of home win (devigged)
        p_over: True probability of over (devigged)
        total_line: Total line (used for CF bucket selection)

    Returns:
        Dict with joint_prob, fair_american, fair_decimal, cf for each combo
    """
    bucket = get_total_bucket(total_line)
    cfs = EMPIRICAL_CFS[bucket]

    p_away = 1 - p_home
    p_under = 1 - p_over

    combos = {
        "home_over": (p_home, p_over, cfs["home_over"]),
        "home_under": (p_home, p_under, cfs["home_under"]),
        "away_over": (p_away, p_over, cfs["away_over"]),
        "away_under": (p_away, p_under, cfs["away_under"]),
    }

    results = {}
    for name, (p_ml, p_total, cf) in combos.items():
        joint_prob = p_ml * p_total * cf
        # Clamp to valid probability range
        joint_prob = min(max(joint_prob, 0.001), 0.999)
        results[name] = {
            "joint_prob": joint_prob,
            "independent_prob": p_ml * p_total,
            "correlation_factor": cf,
            "fair_american": prob_to_american(joint_prob),
            "fair_decimal": 1 / joint_prob,
        }

    results["marginals"] = {
        "p_home": p_home,
        "p_away": p_away,
        "p_over": p_over,
        "p_under": p_under,
        "total_bucket": bucket,
    }

    return results


def compute_edge(fair_prob: float, book_american: int) -> dict:
    """Compute edge of a parlay bet.

    Args:
        fair_prob: True joint probability from our model
        book_american: American odds offered by the book

    Returns:
        Dict with edge_pct, fair_american, book_decimal, fair_decimal
    """
    book_decimal = american_to_decimal(book_american)
    fair_decimal = 1 / fair_prob
    edge_pct = (book_decimal / fair_decimal - 1) * 100

    return {
        "edge_pct": edge_pct,
        "fair_american": prob_to_american(fair_prob),
        "book_american": book_american,
        "book_decimal": book_decimal,
        "fair_decimal": fair_decimal,
        "is_plus_ev": edge_pct > 0,
    }


# =============================================================================
# VALIDATION
# =============================================================================

def validate():
    """Backtest inning-level model on games with odds.

    For each game:
      1. Derive per-inning rates from devigged ML + total line
      2. Simulate inning-by-inning with walk-off rules
      3. Price 4 ML+Total parlay combos
      4. Compare to actual outcomes

    Score with Brier score, calibration, and average CFs by total bucket.
    """
    params = load_params()

    con = duckdb.connect(str(DB_PATH), read_only=True)
    matched = con.execute("""
        WITH consensus AS (
            SELECT home_team, away_team,
                   CAST(commence_time AS DATE) as game_date,
                   MEDIAN(home_ml) as home_ml,
                   MEDIAN(away_ml) as away_ml,
                   MEDIAN(total_line) as total_line
            FROM odds
            WHERE home_ml IS NOT NULL AND away_ml IS NOT NULL AND total_line IS NOT NULL
            GROUP BY home_team, away_team, CAST(commence_time AS DATE)
        )
        SELECT g.home_score, g.away_score, g.total_runs,
               CASE WHEN g.home_won THEN 1 ELSE 0 END as home_won,
               c.home_ml, c.away_ml, c.total_line
        FROM games g
        INNER JOIN consensus c ON g.home_team = c.home_team AND g.away_team = c.away_team
            AND g.date = c.game_date
        WHERE NOT g.went_extras AND NOT g.mercy_rule
    """).fetchall()
    con.close()

    print(f"=== INNING-LEVEL MODEL VALIDATION ({len(matched)} games) ===\n")

    combo_names = ["home_over", "home_under", "away_over", "away_under"]
    preds = {c: [] for c in combo_names}
    actuals = {c: [] for c in combo_names}
    marginal_preds = {"home": [], "over": []}
    marginal_actuals = {"home": [], "over": []}
    cfs = {c: [] for c in combo_names}
    total_lines = []

    for i, (h_score, a_score, total_runs, h_won, h_ml, a_ml, total_line) in enumerate(matched):
        p_home, p_away = devig_ml(h_ml, a_ml)

        # Actual outcomes
        actual_over = 1 if total_runs > total_line else 0
        actual_combos = {
            "home_over": int(h_won == 1 and actual_over == 1),
            "home_under": int(h_won == 1 and actual_over == 0),
            "away_over": int(h_won == 0 and actual_over == 1),
            "away_under": int(h_won == 0 and actual_over == 0),
        }

        try:
            mu_h, mu_a = solve_team_rates(p_home, total_line, params)
            sims = simulate_games(mu_h, mu_a, params, total_line,
                                  n_sims=50_000, seed=i)
            prices = price_parlay(sims, total_line)
        except Exception as e:
            if i < 5:
                print(f"  Warning: game {i} failed: {e}")
            continue

        for combo in combo_names:
            preds[combo].append(prices[combo]["joint_prob"])
            actuals[combo].append(actual_combos[combo])
            cfs[combo].append(prices[combo]["correlation_factor"])

        marginal_preds["home"].append(prices["marginals"]["p_home"])
        marginal_actuals["home"].append(h_won)
        marginal_preds["over"].append(prices["marginals"]["p_over"])
        marginal_actuals["over"].append(actual_over)
        total_lines.append(total_line)

        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{len(matched)} games...")

    n_valid = len(preds["home_over"])
    print(f"\n  Successfully processed: {n_valid}/{len(matched)} games")

    # --- Brier Scores ---
    print(f"\n=== BRIER SCORES (lower = better) ===")
    print(f"{'Model':<20} {'home_over':>12} {'home_under':>12} {'away_over':>12} {'away_under':>12} {'MEAN':>10}")
    print("-" * 82)

    # Inning-level model
    brier_scores = []
    print(f"{'inning_sim':<20}", end="")
    for combo in combo_names:
        p = np.array(preds[combo])
        a = np.array(actuals[combo])
        brier = np.mean((p - a) ** 2)
        brier_scores.append(brier)
        print(f" {brier:>12.6f}", end="")
    print(f" {np.mean(brier_scores):>10.6f}")

    # Independence baseline
    indep_scores = []
    print(f"{'independence':<20}", end="")
    home_p = np.array(marginal_preds["home"])
    over_p = np.array(marginal_preds["over"])
    for combo in combo_names:
        a = np.array(actuals[combo])
        if "home" in combo and "over" in combo:
            indep = home_p * over_p
        elif "home" in combo and "under" in combo:
            indep = home_p * (1 - over_p)
        elif "away" in combo and "over" in combo:
            indep = (1 - home_p) * over_p
        else:
            indep = (1 - home_p) * (1 - over_p)
        brier = np.mean((indep - a) ** 2)
        indep_scores.append(brier)
        print(f" {brier:>12.6f}", end="")
    print(f" {np.mean(indep_scores):>10.6f}")

    improvement = np.mean(indep_scores) - np.mean(brier_scores)
    print(f"\nImprovement over independence: {improvement:+.6f} ({improvement/np.mean(indep_scores)*100:+.2f}%)")

    # --- Marginal calibration ---
    print(f"\n=== MARGINAL BRIER SCORES ===")
    for mkt in ["home", "over"]:
        p = np.array(marginal_preds[mkt])
        a = np.array(marginal_actuals[mkt])
        brier = np.mean((p - a) ** 2)
        print(f"  {mkt:>6}: {brier:.6f}")

    # --- Average CFs by total bucket ---
    print(f"\n=== AVERAGE CORRELATION FACTORS BY TOTAL LINE ===")
    total_arr = np.array(total_lines)

    buckets = {
        "low (<12)": total_arr < 12,
        "medium (12-13.5)": (total_arr >= 12) & (total_arr <= 13.5),
        "high (>13.5)": total_arr > 13.5,
    }

    print(f"{'Bucket':<20} {'N':>5} {'CF(AO)':>10} {'CF(AU)':>10} {'CF(HO)':>10} {'CF(HU)':>10}")
    print("-" * 70)
    for bname, mask in buckets.items():
        n_bucket = mask.sum()
        if n_bucket == 0:
            continue
        print(f"{bname:<20} {n_bucket:>5}", end="")
        for combo in combo_names:
            cf_arr = np.array(cfs[combo])
            print(f" {cf_arr[mask].mean():>10.3f}", end="")
        print()

    # --- Empirical CFs for comparison ---
    print(f"\n=== EMPIRICAL CORRELATION FACTORS (from actual outcomes) ===")
    for bname, mask in buckets.items():
        n_bucket = mask.sum()
        if n_bucket < 10:
            continue
        ho = np.array(actuals["home_over"])[mask]
        hu = np.array(actuals["home_under"])[mask]
        ao = np.array(actuals["away_over"])[mask]
        au = np.array(actuals["away_under"])[mask]

        p_home_b = np.mean(ho) + np.mean(hu)  # P(home)
        p_over_b = np.mean(ho) + np.mean(ao)  # P(over)
        p_away_b = 1 - p_home_b
        p_under_b = 1 - p_over_b

        print(f"{bname:<20} N={n_bucket:>4}", end="")
        for label, joint, m1, m2 in [
            ("AO", np.mean(ao), p_away_b, p_over_b),
            ("AU", np.mean(au), p_away_b, p_under_b),
            ("HO", np.mean(ho), p_home_b, p_over_b),
            ("HU", np.mean(hu), p_home_b, p_under_b),
        ]:
            cf_emp = joint / (m1 * m2) if m1 * m2 > 0 else 0
            print(f"  {label}={cf_emp:.3f}", end="")
        print()

    # --- Calibration plot data ---
    print(f"\n=== CALIBRATION (away_over) ===")
    print(f"{'Pred Bucket':>12} {'N':>6} {'Pred Avg':>10} {'Actual Rate':>12} {'Diff':>8}")
    print("-" * 52)
    p_ao = np.array(preds["away_over"])
    a_ao = np.array(actuals["away_over"])
    for lo in np.arange(0.05, 0.50, 0.05):
        hi = lo + 0.05
        m = (p_ao >= lo) & (p_ao < hi)
        if m.sum() < 5:
            continue
        pred_avg = p_ao[m].mean()
        actual_rate = a_ao[m].mean()
        diff = actual_rate - pred_avg
        print(f"  {lo:.2f}-{hi:.2f}   {m.sum():>6} {pred_avg:>10.4f} {actual_rate:>12.4f} {diff:>+8.4f}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="College Baseball Inning-Level Simulation Engine")
    parser.add_argument("--fit", action="store_true", help="Fit inning-level distributions from linescores")
    parser.add_argument("--validate", action="store_true", help="Backtest on games with odds")
    parser.add_argument("--price", action="store_true", help="Price today's games (production)")

    args = parser.parse_args()

    if not any([args.fit, args.validate, args.price]):
        parser.print_help()
        return

    if args.fit:
        fit_distributions()

    if args.validate:
        validate()

    if args.price:
        print("Production pricing not yet implemented. Run --fit and --validate first.")


if __name__ == "__main__":
    main()
