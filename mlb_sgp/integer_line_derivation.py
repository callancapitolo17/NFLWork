"""
Derive joint fair probabilities for SGP combos at integer total lines
from two adjacent half-point alt SGP prices.

See docs/superpowers/specs/2026-05-08-sgp-integer-line-derivation-design.md
"""

import logging
from typing import Optional

# Floating-point tolerance for matching integer lines (WZ data sometimes
# has 8.000001 / 7.999999 instead of exact 8.0)
_INT_TOLERANCE = 1e-3


def is_integer_line(value: float) -> bool:
    """True if value is at an integer (within FP tolerance), False if half-point."""
    return abs(value - round(value)) < _INT_TOLERANCE


def devig_alt_set(decimals: dict[str, float]) -> tuple[dict[str, float], float]:
    """
    Multiplicative devig of one alt-set's 4 combo decimal odds.

    Returns (devigged_probs, vig_sum). Each devigged_probs[k] is the joint
    probability of combo k under the multiplicative method:
        p_devigged = (1/decimal_k) / sum(1/decimal_i for all i)

    The 4 combos must form a partition of the joint sample space (Home Spread
    + Over, Home Spread + Under, Away Spread + Over, Away Spread + Under at
    one alt total). The returned devigged probs sum to exactly 1.0.

    The vig_sum (= sum of raw implied probs) is also returned for the
    per-alt vig bounds check.
    """
    implied = {k: 1.0 / d for k, d in decimals.items()}
    vig_sum = sum(implied.values())
    devigged = {k: p / vig_sum for k, p in implied.items()}
    return devigged, vig_sum


# ---------------------------------------------------------------------------
# Bounds-check thresholds (config constants — tunable in follow-up commits
# once production logs reveal real-world distributions).
# ---------------------------------------------------------------------------

# Per-alt vig sum bounds. DK ~1.125, FD bimodal 1.13/1.22, PX/NV unknown.
VIG_MIN = 1.05
VIG_MAX = 1.30

# Push-mass cross-consistency: |Δ_from_over - Δ_from_under| / max(...) tolerance
# 0.10 catches systematic vig asymmetry that signals scrape/parse errors at adjacent alts.
PUSH_MASS_REL_TOL = 0.10

# Total marginal push mass plausibility bounds (fraction of games landing on X)
DELTA_TOTAL_MIN = 0.03
DELTA_TOTAL_MAX = 0.18

# Sum of 4 derived fair_probs should ≈ 1.0 (mathematical invariant)
SUM_MIN = 0.97
SUM_MAX = 1.03


def validate_per_alt_vig(vig_sum: float) -> bool:
    """True if per-alt implied-prob sum is within plausible bounds [VIG_MIN, VIG_MAX]."""
    return VIG_MIN <= vig_sum <= VIG_MAX


def validate_push_mass_consistency(delta_a: float, delta_b: float) -> bool:
    """The push mass derived from the Over side must approximately equal the
    push mass derived from the Under side (they're the same joint event).
    Tolerance: PUSH_MASS_REL_TOL relative to the larger of the two.
    """
    max_delta = max(abs(delta_a), abs(delta_b))
    if max_delta == 0:
        return delta_a == delta_b == 0
    return abs(delta_a - delta_b) / max_delta < PUSH_MASS_REL_TOL


def validate_delta_total(delta_total: float) -> bool:
    """True if marginal joint push mass is within plausible bounds [DELTA_TOTAL_MIN, DELTA_TOTAL_MAX]."""
    return DELTA_TOTAL_MIN <= delta_total <= DELTA_TOTAL_MAX


def validate_sum_to_one(fair_probs: list[float]) -> bool:
    """True if sum of derived fair_probs falls within [SUM_MIN, SUM_MAX]. Typically 4 combos summing to ~1.0."""
    s = sum(fair_probs)
    return SUM_MIN <= s <= SUM_MAX


def validate_per_combo_bounds(fair_prob: float) -> bool:
    """True if fair_prob is strictly between 0 and 1. Strict bounds reject degenerate certainty (likely a data error)."""
    return 0.0 < fair_prob < 1.0


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

logger = logging.getLogger(__name__)

# Combo keys — must match exactly across all callers
COMBO_KEYS = ("home_over", "home_under", "away_over", "away_under")


def derive_fair_probs(
    decimals_lo: dict[str, float],
    decimals_hi: dict[str, float],
) -> Optional[dict]:
    """
    Derive 4 joint fair probabilities at an integer total line X from
    decimal odds at the two adjacent half-point alts (X-0.5 and X+0.5).

    Args:
        decimals_lo: 4 decimal odds at the lower alt (X-0.5). Keys: COMBO_KEYS.
        decimals_hi: 4 decimal odds at the higher alt (X+0.5). Keys: COMBO_KEYS.

    Returns:
        dict with keys {'fair_probs': {combo: prob, ...}, 'delta_total': float}
        or None if any bounds check fails.

    On rejection, logs a structured WARN line.
    """
    # Validate inputs have the expected keys
    if set(decimals_lo.keys()) != set(COMBO_KEYS) or set(decimals_hi.keys()) != set(COMBO_KEYS):
        logger.warning(
            "derive_fair_probs: invalid combo keys lo=%s hi=%s",
            sorted(decimals_lo.keys()), sorted(decimals_hi.keys())
        )
        return None

    # Per-alt devig
    devig_lo, vig_lo = devig_alt_set(decimals_lo)
    devig_hi, vig_hi = devig_alt_set(decimals_hi)

    # Bounds check (1): per-alt vig sum
    if not validate_per_alt_vig(vig_lo):
        logger.warning("derive_fair_probs: vig_lo=%.4f out of [%g, %g]", vig_lo, VIG_MIN, VIG_MAX)
        return None
    if not validate_per_alt_vig(vig_hi):
        logger.warning("derive_fair_probs: vig_hi=%.4f out of [%g, %g]", vig_hi, VIG_MIN, VIG_MAX)
        return None

    # Compute push masses two ways for each side
    delta_cover_from_over = devig_lo["home_over"] - devig_hi["home_over"]
    delta_cover_from_under = devig_hi["home_under"] - devig_lo["home_under"]
    delta_uncover_from_over = devig_lo["away_over"] - devig_hi["away_over"]
    delta_uncover_from_under = devig_hi["away_under"] - devig_lo["away_under"]

    # Bounds check (2): push-mass cross-consistency
    if not validate_push_mass_consistency(delta_cover_from_over, delta_cover_from_under):
        logger.warning(
            "derive_fair_probs: cover-side push mass inconsistent: from_over=%.4f from_under=%.4f",
            delta_cover_from_over, delta_cover_from_under
        )
        return None
    if not validate_push_mass_consistency(delta_uncover_from_over, delta_uncover_from_under):
        logger.warning(
            "derive_fair_probs: uncover-side push mass inconsistent: from_over=%.4f from_under=%.4f",
            delta_uncover_from_over, delta_uncover_from_under
        )
        return None

    # Average the consistent estimates per side (more robust than picking one)
    delta_cover = (delta_cover_from_over + delta_cover_from_under) / 2.0
    delta_uncover = (delta_uncover_from_over + delta_uncover_from_under) / 2.0
    delta_total = delta_cover + delta_uncover

    # Bounds check (3): Δ_total reasonableness
    if not validate_delta_total(delta_total):
        logger.warning(
            "derive_fair_probs: delta_total=%.4f out of [%g, %g]",
            delta_total, DELTA_TOTAL_MIN, DELTA_TOTAL_MAX
        )
        return None

    # Apply the formula:
    #   For Over combos: fair_prob = devig_hi["...over"] / (1 - delta_total)
    #   For Under combos: fair_prob = devig_lo["...under"] / (1 - delta_total)
    denom = 1.0 - delta_total
    fair_probs = {
        "home_over":  devig_hi["home_over"]  / denom,
        "home_under": devig_lo["home_under"] / denom,
        "away_over":  devig_hi["away_over"]  / denom,
        "away_under": devig_lo["away_under"] / denom,
    }

    # Bounds check (4): sum to 1
    if not validate_sum_to_one(list(fair_probs.values())):
        logger.warning(
            "derive_fair_probs: sum=%.4f out of [%g, %g] (probs=%s)",
            sum(fair_probs.values()), SUM_MIN, SUM_MAX, fair_probs
        )
        return None

    # Bounds check (5): per-combo (0, 1)
    for k, v in fair_probs.items():
        if not validate_per_combo_bounds(v):
            logger.warning("derive_fair_probs: %s=%.4f out of (0, 1)", k, v)
            return None

    return {"fair_probs": fair_probs, "delta_total": delta_total}
