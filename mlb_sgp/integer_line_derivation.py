"""
Derive joint fair probabilities for SGP combos at integer total lines
from two adjacent half-point alt SGP prices.

See docs/superpowers/specs/2026-05-08-sgp-integer-line-derivation-design.md
"""

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
PUSH_MASS_REL_TOL = 0.10

# Total marginal push mass plausibility bounds (fraction of games landing on X)
DELTA_TOTAL_MIN = 0.03
DELTA_TOTAL_MAX = 0.18

# Sum of 4 derived fair_probs should ≈ 1.0 (mathematical invariant)
SUM_MIN = 0.97
SUM_MAX = 1.03


def validate_per_alt_vig(vig_sum: float) -> bool:
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
    return DELTA_TOTAL_MIN <= delta_total <= DELTA_TOTAL_MAX


def validate_sum_to_one(fair_probs: list[float]) -> bool:
    s = sum(fair_probs)
    return SUM_MIN <= s <= SUM_MAX


def validate_per_combo_bounds(fair_prob: float) -> bool:
    return 0.0 < fair_prob < 1.0
