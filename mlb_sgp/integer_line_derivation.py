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
