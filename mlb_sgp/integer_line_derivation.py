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
