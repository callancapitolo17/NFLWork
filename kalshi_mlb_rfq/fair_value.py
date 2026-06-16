"""Backward-compat shim. Logic lives in kalshi_common.fair_value."""
from kalshi_common.fair_value import *  # noqa: F401,F403
from kalshi_common.fair_value import (
    SpreadLeg, TotalLeg, Leg, model_fair, devig_book, blend, _hit_mask, _probit_devig_n,
)
