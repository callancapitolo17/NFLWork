"""Backward-compat shim. Logic lives in kalshi_common.ev_calc."""
from kalshi_common.ev_calc import *  # noqa: F401,F403
from kalshi_common.ev_calc import (
    fee_per_contract, post_fee_ev_buy_yes, post_fee_ev_buy_no, maker_fee_per_contract,
)
