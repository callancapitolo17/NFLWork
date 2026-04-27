"""Pure-Python helpers for combined parlay pricing in the MLB Dashboard.

Mirrors compute_combined_parlay_pricing() in Tools.R but in Python so the
Flask server doesn't need to call out to R for the live banner pricing.
"""
from __future__ import annotations


def joint_pricing(fair_dec_a: float, fair_dec_b: float, wz_dec: float,
                  bankroll: float | None = None,
                  kelly_mult: float | None = None) -> dict:
    """Joint pricing for two cross-game parlays combined into one ticket.

    Independence across games means joint fair = product of decimals.
    Returns joint fair odds, fair prob, edge, and (if bankroll/kelly_mult
    provided) the recommended Kelly stake in dollars.
    """
    joint_fair_dec = fair_dec_a * fair_dec_b
    joint_fair_prob = 1.0 / joint_fair_dec
    joint_edge = joint_fair_prob * wz_dec - 1.0

    kelly_stake = None
    if bankroll is not None and kelly_mult is not None:
        edge_fraction = joint_edge / (wz_dec - 1.0)
        kelly_stake = round(max(edge_fraction * kelly_mult * bankroll, 0.0), 2)

    return {
        "joint_fair_dec": round(joint_fair_dec, 4),
        "joint_fair_prob": round(joint_fair_prob, 6),
        "joint_edge": round(joint_edge, 4),
        "kelly_stake": kelly_stake,
    }
