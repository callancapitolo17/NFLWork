import math


def kelly_contracts(fair_prob, yes_ask, bankroll, kelly_fraction, per_match_cap_dollars) -> int:
    if not (0 < yes_ask < 1):
        return 0
    b = (1 - yes_ask) / yes_ask
    f = (b * fair_prob - (1 - fair_prob)) / b
    if f <= 0:
        return 0
    stake = min(bankroll * kelly_fraction * f, per_match_cap_dollars)
    return max(0, math.floor(stake / yes_ask))
