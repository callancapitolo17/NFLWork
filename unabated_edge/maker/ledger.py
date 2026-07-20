"""Per-match risk ledger: exact worst-case P&L over the total-goals grid.

Every rung of a match settles on one integer (regulation-time total goals),
so the match portfolio's P&L is a function of that single number. Evaluating
it exactly on g=0..G_MAX replaces per-market caps and prices offsets, hedges,
and middles correctly with no special cases.

A fill is (line, side, contracts, price): side is relative to the Over
market — "yes" wins when g > line, "no" wins when g <= line. Prices in
dollars; P&L per contract is (1 - price) on a win, -price on a loss.
"""

G_MAX = 10   # "10 or more" bucket; P&L is constant above the top rung (tested)


def pnl(fills, g: int) -> float:
    total = 0.0
    for line, side, contracts, price in fills:
        win = (g > line) if side == "yes" else (g <= line)
        total += contracts * ((1.0 - price) if win else -price)
    return total


def pnl_grid(fills) -> list[float]:
    return [pnl(fills, g) for g in range(G_MAX + 1)]


def worst_case(fills) -> float:
    if not fills:
        return 0.0
    return min(pnl_grid(fills))


def max_contracts(fills, line: float, side: str, price: float, budget: float) -> int:
    """Largest integer n >= 0 such that adding (line, side, n, price) keeps
    worst_case >= -budget.

    Closed form: each goal outcome g imposes base_g + n*unit_g >= -budget and
    only outcomes where the candidate loses (unit_g < 0) bind. If the book is
    already beyond budget we never add exposure (fail closed)."""
    if budget <= 0:
        return 0
    base = pnl_grid(fills) if fills else [0.0] * (G_MAX + 1)
    if min(base) < -budget - 1e-9:
        return 0
    bound = None
    for g in range(G_MAX + 1):
        win = (g > line) if side == "yes" else (g <= line)
        unit = (1.0 - price) if win else -price
        if unit < 0:
            allowed = (budget + base[g]) / (-unit)
            bound = allowed if bound is None else min(bound, allowed)
    return int(bound + 1e-9) if bound is not None else 0


def mark_to_fair(fills, fair_by_line) -> float:
    """Unrealized P&L of open fills marked at the current devigged anchor fair.
    fair_by_line: {line: p_over}. A fill's current value ≈ p_over (yes) or
    1−p_over (no); unrealized = contracts × (value − price). Lines with no
    current fair are left unmarked (contribute 0)."""
    total = 0.0
    for line, side, contracts, price in fills:
        p_over = fair_by_line.get(line)
        if p_over is None:
            continue
        val = p_over if side == "yes" else (1.0 - p_over)
        total += contracts * (val - price)
    return total
