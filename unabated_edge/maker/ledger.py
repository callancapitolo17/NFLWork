"""Per-game risk ledger: exact worst-case P&L over a totals ladder.

Sport-agnostic. Every rung settles on one number (the game's final total —
goals, runs, points); portfolio P&L is a piecewise-constant function of that
number whose only breakpoints are the fills' lines. Evaluating one
representative outcome per interval is therefore exact — no fixed grid, no
truncation above a top rung.

A fill is (line, side, contracts, price): side is relative to the Over
market — "yes" wins when t > line, "no" wins when t <= line. Prices in
dollars; P&L per contract is (1 - price) on a win, -price on a loss.
"""


def pnl(fills, t: float) -> float:
    total = 0.0
    for line, side, contracts, price in fills:
        win = (t > line) if side == "yes" else (t <= line)
        total += contracts * ((1.0 - price) if win else -price)
    return total


def outcome_points(lines, extra_line=None) -> list[float]:
    """One representative final-total per payoff interval of the given lines
    (plus an optional candidate line): below the lowest strike, between each
    consecutive pair, above the highest."""
    ls = sorted(set(lines) | ({extra_line} if extra_line is not None else set()))
    if not ls:
        return [0.0]
    return [ls[0] - 0.5] + [line + 0.5 for line in ls]


def pnl_grid(fills) -> list[float]:
    """P&L per payoff interval (variable length — one entry per interval)."""
    return [pnl(fills, t) for t in outcome_points([f[0] for f in fills])]


def worst_case(fills) -> float:
    if not fills:
        return 0.0
    return min(pnl_grid(fills))


def max_contracts(fills, line: float, side: str, price: float, budget: float) -> int:
    """Largest integer n >= 0 such that adding (line, side, n, price) keeps
    worst_case >= -budget.

    Closed form per interval: base_t + n*unit_t >= -budget binds only where
    the candidate loses (unit_t < 0). Intervals derive from existing fills'
    lines PLUS the candidate's line. Fail closed if already beyond budget."""
    if budget <= 0:
        return 0
    points = outcome_points([f[0] for f in fills], extra_line=line)
    base = [pnl(fills, t) for t in points]
    if min(base) < -budget - 1e-9:
        return 0
    bound = None
    for t, base_t in zip(points, base):
        win = (t > line) if side == "yes" else (t <= line)
        unit = (1.0 - price) if win else -price
        if unit < 0:
            allowed = (budget + base_t) / (-unit)
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
