"""Whole-number 1H lines: the push chances, read off one book's SGP prices at
the half-point lines either side (Cal, 2026-11-12: books only), and what
Wagerzon pays on a push.

Pure functions, no network, no DB.

A whole-number away spread `a` splits the away team's 1H result three ways:
away covers (A), push (P), home covers (H). A whole-number total `t` splits
over (O) / push (P) / under (U). One book's 2 x 2 partition at a half-point
line either side covers a band of those outcomes, e.g. at away a - 0.5 the
away side wins only on A and the home side wins on P or H. The bracketing
partitions (two for one whole-number leg, four for two) over-determine the
grid; least squares fits it, then the grid is refused when a cell is
negative past rounding, a leg's push chance is implausible, or the
partitions disagree.

Wagerzon grades a pushed leg by dropping it: the parlay pays the other
leg alone at its straight price; both legs pushing refunds the stake.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np

from nfl_specials.pricing import american_to_decimal
from nfl_specials.wz_halves import CORNERS, WzHalfGame

AWAY, PUSH, HOME = "away", "push", "home"
OVER, UNDER = "over", "under"
HALF_POINT = 0.5

# A fitted cell this far below zero is book rounding, clamped to 0; further
# below, the book's nearby prices contradict each other and the game is skipped.
NEGATIVE_CELL_TOLERANCE = 0.005
# NFL 1H, 2019-2025 (nflverse, 1,960 games): the likeliest margin (3) is 12.3%,
# the likeliest total (20) 8.6%. Past these, the read is a data error.
MAX_SPREAD_PUSH_PROB = 0.20
MAX_TOTAL_PUSH_PROB = 0.15
# Largest gap between a book's devigged cell and the fitted grid's band.
MAX_FIT_RESIDUAL = 0.03


def is_whole_number(line: float) -> bool:
    return line == int(line)


def bracket_lines(line: float) -> tuple[float, ...]:
    """The half-point lines a partition is priced at: the line itself, or the two either side of a whole number."""
    if is_whole_number(line):
        return (line - HALF_POINT, line + HALF_POINT)
    return (line,)


def spread_outcomes(away_spread: float) -> tuple[str, ...]:
    return (AWAY, PUSH, HOME) if is_whole_number(away_spread) else (AWAY, HOME)


def total_outcomes(total: float) -> tuple[str, ...]:
    return (OVER, PUSH, UNDER) if is_whole_number(total) else (OVER, UNDER)


def spread_band(side: str, priced_line: float, wz_line: float) -> frozenset[str]:
    """Grid outcomes on which `side` wins at the book's `priced_line`, for Wagerzon's `wz_line` (both away lines)."""
    if priced_line < wz_line:       # the away side needs one more point than at Wagerzon
        return frozenset({AWAY}) if side == AWAY else frozenset({PUSH, HOME})
    if priced_line > wz_line:
        return frozenset({AWAY, PUSH}) if side == AWAY else frozenset({HOME})
    return frozenset({side})


def total_band(side: str, priced_line: float, wz_line: float) -> frozenset[str]:
    if priced_line < wz_line:       # landing on Wagerzon's total is over the book's line
        return frozenset({OVER, PUSH}) if side == OVER else frozenset({UNDER})
    if priced_line > wz_line:
        return frozenset({OVER}) if side == OVER else frozenset({PUSH, UNDER})
    return frozenset({side})


def double_push_possible(away_spread: float, total: float) -> bool:
    """A 1H margin and total are both even or both odd, so the spread and the
    total can only land exactly together when the two lines share parity."""
    return (int(away_spread) - int(total)) % 2 == 0


@dataclass(frozen=True)
class PushGrid:
    cells: dict[tuple[str, str], float]     # (spread outcome, total outcome) -> probability

    def prob(self, spread: str, total: str) -> float:
        return self.cells.get((spread, total), 0.0)


@dataclass(frozen=True)
class GridFit:
    grid: PushGrid | None
    reason: str | None = None


def fit_push_grid(partitions: dict[tuple[float, float], dict[str, float]],
                  away_spread: float, total: float) -> GridFit:
    """The cover/push x over/push/under grid from one book's devigged corner
    fairs, keyed by the (away line, total) each partition was priced at
    (every bracket_lines pair)."""
    spreads, totals = spread_outcomes(away_spread), total_outcomes(total)
    unknowns = [(s, t) for s in spreads for t in totals
                if not (s == PUSH and t == PUSH and not double_push_possible(away_spread, total))]
    rows, targets = [], []
    for (spread_line, total_line), corner_fairs in partitions.items():
        for corner in CORNERS:
            spread_side, total_side = corner.split("_")
            band = (spread_band(spread_side, spread_line, away_spread), total_band(total_side, total_line, total))
            rows.append([1.0 if s in band[0] and t in band[1] else 0.0 for s, t in unknowns])
            targets.append(corner_fairs[corner])
    rows.append([1.0] * len(unknowns))
    targets.append(1.0)
    fitted, *_ = np.linalg.lstsq(np.array(rows), np.array(targets), rcond=None)
    residual = float(np.max(np.abs(np.array(rows) @ fitted - np.array(targets))))
    if residual > MAX_FIT_RESIDUAL:
        return GridFit(None, f"nearby lines disagree by {residual:.3f}")
    cells = dict(zip(unknowns, (float(p) for p in fitted)))
    worst_cell, worst_prob = min(cells.items(), key=lambda item: item[1])
    if worst_prob < -NEGATIVE_CELL_TOLERANCE:
        return GridFit(None, f"implied {'/'.join(worst_cell)} chance {worst_prob:.3f} is negative")
    cells = {cell: max(prob, 0.0) for cell, prob in cells.items()}
    scale = sum(cells.values())
    cells = {cell: prob / scale for cell, prob in cells.items()}
    grid = PushGrid(cells)
    spread_push = sum(grid.prob(PUSH, t) for t in totals)
    total_push = sum(grid.prob(s, PUSH) for s in spreads)
    if spread_push > MAX_SPREAD_PUSH_PROB:
        return GridFit(None, f"implied spread push {spread_push:.3f} above {MAX_SPREAD_PUSH_PROB}")
    if total_push > MAX_TOTAL_PUSH_PROB:
        return GridFit(None, f"implied total push {total_push:.3f} above {MAX_TOTAL_PUSH_PROB}")
    return GridFit(grid)


def corner_outcomes(grid: PushGrid, corner: str) -> tuple[float, float, float, float]:
    """(both win, spread wins + total pushes, spread pushes + total wins, both push) for one corner."""
    spread_side, total_side = corner.split("_")
    return (grid.prob(spread_side, total_side), grid.prob(spread_side, PUSH),
            grid.prob(PUSH, total_side), grid.prob(PUSH, PUSH))


def push_return(grid: PushGrid, corner: str, game: WzHalfGame) -> float:
    """Expected dollars back per $1 from the corner's push outcomes at
    Wagerzon: the surviving leg at its straight decimal, or the stake."""
    spread_side, total_side = corner.split("_")
    spread_american = game.away_spread_american if spread_side == AWAY else game.home_spread_american
    total_american = game.over_american if total_side == OVER else game.under_american
    _, spread_only, total_only, both_push = corner_outcomes(grid, corner)
    return (spread_only * american_to_decimal(spread_american)
            + total_only * american_to_decimal(total_american)
            + both_push)


def corner_bracket_pairs(game: WzHalfGame) -> list[tuple[float, float]]:
    """Every (away line, total) partition a book prices for `game`."""
    return list(itertools.product(bracket_lines(game.away_spread), bracket_lines(game.total)))
