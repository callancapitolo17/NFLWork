"""Phase 2 on-demand pure math: devig_partition (Route A) and
fair_by_correlation_transfer (Route B) + route selector."""
import math

import pandas as pd
import pytest

from kalshi_common import fair_value as fv


# ---------------------------------------------------------------- #
# Task 2: devig_partition (Route A)                                #
# ---------------------------------------------------------------- #

# A plausible 4-cell SGP grid (decimal odds, ~14% overround for 2 legs).
GRID4 = [3.10, 4.20, 3.60, 3.05]


def test_devig_partition_matches_devig_book_on_4_cells():
    """Regression bridge to Phase 1: same math as devig_book's 4-cell path."""
    rows = pd.DataFrame({
        "combo": ["Home Spread + Over", "Home Spread + Under",
                  "Away Spread + Over", "Away Spread + Under"],
        "sgp_decimal": GRID4,
    })
    expected = fv.devig_book(rows, combo="Home Spread + Over")
    got = fv.devig_partition(GRID4, n_legs=2)
    assert got is not None and expected is not None
    assert abs(got - expected) < 1e-12


def test_devig_partition_result_is_a_probability_and_target_first():
    fair = fv.devig_partition(GRID4, n_legs=2)
    assert 0.0 < fair < 1.0
    # devigging removes margin: fair < raw implied of the target cell
    assert fair < 1.0 / GRID4[0]


def test_devig_partition_overround_above_n_aware_bound_rejected():
    # 2 legs -> bound 1.50; craft cells summing to 1.60 implied
    cells = [1.0 / 0.50, 1.0 / 0.40, 1.0 / 0.40, 1.0 / 0.30]  # sum=1.60
    assert fv.devig_partition(cells, n_legs=2) is None
    # but the same overround is fine for 3 legs (bound 1.75) — pad to 8 cells
    cells8 = [1 / 0.20, 1 / 0.20, 1 / 0.20, 1 / 0.20,
              1 / 0.20, 1 / 0.20, 1 / 0.20, 1 / 0.20]        # sum=1.60
    assert fv.devig_partition(cells8, n_legs=3) is not None


def test_devig_partition_accepts_real_fd_margin():
    """LIVE calibration case (2026-07-10, FD ml+total @ MIN/LAA): a real
    healthy 4-cell SGP partition sums to 1.279 implied — the gate must not
    reject real book margin (the original 0.12/leg bound did, silently
    forcing Route B at high-margin books)."""
    cells = [2.813302325581395, 3.793224, 2.446418604651163, 3.983544]
    fair = fv.devig_partition(cells, n_legs=2)
    assert fair is not None
    assert 0.0 < fair < 1.0 / cells[0]      # devigged below raw implied


def test_devig_partition_underround_rejected():
    cells = [5.0, 5.0, 5.0, 5.0]     # sum of implied = 0.8 < 1.0
    assert fv.devig_partition(cells, n_legs=2) is None


def test_devig_partition_bad_cells_rejected():
    assert fv.devig_partition(None, n_legs=2) is None
    assert fv.devig_partition([], n_legs=2) is None
    assert fv.devig_partition([2.0, 2.0, 2.0], n_legs=2) is None          # wrong count
    assert fv.devig_partition([2.0, 2.0, 2.0, None], n_legs=2) is None    # missing cell
    assert fv.devig_partition([2.0, 2.0, 2.0, 0.9], n_legs=2) is None     # dec <= 1
    assert fv.devig_partition([2.0, 2.0, 2.0, float("nan")], n_legs=2) is None


# ---------------------------------------------------------------- #
# Task 3: fair_by_correlation_transfer (Route B) + route selector  #
# ---------------------------------------------------------------- #

def test_devig_two_way_exact():
    pair = fv.devig_two_way(1.91, 1.91)     # symmetric -110/-110
    assert pair is not None
    a, b = pair
    assert abs(a - 0.5) < 1e-9 and abs(b - 0.5) < 1e-9
    assert fv.devig_two_way(1.91, None) is None
    assert fv.devig_two_way(0.9, 1.91) is None


def test_transfer_recovers_fair_when_vig_cancels():
    # Construct: fair singles .49/.48/.57, independent-ish joint with a known
    # correlation multiplier 1.05; book applies margin m to singles and the
    # SAME compounded margin to the SGP -> transfer must recover fair exactly.
    fairs = [0.49, 0.48, 0.57]
    mult = 1.05
    fair_joint = fairs[0] * fairs[1] * fairs[2] * mult
    m = 1.06                                        # 6% margin per leg
    vigged = [f * m for f in fairs]
    sgp_implied = fair_joint * (m ** 3)             # compounded leg vig
    sgp_decimal = 1.0 / sgp_implied
    singles = list(zip(vigged, fairs))
    got = fv.fair_by_correlation_transfer(sgp_decimal, singles)
    assert got is not None
    assert abs(got - fair_joint) < 1e-9


def test_transfer_frechet_upper_rejects_impossible_joint():
    # joint fair cannot exceed min(single fairs); nested totals give this.
    singles = [(0.58, 0.55), (0.33, 0.30)]
    # SGP price implying joint fair ~0.35 > min(0.30) after transfer:
    # implied = 0.35 * (0.58/0.55) * (0.33/0.30) ≈ 0.406
    sgp_decimal = 1.0 / 0.406
    assert fv.fair_by_correlation_transfer(sgp_decimal, singles) is None


def test_transfer_frechet_lower_rejects_impossible_joint():
    # two near-certain legs: joint fair must be >= p1+p2-1 = 0.9
    singles = [(0.97, 0.95), (0.97, 0.95)]
    # craft an SGP decimal implying transfer fair ~0.5 (< 0.9) -> reject
    sgp_decimal = 1.0 / (0.5 * (0.97 / 0.95) ** 2)
    assert fv.fair_by_correlation_transfer(sgp_decimal, singles) is None


def test_transfer_accepts_high_correlation_multiplier():
    # RL+ML stack: multiplier ≈ 1/p(ml) ≈ 2.04 is legitimate (within Fréchet).
    singles = [(0.52, 0.49), (0.51, 0.48)]
    fair_joint = 0.47                       # just under min(fairs)=0.48
    sgp_decimal = 1.0 / (fair_joint * (0.52 / 0.49) * (0.51 / 0.48))
    got = fv.fair_by_correlation_transfer(sgp_decimal, singles)
    assert got is not None and abs(got - fair_joint) < 1e-9
    # implied multiplier here is fair_joint/(.49*.48) ≈ 2.0 — must NOT be capped
    assert fair_joint / (0.49 * 0.48) > 1.9


def test_transfer_bad_inputs():
    assert fv.fair_by_correlation_transfer(None, [(0.5, 0.48)]) is None
    assert fv.fair_by_correlation_transfer(0.9, [(0.5, 0.48)]) is None     # dec <= 1
    assert fv.fair_by_correlation_transfer(3.0, []) is None
    assert fv.fair_by_correlation_transfer(3.0, [(0.5, None)]) is None
    assert fv.fair_by_correlation_transfer(3.0, [(None, 0.48)]) is None


def test_book_on_demand_fair_prefers_partition():
    fair_a = fv.devig_partition(GRID4, n_legs=2)
    got = fv.book_on_demand_fair(GRID4, sgp_decimal=3.10,
                                 singles=[(0.52, 0.49), (0.51, 0.48)], n_legs=2)
    assert got == (fair_a, "partition")


def test_book_on_demand_fair_falls_back_to_transfer():
    singles = [(0.52, 0.49), (0.51, 0.48)]
    got = fv.book_on_demand_fair(None, sgp_decimal=1.0 / (0.30 * (0.52 / 0.49) * (0.51 / 0.48)),
                                 singles=singles, n_legs=2)
    assert got is not None
    fair, route = got
    assert route == "transfer" and abs(fair - 0.30) < 1e-9
    # incomplete cells also fall back
    got2 = fv.book_on_demand_fair([3.1, None, 3.6, 3.05],
                                  sgp_decimal=1.0 / (0.30 * (0.52 / 0.49) * (0.51 / 0.48)),
                                  singles=singles, n_legs=2)
    assert got2 is not None and got2[1] == "transfer"


def test_book_on_demand_fair_none_when_neither():
    assert fv.book_on_demand_fair(None, None, None, n_legs=2) is None
    assert fv.book_on_demand_fair(None, 3.0, [], n_legs=2) is None


def test_transfer_multiplier_sanity_rejects_px_style_pricer_bugs():
    """Adversarial review #4: PX's known F5-Over bug inflates parlay odds
    5-7x (SANITY_MULT_RATIO defense in the grid path). A 5x-inflated SGP
    decimal deflates the transfer fair; Fréchet's lower bound (often 0)
    doesn't catch it — the implied correlation multiplier does. Bound is
    deliberately LOOSE ([1/3, 3]) so legitimate high-correlation stacks
    (RL+ML ~2x, accepted above) still pass."""
    singles = [(0.52, 0.49), (0.51, 0.48)]
    fair_joint = 0.25
    good_dec = 1.0 / (fair_joint * (0.52 / 0.49) * (0.51 / 0.48))
    assert fv.fair_by_correlation_transfer(good_dec, singles) is not None
    assert fv.fair_by_correlation_transfer(good_dec * 5.0, singles) is None  # bugged 5x
    # inflated the other way (absurdly rich joint) also rejected
    rich_dec = 1.0 / (0.9 * (0.52 / 0.49) * (0.51 / 0.48))
    assert fv.fair_by_correlation_transfer(rich_dec, singles) is None
