import pandas as pd
import pytest

from kalshi_mlb_rfq import fair_value


@pytest.fixture
def synthetic_samples():
    return pd.DataFrame({
        "home_margin": [3, -1, 5, 2, 1, 0, -2, 4, 6, -1],
        "total_final_score": [10, 6, 12, 9, 8, 7, 4, 11, 13, 5],
    })


def test_model_fair_home_minus_1_5_and_over_7_5(synthetic_samples):
    # Home -1.5 = home_margin >= 2; Over 7.5 = total >= 8.
    # Hits: rows where both: (3,10),(5,12),(2,9),(4,11),(6,13) → 5 of 10
    legs = [
        fair_value.SpreadLeg(team_is_home=True, line_n=2, side="yes"),
        fair_value.TotalLeg(line_n=8, side="yes"),
    ]
    fair = fair_value.model_fair(synthetic_samples, legs)
    assert fair == pytest.approx(0.5)


def test_model_fair_no_side_inverts(synthetic_samples):
    # NO side of "Over 7.5" = total_final_score < 8 = 4 of 10 rows (6, 7, 4, 5)
    legs = [fair_value.TotalLeg(line_n=8, side="no")]
    assert fair_value.model_fair(synthetic_samples, legs) == pytest.approx(0.4)


@pytest.fixture
def four_side_dk():
    return pd.DataFrame({
        "combo": ["Home Spread + Over", "Home Spread + Under",
                  "Away Spread + Over", "Away Spread + Under"],
        "bookmaker": ["draftkings"] * 4,
        "sgp_decimal": [3.50, 4.20, 5.00, 4.80],
        "spread_line": [-1.5] * 4,
        "total_line": [8.5] * 4,
    })


def test_devig_with_4_sides(four_side_dk):
    fair = fair_value.devig_book(four_side_dk, combo="Home Spread + Over")
    expected_vig = sum(1 / d for d in [3.50, 4.20, 5.00, 4.80])
    expected_fair = (1 / 3.50) / expected_vig
    assert fair == pytest.approx(expected_fair, abs=1e-4)


def test_devig_fallback_when_fewer_sides():
    df = pd.DataFrame({
        "combo": ["Home Spread + Over"],
        "bookmaker": ["draftkings"],
        "sgp_decimal": [3.50],
        "spread_line": [-1.5],
        "total_line": [8.5],
    })
    fair = fair_value.devig_book(df, combo="Home Spread + Over",
                                 vig_fallback=0.125)
    expected = (1 / 3.50) / 1.125
    assert fair == pytest.approx(expected, abs=1e-4)


def test_blended_requires_two_sources():
    blended = fair_value.blend(model_fair_value=0.30, book_fairs={})
    assert blended is None

    blended = fair_value.blend(model_fair_value=0.30, book_fairs={"dk": 0.34})
    assert blended == pytest.approx(0.32)

    blended = fair_value.blend(model_fair_value=0.30,
                               book_fairs={"dk": 0.34, "fd": 0.28})
    assert blended == pytest.approx((0.30 + 0.34 + 0.28) / 3)
