"""Phase 0: moneyline leg type + 2-way singles devig (pure functions)."""
import pandas as pd

from kalshi_common import fair_value as fv
from kalshi_common.leg_types import _leg_dict_to_typed


# ---- MoneylineLeg model hit-mask -------------------------------------------

def _samples(home_margins):
    # total_final_score is irrelevant for moneyline; fill with a constant.
    return pd.DataFrame({"home_margin": home_margins,
                         "total_final_score": [9] * len(home_margins)})


def test_moneyline_home_yes_counts_home_wins():
    s = _samples([3, -2, 1, -5, 4])  # 3 home wins
    leg = fv.MoneylineLeg(team_is_home=True, side="yes")
    assert abs(fv.model_fair(s, [leg]) - 3 / 5) < 1e-9


def test_moneyline_away_yes_counts_away_wins():
    s = _samples([3, -2, 1, -5, 4])  # 2 away wins
    leg = fv.MoneylineLeg(team_is_home=False, side="yes")
    assert abs(fv.model_fair(s, [leg]) - 2 / 5) < 1e-9


def test_moneyline_no_is_complement():
    s = _samples([3, -2, 1, -5, 4])
    yes = fv.MoneylineLeg(team_is_home=True, side="yes")
    no = fv.MoneylineLeg(team_is_home=True, side="no")
    assert abs(fv.model_fair(s, [yes]) + fv.model_fair(s, [no]) - 1.0) < 1e-9


# ---- leg-typing KXMLBGAME --------------------------------------------------

def test_leg_dict_to_typed_moneyline_home():
    # Event suffix encodes away=TEX, home=LAA. KXMLBGAME market is for LAA (home).
    leg = {"event_ticker": "KXMLBGAME-26MAY232205TEXLAA",
           "market_ticker": "KXMLBGAME-26MAY232205TEXLAA-LAA", "side": "yes"}
    typed = _leg_dict_to_typed(leg, game_id="g")
    assert isinstance(typed, fv.MoneylineLeg)
    assert typed.team_is_home is True
    assert typed.side == "yes"


def test_leg_dict_to_typed_moneyline_away():
    leg = {"event_ticker": "KXMLBGAME-26MAY232205TEXLAA",
           "market_ticker": "KXMLBGAME-26MAY232205TEXLAA-TEX", "side": "yes"}
    typed = _leg_dict_to_typed(leg, game_id="g")
    assert isinstance(typed, fv.MoneylineLeg)
    assert typed.team_is_home is False


# ---- 2-way singles devig ---------------------------------------------------

def test_devig_two_way_sums_to_one():
    p_yes, p_no = fv.devig_two_way(1.90, 1.90)
    assert abs(p_yes + p_no - 1.0) < 1e-9
    assert abs(p_yes - 0.5) < 1e-9


def test_devig_two_way_favorite_above_half():
    # -200 favorite (decimal 1.50) vs +170 dog (decimal 2.70)
    p_fav, p_dog = fv.devig_two_way(1.50, 2.70)
    assert p_fav > 0.5 > p_dog
    assert abs(p_fav + p_dog - 1.0) < 1e-9
