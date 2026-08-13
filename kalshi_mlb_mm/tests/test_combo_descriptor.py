"""Read-side combo descriptor: classify a Kalshi 2-leg combo into a grid lookup
(spread+total or ml+total), with the devig cell derived from the actual sides."""
from kalshi_common.leg_types import combo_descriptor

SUF = "26MAY232205TEXLAA"   # away TEX, home LAA


def _spread(team, n, side):
    return {"event_ticker": f"KXMLBSPREAD-{SUF}",
            "market_ticker": f"KXMLBSPREAD-{SUF}-{team}{n}", "side": side}


def _total(n, side):
    return {"event_ticker": f"KXMLBTOTAL-{SUF}",
            "market_ticker": f"KXMLBTOTAL-{SUF}-{n}", "side": side}


def _ml(team, side):
    return {"event_ticker": f"KXMLBGAME-{SUF}",
            "market_ticker": f"KXMLBGAME-{SUF}-{team}", "side": side}


def test_spread_total_home_over():
    d = combo_descriptor([_spread("LAA", 2, "yes"), _total(9, "yes")])
    assert d.kind == "spread_total"
    assert d.target_combo == "Home Spread + Over"
    assert d.spread_line == -1.5 and d.total_line == 8.5


def test_spread_total_away_under_via_sides():
    # away team LAA? no — away is TEX. Away spread leg + under.
    d = combo_descriptor([_spread("TEX", 2, "yes"), _total(9, "no")])
    assert d.target_combo == "Away Spread + Under"   # NOT hardcoded Home+Over
    # Issue #70: TEX is the away team, so its margin market is the
    # away-favourite grid — home-perspective +1.5, not -1.5.
    assert d.spread_line == 1.5


def test_spread_line_sign_four_way():
    """Issue #70: the four (ticker team x side) combos must land on four
    distinct (spread_line, cell) pairs."""
    got = {}
    for team, side in (("LAA", "yes"), ("LAA", "no"),
                       ("TEX", "yes"), ("TEX", "no")):
        d = combo_descriptor([_spread(team, 2, side), _total(9, "yes")])
        got[(team, side)] = (d.spread_line, d.target_combo)
    assert got == {
        ("LAA", "yes"): (-1.5, "Home Spread + Over"),   # home -1.5
        ("LAA", "no"):  (-1.5, "Away Spread + Over"),   # away +1.5
        ("TEX", "yes"): (1.5, "Away Spread + Over"),    # away -1.5
        ("TEX", "no"):  (1.5, "Home Spread + Over"),    # home +1.5
    }
    assert len(set(got.values())) == 4


def test_spread_no_side_is_other_team():
    # home LAA spread, side 'no' == away covers -> "Away Spread"
    d = combo_descriptor([_spread("LAA", 2, "no"), _total(9, "yes")])
    assert d.target_combo == "Away Spread + Over"


def test_ml_total_home_over():
    d = combo_descriptor([_ml("LAA", "yes"), _total(9, "yes")])
    assert d.kind == "ml_total"
    assert d.target_combo == "Home ML + Over"
    assert d.spread_line is None and d.total_line == 8.5
    assert d.combo_family[0] == "Home ML + Over"


def test_ml_away_team_yes_is_away_ml():
    d = combo_descriptor([_ml("TEX", "yes"), _total(9, "no")])
    assert d.target_combo == "Away ML + Under"


def test_ml_away_team_no_is_home_ml():
    # away TEX 'no' == home wins -> "Home ML"
    d = combo_descriptor([_ml("TEX", "no"), _total(9, "yes")])
    assert d.target_combo == "Home ML + Over"


def test_codes_resolved():
    d = combo_descriptor([_ml("LAA", "yes"), _total(9, "yes")])
    assert d.away_code == "TEX" and d.home_code == "LAA"


def test_unsupported_shapes_return_none():
    assert combo_descriptor([_spread("LAA", 2, "yes")]) is None          # 1 leg
    assert combo_descriptor([_ml("LAA", "yes"), _ml("TEX", "no")]) is None  # ml+ml
    assert combo_descriptor([_spread("LAA", 2, "yes"),
                             _spread("TEX", 2, "no")]) is None           # 2 spreads
