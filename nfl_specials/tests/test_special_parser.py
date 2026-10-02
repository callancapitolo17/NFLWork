from nfl_specials.special_parser import Fecta, Leg, ParseFailure, parse_fecta


def test_plain_trifecta_is_three_period_wins():
    fecta = parse_fecta("SEAHAWKS TRIFECTA (1Q, 1H & GM)")
    assert fecta == Fecta(team="SEA", prop_type="TRIFECTA", legs=(
        Leg("win", "Q1"), Leg("win", "H1"), Leg("win", "GM")))


def test_spread_trifecta_reads_half_points():
    fecta = parse_fecta("PANTHERS TRIFECTA (1Q +½, 1H +2½ & GM +3½)")
    assert fecta.team == "CAR"
    assert fecta.legs == (Leg("spread", "Q1", 0.5), Leg("spread", "H1", 2.5), Leg("spread", "GM", 3.5))


def test_minus_half_without_leading_digit():
    fecta = parse_fecta("SEAHAWKS TRIFECTA (1Q -½, 1H -4½ & GM -7½)")
    assert fecta.legs[0] == Leg("spread", "Q1", -0.5)


def test_superfecta_leads_with_scores_first():
    fecta = parse_fecta("SEAHAWKS SUPERFECTA (SCR 1ST, 1Q, 1H & GM -7½)")
    assert fecta.prop_type == "SUPERFECTA"
    assert fecta.legs == (Leg("scores_first"), Leg("win", "Q1"), Leg("win", "H1"), Leg("spread", "GM", -7.5))


def test_non_fecta_special_is_refused():
    assert isinstance(parse_fecta("CHIEFS, PANTHERS & FALCONS ALL TO WIN"), ParseFailure)


def test_wrong_leg_count_is_refused():
    result = parse_fecta("SEAHAWKS TRIFECTA (1Q & GM)")
    assert isinstance(result, ParseFailure)
    assert "expected 3" in result.reason


def test_unknown_leg_is_refused_not_guessed():
    assert isinstance(parse_fecta("SEAHAWKS TRIFECTA (1Q, 2H & GM)"), ParseFailure)


def test_unknown_team_is_refused():
    assert isinstance(parse_fecta("SONICS TRIFECTA (1Q, 1H & GM)"), ParseFailure)
