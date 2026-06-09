from kalshi_mlb_mm import scope

_LEGS_OK = [
    {"event_ticker": "KXMLBSPREAD-26MAY232205TEXLAA", "market_ticker": "KXMLBSPREAD-26MAY232205TEXLAA-TEX2", "side": "no"},
    {"event_ticker": "KXMLBTOTAL-26MAY232205TEXLAA", "market_ticker": "KXMLBTOTAL-26MAY232205TEXLAA-10", "side": "no"},
]

def test_decode_legs_reads_mve_selected_legs():
    market = {"mve_selected_legs": _LEGS_OK}
    assert scope.decode_legs(market) == _LEGS_OK

def test_decode_legs_none_when_absent():
    assert scope.decode_legs({"ticker": "x"}) is None

def test_in_scope_true_for_spread_total_2leg():
    assert scope.is_spread_total_2leg(_LEGS_OK) is True

def test_in_scope_false_for_three_legs():
    assert scope.is_spread_total_2leg(_LEGS_OK + [_LEGS_OK[0]]) is False

def test_in_scope_false_for_two_spreads():
    legs = [_LEGS_OK[0], dict(_LEGS_OK[0])]
    assert scope.is_spread_total_2leg(legs) is False
