from kalshi_common import legset

EVT = "KXMLBGAME-25JUN271905NYYBOS"  # away=NYY, home=BOS

def _spread(team, n, side):   # team code, ticker N (2 => ±1.5), yes/no
    return {"market_ticker": f"KXMLBSPREAD-25JUN271905NYYBOS-{team}{n}",
            "event_ticker": EVT, "side": side}

def _total(n, side):
    return {"market_ticker": f"KXMLBTOTAL-25JUN271905NYYBOS-{n}",
            "event_ticker": EVT, "side": side}

def _ml(team, side):
    return {"market_ticker": f"KXMLBGAME-25JUN271905NYYBOS-{team}",
            "event_ticker": EVT, "side": side}

def test_parse_spread_home_minus_1_5():
    leg = legset.parse_leg(_spread("BOS", 2, "yes"))   # home -1.5, YES
    assert leg == legset.CanonicalLeg(EVT, "spread", -1.5, "home")

def test_parse_total_over_8_5():
    leg = legset.parse_leg(_total(9, "yes"))            # Over 8.5
    assert leg == legset.CanonicalLeg(EVT, "total", 8.5, "over")

def test_parse_ml_home():
    leg = legset.parse_leg(_ml("BOS", "yes"))           # home ML
    assert leg == legset.CanonicalLeg(EVT, "ml", None, "home")

def test_parse_spread_no_side_flips_to_away():
    # home -1.5 NO == away covers +1.5
    leg = legset.parse_leg(_spread("BOS", 2, "no"))
    assert leg == legset.CanonicalLeg(EVT, "spread", -1.5, "away")

def test_unparseable_leg_returns_none():
    bad = {"market_ticker": "KXMLBPLAYER-foo", "event_ticker": EVT, "side": "yes"}
    assert legset.parse_leg(bad) is None

def test_parse_legs_all_or_nothing():
    good = [_spread("BOS", 2, "yes"), _total(9, "yes")]
    assert legset.parse_legs(good) is not None
    assert len(legset.parse_legs(good)) == 2
    bad = good + [{"market_ticker": "KXMLBPLAYER-foo", "event_ticker": EVT, "side": "yes"}]
    assert legset.parse_legs(bad) is None
