"""Issue #84 maker-side scope tests: F5 combos enter scope (route on_demand),
F5-winner combos decline with a distinct labeled reason."""
from kalshi_common import legset

GAME = "26AUG122210KCLAD"            # away=KC, home=LAD


def _f5_spread(team, n, side):
    return {"market_ticker": f"KXMLBF5SPREAD-{GAME}-{team}{n}",
            "event_ticker": f"KXMLBF5SPREAD-{GAME}", "side": side}


def _f5_total(n, side):
    return {"market_ticker": f"KXMLBF5TOTAL-{GAME}-{n}",
            "event_ticker": f"KXMLBF5TOTAL-{GAME}", "side": side}


def _f5_winner(team, side):
    return {"market_ticker": f"KXMLBF5-{GAME}-{team}",
            "event_ticker": f"KXMLBF5-{GAME}", "side": side}


def test_f5_spread_total_rfq_is_priceable():
    """Acceptance: an RFQ of KXMLBF5SPREAD + KXMLBF5TOTAL legs parses and is
    in scope (classifies on_demand), so it reaches the on-demand engine and
    can only decline for lack of book fairs pre-#85 — never out_of_scope*."""
    from kalshi_mlb_mm.main import _priceable
    canon = legset.parse_legs([_f5_spread("LAD", 2, "yes"), _f5_total(7, "yes")])
    assert canon is not None
    assert _priceable(canon)


def test_f5_winner_rfq_is_not_priceable():
    from kalshi_mlb_mm.main import _priceable
    canon = legset.parse_legs([_f5_winner("LAD", "yes"), _f5_total(7, "yes")])
    assert canon is not None
    assert not _priceable(canon)


def test_f5_winner_combo_gets_distinct_scope_reason():
    """The #86 guard must be measurable in the firehose: not the generic
    out_of_scope, and never the non_mlb mislabel."""
    from kalshi_mlb_mm.main import _out_of_scope_reason
    legs = [_f5_winner("LAD", "yes"), _f5_total(7, "yes")]
    canon = legset.parse_legs(legs)
    assert _out_of_scope_reason(legs, canon) == "out_of_scope_f5_winner"


def test_generic_out_of_scope_reason_unchanged_for_fg():
    from kalshi_mlb_mm.main import _out_of_scope_reason
    over = {"market_ticker": f"KXMLBTOTAL-{GAME}-9",
            "event_ticker": f"KXMLBTOTAL-{GAME}", "side": "yes"}
    under = {"market_ticker": f"KXMLBTOTAL-{GAME}-9",
             "event_ticker": f"KXMLBTOTAL-{GAME}", "side": "no"}
    legs = [over, under]               # contradictory pair -> unpriceable
    canon = legset.parse_legs(legs)
    assert _out_of_scope_reason(legs, canon) == "out_of_scope"
