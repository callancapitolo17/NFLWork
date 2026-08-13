"""Issue #87 maker-side scope tests: RFI-containing RFQs enter scope.

KXMLBRFI legs re-encode at parse as I1 totals at 0.5, so once parse_legs
types them the RFQ is in scope automatically — these tests pin that (a
parse regression would silently demote RFI RFQs to out_of_scope_non_mlb,
which #83 showed is a measurement blind spot).
"""
from kalshi_common import legset

GAME = "26AUG132210MILLAD"           # away=MIL, home=LAD


def _rfi(side):
    return {"market_ticker": f"KXMLBRFI-{GAME}",
            "event_ticker": f"KXMLBRFI-{GAME}", "side": side}


def _fg_ml(team, side):
    return {"market_ticker": f"KXMLBGAME-{GAME}-{team}",
            "event_ticker": f"KXMLBGAME-{GAME}", "side": side}


def _fg_total(n, side):
    return {"market_ticker": f"KXMLBTOTAL-{GAME}-{n}",
            "event_ticker": f"KXMLBTOTAL-{GAME}", "side": side}


def test_rfi_plus_ml_rfq_is_priceable():
    from kalshi_mlb_mm.main import _priceable
    canon = legset.parse_legs([_rfi("yes"), _fg_ml("LAD", "yes")])
    assert canon is not None
    assert _priceable(canon)


def test_rfi_plus_total_rfq_is_priceable():
    from kalshi_mlb_mm.main import _priceable
    canon = legset.parse_legs([_rfi("no"), _fg_total(9, "yes")])
    assert canon is not None
    assert _priceable(canon)


def test_rfi_rfq_never_labels_non_mlb():
    """The pre-#87 failure mode: parse_legs returned None on the RFI leg and
    the whole RFQ was misfiled out_of_scope_non_mlb."""
    from kalshi_mlb_mm.main import _out_of_scope_reason
    legs = [_rfi("yes"), _fg_ml("LAD", "yes")]
    canon = legset.parse_legs(legs)
    assert canon is not None
    assert _out_of_scope_reason(legs, canon) != "out_of_scope_non_mlb"


def test_contradictory_rfi_pair_stays_unpriceable():
    from kalshi_mlb_mm.main import _priceable
    canon = legset.parse_legs([_rfi("yes"), _rfi("no")])
    assert canon is not None
    assert not _priceable(canon)
