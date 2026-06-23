"""Phase 3: generalized scope gate (spread/total/moneyline; props out)."""
from kalshi_mlb_mm import scope


def _leg(mt):
    return {"event_ticker": mt.rsplit("-", 1)[0], "market_ticker": mt, "side": "yes"}


def test_moneyline_only_in_scope():
    assert scope.is_in_scope([_leg("KXMLBGAME-26MAY232205TEXLAA-LAA")]) is True


def test_moneyline_parlay_in_scope():
    assert scope.is_in_scope([
        _leg("KXMLBGAME-26MAY232205TEXLAA-LAA"),
        _leg("KXMLBGAME-26MAY232205NYYBOS-BOS"),
        _leg("KXMLBGAME-26MAY232205CHCMIL-MIL"),
    ]) is True


def test_mixed_spread_total_moneyline_in_scope():
    assert scope.is_in_scope([
        _leg("KXMLBSPREAD-26MAY232205TEXLAA-LAA2"),
        _leg("KXMLBTOTAL-26MAY232205TEXLAA-9"),
        _leg("KXMLBGAME-26MAY232205NYYBOS-BOS"),
    ]) is True


def test_prop_legs_in_scope():
    # Props are now supported (priced off Odds API per-event prop markets).
    assert scope.is_in_scope([
        _leg("KXMLBGAME-26MAY232205TEXLAA-LAA"),
        _leg("KXMLBHR-26MAY232205TEXLAA-LAASOMEONE1"),
    ]) is True
    assert scope.is_in_scope([_leg("KXMLBKS-26MAY232205TEXLAA-TEXPITCHER5")]) is True


def test_non_mlb_leg_out_of_scope():
    assert scope.is_in_scope([_leg("KXWCGAME-SOMETHING-X")]) is False


def test_empty_out_of_scope():
    assert scope.is_in_scope([]) is False
