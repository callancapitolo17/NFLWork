"""Pin _spread_line_from_legs' signed contract (issue #70).

The helper's only production caller is the taker's research telemetry, so
nothing else exercises it — without these tests a refactor could silently
revert it to the pre-#70 always-negative collapse.
"""
from kalshi_common.leg_types import _spread_line_from_legs

GAME = "25JUN271905NYYBOS"          # away=NYY, home=BOS


def _leg(market_ticker, event_ticker, side="yes"):
    return {"market_ticker": market_ticker,
            "event_ticker": event_ticker, "side": side}


def _spread(team, n):
    return _leg(f"KXMLBSPREAD-{GAME}-{team}{n}", f"KXMLBSPREAD-{GAME}")


def test_home_margin_ticker_is_negative():
    assert _spread_line_from_legs([_spread("BOS", 2)]) == -1.5


def test_away_margin_ticker_is_positive():
    assert _spread_line_from_legs([_spread("NYY", 2)]) == 1.5


def test_no_spread_leg_keeps_legacy_zero():
    total = _leg(f"KXMLBTOTAL-{GAME}-9", f"KXMLBTOTAL-{GAME}")
    assert _spread_line_from_legs([total]) == 0.0


def test_unresolvable_event_ticker_returns_none():
    # A silently wrong sign is worse than a NULL in telemetry.
    bad = _leg(f"KXMLBSPREAD-{GAME}-BOS2", "KXMLBSPREAD-GARBAGE")
    assert _spread_line_from_legs([bad]) is None


def test_malformed_leg_returns_none():
    no_side = {"market_ticker": f"KXMLBSPREAD-{GAME}-BOS2",
               "event_ticker": f"KXMLBSPREAD-{GAME}"}
    assert _spread_line_from_legs([no_side]) is None
