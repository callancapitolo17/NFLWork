"""Novig moneyline-in-SGP: the ML×total emission path.

Novig names its moneyline market MONEY (FG) / MONEY_1H (F5); the outcomes carry
competitor.symbol exactly like SPREAD, so fetch_event_legs reuses the spread
symbol-matcher to resolve home/away ML legs. price_sgps then prices the 4
FG ML×total combos inline (reusing the per-target over/under legs), tagging
spread_line=None — the NULL marker for "no spread leg".
"""
import sys
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock

sys.path.insert(0, str(Path(__file__).parent.parent))

from mlb_sgp._shared import TargetLine  # noqa: E402


def test_money_type_constant():
    from scraper_novig_sgp import MONEY_TYPE
    assert MONEY_TYPE == {"fg": "MONEY", "f5": "MONEY_1H"}


def _markets(home_sym="DET", away_sym="HOU"):
    """Minimal FG market tree: MONEY + SPREAD(-1.5) + TOTAL(8.5)."""
    return [
        {"type": "MONEY", "strike": 0.0, "is_consensus": True, "outcomes": [
            {"id": "ml_h", "description": home_sym, "available": 0.52,
             "competitor": {"symbol": home_sym}},
            {"id": "ml_a", "description": away_sym, "available": 0.49,
             "competitor": {"symbol": away_sym}},
        ]},
        {"type": "SPREAD", "strike": -1.5, "is_consensus": True, "outcomes": [
            {"id": "sp_h", "description": f"{home_sym} -1.5", "available": 0.40,
             "competitor": {"symbol": home_sym}},
            {"id": "sp_a", "description": f"{away_sym} +1.5", "available": 0.62,
             "competitor": {"symbol": away_sym}},
        ]},
        {"type": "TOTAL", "strike": 8.5, "is_consensus": True, "outcomes": [
            {"id": "ov", "description": "Over 8.5", "available": 0.50},
            {"id": "un", "description": "Under 8.5", "available": 0.50},
        ]},
    ]


def test_fetch_event_legs_extracts_ml(monkeypatch):
    """fetch_event_legs pulls home_ml/away_ml from the MONEY market by symbol."""
    import scraper_novig_sgp as nv
    # NB: do NOT bind a local named `_markets` here — it would shadow the
    # module-level `_markets()` helper that the mock lambda closes over, and
    # the resulting NameError would be swallowed by fetch_event_legs' try/except.
    monkeypatch.setattr(nv, "_gql",
                        lambda *_a, **_k: {"data": {"event": [{"markets": _markets()}]}})
    game = {
        "game_id": "g1", "nv_event_id": "e1",
        "nv_home_sym": "DET", "nv_away_sym": "HOU",
        "fg_spread_line": -1.5, "fg_total_line": 8.5,
        "f5_spread_line": None, "f5_total_line": None,
    }
    legs, _mkts = nv.fetch_event_legs(MagicMock(), game, verbose=False)
    assert legs["fg"]["home_ml"] == {"id": "ml_h", "available": 0.52}
    assert legs["fg"]["away_ml"] == {"id": "ml_a", "available": 0.49}


def test_price_sgps_emits_ml_rows(monkeypatch):
    """End-to-end: price_sgps emits the 4 FG ML×total combos, spread_line=None."""
    import scraper_novig_sgp as nv
    from mlb_sgp import novig

    markets = _markets()
    legs = {
        "fg": {"home_spread": {"id": "sp_h", "available": 0.40},
               "away_spread": {"id": "sp_a", "available": 0.62},
               "over": {"id": "ov", "available": 0.50},
               "under": {"id": "un", "available": 0.50},
               "home_ml": {"id": "ml_h", "available": 0.52},
               "away_ml": {"id": "ml_a", "available": 0.49}},
        "f5": {"home_spread": None, "away_spread": None, "over": None,
               "under": None, "home_ml": None, "away_ml": None},
    }

    monkeypatch.setattr(nv, "match_events", lambda nv_events, parlay_lines: [{
        "game_id": "g1", "nv_event_id": "e1",
        "nv_home_sym": "DET", "nv_away_sym": "HOU",
        "fg_spread_line": -1.5, "fg_total_line": 8.5,
        "f5_spread_line": None, "f5_total_line": None,
    }])
    monkeypatch.setattr(nv, "fetch_event_legs", lambda *_a, **_k: (legs, markets))
    monkeypatch.setattr(nv, "submit_parlay",
                        lambda _s, _ids, **_k: ({"decimal": 3.5, "american": 250}, False))

    client = MagicMock()
    client.list_events.return_value = [
        MagicMock(event_id="e1", home_team="Tigers", away_team="Astros",
                  home_sym="DET", away_sym="HOU", start_time="")
    ]

    targets = [TargetLine("g1", "Tigers", "Astros",
                          datetime(2026, 6, 25, 23, 0, tzinfo=timezone.utc),
                          "FG", -1.5, 8.5)]
    rows = novig.price_sgps(targets, client=client, verbose=False)

    ml_rows = [r for r in rows if "ML" in r.combo]
    assert {r.combo for r in ml_rows} == {
        "Home ML + Over", "Home ML + Under", "Away ML + Over", "Away ML + Under"}
    for r in ml_rows:
        assert r.spread_line is None
        assert r.total_line == 8.5
        assert r.bookmaker == "novig"
        assert abs(r.sgp_decimal - 3.5) < 1e-6
