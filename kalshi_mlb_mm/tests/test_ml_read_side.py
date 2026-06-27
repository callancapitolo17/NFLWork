"""Read-side end-to-end: the maker looks up an ML×total combo in the grid
(spread_line NULL) and devigs the correct cell; and the spread+total path
devigs the cell matching the actual sides (not the old hardcoded Home+Over)."""
import numpy as np
import pandas as pd

from kalshi_mlb_mm import main
from kalshi_common.leg_types import combo_descriptor

SUF = "26MAY232205TEXLAA"   # away TEX, home LAA
_BOOKS = ["dk", "fd", "px"]


def _grid(rows):
    cols = ["game_id", "combo", "period", "bookmaker",
            "sgp_decimal", "fetch_time", "spread_line", "total_line"]
    return pd.DataFrame(rows, columns=cols)


def _ml_grid():
    # 4-cell ML×total per book; spread_line = NaN (NULL). Implied probs
    # [0.30,0.25,0.20,0.30] (sum 1.05 vig) for Home/Away ML × Over/Under.
    decs = {"Home ML + Over": 1 / 0.30, "Home ML + Under": 1 / 0.25,
            "Away ML + Over": 1 / 0.20, "Away ML + Under": 1 / 0.30}
    rows = []
    for b in _BOOKS:
        for combo, dec in decs.items():
            rows.append(["g1", combo, "FG", b, dec, None, np.nan, 8.5])
    return _grid(rows)


def test_ml_total_book_fairs_from_grid(monkeypatch):
    monkeypatch.setattr(main, "_SGP_ODDS", _ml_grid())
    legs = [{"event_ticker": f"KXMLBGAME-{SUF}",
             "market_ticker": f"KXMLBGAME-{SUF}-LAA", "side": "yes"},
            {"event_ticker": f"KXMLBTOTAL-{SUF}",
             "market_ticker": f"KXMLBTOTAL-{SUF}-9", "side": "yes"}]
    desc = combo_descriptor(legs)
    assert desc.kind == "ml_total" and desc.target_combo == "Home ML + Over"
    fairs = main._book_fairs("g1", desc)
    assert set(fairs) == {"dk", "fd", "px"}          # all 3 books priced
    # Devigged Home ML + Over ~ 0.30/1.05 area; just assert sane + consistent.
    for f in fairs.values():
        assert 0.25 < f < 0.32


def test_ml_lookup_ignores_spread_total_rows(monkeypatch):
    # Grid has BOTH an ML×total family and a spread+total family for the game.
    df_ml = _ml_grid()
    sp_rows = []
    for b in _BOOKS:
        for combo in ("Home Spread + Over", "Home Spread + Under",
                      "Away Spread + Over", "Away Spread + Under"):
            sp_rows.append(["g1", combo, "FG", b, 2.0, None, -1.5, 8.5])
    df = pd.concat([df_ml, _grid(sp_rows)], ignore_index=True)
    monkeypatch.setattr(main, "_SGP_ODDS", df)
    legs = [{"event_ticker": f"KXMLBGAME-{SUF}",
             "market_ticker": f"KXMLBGAME-{SUF}-LAA", "side": "yes"},
            {"event_ticker": f"KXMLBTOTAL-{SUF}",
             "market_ticker": f"KXMLBTOTAL-{SUF}-9", "side": "yes"}]
    fairs = main._book_fairs("g1", combo_descriptor(legs))
    # 4 ML cells per book only (spread rows excluded by combo family + NULL filter)
    assert set(fairs) == {"dk", "fd", "px"}
    for f in fairs.values():
        assert 0.25 < f < 0.32


def test_spread_total_devigs_correct_cell(monkeypatch):
    # Asymmetric 4-cell so the chosen cell matters. Away Spread + Under should
    # NOT return the Home Spread + Over value (the old hardcode bug).
    decs = {"Home Spread + Over": 1 / 0.40, "Home Spread + Under": 1 / 0.20,
            "Away Spread + Over": 1 / 0.20, "Away Spread + Under": 1 / 0.25}
    rows = []
    for b in _BOOKS:
        for combo, dec in decs.items():
            rows.append(["g1", combo, "FG", b, dec, None, -1.5, 8.5])
    monkeypatch.setattr(main, "_SGP_ODDS", _grid(rows))
    # Away team (TEX) spread yes + Under -> "Away Spread + Under"
    legs = [{"event_ticker": f"KXMLBSPREAD-{SUF}",
             "market_ticker": f"KXMLBSPREAD-{SUF}-TEX2", "side": "yes"},
            {"event_ticker": f"KXMLBTOTAL-{SUF}",
             "market_ticker": f"KXMLBTOTAL-{SUF}-9", "side": "no"}]
    desc = combo_descriptor(legs)
    assert desc.target_combo == "Away Spread + Under"
    fairs = main._book_fairs("g1", desc)
    # devigged Away Spread + Under ~ 0.25/1.05 ≈ 0.238, clearly NOT the
    # Home Spread + Over cell (~0.40/1.05 ≈ 0.38).
    for f in fairs.values():
        assert 0.21 < f < 0.27
