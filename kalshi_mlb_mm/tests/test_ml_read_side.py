"""Read-side grid semantics on the ROUTER's grid path (#81 ported these off
the deleted main._book_fairs oracle — router.grid_cell_fairs is the one
surviving implementation): an ML×total combo matches grid rows with
spread_line NULL and devigs the correct cell; the spread+total path devigs
the cell matching the actual sides (not the old hardcoded Home+Over)."""
import numpy as np
import pandas as pd

from kalshi_mlb_mm import router
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


def _cell_fairs(desc, df):
    return router.grid_cell_fairs("g1", desc.combo_family, desc.spread_line,
                                  desc.total_line, desc.target_combo, df)


def test_ml_total_book_fairs_from_grid():
    legs = [{"event_ticker": f"KXMLBGAME-{SUF}",
             "market_ticker": f"KXMLBGAME-{SUF}-LAA", "side": "yes"},
            {"event_ticker": f"KXMLBTOTAL-{SUF}",
             "market_ticker": f"KXMLBTOTAL-{SUF}-9", "side": "yes"}]
    desc = combo_descriptor(legs)
    assert desc.kind == "ml_total" and desc.target_combo == "Home ML + Over"
    fairs = _cell_fairs(desc, _ml_grid())
    assert set(fairs) == {"dk", "fd", "px"}          # all 3 books priced
    # Devigged Home ML + Over ~ 0.30/1.05 area; just assert sane + consistent.
    for f in fairs.values():
        assert 0.25 < f < 0.32


def test_ml_lookup_ignores_spread_total_rows():
    # Grid has BOTH an ML×total family and a spread+total family for the game.
    df_ml = _ml_grid()
    sp_rows = []
    for b in _BOOKS:
        for combo in ("Home Spread + Over", "Home Spread + Under",
                      "Away Spread + Over", "Away Spread + Under"):
            sp_rows.append(["g1", combo, "FG", b, 2.0, None, -1.5, 8.5])
    df = pd.concat([df_ml, _grid(sp_rows)], ignore_index=True)
    legs = [{"event_ticker": f"KXMLBGAME-{SUF}",
             "market_ticker": f"KXMLBGAME-{SUF}-LAA", "side": "yes"},
            {"event_ticker": f"KXMLBTOTAL-{SUF}",
             "market_ticker": f"KXMLBTOTAL-{SUF}-9", "side": "yes"}]
    fairs = _cell_fairs(combo_descriptor(legs), df)
    # 4 ML cells per book only (spread rows excluded by combo family + NULL filter)
    assert set(fairs) == {"dk", "fd", "px"}
    for f in fairs.values():
        assert 0.25 < f < 0.32


def test_spread_total_devigs_correct_cell():
    # Asymmetric 4-cell so the chosen cell matters. Away Spread + Under should
    # NOT return the Home Spread + Over value (the old hardcode bug).
    # spread_line is +1.5: an away-team margin ticker canonicalizes to
    # +(N-0.5) home-perspective (#70 signed convention). The pre-#81 version
    # of this test seeded -1.5 rows and passed VACUOUSLY (empty fairs dict,
    # value-loop over nothing) — the set assertion below prevents a repeat.
    decs = {"Home Spread + Over": 1 / 0.40, "Home Spread + Under": 1 / 0.20,
            "Away Spread + Over": 1 / 0.20, "Away Spread + Under": 1 / 0.25}
    rows = []
    for b in _BOOKS:
        for combo, dec in decs.items():
            rows.append(["g1", combo, "FG", b, dec, None, 1.5, 8.5])
    # Away team (TEX) spread yes + Under -> "Away Spread + Under"
    legs = [{"event_ticker": f"KXMLBSPREAD-{SUF}",
             "market_ticker": f"KXMLBSPREAD-{SUF}-TEX2", "side": "yes"},
            {"event_ticker": f"KXMLBTOTAL-{SUF}",
             "market_ticker": f"KXMLBTOTAL-{SUF}-9", "side": "no"}]
    desc = combo_descriptor(legs)
    assert desc.target_combo == "Away Spread + Under"
    fairs = _cell_fairs(desc, _grid(rows))
    # devigged Away Spread + Under ~ 0.25/1.05 ≈ 0.238, clearly NOT the
    # Home Spread + Over cell (~0.40/1.05 ≈ 0.38).
    assert set(fairs) == {"dk", "fd", "px"}
    for f in fairs.values():
        assert 0.21 < f < 0.27
