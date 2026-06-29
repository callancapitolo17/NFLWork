"""Phase-1 regression: existing 2-leg combos price identically through the
router as through the legacy _book_fairs + fairs.blended_fair path.

The brief specified a harness using combo_descriptor(legs), but that function
was removed from leg_types.py before this task landed.  We achieve the same
goal by calling main._book_fairs(EVT, spread_line, total_line) directly (the
actual legacy signature) and comparing its median to router.combo_fair.

Three books are included in the fixture so that config.MIN_AGREEING_BOOKS (=3
by default) is satisfied by both paths.
"""
import statistics

import pandas as pd

from kalshi_mlb_mm import config, fairs, main, router

EVT = "KXMLBGAME-25JUN271905NYYBOS"

ST_CELLS = {
    "Home Spread + Over": 4.2,
    "Home Spread + Under": 4.2,
    "Away Spread + Over": 3.8,
    "Away Spread + Under": 3.8,
}


def _grid_rows(book, game_id, spread_line, total_line, cells):
    return [
        {
            "game_id": game_id,
            "combo": c,
            "period": "FG",
            "bookmaker": book,
            "sgp_decimal": d,
            "fetch_time": None,
            "spread_line": spread_line,
            "total_line": total_line,
        }
        for c, d in cells.items()
    ]


# KXMLBSPREAD-25JUN271905NYYBOS-BOS2: BOS is home, '2' -> line_n=2 -> spread_line=-(2-0.5)=-1.5
# KXMLBTOTAL-25JUN271905NYYBOS-9: line_n=9 -> total_line=9-0.5=8.5; side=yes -> Over
LEGS = [
    {
        "market_ticker": "KXMLBSPREAD-25JUN271905NYYBOS-BOS2",
        "event_ticker": EVT,
        "side": "yes",
    },
    {
        "market_ticker": "KXMLBTOTAL-25JUN271905NYYBOS-9",
        "event_ticker": EVT,
        "side": "yes",
    },
]


def test_router_matches_legacy_two_leg_fair():
    """router.combo_fair must equal median(main._book_fairs(...)) for a 2-leg
    spread×total combo — this is the regression lock that Phase-1 must preserve."""
    # three books to satisfy MIN_AGREEING_BOOKS=3
    df = pd.DataFrame(
        _grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS)
        + _grid_rows("fd", EVT, -1.5, 8.5, ST_CELLS)
        + _grid_rows("px", EVT, -1.5, 8.5, ST_CELLS)
    )

    # --- legacy path ---
    old_odds = main._SGP_ODDS
    main._SGP_ODDS = df
    try:
        legacy_book_fairs = main._book_fairs(EVT, -1.5, 8.5)
    finally:
        main._SGP_ODDS = old_odds

    assert legacy_book_fairs, (
        "legacy _book_fairs returned no fairs — fixture must yield >=3 agreeing books"
    )
    # fairs.blended_fair just medians the values; replicate that here
    _, legacy_fair = fairs.blended_fair(LEGS, EVT, legacy_book_fairs)
    assert legacy_fair is not None

    # --- router path ---
    router_fair = router.combo_fair(
        LEGS,
        df,
        lambda gl: EVT,  # identity resolver: event_ticker IS the game_id in the fixture
        config.MIN_AGREEING_BOOKS,
        config.BOOK_CONSENSUS_BAND,
    )
    assert router_fair is not None, "router.combo_fair returned None for a valid 2-leg combo"

    assert abs(router_fair - legacy_fair) < 1e-9, (
        f"router={router_fair!r} legacy={legacy_fair!r} — "
        "2-leg pricing must be identical through both paths"
    )


def test_router_returns_none_for_untypeable_leg():
    """Sanity: router fails safely on a leg it can't classify."""
    legs = [{"market_ticker": "KXMLBPLAYER-foo", "event_ticker": EVT, "side": "yes"}]
    assert router.combo_fair(legs, pd.DataFrame(), lambda gl: EVT, 2, 0.02) is None


def test_router_returns_none_when_resolver_returns_none():
    """Sanity: router fails safely when game resolution fails."""
    df = pd.DataFrame(
        _grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS)
        + _grid_rows("fd", EVT, -1.5, 8.5, ST_CELLS)
        + _grid_rows("px", EVT, -1.5, 8.5, ST_CELLS)
    )
    assert router.combo_fair(LEGS, df, lambda gl: None, 3, 0.02) is None
