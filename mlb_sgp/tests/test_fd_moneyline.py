"""FanDuel moneyline-in-SGP: market-map entry + the ML×total combo emission.
FD's implyBets prices any two (marketId, selectionId) legs, so ML+total is the
same call as spread+total (verified live: market 'Moneyline', team-named
runners)."""
import sys
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock

sys.path.insert(0, str(Path(__file__).parent.parent))

from mlb_sgp._shared import TargetLine  # noqa: E402


def test_market_map_has_moneyline():
    from scraper_fanduel_sgp import _MARKET_MAP
    assert _MARKET_MAP["Moneyline"] == ("fg", "moneyline", "main")
    assert _MARKET_MAP["First 5 Innings Money Line"] == ("f5", "moneyline", "main")


def _runners_cache_with_ml():
    return {
        "g1": {
            "fg": {
                "moneyline": {"home": ("734.ml", 29164), "away": ("734.ml", 29163)},
                "totals": {
                    ("O", 8.5): ("734.tot", 1001),
                    ("U", 8.5): ("734.tot", 1002),
                },
            }
        }
    }


def test_fd_ml_total_emits_four_cell_grid():
    from mlb_sgp.fanduel import _price_ml_total_for_games

    targets = [TargetLine("g1", "Toronto Blue Jays", "Texas Rangers",
                          datetime(2026, 6, 23, 23, 0, tzinfo=timezone.utc),
                          "FG", -1.5, 8.5)]
    price_combo = MagicMock(return_value={"decimal": 2.6, "american": 160})
    rows = _price_ml_total_for_games(
        MagicMock(), _runners_cache_with_ml(), targets, price_combo,
        datetime.now(timezone.utc), 2, False)

    combos = {r.combo for r in rows}
    assert combos == {"Home ML + Over", "Home ML + Under",
                      "Away ML + Over", "Away ML + Under"}
    for r in rows:
        assert r.spread_line is None
        assert r.total_line == 8.5
        assert r.bookmaker == "fanduel"
        assert r.source == "fanduel_direct"
        assert r.sgp_decimal == 2.6
    # home ML selection 29164 and away 29163 both used as the FIRST leg.
    first_leg_sels = {call.args[2] for call in price_combo.call_args_list}
    assert first_leg_sels == {29164, 29163}


def test_fd_ml_total_skips_without_moneyline():
    from mlb_sgp.fanduel import _price_ml_total_for_games
    cache = _runners_cache_with_ml()
    cache["g1"]["fg"]["moneyline"] = {}
    targets = [TargetLine("g1", "TOR", "TEX",
                          datetime(2026, 6, 23, 23, 0, tzinfo=timezone.utc),
                          "FG", -1.5, 8.5)]
    rows = _price_ml_total_for_games(
        MagicMock(), cache, targets, MagicMock(return_value={"decimal": 2.0, "american": 100}),
        datetime.now(timezone.utc), 2, False)
    assert rows == []


def test_fd_ml_total_skips_when_price_none():
    from mlb_sgp.fanduel import _price_ml_total_for_games
    targets = [TargetLine("g1", "TOR", "TEX",
                          datetime(2026, 6, 23, 23, 0, tzinfo=timezone.utc),
                          "FG", -1.5, 8.5)]
    rows = _price_ml_total_for_games(
        MagicMock(), _runners_cache_with_ml(), targets, MagicMock(return_value=None),
        datetime.now(timezone.utc), 2, False)
    assert rows == []
