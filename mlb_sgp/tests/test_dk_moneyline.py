"""DraftKings moneyline-in-SGP: selection extraction, market-num capture, and
the ML×total combo emission. Mirrors the inline-fixture style of the other DK
tests (no captured recon blobs)."""
import sys
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock

sys.path.insert(0, str(Path(__file__).parent.parent))

from mlb_sgp._shared import TargetLine, NO_LINE  # noqa: E402


# ---------------------------------------------------------------------------
# _extract_moneyline_selections — pure parse of 0ML{mnum}_{participant}
# ---------------------------------------------------------------------------
def test_extract_moneyline_home_away():
    from scraper_draftkings_sgp import _extract_moneyline_selections
    # Real DK shape: 0ML{mnum}_1 (home), 0ML{mnum}_3 (away). Embed in noise.
    text = 'x"0ML85265198_1"y "0ML85265198_3" "0HC85265198N150_3" "0OU111O850_1"'
    out = _extract_moneyline_selections(text, "85265198")
    assert out["1"] == ["0ML85265198_1"]
    assert out["3"] == ["0ML85265198_3"]


def test_extract_moneyline_filters_other_market_nums():
    from scraper_draftkings_sgp import _extract_moneyline_selections
    # A different ML market_num (e.g. an F5 moneyline) must NOT leak in.
    text = '"0ML85265198_1" "0ML99999999_1" "0ML99999999_3"'
    out = _extract_moneyline_selections(text, "85265198")
    assert out["1"] == ["0ML85265198_1"]
    assert out["3"] == []


def test_extract_moneyline_empty_when_no_mnum():
    from scraper_draftkings_sgp import _extract_moneyline_selections
    assert _extract_moneyline_selections("0ML123_1", None) == {"1": [], "3": []}


# ---------------------------------------------------------------------------
# fetch_main_market_nums — captures the Moneyline market number
# ---------------------------------------------------------------------------
def test_fetch_main_market_nums_captures_moneyline(monkeypatch):
    import scraper_draftkings_sgp as dk

    def fake_subcat(session, eid, subcat):
        if subcat == "4519":
            return [("1_85265198", "Moneyline"), ("2_85265198", "Run Line"),
                    ("3_85265198", "Total")]
        return []  # no F5
    monkeypatch.setattr(dk, "_fetch_subcat_markets", fake_subcat)
    out = dk.fetch_main_market_nums(MagicMock(), "evt")
    assert out["fg"]["moneyline"] == "85265198"
    assert out["fg"]["run_line"] == "85265198"
    assert out["fg"]["total"] == "85265198"


# ---------------------------------------------------------------------------
# _price_ml_total_for_games — emits the 4-cell ML×total grid
# ---------------------------------------------------------------------------
def _cache_with_ml():
    return {
        "g1": {
            "sel_ids_all": {
                "fg": {
                    "moneyline": {"1": ["0ML200_1"], "3": ["0ML200_3"]},
                    "totals": {
                        ("O", 8.5): ["0OU201O850_1"],
                        ("U", 8.5): ["0OU201U850_1"],
                    },
                    "canonical": {"201"},
                }
            }
        }
    }


def test_ml_total_emits_four_cell_grid():
    from mlb_sgp.draftkings import _price_ml_total_for_games
    from mlb_sgp.scraper_draftkings_sgp import _market_num

    targets = [TargetLine("g1", "NYY", "BOS",
                          datetime(2026, 6, 23, 23, 0, tzinfo=timezone.utc),
                          "FG", -1.5, 8.5)]
    calc = MagicMock(return_value={"trueOdds": 2.5, "displayOdds": "+150"})
    rows = _price_ml_total_for_games(
        MagicMock(), _cache_with_ml(), targets, calc, _market_num,
        datetime.now(timezone.utc), False, 2)

    combos = {r.combo for r in rows}
    assert combos == {"Home ML + Over", "Home ML + Under",
                      "Away ML + Over", "Away ML + Under"}
    for r in rows:
        assert r.spread_line == NO_LINE      # no spread leg
        assert r.total_line == 8.5
        assert r.bookmaker == "draftkings"
        assert r.source == "draftkings_direct"
        assert r.sgp_decimal == 2.5
    # Home ML paired with over/under sels; calc called with the ML sel first.
    ml_args = {call.args[1] for call in calc.call_args_list}
    assert "0ML200_1" in ml_args and "0ML200_3" in ml_args


def test_ml_total_skips_game_without_moneyline():
    from mlb_sgp.draftkings import _price_ml_total_for_games
    from mlb_sgp.scraper_draftkings_sgp import _market_num
    cache = _cache_with_ml()
    cache["g1"]["sel_ids_all"]["fg"]["moneyline"] = {"1": [], "3": []}
    targets = [TargetLine("g1", "NYY", "BOS",
                          datetime(2026, 6, 23, 23, 0, tzinfo=timezone.utc),
                          "FG", -1.5, 8.5)]
    rows = _price_ml_total_for_games(
        MagicMock(), cache, targets, MagicMock(return_value={"trueOdds": 2.0}),
        _market_num, datetime.now(timezone.utc), False, 2)
    assert rows == []


def test_ml_total_skips_when_calc_returns_none():
    from mlb_sgp.draftkings import _price_ml_total_for_games
    from mlb_sgp.scraper_draftkings_sgp import _market_num
    targets = [TargetLine("g1", "NYY", "BOS",
                          datetime(2026, 6, 23, 23, 0, tzinfo=timezone.utc),
                          "FG", -1.5, 8.5)]
    rows = _price_ml_total_for_games(
        MagicMock(), _cache_with_ml(), targets, MagicMock(return_value=None),
        _market_num, datetime.now(timezone.utc), False, 2)
    assert rows == []
