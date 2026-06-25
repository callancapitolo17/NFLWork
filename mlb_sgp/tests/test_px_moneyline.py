"""ProphetX moneyline-in-SGP: the ML×total emission via the RFQ pricer.

PX prices a parlay by submitting an RFQ (submit_parlay_rfq) with legs built from
(market_id, selection). ML outcomes live in the market's top-level `outcomes`
list keyed by competitorId. We reuse _price_combos_parallel with the moneyline
market id as the first leg's market — exactly an ML+total parlay."""
import sys
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock

sys.path.insert(0, str(Path(__file__).parent.parent))

from mlb_sgp._shared import TargetLine  # noqa: E402


def test_px_market_names_has_moneyline():
    from scraper_prophetx_sgp import MARKET_NAMES
    assert MARKET_NAMES["fg"]["moneyline"] == "Moneyline"


def _markets():
    # Moneyline market: flat top-level outcomes keyed by competitorId.
    ml = {
        "id": "ml1", "name": "Moneyline", "marketLines": [],
        "outcomes": [
            {"id": "o_home", "competitorId": 100, "odds": -150, "line": None},
            {"id": "o_away", "competitorId": 200, "odds": 130, "line": None},
        ],
    }
    # Total market: over/under under marketLines[].selections at line 8.5.
    total = {
        "id": "tot1", "name": "Total Runs",
        "marketLines": [{"selections": [
            [{"id": "over", "name": "Over", "line": 8.5, "odds": -110}],
            [{"id": "under", "name": "Under", "line": 8.5, "odds": -110}],
        ]}],
        "outcomes": [],
    }
    return [ml, total]


def _matched():
    return {"g1": {"px_event_id": "e1", "px_home_competitor_id": 100,
                   "px_away_competitor_id": 200, "game_id": "g1"}}


def test_px_ml_total_emits_four_cell_grid():
    from mlb_sgp.prophetx import _price_ml_total_for_games

    targets = [TargetLine("g1", "Reds", "Brewers",
                          datetime(2026, 6, 24, 23, 0, tzinfo=timezone.utc),
                          "FG", -1.5, 8.5)]
    client = MagicMock()
    # offer ladder: american +250 -> decimal 3.5; passes sanity vs ~1.67*1.91.
    client.submit_parlay_rfq.return_value = ({"odds": 250}, False)

    rows = _price_ml_total_for_games(
        client, {"g1": _markets()}, _matched(), targets,
        datetime.now(timezone.utc), 2, False)

    combos = {r.combo for r in rows}
    assert combos == {"Home ML + Over", "Home ML + Under",
                      "Away ML + Over", "Away ML + Under"}
    for r in rows:
        assert r.spread_line is None
        assert r.total_line == 8.5
        assert r.bookmaker == "prophetx"
        assert r.source == "prophetx_direct"
        assert abs(r.sgp_decimal - 3.5) < 1e-6


def test_px_ml_total_skips_without_moneyline_outcomes():
    from mlb_sgp.prophetx import _price_ml_total_for_games
    markets = _markets()
    markets[0]["outcomes"] = []   # no ML outcomes
    targets = [TargetLine("g1", "Reds", "Brewers",
                          datetime(2026, 6, 24, 23, 0, tzinfo=timezone.utc),
                          "FG", -1.5, 8.5)]
    client = MagicMock()
    client.submit_parlay_rfq.return_value = ({"odds": 250}, False)
    rows = _price_ml_total_for_games(
        client, {"g1": markets}, _matched(), targets,
        datetime.now(timezone.utc), 2, False)
    assert rows == []
