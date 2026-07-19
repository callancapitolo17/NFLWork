"""Phase 2 shared on-demand types."""
import pytest
from mlb_sgp._shared import ResolvedLeg, OnDemandBookResult, GameRef


def test_resolved_leg_defaults_and_frozen():
    rl = ResolvedLeg(ref="0HC123_1")
    assert rl.opposite_ref is None and rl.single_decimal is None \
        and rl.opposite_decimal is None
    with pytest.raises(Exception):
        rl.ref = "x"


def test_on_demand_book_result_frozen():
    r = OnDemandBookResult(book="draftkings", fair=0.141, route="partition",
                           n_cells_priced=8, latency_sec=3.2)
    assert r.route in ("partition", "transfer")
    with pytest.raises(Exception):
        r.fair = 0.5


def test_game_ref_fields():
    g = GameRef(game_id="abc", home_team="Boston Red Sox",
                away_team="New York Yankees", commence_time=None)
    assert g.home_team.startswith("Boston")
