import pytest

from kalshi_mlb_rfq import combo_enumerator


def test_enumerate_yields_n_times_m_combos():
    spread_legs = [(-1.5, "home"), (-1.5, "away"), (-2.5, "home"), (-2.5, "away")]
    total_legs = [7.5, 8.5, 9.5]

    combos = list(combo_enumerator.enumerate_2leg(
        game_id="gid",
        event_suffix="26APR282005NYYTEX",
        home_code="TEX", away_code="NYY",
        available_spreads=spread_legs,
        available_totals=total_legs,
    ))
    # 4 spread legs × 2 sides = 8 spread leg specs
    # 3 total legs × 2 sides = 6 total leg specs
    # combos = 8 × 6 = 48
    assert len(combos) == 48
    assert all(len(c.legs) == 2 for c in combos)


def test_canonical_leg_set_hash_is_order_invariant():
    legs_a = [{"market_ticker": "A", "side": "yes"}, {"market_ticker": "B", "side": "no"}]
    legs_b = [{"market_ticker": "B", "side": "no"}, {"market_ticker": "A", "side": "yes"}]
    assert combo_enumerator.canonical_leg_set_hash(legs_a) == \
           combo_enumerator.canonical_leg_set_hash(legs_b)


def test_priority_queue_orders_by_edge_magnitude():
    candidates_with_scores = [
        (combo_enumerator.ComboCandidate(game_id="g", legs=(), leg_set_hash="a",
                                         descriptor="A"), 0.32, 0.30),  # |Δ|=0.02
        (combo_enumerator.ComboCandidate(game_id="g", legs=(), leg_set_hash="b",
                                         descriptor="B"), 0.40, 0.20),  # |Δ|=0.20  ← top
        (combo_enumerator.ComboCandidate(game_id="g", legs=(), leg_set_hash="c",
                                         descriptor="C"), 0.50, 0.45),  # |Δ|=0.05
    ]
    ranked = combo_enumerator.rank_by_edge(candidates_with_scores)
    assert [c.leg_set_hash for c in ranked] == ["b", "c", "a"]


def test_score_falls_back_to_distance_from_half_when_ref_zero():
    score = combo_enumerator.edge_score(blended_fair=0.30, kalshi_ref=0.0)
    assert score == pytest.approx(0.20)
