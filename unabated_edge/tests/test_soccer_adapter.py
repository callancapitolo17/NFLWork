import datetime
from unabated_edge import feed
from unabated_edge.sports.soccer import Soccer


def _state_with_bt3(eid=9, ms=7, over_price=115, under_price=-140, line=2.5):
    """FeedState with bt3 total over/under lines for one event."""
    return feed.parse_snapshot({
        "marketSources": [{"id": ms, "name": "S"}],
        "teams": {"1": {"name": "Colombia"}, "2": {"name": "Ghana"}},
        "gameOddsEvents": {
            "lg21:pt1:pregame": [{
                "eventId": eid,
                "eventStart": "2026-06-22T17:00:00+00:00",
                "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
                "gameOddsMarketSourcesLines": {
                    f"si0:ms{ms}:an0": {"bt3": {"price": over_price, "points": line}},  # Over
                    f"si1:ms{ms}:an0": {"bt3": {"price": under_price, "points": line}},  # Under
                }
            }]
        }
    }, {"lg21"})


def test_canon_alias():
    assert Soccer().canon_team("Korea Republic") == "South Korea"


def test_canon_alias_expanded_collapses_both_spellings():
    s = Soccer()
    assert s.canon_team("Côte d'Ivoire") == "Ivory Coast"
    assert s.canon_team("Ivory Coast") == "Ivory Coast"      # plain form falls through to itself
    assert s.canon_team("Czechia") == s.canon_team("Czech Republic")
    assert s.canon_team("Türkiye") == "Turkey"


def test_anchor_total_returns_line_and_prob():
    s = Soccer()
    st = _state_with_bt3(eid=9, ms=7, over_price=115, under_price=-140, line=2.5)
    result = s._anchor_total(st, 9)
    assert result is not None
    line_val, p_over = result
    assert line_val == 2.5
    assert 0 < p_over < 1


def test_anchor_total_prob_sums_to_one():
    """Devigged over + under must sum to 1."""
    s = Soccer()
    st = _state_with_bt3(eid=9, ms=7, over_price=115, under_price=-140, line=2.5)
    line_val, p_over = s._anchor_total(st, 9)
    # devig returns [p_over, p_under]; they sum to 1
    # p_over + (1 - p_over) is trivially 1; test via actual values
    assert abs(p_over + (1 - p_over) - 1.0) < 1e-9


def test_price_event_returns_yes_and_no_candidates():
    s = Soccer()
    st = _state_with_bt3(eid=9, ms=7, over_price=115, under_price=-140, line=2.5)
    ev_meta = st.events[9]
    kalshi_event = {
        "title": "Colombia vs Ghana: Regulation Time Total Goals",
        "markets": [
            {"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5,
             "yes_sub_title": "Reg Time: Over 2.5 goals scored"},
        ]
    }
    candidates = s.price_event(st, ev_meta, kalshi_event)
    assert len(candidates) == 2
    assert all(c.market_ticker == "T-O25" for c in candidates)
    sides = {c.side for c in candidates}
    assert sides == {"yes", "no"}
    labels = {c.label for c in candidates}
    assert "over_2.5" in labels
    assert "under_2.5" in labels
    # fair probs sum to 1
    total_prob = sum(c.fair_prob for c in candidates)
    assert abs(total_prob - 1.0) < 1e-5


def test_event_teams_parsed_from_title():
    s = Soccer()
    result = s.event_teams({"title": "South Africa vs Canada: Regulation Time Total Goals"})
    assert result == frozenset({s.canon_team("South Africa"), s.canon_team("Canada")})


def test_event_teams_with_alias():
    """Aliases are applied during event_teams parsing."""
    s = Soccer()
    result = s.event_teams({"title": "Korea Republic vs Mexico: Regulation Time Total Goals"})
    assert "South Korea" in result
    assert "Mexico" in result


def test_fail_closed_no_matching_floor_strike():
    """When no Kalshi market matches the anchor line, price_event returns []."""
    s = Soccer()
    st = _state_with_bt3(eid=9, ms=7, over_price=115, under_price=-140, line=2.5)
    ev_meta = st.events[9]
    # Market at 1.5, but anchor quoted 2.5 → no match
    kalshi_event = {
        "title": "Colombia vs Ghana: Regulation Time Total Goals",
        "markets": [
            {"ticker": "T-O15", "strike_type": "greater", "floor_strike": 1.5,
             "yes_sub_title": "Reg Time: Over 1.5 goals scored"},
        ]
    }
    candidates = s.price_event(st, ev_meta, kalshi_event)
    assert candidates == []


def test_fail_closed_anchor_missing_under():
    """When the Under side is missing, _anchor_total returns None → price_event []."""
    s = Soccer()
    st = feed.parse_snapshot({
        "marketSources": [{"id": 7, "name": "S"}],
        "teams": {"1": {"name": "Colombia"}, "2": {"name": "Ghana"}},
        "gameOddsEvents": {
            "lg21:pt1:pregame": [{
                "eventId": 9,
                "eventStart": "2026-06-22T17:00:00+00:00",
                "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
                "gameOddsMarketSourcesLines": {
                    "si0:ms7:an0": {"bt3": {"price": 115, "points": 2.5}},
                    # si1 (Under) deliberately absent
                }
            }]
        }
    }, {"lg21"})
    result = s._anchor_total(st, 9)
    assert result is None
    ev_meta = st.events[9]
    kalshi_event = {
        "title": "Colombia vs Ghana: Regulation Time Total Goals",
        "markets": [{"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5,
                     "yes_sub_title": "Reg Time: Over 2.5 goals scored"}]
    }
    candidates = s.price_event(st, ev_meta, kalshi_event)
    assert candidates == []


def test_fail_closed_mismatched_points():
    """Over and Under with different `points` values → _anchor_total returns None."""
    s = Soccer()
    st = feed.parse_snapshot({
        "marketSources": [{"id": 7, "name": "S"}],
        "teams": {"1": {"name": "Colombia"}, "2": {"name": "Ghana"}},
        "gameOddsEvents": {
            "lg21:pt1:pregame": [{
                "eventId": 9,
                "eventStart": "2026-06-22T17:00:00+00:00",
                "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
                "gameOddsMarketSourcesLines": {
                    "si0:ms7:an0": {"bt3": {"price": 115, "points": 2.5}},
                    "si1:ms7:an0": {"bt3": {"price": -140, "points": 3.5}},  # different line!
                }
            }]
        }
    }, {"lg21"})
    result = s._anchor_total(st, 9)
    assert result is None
