import datetime
from unabated_edge import feed
from unabated_edge.sports.soccer import Soccer

def _state_with_event():
    return feed.parse_snapshot({"marketSources":[{"id":7,"name":"S"}],
        "teams":{"1":{"name":"Korea Republic"},"2":{"name":"Mexico"}},
        "gameOddsEvents":{"lg21:pt1:pregame":[{"eventId":9,"eventStart":"2026-06-22T17:00:00+00:00",
            "eventTeams":{"1":{"id":1},"0":{"id":2}},
            "gameOddsMarketSourcesLines":{
                "si1:ms7:an0":{"bt1":{"price":-150},"bt4":{"price":300}},
                "si0:ms7:an0":{"bt1":{"price":400}}}}]}}, {"lg21"})

def test_canon_alias():
    assert Soccer().canon_team("Korea Republic") == "South Korea"


def test_canon_alias_expanded_collapses_both_spellings():
    """A divergent spelling and its plain-English form must collapse to one string."""
    s = Soccer()
    assert s.canon_team("Côte d'Ivoire") == "Ivory Coast"
    assert s.canon_team("Ivory Coast") == "Ivory Coast"      # plain form falls through to itself
    assert s.canon_team("Czechia") == s.canon_team("Czech Republic")
    assert s.canon_team("Türkiye") == "Turkey"

def test_fair_three_way_sums_to_one():
    s = Soccer(); st = _state_with_event()
    fair = s.fair(st, st.events[9])
    assert set(fair) == {"home","draw","away"}
    assert abs(sum(fair.values()) - 1.0) < 1e-6
