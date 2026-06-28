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


def _state_split_books(draw_book):
    """Home+away moneylines on book 7; draw bt4 on `draw_book`.
    bt4 merges into book 7's home dict when draw_book==7 (one dict per side+book,
    as the real feed nests it), else lives in a separate book's entry."""
    home7 = {"bt1": {"price": -150}}
    lines = {"si1:ms7:an0": home7, "si0:ms7:an0": {"bt1": {"price": 400}}}
    if draw_book == 7:
        home7["bt4"] = {"price": 300}
    else:
        lines[f"si1:ms{draw_book}:an0"] = {"bt4": {"price": 300}}
    return feed.parse_snapshot({"marketSources":[{"id":7,"name":"S7"},{"id":6,"name":"S6"}],
        "teams":{"1":{"name":"A"},"2":{"name":"B"}},
        "gameOddsEvents":{"lg21:pt1:pregame":[{"eventId":9,"eventStart":"2026-06-22T17:00:00+00:00",
            "eventTeams":{"1":{"id":1},"0":{"id":2}},
            "gameOddsMarketSourcesLines":lines}]}}, {"lg21"})


def test_fair_rejects_cross_book_three_way():
    """No single book has home+away+draw -> reject (no Frankenstein devig)."""
    s = Soccer(); st = _state_split_books(draw_book=6)   # draw on book 6, ML on book 7
    assert s.fair(st, st.events[9]) is None


def test_fair_accepts_single_book_three_way():
    """All three legs on the same book -> fair is produced."""
    s = Soccer(); st = _state_split_books(draw_book=7)   # everything on book 7
    fair = s.fair(st, st.events[9])
    assert fair is not None and abs(sum(fair.values()) - 1.0) < 1e-6


def test_fair_rejects_spread_as_moneyline():
    """A bt1 line carrying `points` is a spread, not a moneyline -> not used as anchor."""
    st = feed.parse_snapshot({"marketSources":[{"id":7,"name":"S"}],
        "teams":{"1":{"name":"A"},"2":{"name":"B"}},
        "gameOddsEvents":{"lg21:pt1:pregame":[{"eventId":9,"eventStart":"2026-06-22T17:00:00+00:00",
            "eventTeams":{"1":{"id":1},"0":{"id":2}},
            "gameOddsMarketSourcesLines":{
                "si1:ms7:an0":{"bt1":{"price":-150,"points":-1.5},"bt4":{"price":300}},  # spread, not ML
                "si0:ms7:an0":{"bt1":{"price":400}}}}]}}, {"lg21"})
    assert Soccer().fair(st, st.events[9]) is None
