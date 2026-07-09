from unabated_edge import feed
def test_apply_deltas_updates_line():
    st = feed.parse_snapshot({"marketSources":[{"id":7,"name":"S"}],"teams":{"1":{"name":"A"},"2":{"name":"B"}},
        "gameOddsEvents":{"lg21:pt1:pregame":[{"eventId":1,"eventStart":"2026-06-22T17:00:00+00:00",
            "eventTeams":{"1":{"id":1},"0":{"id":2}},
            "gameOddsMarketSourcesLines":{"si1:ms7:an0":{"bt1":{"price":-150}}}}]}}, {"lg21"})
    feed.apply_deltas(st, [{"lg21:pt1:pregame":[{"eventId":1,"eventStart":"2026-06-22T17:00:00+00:00",
        "eventTeams":{"1":{"id":1},"0":{"id":2}},
        "gameOddsMarketSourcesLines":{"si1:ms7:an0":{"bt1":{"price":-135}}}}]}], {"lg21"})
    assert st.lines["1|1|7|bt1"]["price"] == -135
