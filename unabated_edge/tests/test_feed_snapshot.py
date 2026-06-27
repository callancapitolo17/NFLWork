import datetime
from unabated_edge import feed

def _snap():
    return {"marketSources":[{"id":7,"name":"Sharp Book Price"},{"id":2,"name":"FanDuel"}],
            "teams":{"2063":{"name":"Argentina"},"2038":{"name":"Austria"}},
            "gameOddsEvents":{"lg21:pt1:pregame":[{
                "eventId":126539,"eventStart":"2026-06-22T17:00:00+00:00",
                "eventTeams":{"1":{"id":2063},"0":{"id":2038}},
                "gameOddsMarketSourcesLines":{
                    "si1:ms7:an0":{"bt1":{"marketSourceId":7,"points":None,"price":-150}},
                    "si0:ms7:an0":{"bt1":{"marketSourceId":7,"points":None,"price":130}}}}],
              "lg5:pt1:pregame":[{"eventId":1,"eventStart":"2026-06-22T17:00:00+00:00",
                "eventTeams":{"1":{"id":2063},"0":{"id":2038}},"gameOddsMarketSourcesLines":{}}]}}

def test_parse_filters_to_registered_leagues():
    st = feed.parse_snapshot(_snap(), {"lg21"})
    assert st.books[7] == "Sharp Book Price"
    assert 126539 in st.events and 1 not in st.events      # lg5 filtered out
    e = st.events[126539]
    assert e.home == "Argentina" and e.away == "Austria"
    assert e.start_utc == datetime.datetime(2026,6,22,17,0,tzinfo=datetime.timezone.utc)
    assert st.lines["126539|1|7|bt1"]["price"] == -150
