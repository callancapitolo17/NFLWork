import datetime
from unabated_edge import feed


def _v2_raw():
    """Miniature of the live v2 per-league file (verified shape, 2026-07-08)."""
    return {
        "marketSources": [{"id": 7, "name": "Sharp Book Price"}, {"id": 2, "name": "FanDuel"}],
        "teams": {"2024": {"name": "France"}, "2058": {"name": "Morocco"}},
        "odds": {
            "lg21:pt1:pregame": [
                {
                    "eventId": 128061, "eventStart": "2026-07-09T20:00:00",
                    "eventName": "Morocco @ France", "betTypeId": 3,
                    "key": "pt1:pregame:bt3:e128061",
                    "sides": {
                        "si0:tid2058": {
                            "ms7": {"points": 2.5, "americanPrice": 104, "isBlurred": False,
                                    "alternateLines": [{"points": 3.5, "americanPrice": 254}]},
                            "ms2": {"points": 2.5, "americanPrice": -110, "isBlurred": True},  # blurred
                        },
                        "si1:tid2024": {
                            "ms7": {"points": 2.5, "americanPrice": -117, "isBlurred": False,
                                    "alternateLines": [{"points": 3.5, "americanPrice": -347}]},
                        },
                    },
                },
                {   # betType outside _V2_BET_TYPES (e.g. corners) must be skipped
                    "eventId": 128061, "eventStart": "2026-07-09T20:00:00",
                    "eventName": "Morocco @ France", "betTypeId": 1456,
                    "sides": {"si0:tid2058": {"ms7": {"points": 9.5, "americanPrice": -105}}},
                },
                {   # malformed row (no eventId) must be isolated, not kill the batch
                    "betTypeId": 3, "sides": {},
                },
            ],
            "lg21:pt2:pregame": [   # non-pt1 period must be ignored
                {"eventId": 128061, "eventStart": "2026-07-09T20:00:00", "betTypeId": 3,
                 "sides": {"si0:tid2058": {"ms7": {"points": 1.0, "americanPrice": 100}}}},
            ],
        },
    }


def test_parse_v2_builds_state():
    st = feed.parse_v2(_v2_raw(), "lg21")
    e = st.events[128061]
    assert e.home == "France" and e.away == "Morocco"     # si1=home, si0=away
    assert e.start_utc == datetime.datetime(2026, 7, 9, 20, 0, tzinfo=datetime.timezone.utc)
    # line keys use the same "eid|si|ms|bt" convention as the legacy feed
    over = st.lines["128061|0|7|bt3"]
    assert feed.line_american_price(over) == 104 and over["points"] == 2.5
    assert over["alternateLines"][0]["points"] == 3.5     # alt ladder preserved
    assert "128061|1|7|bt3" in st.lines


def test_parse_v2_blur_guard_and_filters():
    st = feed.parse_v2(_v2_raw(), "lg21")
    assert "128061|0|2|bt3" not in st.lines               # blurred line dropped
    assert not any("bt1456" in k for k in st.lines)       # non-core betType skipped
    assert not any(v.get("points") == 1.0 for v in st.lines.values())  # pt2 period ignored
    # malformed row was isolated: the good rows still landed
    assert len(st.lines) == 2
