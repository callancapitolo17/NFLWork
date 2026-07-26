import logging

import pytest
from unabated_edge import feed, storage


def build_bt3_state(league_prefix, home_name, away_name, eid=9, ms=7,
                     over_price=115, under_price=-140, line=2.5,
                     event_start="2026-06-22T17:00:00+00:00"):
    """FeedState with bt3 total over/under lines for one event on `league_prefix`.

    Shared feed-state builder for the totals-ladder sport adapters
    (test_soccer_adapter.py, test_mlb_adapter.py) — moved here so both test
    files can build a minimal Unabated snapshot without duplicating the
    gameOddsEvents scaffolding."""
    return feed.parse_snapshot({
        "marketSources": [{"id": ms, "name": "S"}],
        "teams": {"1": {"name": home_name}, "2": {"name": away_name}},
        "gameOddsEvents": {
            f"{league_prefix}:pt1:pregame": [{
                "eventId": eid,
                "eventStart": event_start,
                "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
                "gameOddsMarketSourcesLines": {
                    f"si0:ms{ms}:an0": {"bt3": {"price": over_price, "points": line}},  # Over
                    f"si1:ms{ms}:an0": {"bt3": {"price": under_price, "points": line}},  # Under
                }
            }]
        }
    }, {league_prefix})


@pytest.fixture(autouse=True)
def _clean_storage_buffer():
    storage._BUFFER.clear()
    yield
    storage._BUFFER.clear()


@pytest.fixture(autouse=True, scope="session")
def _no_live_log_pollution():
    """Detach the rotating bot.log handler during tests: test EDGE/WARNING lines
    in the live log masquerade as real engine output (found during monitoring).
    caplog still works — it captures via propagation, not via these handlers."""
    log = logging.getLogger("unabated_edge")
    saved = log.handlers[:]
    log.handlers = []
    yield
    log.handlers = saved
