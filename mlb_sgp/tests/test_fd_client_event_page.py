"""Tests for FanDuelClient.fetch_event_page (both-tab combined fetch)."""
from mlb_sgp.fd_client import FanDuelClient, Market, Runner


class _FakeResp:
    def __init__(self, payload):
        self._payload = payload
        self.status_code = 200

    def json(self):
        return self._payload


class _FakeSession:
    """Records the last URL requested; returns a fixed payload.

    Mirrors the unit-test FakeSession contract the client already supports:
    its get() takes no headers/timeout kwargs, so the client's
    `except TypeError` fallback path is exercised.
    """
    def __init__(self, payload):
        self._payload = payload
        self.last_url = None

    def get(self, url):
        self.last_url = url
        return _FakeResp(self._payload)


def _client_with(payload):
    # Bypass __init__ (which builds a real curl_cffi session) and inject a fake.
    c = FanDuelClient.__new__(FanDuelClient)
    c.session = _FakeSession(payload)
    c.verbose = False
    return c


_PAYLOAD = {
    "attachments": {
        "markets": {
            "m1": {
                "marketId": "m1",
                "marketName": "First 7 Innings Total Runs",
                "runners": [
                    {"selectionId": "r1", "runnerName": "Over",
                     "runnerStatus": "ACTIVE", "handicap": 7.5,
                     "winRunnerOdds": {"americanDisplayOdds": {"americanOdds": -128}}},
                    {"selectionId": "r2", "runnerName": "Under",
                     "runnerStatus": "ACTIVE", "handicap": 7.5,
                     "winRunnerOdds": {"americanDisplayOdds": {"americanOdds": 104}}},
                ],
            }
        }
    }
}


def test_fetch_event_page_returns_markets_and_runners():
    c = _client_with(_PAYLOAD)
    markets, runners = c.fetch_event_page("99", tab="")
    assert markets == [Market(market_id="m1", name="First 7 Innings Total Runs")]
    assert {r.runner_id for r in runners} == {"r1", "r2"}
    over = next(r for r in runners if r.runner_id == "r1")
    assert over.american_odds == -128 and over.line == 7.5


def test_fetch_event_page_uses_requested_tab():
    c = _client_with(_PAYLOAD)
    c.fetch_event_page("99", tab="")
    assert "tab=&" in c.session.last_url  # tab="" is passed through literally
    c.fetch_event_page("99", tab="same-game-parlay-")
    assert "tab=same-game-parlay-" in c.session.last_url


def test_wrappers_delegate_to_fetch_event_page():
    c = _client_with(_PAYLOAD)
    assert c.fetch_event_markets("99", tab="") == [
        Market(market_id="m1", name="First 7 Innings Total Runs")
    ]
    assert {r.runner_id for r in c.fetch_event_runners("99", tab="")} == {"r1", "r2"}
