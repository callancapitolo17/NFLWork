"""ProphetX's live prices for the Edges overlay (prophetx_live.py, 2026-11-05):
the reader's game reads, its order and expiry of wanted games, and the bets
service's POST /prophetx_live.json. Scripted client, no network."""
import json
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import prophetx_live, service
from unabated_ticket.bets_service.prophetx_client import ProphetXError
from unabated_ticket.bets_service.store import BetsStore
from unabated_ticket.bets_service.tests.test_prophetx import EVENT, MARKETS, PX_EVENT_ID, team_keys_reader

LINE = {"league": "nfl", "awayTeam": "Buffalo Bills", "homeTeam": "Kansas City Chiefs", "awayTeamId": "1", "homeTeamId": "2"}
GAME = {"gameKey": "123742", "eventStart": "2026-11-01T21:25:00.000Z", "strikeIds": ["e" * 32], "line": LINE}


class Clock:
    def __init__(self):
        self.now = 1_790_000_000.0

    def __call__(self):
        return self.now


class LiveClient:
    """The reads the live reader makes, counted."""

    def __init__(self, event=None, markets=None):
        self.access_key = "live1234abcd"
        self.blocked_until = 0.0
        self.calls = []
        self.event = EVENT if event is None else event
        self.markets_list = MARKETS if markets is None else markets

    def strikes(self, strike_ids):
        self.calls.append("strikes")
        return [{"event_id": PX_EVENT_ID, "strike_id": strike_id} for strike_id in strike_ids]

    def sport_events(self, event_ids=None):
        self.calls.append("sport_events")
        return [self.event]

    def markets(self, event_id):
        self.calls.append("markets")
        if isinstance(self.markets_list, Exception):
            raise self.markets_list
        return self.markets_list


def reader(client=None, clock=None, team_keys=None):
    return prophetx_live.ProphetXLiveBooks(client or LiveClient(), team_keys or team_keys_reader(), clock=clock or Clock())


def side(read, bet_type, period, which, points):
    found = [line for line in read["lines"] if (line["betType"], line["period"], line["side"], line["points"]) == (bet_type, period, which, points)]
    assert len(found) == 1, found
    return found[0]


def test_a_game_read_prices_every_side_at_its_best_offer_and_reuses_the_event():
    client = LiveClient()
    books = reader(client)
    read = books.read_game(GAME)
    assert read["prophetxEventId"] == PX_EVENT_ID and read["totalsOnly"] is False
    assert side(read, "spread", "FG", "away", 3.5) == {"betType": "spread", "period": "FG", "side": "away", "points": 3.5,
                                                      "strikeId": "c" * 32, "price": -110, "depth": 200.0}
    assert side(read, "total", "FG", "under", 47.5)["price"] == -115
    assert side(read, "moneyline", "FG", "home", None)["price"] == -160
    assert client.calls == ["strikes", "sport_events", "markets"]
    books.read_game(GAME)
    assert client.calls[3:] == ["markets"]  # the event is reused: one ProphetX call per later read


def test_teams_the_other_way_round_keep_totals_only_and_other_teams_read_as_an_error():
    swapped = {**EVENT, "competitors": [{"id": 11, "name": "Kansas City Chiefs", "side": "away"},
                                        {"id": 12, "name": "Buffalo Bills", "side": "home"}]}
    assert reader(LiveClient(event=swapped)).read_game(GAME)["totalsOnly"] is True
    other = {**EVENT, "competitors": [{"id": 11, "name": "Kansas City Chiefs", "side": "home"},
                                      {"id": 13, "name": "Miami Dolphins", "side": "away"}]}
    assert "does not read as the board's" in reader(LiveClient(event=other)).read_game(GAME)["error"]


def test_a_started_game_or_a_failed_read_is_an_error_and_the_event_is_found_again_after_a_failure():
    assert reader(LiveClient(event={**EVENT, "status": "live"})).read_game(GAME)["error"] == f"ProphetX game {PX_EVENT_ID} is live"
    client = LiveClient(markets=ProphetXError("GET /mm/get_markets: HTTP 500"))
    books = reader(client)
    assert "could not be read" in books.read_game(GAME)["error"]
    books.read_game(GAME)
    assert client.calls == ["strikes", "sport_events", "markets", "strikes", "sport_events", "markets"]


def test_the_ask_never_waits_and_the_thread_reads_the_game_read_longest_ago_first():
    clock = Clock()
    books = reader(clock=clock)
    other = {**GAME, "gameKey": "999"}
    answer = books.want([GAME, other])
    assert answer["games"] == {"123742": None, "999": None} and answer["reader"]["wanted"] == 2
    first = books._next_game()
    books._reads[first["gameKey"]] = books.read_game(first)
    clock.now += 2
    assert books._next_game()["gameKey"] != first["gameKey"]
    assert books.want([first])["games"][first["gameKey"]]["prophetxEventId"] == PX_EVENT_ID


def test_a_game_no_longer_asked_for_is_dropped_with_its_read():
    clock = Clock()
    books = reader(clock=clock)
    books.want([GAME])
    books._reads["123742"] = books.read_game(GAME)
    clock.now += prophetx_live.WANTED_FOR_SEC + 1
    assert books._next_game() is None
    assert books.want([GAME])["games"]["123742"] is None


def test_the_thread_stops_on_its_event():
    books = reader()
    books.want([GAME])
    stop = threading.Event()
    original = books.read_game

    def read_once(game):
        stop.set()
        return original(game)
    books.read_game = read_once
    books.run(stop)
    assert books.status()["readsDone"] == 1 and books.want([GAME])["games"]["123742"]["prophetxEventId"] == PX_EVENT_ID


def test_validate_names_what_is_wrong():
    assert prophetx_live.validate_prophetx_live_request({"games": [GAME]}) == [GAME]
    assert "body must be" in prophetx_live.validate_prophetx_live_request({"games": "x"})
    assert "gameKey" in prophetx_live.validate_prophetx_live_request({"games": [{"line": {}}]})
    assert "line object" in prophetx_live.validate_prophetx_live_request({"games": [{"gameKey": "1"}]})
    too_many = {"games": [GAME] * (prophetx_live.MAX_GAMES_PER_ASK + 1)}
    assert "at most" in prophetx_live.validate_prophetx_live_request(too_many)


# ---- the route --------------------------------------------------------------------------

@pytest.fixture
def served(tmp_path):
    def serve(live):
        store = BetsStore(tmp_path / "bets.duckdb", 0)
        handler = service.make_handler(store, 0.0, [], runner_url="http://127.0.0.1:1", extra_hosts=(), prophetx_live=live)
        server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
        threading.Thread(target=server.serve_forever, daemon=True).start()
        servers.append((server, store))
        return f"http://127.0.0.1:{server.server_address[1]}"
    servers = []
    yield serve
    for server, store in servers:
        server.shutdown()
        store.close()


def post(base, body):
    request = urllib.request.Request(f"{base}/prophetx_live.json", data=json.dumps(body).encode(), method="POST",
                                     headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=5) as reply:
            return reply.status, json.loads(reply.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


def test_the_route_answers_the_reads_and_503s_without_the_live_key(served):
    status, body = post(served(reader()), {"games": [GAME]})
    assert status == 200 and body["games"] == {"123742": None} and body["reader"]["accessKey"] == "live1234"
    status, body = post(served(None), {"games": [GAME]})
    assert status == 503 and "PROPHETX_LIVE_ACCESS_KEY" in body["error"]
    assert post(served(reader()), {"games": "x"})[0] == 400
