"""ProphetX's live prices off ProphetX's market websocket (prophetx_live.py,
2026-11-06): one seed read per game, frames applied to it, reconnects and
unwanted games dropped, what Place reads, and the bets service's POST
/prophetx_live.json. Scripted client and socket, no network."""
import base64
import copy
import json
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import prophetx_live, service
from unabated_ticket.bets_service.prophetx_client import ProphetXError
from unabated_ticket.bets_service.store import BetsStore
from unabated_ticket.bets_service.tests.test_prophetx import EVENT, MARKETS, PX_EVENT_ID, SPREAD, offer, team_keys_reader

LINE = {"league": "nfl", "awayTeam": "Buffalo Bills", "homeTeam": "Kansas City Chiefs", "awayTeamId": "1", "homeTeamId": "2"}
GAME = {"gameKey": "123742", "eventStart": "2026-11-01T21:25:00.000Z", "strikeIds": ["e" * 32], "line": LINE}
PRIVATE_CHANNEL = "private-service=6-device_type=5-user=abc"
BEFORE_KICKOFF = 1_790_000_000.0  # 2026-09-21, well before EVENT's 2026-11-01 start


class Clock:
    def __init__(self):
        self.now = BEFORE_KICKOFF

    def __call__(self):
        return self.now


class FakeSocket:
    """Hands out queued frames one per recv (None when none is queued); records what is sent."""

    def __init__(self, clock):
        self.clock = clock
        self.queued = [pusher("pusher:connection_established", {"socket_id": "1.2"}),
                       pusher("pusher:signin_success", {}),
                       pusher("pusher_internal:subscription_succeeded", {}, PRIVATE_CHANNEL)]
        self.sent = []
        self.closed = False

    def send(self, text):
        self.sent.append(json.loads(text))

    def recv(self, timeout_sec):
        self.clock.now += 0.1
        return self.queued.pop(0) if self.queued else None

    def close(self):
        self.closed = True


class LiveClient:
    """The live key's REST calls, counted, with the register answer ProphetX gave in the recon."""

    def __init__(self, event=None, markets=None, channel_limit=500):
        self.access_key = "live1234abcd"
        self.blocked_until = 0.0
        self.calls = []
        self.registered = []
        self.event = EVENT if event is None else event
        self.markets_list = MARKETS if markets is None else markets
        self.channel_limit = channel_limit

    def websocket_config(self):
        self.calls.append("config")
        return {"app_id": "1", "key": "k", "cluster": "mt1", "service": "6", "ws_host": "ws.example"}

    def request(self, method, path, params=None, body=None, extra_headers=None):
        self.calls.append("register")
        ids = [int(i) for sub in body["subscriptions"] for i in sub["ids"]]
        self.registered.append(ids)
        channels = [{"channel_name": PRIVATE_CHANNEL, "auth": "p", "scope": None}]
        channels += [{"channel_name": f"private-broadcast-service=6-device_type=5-event={i}", "auth": f"a{i}",
                      "scope": {"event_id": i}} for i in ids]
        return 200, {"authenticated": {"auth": "s", "user_data": "{}"}, "authorized_channel": channels,
                     "channel_count": len(ids), "channel_limit": self.channel_limit}

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
        return copy.deepcopy(self.markets_list)


def pusher(event, data, channel=None):
    frame = {"event": event, "data": json.dumps(data)}
    if channel:
        frame["channel"] = channel
    return json.dumps(frame)


def selections_frame(info, event_id=PX_EVENT_ID):
    payload = {"sport_event_id": event_id, "market_id": info.get("id", 901), "info": info}
    return pusher("v2.market_selections", {"change_type": "market_selections", "op": "u",
                                           "payload": base64.b64encode(json.dumps(payload).encode()).decode()},
                  f"private-broadcast-service=6-device_type=5-event={event_id}")


def spread_info(away_side, home_side, sequence=10, status="active"):
    return {"id": 901, "name": "Fixed away +3.5", "type": "spread", "status": status, "sequence_number": sequence,
            "selections": [side for side in (away_side, home_side) if side is not None]}


class Live:
    """A reader with its socket opened by hand, driven one step at a time."""

    def __init__(self, client=None, team_keys=None):
        self.clock = Clock()
        self.client = client or LiveClient()
        self.socket = FakeSocket(self.clock)
        self.books = prophetx_live.ProphetXLiveBooks(self.client, team_keys or team_keys_reader(),
                                                    connect=lambda url, headers: self.socket, clock=self.clock)
        self.connection, self.socket_id, self.service = self.books._open()

    def follow(self):
        self.clock.now += prophetx_live.REGISTER_GAP_SEC
        self.books._follow_wanted(self.connection, self.socket_id, self.service)

    def deliver(self, *frames):
        """Frames arriving on the socket, read the way the pump reads them."""
        for text in frames:
            event, channel, data = prophetx_live.parse_frame(text)
            self.books._heard()
            self.books._handle(self.connection, event, channel, data)

    def ready(self, game=GAME):
        """Wanted, event found, subscribed and seeded."""
        self.books.want([game])
        assert self.books.seed_one()  # the event
        self.follow()
        assert self.books.seed_one()  # the seed read

    def read(self, game=GAME):
        return self.books.want([game])["games"][game["gameKey"]]


def side(read, bet_type, period, which, points):
    found = [line for line in read["lines"] if (line["betType"], line["period"], line["side"], line["points"]) == (bet_type, period, which, points)]
    assert len(found) == 1, found
    return found[0]


def test_a_game_is_read_once_then_served_off_the_stream():
    live = Live()
    assert live.read() is None  # asked, not read yet
    live.ready()
    read = live.read()
    assert read["prophetxEventId"] == PX_EVENT_ID and read["totalsOnly"] is False
    assert side(read, "spread", "FG", "away", 3.5) == {"betType": "spread", "period": "FG", "side": "away", "points": 3.5,
                                                      "strikeId": "c" * 32, "price": -110, "depth": 200.0}
    assert side(read, "moneyline", "FG", "home", None)["price"] == -160
    subscribes = [sent["data"]["channel"] for sent in live.socket.sent if sent["event"] == "pusher:subscribe"]
    assert subscribes == [PRIVATE_CHANNEL, f"private-broadcast-service=6-device_type=5-event={PX_EVENT_ID}"]
    assert live.client.calls == ["config", "register", "strikes", "sport_events", "register", "markets"]
    for _ in range(20):
        live.read()
        live.books.seed_one()
    assert live.client.calls.count("markets") == 1, "no ProphetX read per ask"


def test_a_frame_replaces_its_strike_pairs_book_and_leaves_the_rest():
    live = Live()
    live.ready()
    live.deliver(selections_frame(spread_info([offer("c" * 32, -102, 75, 12, "Buffalo Bills", 3.5, "+3.5")],
                                              [offer("d" * 32, -118, 40, 11, "Kansas City Chiefs", -3.5, "-3.5")])))
    read = live.read()
    assert side(read, "spread", "FG", "away", 3.5)["price"] == -102 and side(read, "spread", "FG", "away", 3.5)["depth"] == 75
    assert side(read, "spread", "FG", "home", -3.5)["price"] == -118
    assert side(read, "spread", "FG", "away", 2.5)["price"] == 120, "the alt number is another pair"
    assert live.books.status()["frames"] == 1


def test_a_side_the_frame_leaves_out_or_an_inactive_strike_has_nothing_resting():
    live = Live()
    live.ready()
    live.deliver(selections_frame(spread_info([offer("c" * 32, -102, 75, 12, "Buffalo Bills", 3.5, "+3.5")],
                                              [{"price": None, "quantity": 0}], sequence=10)))
    read = live.read()
    assert side(read, "spread", "FG", "home", -3.5)["price"] is None and side(read, "spread", "FG", "home", -3.5)["depth"] == 0.0
    live.deliver(selections_frame(spread_info([offer("c" * 32, -102, 75, 12, "Buffalo Bills", 3.5, "+3.5")],
                                              [offer("d" * 32, -118, 40, 11, "Kansas City Chiefs", -3.5, "-3.5")],
                                              sequence=11, status="suspended")))
    read = live.read()
    assert side(read, "spread", "FG", "away", 3.5)["price"] is None and side(read, "spread", "FG", "home", -3.5)["price"] is None


def test_a_frame_older_than_the_last_applied_is_skipped():
    live = Live()
    live.ready()
    live.deliver(selections_frame(spread_info([offer("c" * 32, -102, 75, 12, "Buffalo Bills", 3.5, "+3.5")],
                                              SPREAD["selections"][1], sequence=20)),
                 selections_frame(spread_info([offer("c" * 32, -130, 75, 12, "Buffalo Bills", 3.5, "+3.5")],
                                              SPREAD["selections"][1], sequence=19)))
    assert side(live.read(), "spread", "FG", "away", 3.5)["price"] == -102


def test_frames_between_the_subscribe_and_the_seed_apply_only_when_newer_than_the_seed():
    seeded = copy.deepcopy(MARKETS)
    for market in seeded:
        for sub in [market] + market.get("market_strikes", []):
            for ladder in sub["selections"]:
                for entry in ladder:
                    entry["updated_at"] = 100
    live = Live(LiveClient(markets=seeded))
    live.books.want([GAME])
    live.books.seed_one()
    live.follow()
    older = spread_info([offer("c" * 32, -150, 75, 12, "Buffalo Bills", 3.5, "+3.5")], SPREAD["selections"][1], sequence=90)
    newer_alt = {"id": 902, "status": "active", "sequence_number": 150, "selections": [
        [offer("e" * 32, 135, 20, 12, "Buffalo Bills", 2.5, "+2.5")], [offer("f" * 32, -160, 20, 11, "Kansas City Chiefs", -2.5, "-2.5")]]}
    live.deliver(selections_frame(older), selections_frame(newer_alt))
    live.books.seed_one()
    read = live.read()
    assert side(read, "spread", "FG", "away", 3.5)["price"] == -110, "the seed already shows the older frame"
    assert side(read, "spread", "FG", "away", 2.5)["price"] == 135


def test_a_silent_stream_reads_as_down_and_a_reconnect_seeds_every_game_again():
    live = Live()
    live.ready()
    live.clock.now += prophetx_live.SILENT_RECONNECT_SEC + 1
    assert live.read()["error"] == prophetx_live.STREAM_DOWN_WORDS
    assert live.books.streamed_game(GAME["gameKey"]) is None
    assert not live.books.seed_one() or live.client.calls[-1] != "markets", "nothing is seeded while the stream is down"
    live.books._drop(live.connection, "TimeoutError: silent")
    live.socket = FakeSocket(live.clock)
    live.connection, live.socket_id, live.service = live.books._open()
    assert live.read() is None, "a book from the old socket is never served"
    live.follow()
    live.books.seed_one()
    assert live.read()["prophetxEventId"] == PX_EVENT_ID and live.client.calls.count("markets") == 2


def test_a_game_no_longer_asked_for_is_dropped_and_unsubscribed():
    live = Live()
    live.ready()
    live.clock.now += prophetx_live.WANTED_FOR_SEC + 1
    assert not live.books.seed_one()
    live.follow()
    assert live.socket.sent[-1] == {"event": "pusher:unsubscribe",
                                    "data": {"channel": f"private-broadcast-service=6-device_type=5-event={PX_EVENT_ID}"}}
    assert live.books.status()["subscribed"] == 0 and live.read() is None


def test_a_frame_for_a_strike_added_after_the_seed_asks_for_a_new_seed():
    live = Live()
    live.ready()
    added = {"id": 950, "status": "active", "sequence_number": 5, "selections": [
        [offer("y" * 32, 110, 20, 12, "Buffalo Bills", 1.5, "+1.5")], [offer("z" * 32, -130, 20, 11, "Kansas City Chiefs", -1.5, "-1.5")]]}
    live.deliver(selections_frame(added))
    assert live.books.status()["unknownStrikeFrames"] == 1
    assert not live.books.seed_one(), "not before NEW_STRIKE_RESEED_SEC"
    live.clock.now += prophetx_live.NEW_STRIKE_RESEED_SEC + 1
    live.books._heard()  # the private channel's health_check every 5 s
    live.books.seed_one()
    assert live.client.calls.count("markets") == 2


def test_teams_the_other_way_round_keep_totals_only_and_other_teams_or_a_started_game_read_as_an_error():
    swapped = {**EVENT, "competitors": [{"id": 11, "name": "Kansas City Chiefs", "side": "away"},
                                        {"id": 12, "name": "Buffalo Bills", "side": "home"}]}
    live = Live(LiveClient(event=swapped))
    live.ready()
    assert live.read()["totalsOnly"] is True
    other = {**EVENT, "competitors": [{"id": 11, "name": "Kansas City Chiefs", "side": "home"},
                                      {"id": 13, "name": "Miami Dolphins", "side": "away"}]}
    live = Live(LiveClient(event=other))
    live.books.want([GAME])
    live.books.seed_one()
    assert "does not read as the board's" in live.read()["error"]
    live = Live(LiveClient(event={**EVENT, "status": "live"}))
    live.books.want([GAME])
    live.books.seed_one()
    assert live.read()["error"] == f"ProphetX game {PX_EVENT_ID} is live"
    live = Live()
    live.ready()
    live.clock.now = 1_794_000_000.0  # 2026-11-07, past the 2026-11-01 start
    live.books._heard()
    assert live.read()["error"] == f"ProphetX game {PX_EVENT_ID} has started"


def test_a_failed_seed_finds_the_event_again():
    live = Live(LiveClient(markets=ProphetXError("GET /v4/mm/get_markets: HTTP 429")))
    live.books.want([GAME])
    live.books.seed_one()
    live.follow()
    live.books.seed_one()
    assert "could not be seeded" in live.books.status()["lastError"]
    live.books.seed_one()
    assert live.client.calls[-2:] == ["strikes", "sport_events"]


def test_place_reads_the_streamed_game_only_while_the_book_is_current():
    live = Live()
    assert live.books.streamed_game(GAME["gameKey"]) is None
    live.ready()
    event, markets = live.books.streamed_game(GAME["gameKey"])
    assert event["event_id"] == PX_EVENT_ID and markets[0]["id"] == MARKETS[0]["id"]
    markets[0]["selections"] = []
    assert live.books.streamed_game(GAME["gameKey"])[1][0]["selections"], "Place gets a copy"


def test_a_socket_near_its_channel_limit_is_replaced():
    live = Live(LiveClient(channel_limit=prophetx_live.CHANNEL_HEADROOM))
    live.books.want([GAME])
    live.books.seed_one()
    with pytest.raises(ConnectionError, match="a new socket starts clean"):
        live.follow()


def test_the_thread_stops_on_its_event():
    clock = Clock()
    client = LiveClient()
    books = prophetx_live.ProphetXLiveBooks(client, team_keys_reader(), connect=lambda url, headers: FakeSocket(clock),
                                           clock=clock, read_gap_sec=0.0)
    stop = threading.Event()
    stop.set()
    books.run(stop)
    assert books.status()["streamLive"] is False


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
    status, body = post(served(Live().books), {"games": [GAME]})
    assert status == 200 and body["games"] == {"123742": None} and body["reader"]["accessKey"] == "live1234"
    status, body = post(served(None), {"games": [GAME]})
    assert status == 503 and "PROPHETX_LIVE_ACCESS_KEY" in body["error"]
    assert post(served(Live().books), {"games": "x"})[0] == 400
