"""The Ticket's Place button on Kalshi (manual_orders.py) against a scripted
Kalshi and a scripted runner quote: the request shape, the refusals before
anything is sent (started, price over 0.1c worse, other side held, an earlier
fill the stake has not seen, a repeat click), the IOC body at the live ask, and
what each Kalshi answer reads as. No network."""
from datetime import datetime, timedelta, timezone

import pytest

from unabated_ticket.bets_service.manual_orders import ManualKalshiDesk, validate_place_request
from unabated_ticket.bets_service.store import BetsStore

NOW = datetime(2026, 10, 16, 21, 0, tzinfo=timezone.utc)
TICKER = "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36"
EVENT_ID = 123742
ALL_IN = 0.587157
BODY = {"clickId": "c1a2b3c4-0001", "eventId": EVENT_ID, "periodTypeId": 1, "betTypeId": 2, "sideIndex": 1,
        "points": -35.5, "stake": 100, "ticketAllIn": 0.5872, "kalshiReadAt": "2026-10-16T20:59:30Z"}
QUOTE = {"ok": True, "ticker": TICKER, "side": "yes", "ask": 0.57, "allIn": ALL_IN, "depth": 40,
         "eventId": EVENT_ID, "eventStart": "2026-10-17T19:30:00.000Z", "label": "Washington St -35.5 full game"}


class ScriptedKalshi:
    def __init__(self, post_status=201, filled="170.00", remaining="0.00", fills_count="170", raises=False):
        self.calls = []
        self.post_status = post_status
        self.filled = filled
        self.remaining = remaining
        self.fills_count = fills_count
        self.raises = raises

    def request(self, method, path, body=None):
        self.calls.append((method, path, body))
        if method == "POST":
            if self.raises:
                raise OSError("timed out")
            if self.post_status != 201:
                return self.post_status, {"error": {"code": "insufficient_balance", "message": "insufficient balance"}}
            order = {"order_id": "k-9", "fill_count": self.filled}
            if self.remaining is not None:
                order["remaining_count"] = self.remaining
            return 201, order
        if path.startswith("/portfolio/fills"):
            return 200, {"fills": [{"count_fp": self.fills_count, "yes_price_dollars": "0.5700",
                                    "no_price_dollars": "0.4300", "fee_cost": "2.92"}], "cursor": ""}
        raise AssertionError(f"unexpected {method} {path}")

    def posts(self):
        return [call for call in self.calls if call[0] == "POST"]


@pytest.fixture
def store(tmp_path):
    held = BetsStore(tmp_path / "bets.duckdb", 0)
    yield held
    held.close()


def desk_with(store, kalshi, quote=None, now=NOW):
    quotes = []

    def reader(line):
        quotes.append(line)
        return dict(QUOTE if quote is None else quote)

    desk = ManualKalshiDesk(kalshi, store, reader, now=lambda: now)
    desk.quotes = quotes
    return desk


def request(**changes):
    return validate_place_request({**BODY, **changes})


def test_validate_names_what_is_wrong():
    assert validate_place_request({**BODY, "clickId": "x"}).startswith("clickId")
    assert validate_place_request({**BODY, "sideIndex": 2}).startswith("sideIndex")
    assert validate_place_request({**BODY, "stake": 0}).startswith("stake")
    assert validate_place_request({**BODY, "ticketAllIn": 1.2}).startswith("ticketAllIn")
    assert validate_place_request({**BODY, "kalshiReadAt": None}).startswith("kalshiReadAt")
    assert validate_place_request({**BODY, "points": None})["points"] is None


def test_a_click_sends_one_ioc_buy_at_the_live_ask_and_reads_what_filled(store):
    kalshi = ScriptedKalshi()
    desk = desk_with(store, kalshi)
    result = desk.place(request())
    assert desk.quotes == [{"eventId": EVENT_ID, "periodTypeId": 1, "betTypeId": 2, "sideIndex": 1, "points": -35.5}]
    [(_, path, body)] = kalshi.posts()
    assert path == "/portfolio/events/orders"
    # floor(100 / 0.587157) = 170 contracts, YES bid at the 57c ask, nothing rests.
    assert body["client_order_id"] == "uman-c1a2b3c4-0001"
    assert (body["ticker"], body["side"], body["price"], body["count"], body["time_in_force"]) == (
        TICKER, "bid", "0.5700", "170.00", "immediate_or_cancel")
    assert result["status"] == "filled"
    assert result["cost"] == pytest.approx(170 * 0.57 + 2.92)
    assert result["message"].startswith("Bought 170 contracts at 58.7")


def test_the_same_click_is_never_sent_twice(store):
    kalshi = ScriptedKalshi()
    desk = desk_with(store, kalshi)
    desk.place(request())
    again = desk.place(request(kalshiReadAt="2026-10-16T21:30:00Z"))
    assert again["status"] == "refused" and "already sent" in again["message"]
    assert len(kalshi.posts()) == 1


@pytest.mark.parametrize("quote_changes, words", [
    ({"eventStart": "2026-10-16T20:00:00Z"}, "has started"),
    ({"allIn": 0.5884}, "more than 0.1"),
    ({"ok": False, "reason": "no_line", "detail": "the board holds 0 Kalshi lines"}, "0 Kalshi lines"),
])
def test_refusals_send_nothing(store, quote_changes, words):
    kalshi = ScriptedKalshi()
    result = desk_with(store, kalshi, {**QUOTE, **quote_changes}).place(request())
    assert result["status"] == "refused" and words in result["message"]
    assert kalshi.posts() == []


def test_a_price_within_a_tenth_of_a_cent_worse_or_any_better_places(store):
    kalshi = ScriptedKalshi()
    assert desk_with(store, kalshi, {**QUOTE, "allIn": 0.5881}).place(request(clickId="c1a2b3c4-0002"))["status"] == "filled"
    better = desk_with(store, kalshi, {**QUOTE, "allIn": 0.55}).place(request(clickId="c1a2b3c4-0003",
                                                                              kalshiReadAt="2026-10-16T21:30:00Z"))
    assert better["status"] == "partial" and better["count"] == 181  # floor(100 / 0.55); the script fills 170


def test_the_other_side_held_is_refused(store):
    store.upsert_bets([{"id": f"kalshi:{TICKER}:no", "venue": "kalshi", "status": "open"}], NOW)
    kalshi = ScriptedKalshi()
    result = desk_with(store, kalshi).place(request())
    assert result["status"] == "refused" and "would sell it" in result["message"]
    assert kalshi.posts() == []


def test_an_add_waits_until_the_last_fill_on_the_game_was_read(store):
    kalshi = ScriptedKalshi()
    desk_with(store, kalshi).place(request())
    # The stake was sized on a Kalshi read from before that fill: refused.
    stale = desk_with(store, kalshi, now=NOW + timedelta(seconds=20)).place(request(clickId="c1a2b3c4-0009"))
    assert stale["status"] == "refused" and "not on the Bets tab yet" in stale["message"]
    # A read after it: placed.
    fresh = desk_with(store, kalshi, now=NOW + timedelta(seconds=90)).place(
        request(clickId="c1a2b3c4-0010", kalshiReadAt=(NOW + timedelta(seconds=60)).isoformat()))
    assert fresh["status"] == "filled"


def test_kalshi_refusing_is_rejected_with_its_words(store):
    result = desk_with(store, ScriptedKalshi(post_status=400)).place(request())
    assert result["status"] == "rejected" and "insufficient balance" in result["message"]


def test_no_answer_or_no_closed_count_is_unknown(store):
    assert desk_with(store, ScriptedKalshi(raises=True)).place(request())["status"] == "unknown"
    assert desk_with(store, ScriptedKalshi(remaining=None)).place(
        request(clickId="c1a2b3c4-0004", kalshiReadAt="2026-10-16T21:30:00Z"))["status"] == "unknown"


def test_partial_and_unfilled(store):
    partial = desk_with(store, ScriptedKalshi(filled="40.00", fills_count="40")).place(request())
    assert partial["status"] == "partial" and "40 of 170" in partial["message"]
    nothing = desk_with(store, ScriptedKalshi(filled="0.00")).place(
        request(clickId="c1a2b3c4-0005", kalshiReadAt="2026-10-16T21:30:00Z"))
    assert nothing["status"] == "unfilled" and "Nothing rests" in nothing["message"]


def test_a_fills_feed_that_lags_leaves_the_cost_unknown(store):
    result = desk_with(store, ScriptedKalshi(fills_count="100")).place(request())
    assert result["status"] == "filled" and result["cost"] is None


# ---- the route: POST /place_kalshi.json ----
import json  # noqa: E402
import threading  # noqa: E402
import urllib.error  # noqa: E402
import urllib.request  # noqa: E402
from http.server import ThreadingHTTPServer  # noqa: E402

from unabated_ticket.bets_service import service  # noqa: E402

TAILNET_HOST = "mac.tail1234.ts.net"


class StubDesk:
    def __init__(self):
        self.requests = []

    def place(self, placed):
        self.requests.append(placed)
        return {"status": "filled", "message": "Bought 170 contracts at 58.7¢ all-in for $99.82."}


@pytest.fixture
def serve_with(store):
    servers = []

    def start(desk):
        server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(
            store, 0.0, ["kalshi"], extra_hosts=(TAILNET_HOST,), kalshi_place_desk=desk))
        threading.Thread(target=server.serve_forever, daemon=True).start()
        servers.append(server)
        return f"http://127.0.0.1:{server.server_address[1]}"

    yield start
    for server in servers:
        server.shutdown()
        server.server_close()


def post_place(url, body, host=None):
    headers = {"Content-Type": "application/json", **({"Host": host} if host else {})}
    sent = urllib.request.Request(f"{url}/place_kalshi.json", data=json.dumps(body).encode(), method="POST", headers=headers)
    try:
        with urllib.request.urlopen(sent, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


def test_the_route_places_through_the_tailnet_name_and_refuses_a_bad_body_or_no_key(serve_with):
    desk = StubDesk()
    url = serve_with(desk)
    status, reply = post_place(url, BODY, host=TAILNET_HOST)
    assert status == 200 and reply["ok"] is True and reply["status"] == "filled"
    assert desk.requests == [validate_place_request(BODY)]
    status, reply = post_place(url, {**BODY, "stake": "100"})
    assert status == 400 and "stake" in reply["error"]
    status, reply = post_place(url, BODY, host="evil.example")
    assert status == 403
    status, reply = post_place(serve_with(None), BODY)
    assert status == 503 and "KALSHI_API_KEY_ID" in reply["error"]
