"""The Ticket's Place button on Kalshi (manual_orders.py) against a scripted
Kalshi and a scripted runner quote: the request shape, the refusals before
anything is sent (started, price over 0.1c worse, other side held, an earlier
fill the stake has not seen, a repeat click), the IOC body at the live ask, and
what each Kalshi answer reads as. No network."""
from datetime import datetime, timedelta, timezone

import pytest

from unabated_ticket.bets_service.manual_orders import (ManualKalshiDesk, contracts_detail, filled_bet_message,
                                                       unpriced_fill_message, validate_place_request)
from unabated_ticket.bets_service.store import BetsStore

NOW = datetime(2026, 10, 16, 21, 0, tzinfo=timezone.utc)
TICKER = "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36"
EVENT_ID = 123742
ALL_IN = 0.587157
BODY = {"clickId": "c1a2b3c4-0001", "eventId": EVENT_ID, "periodTypeId": 1, "betTypeId": 2, "sideIndex": 1,
        "points": -35.5, "stake": 100, "ticketAllIn": 0.5872, "kalshiReadAt": "2026-10-16T20:59:30Z"}
QUOTE = {"ok": True, "ticker": TICKER, "side": "yes", "ask": 0.57, "allIn": ALL_IN, "depth": 40,
         "eventId": EVENT_ID, "eventStart": "2026-10-17T19:30:00.000Z", "label": "Washington St -35.5 full game",
         "kalshiBaseUrl": "https://api.elections.kalshi.com/trade-api/v2"}


class ScriptedKalshi:
    base_url = "https://api.elections.kalshi.com/trade-api/v2"

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
        if method == "DELETE":
            return 200, {}
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
    assert validate_place_request({**BODY, "novigReadAt": "yesterday"}).startswith("novigReadAt")
    assert validate_place_request(BODY)["novig_read_at"] is None
    assert validate_place_request({**BODY, "points": None})["points"] is None


def test_a_click_sends_one_ioc_buy_at_the_live_ask_and_reads_what_filled(store):
    kalshi = ScriptedKalshi()
    desk = desk_with(store, kalshi)
    result = desk.place(request())
    assert desk.quotes == [{"eventId": EVENT_ID, "periodTypeId": 1, "betTypeId": 2, "sideIndex": 1, "points": -35.5, "live": False}]
    [(_, path, body)] = kalshi.posts()
    assert path == "/portfolio/events/orders"
    # floor(100 / 0.587157) = 170 contracts, YES bid at the 57c ask, nothing rests.
    assert body["client_order_id"] == "uman-c1a2b3c4-0001"
    assert (body["ticker"], body["side"], body["price"], body["count"], body["time_in_force"]) == (
        TICKER, "bid", "0.5700", "170.00", "immediate_or_cancel")
    assert "subaccount" not in body  # Cal's primary account, never the bot's
    assert result["status"] == "filled"
    assert result["cost"] == pytest.approx(170 * 0.57 + 2.92)
    assert result["message"] == ("Bet $99.82 on Washington St -35.5 full game at -142 to win $70.18 "
                                 "(170 contracts, 58.7¢ each with fee).")


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
    ({"kalshiBaseUrl": "https://demo-api.kalshi.co/trade-api/v2"}, "this desk trades"),
    ({"ask": None}, "is not a price"),
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


def test_the_bots_subaccount_holding_the_other_side_does_not_block(store):
    store.upsert_bets([{"id": f"kalshi:{TICKER}:no:s1", "venue": "kalshi", "status": "open"}], NOW)
    assert desk_with(store, ScriptedKalshi()).place(request())["status"] == "filled"


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


def test_an_add_waits_on_a_novig_place_the_novig_read_has_not_seen(store):
    """The stake counts Novig bets off the Novig read: a Novig Place after that
    read is missing from it, however fresh the Kalshi read is."""
    store.insert_manual_order({
        "client_order_id": "nman-00000000-0000-4000-8000-000000000001", "placed_at": NOW - timedelta(seconds=30),
        "venue": "novig", "event_id": EVENT_ID, "ticker": "01a0a05e-e864-7b73-9b80-9e502e4e98d5", "contract_side": "buy",
        "label": None, "ticket_all_in": 0.675, "ask_price": 0.675, "all_in_price": 0.675, "stake": 50,
        "count_requested": 7407, "status": "filled", "updated_at": NOW})
    kalshi = ScriptedKalshi()
    read_before = request(kalshiReadAt=NOW.isoformat(), novigReadAt=(NOW - timedelta(seconds=60)).isoformat())
    assert desk_with(store, kalshi).place(read_before)["status"] == "refused"
    never_read = request(clickId="c1a2b3c4-0020", kalshiReadAt=NOW.isoformat(), novigReadAt=None)
    assert desk_with(store, kalshi).place(never_read)["status"] == "refused"
    assert kalshi.calls == []
    read_after = request(clickId="c1a2b3c4-0021", kalshiReadAt=NOW.isoformat(), novigReadAt=NOW.isoformat())
    assert desk_with(store, kalshi).place(read_after)["status"] == "filled"


def test_kalshi_refusing_is_rejected_with_its_words(store):
    result = desk_with(store, ScriptedKalshi(post_status=400)).place(request())
    assert result["status"] == "rejected" and "insufficient balance" in result["message"]


def test_no_answer_or_no_closed_count_is_unknown(store):
    assert desk_with(store, ScriptedKalshi(raises=True)).place(request())["status"] == "unknown"
    assert desk_with(store, ScriptedKalshi(remaining=None)).place(
        request(clickId="c1a2b3c4-0004", kalshiReadAt="2026-10-16T21:30:00Z"))["status"] == "unknown"


def test_a_reply_still_open_gets_one_cancel(store):
    kalshi = ScriptedKalshi(remaining="30.00")
    result = desk_with(store, kalshi).place(request())
    assert result["status"] == "unknown" and "rest was cancelled" in result["message"]
    assert [call[1] for call in kalshi.calls if call[0] == "DELETE"] == [f"/portfolio/events/orders/k-9?market_ticker={TICKER}"]


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
        return {"status": "filled", "message": "Bet $99.82 on Washington St -35.5 at -142 to win $70.18."}


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
    assert status == 503 and "KALSHI_PLACE_API_KEY_ID" in reply["error"]



# ---- which key Place signs with ----
def test_place_is_off_without_its_own_trade_key(monkeypatch, store):
    monkeypatch.setattr(service.config, "KALSHI_PLACE_API_KEY_ID", None)
    monkeypatch.setattr(service.config, "KALSHI_PLACE_PRIVATE_KEY_PATH", "/data/kalshi_place.pem")
    monkeypatch.setattr(service.config, "KALSHI_API_KEY_ID", "read-only-id")
    desk, reason = service.open_kalshi_place_desk(store)
    assert desk is None and "KALSHI_PLACE_API_KEY_ID" in reason


# ---- a Place reads its venue before it answers ----
class CountingSource:
    poll_sec = 60

    def __init__(self, name):
        self.name = name
        self.fetches = 0

    def fetch(self):
        self.fetches += 1
        return []


class StatusDesk:
    def __init__(self, status):
        self.status = status

    def place(self, placed):
        return {"status": self.status, "message": self.status}


def serve_reading(store, desk, sources):
    reader = service.SourceReader(sources, store)
    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(
        store, 0.0, [source.name for source in sources], kalshi_place_desk=desk, source_reader=reader))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, f"http://127.0.0.1:{server.server_address[1]}"


@pytest.mark.parametrize("status, reads", [("filled", 1), ("partial", 1), ("unknown", 1),
                                           ("unfilled", 0), ("refused", 0), ("rejected", 0)])
def test_a_kalshi_place_that_may_hold_contracts_reads_kalshi_before_the_reply(store, status, reads):
    kalshi, novig = CountingSource("kalshi"), CountingSource("novig")
    server, url = serve_reading(store, StatusDesk(status), [kalshi, novig])
    try:
        assert post_place(url, BODY)[0] == 200
    finally:
        server.shutdown()
        server.server_close()
    assert (kalshi.fetches, novig.fetches) == (reads, 0)  # only the venue traded at
    if reads:
        assert store.source_status()["kalshi"]["ok"] is True  # the read the next stake counts


def test_read_now_reads_each_named_source_and_skips_one_not_registered(store):
    novig, place = CountingSource("novig"), CountingSource("novig_place")
    reader = service.SourceReader([novig, place], store)
    assert reader.read_now(service.PLACE_READ_SOURCES["novig"], 1.0) == []  # no novig_bot source: skipped
    assert (novig.fetches, place.fetches) == (1, 1)


def test_read_now_gives_up_on_a_source_whose_read_is_still_running(store):
    novig = CountingSource("novig")
    reader = service.SourceReader([novig], store)
    reader._locks["novig"].acquire()  # the poll lane is mid-read
    try:
        assert reader.read_now(("novig",), 0.05) == ["novig"]
    finally:
        reader._locks["novig"].release()
    assert novig.fetches == 0


def test_fill_message_speaks_in_the_bets_terms():
    # Cal's 2026-10-31 Polymarket US fill: Stonehill -14.5 was the market's NO side.
    message = filled_bet_message("Stonehill College -14.5 FG", 231.78, 749, contracts_detail(749, 749, 231.78))
    assert message == "Bet $231.78 on Stonehill College -14.5 at +223 to win $517.22 (749 contracts, 30.9¢ each with fee)"
    assert contracts_detail(60, 100, 18.0) == "60 of 100 contracts, 30.0¢ each with fee"
    assert filled_bet_message("Over 44.5 1H", 50.0, 100.0, "x").startswith("Bet $50.00 on Over 44.5 1H at +100 ")
    assert filled_bet_message(None, 60.0, 100.0, "x").startswith("Bet $60.00 on your line at -150 ")
    assert unpriced_fill_message("Kalshi", "A +1.5 FG", 170, 170, 0.587, ".") == (
        "Bought 170 contracts on A +1.5 at 58.7¢ quoted; Kalshi has not said yet what they cost.")
    # A fill that pays back no more than it cost has no odds to show.
    assert filled_bet_message("A +1.5 FG", 10.0, 10.0, "x") == "Bet $10.00 on A +1.5 to win $0.00 (x)"


# ---- Kalshi cash pots: a market on a short pot is topped up from the main pot first ----
class PottedKalshi(ScriptedKalshi):
    """ScriptedKalshi plus Kalshi's per-pot balance and pot-to-pot transfers."""

    def __init__(self, pot_cents=11386, transfer_status=200, transfer_reply=None, completes_after=1,
                 balance_statuses=(200,), **kwargs):
        super().__init__(**kwargs)
        self.pot_cents = pot_cents
        self.transfer_status = transfer_status
        self.transfer_reply = transfer_reply
        self.completes_after = completes_after
        self.balance_statuses = list(balance_statuses)
        self.status_reads = 0

    def request(self, method, path, body=None):
        if path.startswith("/portfolio/balance"):
            self.calls.append((method, path, body))
            status = self.balance_statuses.pop(0) if self.balance_statuses else 200
            return status, ({"balance": self.pot_cents} if status == 200 else {"error": {"code": "too_many_requests"}})
        if path == "/portfolio/intra_exchange_instance_transfer":
            self.calls.append((method, path, body))
            if self.transfer_reply is not None:
                return self.transfer_status, self.transfer_reply
            return self.transfer_status, {"transfer_id": "t-1"}
        if path.startswith("/portfolio/intra_exchange_instance_transfers/"):
            self.calls.append((method, path, body))
            self.status_reads += 1
            done = self.status_reads >= self.completes_after
            return 200, {"transfer": {"transfer_id": "t-1", "status": "complete" if done else "pending"}}
        return super().request(method, path, body)

    def transfers(self):
        return [call for call in self.calls if call[1] == "/portfolio/intra_exchange_instance_transfer"]

    def orders(self):
        return [call for call in self.posts() if call[1] == "/portfolio/events/orders"]


class FakeClock:
    def __init__(self):
        self.now = 0.0

    def sleep(self, seconds):
        self.now += seconds

    def monotonic(self):
        return self.now


def potted_desk(store, kalshi, quote_changes, clock=None):
    clock = clock or FakeClock()
    return ManualKalshiDesk(kalshi, store, lambda line: {**QUOTE, **quote_changes}, now=lambda: NOW,
                            sleep=clock.sleep, monotonic=clock.monotonic)


def test_the_shortfall_is_moved_in_centicents_and_the_order_waits_for_complete(store):
    kalshi = PottedKalshi(pot_cents=5000, completes_after=2)
    result = potted_desk(store, kalshi, {"exchangeIndex": 3}).place(request())
    [(_, _, transfer)] = kalshi.transfers()
    # cost bound = ceil(170 x 0.587157 = 99.8167) = $99.82 + $0.01 = $99.83; pot $50.00 -> move $49.83.
    assert transfer == {**transfer, "source": "event_contract", "destination": "event_contract",
                        "source_exchange_shard": 0, "destination_exchange_shard": 3, "amount": 498300}
    assert transfer["client_transfer_id"]
    assert kalshi.status_reads == 2
    order_index = kalshi.calls.index(kalshi.orders()[0])
    assert order_index > kalshi.calls.index(kalshi.transfers()[0]), "the order goes after the top-up"
    assert result["status"] == "filled"
    assert result["message"].endswith("(Moved $49.83 into the basketball/baseball/tennis pot first.)")


def test_a_pot_that_covers_the_cost_moves_nothing(store):
    kalshi = PottedKalshi(pot_cents=9983)
    assert potted_desk(store, kalshi, {"exchangeIndex": 3}).place(request())["status"] == "filled"
    assert kalshi.transfers() == []


@pytest.mark.parametrize("quote_changes", [{"exchangeIndex": 0}, {"exchangeIndex": None}, {}])
def test_a_main_pot_or_unknown_pot_market_reads_and_moves_nothing(store, quote_changes):
    kalshi = PottedKalshi(pot_cents=0)
    assert potted_desk(store, kalshi, quote_changes).place(request())["status"] == "filled"
    assert [call[1] for call in kalshi.calls if call[1].startswith("/portfolio/balance")] == []
    assert kalshi.transfers() == []


def test_a_rate_limited_balance_read_retries_once(store):
    kalshi = PottedKalshi(pot_cents=5000, balance_statuses=(429, 200))
    assert potted_desk(store, kalshi, {"exchangeIndex": 3}).place(request())["status"] == "filled"
    assert len(kalshi.transfers()) == 1, "the retry read the pot"


def test_a_pot_kalshi_will_not_report_orders_unchecked(store):
    unread = PottedKalshi(pot_cents=5000, balance_statuses=(429, 429))
    result = potted_desk(store, unread, {"exchangeIndex": 3}).place(request())
    assert unread.transfers() == [] and len(unread.orders()) == 1 and result["status"] == "filled"


def test_a_transfer_kalshi_refuses_sends_no_order(store):
    kalshi = PottedKalshi(pot_cents=5000, transfer_status=400,
                          transfer_reply={"error": {"code": "insufficient_balance", "message": "insufficient balance"}})
    result = potted_desk(store, kalshi, {"exchangeIndex": 3}).place(request())
    assert result["status"] == "refused" and kalshi.orders() == []
    assert "pot has $50.00" in result["message"] and "insufficient balance" in result["message"]
    assert desk_never_recorded(store)


def test_a_key_without_transfer_permission_says_so(store):
    kalshi = PottedKalshi(pot_cents=5000, transfer_status=403,
                          transfer_reply={"error": {"code": "insufficient_scope", "message": "insufficient scope"}})
    result = potted_desk(store, kalshi, {"exchangeIndex": 3}).place(request())
    assert result["status"] == "refused" and "transfer permission" in result["message"]
    assert kalshi.orders() == []


def test_a_transfer_that_does_not_complete_in_time_sends_no_order_and_the_click_can_retry(store):
    kalshi = PottedKalshi(pot_cents=5000, completes_after=10_000)
    result = potted_desk(store, kalshi, {"exchangeIndex": 3}).place(request())
    assert result["status"] == "refused" and "place again in a few seconds" in result["message"]
    assert kalshi.orders() == []
    # The same click, once the money landed, places (nothing was recorded for it).
    kalshi.pot_cents = 9983
    assert potted_desk(store, kalshi, {"exchangeIndex": 3}).place(request())["status"] == "filled"


def desk_never_recorded(store):
    return not store.manual_orders_maybe_filled_since(EVENT_ID, NOW - timedelta(days=1), venue="kalshi")


LIVE_QUOTE = {**QUOTE, "eventStart": "2026-10-16T19:30:00.000Z", "live": True}


def test_validate_live_fields():
    assert request()["live"] is False and request()["live_fair_read_at"] is None
    assert validate_place_request({**BODY, "live": "yes"}).startswith("live must be")
    assert validate_place_request({**BODY, "live": True}).startswith("liveFairReadAt")
    live = request(live=True, liveFairReadAt="2026-10-16T20:59:50Z")
    assert live["live"] is True and live["live_fair_read_at"] == datetime(2026, 10, 16, 20, 59, 50, tzinfo=timezone.utc)


def test_live_place_on_the_live_board_sends_one_order_and_asks_for_the_live_line(store):
    kalshi = ScriptedKalshi()
    desk = desk_with(store, kalshi, quote=LIVE_QUOTE)
    result = desk.place(request(live=True, liveFairReadAt="2026-10-16T20:59:50Z"))
    assert result["status"] == "filled"
    assert len(kalshi.posts()) == 1
    assert desk.quotes[0]["live"] is True


def test_live_place_refuses_a_stale_fair_and_a_line_off_the_live_board(store):
    kalshi = ScriptedKalshi()
    stale = desk_with(store, kalshi, quote=LIVE_QUOTE).place(request(live=True, liveFairReadAt="2026-10-16T20:58:00Z"))
    assert stale["status"] == "refused" and "last seen ready 120 s ago" in stale["message"]
    off_board = desk_with(store, kalshi, quote={**LIVE_QUOTE, "live": False}).place(
        request(live=True, liveFairReadAt="2026-10-16T20:59:50Z", clickId="c1a2b3c4-0002"))
    assert off_board["status"] == "refused" and "not on Unabated's live board" in off_board["message"]
    assert kalshi.posts() == []


def test_a_pregame_click_on_a_started_game_is_still_refused(store):
    kalshi = ScriptedKalshi()
    result = desk_with(store, kalshi, quote=LIVE_QUOTE).place(request())
    assert result["status"] == "refused" and "the game has started" in result["message"]
    assert kalshi.posts() == []
