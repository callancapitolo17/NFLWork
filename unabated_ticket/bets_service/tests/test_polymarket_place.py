"""The Ticket's Place button at Polymarket US (polymarket_place.py) against a
scripted order client and a scripted runner quote: the request shape, the
refusals before anything is sent (started, all-in over 0.1c worse, other side
held, a stale Polymarket US read, an earlier Place the stake has not seen, a
repeat click), the contract count at the all-in price, the YES price a NO
buy sends, what each answer reads as; plus the order client's signing and
body (polymarket_orders.py). No network."""
import base64
import json
import threading
import urllib.error
import urllib.request
from datetime import datetime, timedelta, timezone

import pytest
from cryptography.hazmat.primitives.asymmetric import ed25519
from cryptography.hazmat.primitives.serialization import Encoding, PublicFormat

from http.server import ThreadingHTTPServer

from unabated_ticket.bets_service import manual_orders, polymarket_orders, service
from unabated_ticket.bets_service.polymarket_orders import (PolymarketError, PolymarketOrderClient, ioc_order_body,
                                                            order_fill_totals)
from unabated_ticket.bets_service.polymarket_place import (ManualPolymarketDesk, held_id_of,
                                                           validate_polymarket_place_request)
from unabated_ticket.bets_service.store import BetsStore

NOW = datetime(2026, 10, 22, 15, 0, tzinfo=timezone.utc)
EVENT_ID = 123742
SLUG = "asc-cfb-duq-wst-2026-10-24-pos-35pt5"
ORDER_ID = "ord-1"
GATEWAY = "https://gateway.polymarket.us"
BODY = {"clickId": "abc12345-0001", "eventId": EVENT_ID, "periodTypeId": 1, "betTypeId": 2, "sideIndex": 1,
        "points": -35.5, "stake": 50, "ticketAllIn": 0.65, "polymarketReadAt": "2026-10-22T14:59:40Z",
        "kalshiReadAt": "2026-10-22T14:59:50Z", "novigReadAt": None}
LINE = {"league": "cfb", "betType": "spread", "period": "FG", "side": "home", "points": -35.5,
        "awayTeam": "Duquesne", "homeTeam": "Washington State", "awayTeamId": "927", "homeTeamId": "717"}
# NO at 63c: 63c + 0.06 x 0.63 x 0.37 = 0.643986 all-in.
QUOTE = {"ok": True, "slug": SLUG, "side": "no", "ask": 0.63, "allIn": 0.643986, "depth": 500, "eventId": EVENT_ID,
         "eventStart": "2026-10-24T19:30:00.000Z", "label": "Washington State -35.5 FG", "polymarketGatewayUrl": GATEWAY,
         "line": LINE}


def final_order(quantity, filled, state="ORDER_STATE_FILLED", average="0.3700", commission="0.44"):
    return {"id": ORDER_ID, "marketSlug": SLUG, "state": state, "quantity": quantity, "cumQuantity": filled,
            "leavesQuantity": quantity - filled, "avgPx": {"value": average, "currency": "USD"},
            "commissionNotionalTotalCollected": {"value": commission, "currency": "USD"}}


class ScriptedPolymarket:
    """PolymarketOrderClient's order calls, scripted. By default the IOC is
    taken and reads FILLED for every contract asked for."""
    gateway_url = GATEWAY

    def __init__(self, send=None, order_reads=None, raises=False):
        self.calls = []
        self.cancels = []
        self.send = send
        self.order_reads = list(order_reads or [])
        self.raises = raises

    def place_ioc_order(self, slug, contract_side, side_price, quantity):
        self.calls.append((slug, contract_side, side_price, quantity))
        self.last_quantity = quantity
        if self.raises:
            raise OSError("timed out")
        return self.send or (200, {"id": ORDER_ID, "executions": []})

    def get_order(self, order_id):
        if self.order_reads:
            read = self.order_reads.pop(0)
            if isinstance(read, Exception):
                raise read
            return read
        return 200, {"order": final_order(self.last_quantity, self.last_quantity)}

    def cancel_order(self, order_id, slug):
        self.cancels.append((order_id, slug))
        return 200, {}


@pytest.fixture
def store(tmp_path):
    held = BetsStore(tmp_path / "bets.duckdb", 0)
    held.log_source_run("polymarket_us", NOW - timedelta(seconds=40), NOW - timedelta(seconds=30), True, None, 3)
    yield held
    held.close()


def desk_with(store, client, quote=None, now=NOW):
    quotes = []

    def reader(line):
        quotes.append(line)
        return dict(QUOTE if quote is None else quote)

    desk = ManualPolymarketDesk(client, store, reader, now=lambda: now, sleep=lambda _sec: None)
    desk.quotes = quotes
    return desk


def request(**changes):
    return validate_polymarket_place_request({**BODY, **changes})


def test_validate_names_what_is_wrong():
    assert validate_polymarket_place_request({**BODY, "clickId": "x"}).startswith("clickId")
    assert validate_polymarket_place_request({**BODY, "betTypeId": 4}).startswith("betTypeId")
    assert validate_polymarket_place_request({**BODY, "betTypeId": 1})["points"] is None
    assert validate_polymarket_place_request({**BODY, "points": None}).startswith("points")
    assert validate_polymarket_place_request({**BODY, "polymarketReadAt": None}).startswith("polymarketReadAt")
    assert validate_polymarket_place_request({**BODY, "kalshiReadAt": "yesterday"}).startswith("kalshiReadAt")
    assert validate_polymarket_place_request({**BODY, "ticketAllIn": 1.0}).startswith("ticketAllIn")
    assert validate_polymarket_place_request({**BODY, "stake": 0}).startswith("stake")
    assert validate_polymarket_place_request([BODY]).startswith("body")
    accepted = validate_polymarket_place_request(BODY)
    assert accepted["novig_read_at"] is None and accepted["kalshi_read_at"] is not None


def test_a_click_sends_one_ioc_buy_counted_at_the_all_in_price(store):
    client = ScriptedPolymarket()
    desk = desk_with(store, client)
    result = desk.place(request())
    assert desk.quotes == [{"eventId": EVENT_ID, "periodTypeId": 1, "betTypeId": 2, "sideIndex": 1, "points": -35.5}]
    # floor(50 / 0.643986) = 77 contracts, at NO's own 63c.
    assert client.calls == [(SLUG, "no", 0.63, 77)]
    assert result["status"] == "filled"
    # 77 NO contracts at YES avgPx 0.37 cost 77 x 0.63 = 48.51, plus 0.44 commission.
    assert result["cost"] == pytest.approx(48.95)
    assert result["message"] == "Bought 77 NO contracts at 64.4¢ all-in for $48.95."


def test_the_row_is_held_before_the_send_and_the_same_click_is_never_sent_twice(store):
    client = ScriptedPolymarket()
    desk_with(store, client).place(request())
    again = desk_with(store, client).place(request(polymarketReadAt="2026-10-22T15:30:00Z"))
    assert again["status"] == "refused" and "already sent" in again["message"]
    assert len(client.calls) == 1


@pytest.mark.parametrize("quote_changes, words", [
    ({"eventStart": "2026-10-22T14:00:00Z"}, "has started"),
    ({"eventStart": None}, "start time is unknown"),
    ({"allIn": 0.652}, "more than 0.1"),
    ({"ok": False, "reason": "no_market", "detail": "0 of Polymarket US's 3 events are X at its start"}, "0 of Polymarket"),
    ({"polymarketGatewayUrl": "https://gateway.test"}, "this desk trades"),
    ({"ask": None}, "is not a price"),
    ({"allIn": 0.6}, "is not a price"),
    ({"side": "maybe"}, "the market and its side"),
    ({"slug": ""}, "the market and its side"),
])
def test_refusals_send_nothing(store, quote_changes, words):
    client = ScriptedPolymarket()
    result = desk_with(store, client, {**QUOTE, **quote_changes}).place(request())
    assert result["status"] == "refused" and words in result["message"], result["message"]
    assert client.calls == []


def test_a_price_within_a_tenth_of_a_cent_worse_or_any_better_places(store):
    client = ScriptedPolymarket()
    assert desk_with(store, client, {**QUOTE, "allIn": 0.6509}).place(request())["status"] == "filled"
    assert desk_with(store, client, {**QUOTE, "ask": 0.5, "allIn": 0.515}).place(
        request(clickId="abc12345-0002", polymarketReadAt=(NOW + timedelta(seconds=1)).isoformat()))["status"] == "filled"


def test_holding_the_other_side_refuses(store):
    store.upsert_bets([{
        "id": held_id_of(SLUG, "yes"), "venue": "polymarket_us", "source": "polymarket_us_api", "status": "open",
        "placedAt": "2026-10-21T12:00:00Z", "sourceFetchedAt": NOW.isoformat(),
    }], NOW)
    client = ScriptedPolymarket()
    result = desk_with(store, client).place(request())
    assert result["status"] == "refused" and "you hold the YES side" in result["message"]
    assert client.calls == []


def test_a_failing_or_old_polymarket_read_refuses(store):
    client = ScriptedPolymarket()
    old = desk_with(store, client, now=NOW + timedelta(minutes=6)).place(request())
    assert old["status"] == "refused" and "over 5 minutes old" in old["message"]
    store.log_source_run("polymarket_us", NOW, NOW, False, "HTTP 500", 0)
    failing = desk_with(store, client).place(request())
    assert failing["status"] == "refused" and "read is failing" in failing["message"]
    assert client.calls == []


def test_an_earlier_place_the_stake_has_not_seen_refuses_at_any_venue(store):
    client = ScriptedPolymarket()
    desk_with(store, client).place(request())
    later = desk_with(store, client, now=NOW + timedelta(seconds=5)).place(request(clickId="abc12345-0002"))
    assert later["status"] == "refused" and "not on the Bets tab yet" in later["message"]
    store.insert_manual_order({
        "client_order_id": "uman-k1", "placed_at": NOW, "venue": "kalshi", "event_id": EVENT_ID + 1, "ticker": "T",
        "contract_side": "yes", "label": None, "ticket_all_in": 0.5, "ask_price": 0.5, "all_in_price": 0.51,
        "stake": 10, "count_requested": 19, "status": "filled", "updated_at": NOW})
    kalshi_unseen = desk_with(store, client).place(request(clickId="abc12345-0003", eventId=EVENT_ID + 1,
                                                            kalshiReadAt="2026-10-22T14:59:00Z"))
    assert kalshi_unseen["status"] == "refused" and "not on the Bets tab yet" in kalshi_unseen["message"]
    assert len(client.calls) == 1


def test_a_rejection_places_nothing_and_no_answer_is_unknown(store):
    rejected = desk_with(store, ScriptedPolymarket(send=(400, {"message": "insufficient buying power"}))).place(request())
    assert rejected["status"] == "rejected" and "insufficient buying power" in rejected["message"]
    lost = desk_with(store, ScriptedPolymarket(raises=True)).place(request(clickId="abc12345-0002"))
    assert lost["status"] == "unknown" and "did not answer" in lost["message"]
    # The lost order may hold contracts: the game's next click waits for a read after it.
    blocked = desk_with(store, ScriptedPolymarket()).place(request(clickId="abc12345-0003"))
    assert blocked["status"] == "refused" and "not on the Bets tab yet" in blocked["message"]
    odd = desk_with(store, ScriptedPolymarket(send=(502, "bad gateway"))).place(request(clickId="abc12345-0004", eventId=EVENT_ID + 1))
    assert odd["status"] == "unknown" and "HTTP 502" in odd["message"]


def test_partial_unfilled_and_resting_orders(store):
    partial = ScriptedPolymarket(order_reads=[(200, {"order": final_order(77, 30, state="ORDER_STATE_CANCELED")})])
    result = desk_with(store, partial).place(request())
    assert result["status"] == "partial" and result["filled"] == 30
    unfilled = ScriptedPolymarket(order_reads=[(200, {"order": final_order(77, 0, state="ORDER_STATE_CANCELED")})])
    assert desk_with(store, unfilled).place(request(clickId="abc12345-0002", eventId=EVENT_ID + 2))["status"] == "unfilled"
    resting = ScriptedPolymarket(order_reads=[(200, {"order": final_order(77, 0, state="ORDER_STATE_NEW")})] * 10)
    rest = desk_with(store, resting).place(request(clickId="abc12345-0003", eventId=EVENT_ID + 3))
    assert rest["status"] == "unknown" and resting.cancels == [(ORDER_ID, SLUG)]
    never = ScriptedPolymarket(order_reads=[OSError("reset")] * 10)
    assert desk_with(store, never).place(request(clickId="abc12345-0004", eventId=EVENT_ID + 4))["status"] == "unknown"


# ---- the order client ----

def fake_secret_key():
    private = ed25519.Ed25519PrivateKey.generate()
    seed = private.private_bytes_raw()
    public = private.public_key().public_bytes(Encoding.Raw, PublicFormat.Raw)
    return base64.b64encode(seed + public).decode(), private.public_key()


def test_the_client_signs_timestamp_method_and_path_without_the_query():
    secret, public = fake_secret_key()
    sent = []

    def transport(method, url, headers, body):
        sent.append((method, url, headers, body))
        return 200, b'{"id": "o1", "executions": []}'

    client = PolymarketOrderClient("key-1", secret, transport=transport, clock=lambda: 1_760_000_000.123)
    status, reply = client.place_ioc_order(SLUG, "no", 0.63, 77)
    assert (status, reply) == (200, {"id": "o1", "executions": []})
    method, url, headers, body = sent[0]
    assert (method, url) == ("POST", "https://api.polymarket.us/v1/orders")
    assert headers["X-PM-Access-Key"] == "key-1" and headers["X-PM-Timestamp"] == "1760000000123"
    public.verify(base64.b64decode(headers["X-PM-Signature"]), b"1760000000123POST/v1/orders")
    assert b'"price": {"value": "0.3700", "currency": "USD"}' in body
    assert b'"intent": "ORDER_INTENT_BUY_SHORT"' in body and b'"tif": "TIME_IN_FORCE_IMMEDIATE_OR_CANCEL"' in body


def test_a_bad_key_fails_at_start():
    with pytest.raises(PolymarketError):
        PolymarketOrderClient("key-1", "not base64!")


def test_order_body_sends_the_yes_price_and_refuses_nonsense():
    assert ioc_order_body(SLUG, "yes", 0.42, 10)["price"]["value"] == "0.4200"
    assert ioc_order_body(SLUG, "no", 0.835, 10)["price"]["value"] == "0.1650"
    assert ioc_order_body(SLUG, "yes", 0.42, 10)["intent"] == polymarket_orders.INTENT_BUY_LONG
    with pytest.raises(ValueError):
        ioc_order_body(SLUG, "yes", 0.42, 0)
    with pytest.raises(ValueError):
        ioc_order_body(SLUG, "yes", 1.2, 5)


def test_fill_totals_read_a_no_fill_off_the_yes_average():
    assert order_fill_totals(final_order(77, 77), "no") == (77.0, pytest.approx(48.95))
    assert order_fill_totals(final_order(77, 77, average="0.6300"), "yes") == (77.0, pytest.approx(48.95))
    assert order_fill_totals(final_order(77, 0), "yes") == (0.0, 0.0)
    assert order_fill_totals({"state": "ORDER_STATE_FILLED"}, "yes").startswith("no cumQuantity")
    assert order_fill_totals({**final_order(5, 5), "avgPx": None}, "yes") == (5.0, None)


# ---- the route and the desk's start ----

class StubDesk:
    def __init__(self):
        self.requests = []

    def place(self, placed):
        self.requests.append(placed)
        return {"status": "filled", "message": "Bought 77 NO contracts at 64.4¢ all-in for $48.95."}


def post_place(store, desk):
    options = {"polymarket_place_desk": desk} if desk is not None else {}
    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(store, 0.0, ["polymarket_us"], **options))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{server.server_address[1]}/place_polymarket.json"

    def send(body):
        sent = urllib.request.Request(url, data=json.dumps(body).encode(), method="POST",
                                      headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(sent, timeout=5) as response:
                return response.status, json.loads(response.read())
        except urllib.error.HTTPError as error:
            return error.code, json.loads(error.read())

    return server, send


def test_the_route_places_and_refuses_a_bad_body_or_no_desk(store):
    desk = StubDesk()
    server, send = post_place(store, desk)
    try:
        status, reply = send(BODY)
        assert status == 200 and reply["ok"] is True and reply["status"] == "filled"
        assert desk.requests == [validate_polymarket_place_request(BODY)]
        status, reply = send({**BODY, "betTypeId": 4})
        assert status == 400 and "betTypeId" in reply["error"]
    finally:
        server.shutdown()
        server.server_close()
    off, send_off = post_place(store, None)
    try:
        status, reply = send_off(BODY)
        assert status == 503 and "no Polymarket US key" in reply["error"]
    finally:
        off.shutdown()
        off.server_close()


class BalanceClient:
    gateway_url = GATEWAY

    def __init__(self, status=200):
        self.status = status

    def balance(self):
        return self.status, {"balances": []}


def test_the_service_arms_place_only_on_a_key_that_signs_with_the_read_in_place(store, monkeypatch):
    from unabated_ticket.bets_service import config
    monkeypatch.setattr(config, "POLYMARKET_US_KEY_ID", "")
    assert "no Polymarket US key" in service.open_polymarket_place_desk(store, ["polymarket_us"])[1]
    monkeypatch.setattr(config, "POLYMARKET_US_KEY_ID", "key-1")
    monkeypatch.setattr(config, "POLYMARKET_US_SECRET_KEY", "fake-secret")
    assert "does not read polymarket_us" in service.open_polymarket_place_desk(
        store, ["kalshi"], client_factory=lambda *_: BalanceClient())[1]
    refused = service.open_polymarket_place_desk(store, ["polymarket_us"], client_factory=lambda *_: BalanceClient(401))
    assert refused[0] is None and "HTTP 401" in refused[1]

    def unreadable(*_):
        raise PolymarketError("not a key")

    assert "not a key" in service.open_polymarket_place_desk(store, ["polymarket_us"], client_factory=unreadable)[1]
    armed, _ = service.open_polymarket_place_desk(store, ["polymarket_us"], client_factory=lambda *_: BalanceClient())
    assert isinstance(armed, ManualPolymarketDesk)
    opened = []

    class Opener:
        def open(self, url, timeout):
            opened.append((url, timeout))
            raise urllib.error.URLError("timed out")

    monkeypatch.setattr(manual_orders, "_RUNNER_OPENER", Opener())
    answer = armed._quote_reader({"eventId": EVENT_ID, "periodTypeId": 1, "betTypeId": 2, "sideIndex": 1, "points": -35.5})
    assert answer["reason"] == "runner"
    assert [(url.split("?")[0].endswith("/polymarket_quote.json"), timeout) for url, timeout in opened] == [(True, 30)]
