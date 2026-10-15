"""The bot's order desk (taker_orders.py) against a scripted Kalshi: the stop
file, the $500 ceiling on open exposure, the row written before the send,
limit-at-the-ask bodies, the cancel routed by market and verified on a 404,
the fills read, the unseen-fill rule, lost-order recovery, the untracked
resting sweep, the final re-read and settlement freeing the money. No network."""
import threading
from datetime import datetime, timedelta, timezone

import pytest

from unabated_ticket.bets_service import taker_orders
from unabated_ticket.bets_service.taker_orders import TakerDesk, order_body, validate_order_request
from unabated_ticket.bets_service.taker_store import TakerStore

NOW = datetime(2026, 10, 14, 12, 0, tzinfo=timezone.utc)
TICKER = "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36"
OTHER_TICKER = "KXNCAAFSPREAD-26SEP19DUQWSU-WSU42"
EVENT_ID = 123742
REQUEST_BODY = {"clientOrderId": "ubot-" + "a" * 32, "ticker": TICKER, "side": "yes", "ask": 0.57,
                "allIn": 0.587157, "count": 10, "rowKey": "420942308:ms105:si1:tid717:alt-35.5",
                "eventId": EVENT_ID, "heldIds": [], "eventStart": "2026-10-15T19:30:00.000Z",
                "decidedAt": "2026-10-14T11:59:50.000Z", "kalshiBaseUrl": "https://kalshi.test/trade-api/v2/"}
HELD_FIRST = f"kalshi:{TICKER}:yes"


class ScriptedKalshi:
    """Answers POST/DELETE/GET by path; records every call. `resting` is what
    GET /portfolio/orders?status=resting lists; `listed` what the ticker
    listing does; None for either makes that read fail."""

    base_url = "https://kalshi.test/trade-api/v2"

    def __init__(self, filled: float = 10, market_status: str = "active", post_status: int = 201,
                 delete_status: int = 200, resting=(), listed=(), exchange_index=None, fills_pages=None,
                 remaining_count="0.00"):
        self.calls = []
        self.filled = filled
        # What the create reply says is still open; None leaves the field out (an older reply shape).
        self.remaining_count = remaining_count
        # [(contracts, next cursor), ...] served page by page; None answers one page of `filled`.
        self.fills_pages = list(fills_pages) if fills_pages is not None else None
        self.market_status = market_status
        self.post_status = post_status
        self.delete_status = delete_status
        self.resting = list(resting) if resting is not None else None
        self.listed = list(listed) if listed is not None else None
        self.exchange_index = exchange_index

    def request(self, method, path, body=None):
        self.calls.append((method, path, body))
        if method == "POST":
            if self.post_status != 201:
                return self.post_status, {"error": "no"}
            # Create Order V2 answers flat: order_id, fill_count, remaining_count, ts_ms.
            order = {"order_id": "k-1", "fill_count": f"{float(self.filled):.2f}", "ts_ms": 1760000000000}
            if self.remaining_count is not None:
                order["remaining_count"] = self.remaining_count
            if self.exchange_index is not None:
                order["exchange_index"] = self.exchange_index
            return 201, order
        if method == "DELETE":
            return self.delete_status, {}
        if path.startswith("/portfolio/fills"):
            filled, cursor = (self.fills_pages.pop(0) if self.fills_pages else (self.filled, ""))
            fills = [{"count_fp": str(filled), "yes_price_dollars": "0.5700", "no_price_dollars": "0.4300",
                      "fee_cost": "0.17"}] if filled else []
            return 200, {"fills": fills, "cursor": cursor}
        if path.startswith("/portfolio/orders?status=resting"):
            return (500, "down") if self.resting is None else (200, {"orders": self.resting, "cursor": ""})
        if path.startswith("/portfolio/orders?ticker="):
            return (500, "down") if self.listed is None else (200, {"orders": self.listed, "cursor": ""})
        if path.startswith("/markets/"):
            return 200, {"market": {"status": self.market_status}}
        raise AssertionError(f"unexpected {method} {path}")

    def methods(self):
        return [call[0] for call in self.calls]


@pytest.fixture
def store(tmp_path):
    taker_store = TakerStore(tmp_path / "taker.duckdb")
    yield taker_store
    taker_store.close()


def desk_for(store, kalshi, tmp_path, clock=lambda: NOW):
    return TakerDesk(kalshi, store, tmp_path / "taker.stop", now=clock)


def request(**overrides):
    validated = validate_order_request({**REQUEST_BODY, **overrides})
    assert isinstance(validated, dict), validated
    return validated


def order_row(store, client_order_id):
    return next(order for order in store.open_orders() if order["client_order_id"] == client_order_id)


def test_order_body_buys_yes_as_a_bid_and_no_as_an_ask_at_one_minus():
    assert order_body(TICKER, "yes", 0.57, 10, "ubot-x") == {
        "ticker": TICKER, "side": "bid", "price": "0.5700", "count": "10.00", "time_in_force": "immediate_or_cancel",
        "self_trade_prevention_type": "taker_at_cross", "client_order_id": "ubot-x"}
    assert order_body(TICKER, "no", 0.43, 3, "ubot-y")["side"] == "ask"
    assert order_body(TICKER, "no", 0.43, 3, "ubot-y")["price"] == "0.5700"


def test_a_filled_order_is_sent_immediate_or_cancel_and_its_fills_recorded(store, tmp_path):
    kalshi = ScriptedKalshi()
    result = desk_for(store, kalshi, tmp_path).place(request())
    assert result["status"] == "filled" and result["filled"] == 10
    assert result["cost"] == pytest.approx(10 * 0.57 + 0.17, abs=1e-3)
    # Kalshi cancelled the rest in the same request (remaining_count 0): no cancel call, just the fills.
    assert kalshi.methods() == ["POST", "GET"]
    assert kalshi.calls[0][2]["client_order_id"] == REQUEST_BODY["clientOrderId"]
    assert kalshi.calls[0][2]["time_in_force"] == "immediate_or_cancel"
    row = order_row(store, REQUEST_BODY["clientOrderId"])
    assert (row["exchange_index"], row["cancelled_at"], row["fills_final"]) == (None, NOW, False)
    # Until the final fills read the whole requested cost is held.
    assert store.open_exposure() == pytest.approx(10 * 0.587157 + 0.01, abs=1e-3)


@pytest.mark.parametrize("remaining_count", [None, "4.00"])
def test_a_reply_without_a_zero_remaining_count_is_cancelled_by_market_and_shard(store, tmp_path, remaining_count):
    kalshi = ScriptedKalshi(exchange_index=3, remaining_count=remaining_count)
    result = desk_for(store, kalshi, tmp_path).place(request())
    assert result["status"] == "filled" and result["filled"] == 10
    assert kalshi.methods() == ["POST", "DELETE", "GET"]
    # The cancel names the market and the shard Kalshi's reply gave, never a bare id (shard-0 routing).
    assert kalshi.calls[1][1] == f"/portfolio/events/orders/k-1?market_ticker={TICKER}&exchange_index=3"
    assert order_row(store, REQUEST_BODY["clientOrderId"])["cancelled_at"] == NOW


def test_the_same_client_order_id_is_never_sent_twice(store, tmp_path):
    kalshi = ScriptedKalshi()
    desk = desk_for(store, kalshi, tmp_path)
    desk.place(request())
    assert desk.place(request())["status"] == "duplicate"
    assert kalshi.methods().count("POST") == 1


def test_the_stop_file_sends_nothing(store, tmp_path):
    (tmp_path / "taker.stop").touch()
    kalshi = ScriptedKalshi()
    assert desk_for(store, kalshi, tmp_path).place(request())["status"] == "stopped"
    assert kalshi.calls == []


def test_the_ceiling_cuts_contracts_then_refuses(store, tmp_path, monkeypatch):
    monkeypatch.setattr(taker_orders, "BOT_BANKROLL", 10.0)
    kalshi = ScriptedKalshi(filled=0)
    desk = desk_for(store, kalshi, tmp_path)
    first = desk.place(request(count=100))
    # $10 - 1c slack at 0.587157 a contract: 17 contracts.
    assert first["count"] == 17
    assert kalshi.calls[0][2]["count"] == "17.00"
    refused = desk.place(request(clientOrderId="ubot-" + "b" * 32, eventId=EVENT_ID + 1))
    assert refused["status"] == "ceiling"
    assert kalshi.methods().count("POST") == 1


def test_a_refused_order_holds_nothing(store, tmp_path):
    result = desk_for(store, ScriptedKalshi(post_status=400), tmp_path).place(request())
    assert result["status"] == "rejected"
    assert store.open_exposure() == 0


def test_final_fills_then_settlement_free_the_money(store, tmp_path):
    clock = {"now": NOW}
    kalshi = ScriptedKalshi(filled=4)
    desk = desk_for(store, kalshi, tmp_path, clock=lambda: clock["now"])
    assert desk.place(request())["status"] == "partial"
    # A minute after the cancel the fills are final: only what filled is held.
    clock["now"] = NOW + timedelta(seconds=61)
    desk.place(request(clientOrderId="ubot-" + "c" * 32, count=1, ticker=OTHER_TICKER, heldIds=[HELD_FIRST],
                       decidedAt=clock["now"].isoformat()))
    first_cost = 4 * 0.57 + 0.17
    assert store.open_exposure() == pytest.approx(first_cost + 1 * 0.587157 + 0.01, abs=1e-3)
    # After the game starts and the market settles, both orders are freed.
    kalshi.market_status = "finalized"
    clock["now"] = datetime(2026, 10, 16, tzinfo=timezone.utc)
    desk.refresh()
    assert store.open_exposure() == 0


def test_an_add_on_the_game_waits_until_the_earlier_fill_is_among_the_held_bets(store, tmp_path):
    kalshi = ScriptedKalshi(filled=4)
    desk = desk_for(store, kalshi, tmp_path)
    assert desk.place(request())["status"] == "partial"
    # The same game, another rung: the stake was sized without the first fill.
    blind = desk.place(request(clientOrderId="ubot-" + "b" * 32, ticker=OTHER_TICKER))
    assert blind["status"] == "unseen_fill"
    assert HELD_FIRST in blind["detail"]
    # Another game is not blocked; the same game is, once the fill is among the held ids.
    assert desk.place(request(clientOrderId="ubot-" + "c" * 32, eventId=EVENT_ID + 1))["status"] == "partial"
    seen = desk.place(request(clientOrderId="ubot-" + "d" * 32, ticker=OTHER_TICKER, heldIds=[HELD_FIRST]))
    assert seen["status"] == "partial"
    assert kalshi.methods().count("POST") == 3


def test_an_unfilled_order_blocks_the_game_only_until_its_fills_are_final(store, tmp_path):
    clock = {"now": NOW}
    kalshi = ScriptedKalshi(filled=0)
    desk = desk_for(store, kalshi, tmp_path, clock=lambda: clock["now"])
    assert desk.place(request())["status"] == "unfilled"
    # Right after the cancel a fill could still trail in: blocked.
    assert desk.place(request(clientOrderId="ubot-" + "b" * 32, ticker=OTHER_TICKER))["status"] == "unseen_fill"
    clock["now"] = NOW + timedelta(seconds=61)
    assert desk.place(request(clientOrderId="ubot-" + "c" * 32, ticker=OTHER_TICKER,
                              decidedAt=clock["now"].isoformat()))["status"] == "unfilled"


def test_a_cancel_404_still_resting_holds_the_order_and_the_refresh_retries_it(store, tmp_path):
    clock = {"now": NOW}
    kalshi = ScriptedKalshi(filled=0, delete_status=404, resting=[{"order_id": "k-1", "ticker": TICKER, "exchange_index": 3}],
                            remaining_count=None)
    desk = desk_for(store, kalshi, tmp_path, clock=lambda: clock["now"])
    result = desk.place(request())
    assert result["status"] == "resting"
    assert kalshi.methods() == ["POST", "DELETE", "GET"]
    assert kalshi.calls[2][1].startswith("/portfolio/orders?status=resting")
    assert order_row(store, REQUEST_BODY["clientOrderId"])["status"] == "resting"
    assert store.open_exposure() == pytest.approx(10 * 0.587157 + 0.01, abs=1e-3)
    # The next refresh: the cancel goes through, the fills are read, and a minute on they are final.
    kalshi.delete_status = 200
    kalshi.resting = []
    clock["now"] = NOW + timedelta(seconds=30)
    desk.refresh()
    row = order_row(store, REQUEST_BODY["clientOrderId"])
    assert (row["status"], row["cancelled_at"]) == ("unfilled", clock["now"])
    clock["now"] = NOW + timedelta(seconds=95)
    desk.refresh()
    assert order_row(store, REQUEST_BODY["clientOrderId"])["fills_final"] is True
    assert store.open_exposure() == 0


def test_a_cancel_404_verified_gone_counts_as_cancelled(store, tmp_path):
    kalshi = ScriptedKalshi(filled=10, delete_status=404, resting=[{"order_id": "someone-else", "ticker": TICKER}], remaining_count=None)
    assert desk_for(store, kalshi, tmp_path).place(request())["status"] == "filled"
    assert order_row(store, REQUEST_BODY["clientOrderId"])["cancelled_at"] == NOW


def test_a_cancel_404_with_the_resting_list_down_holds(store, tmp_path):
    kalshi = ScriptedKalshi(filled=10, delete_status=404, resting=None, remaining_count=None)
    assert desk_for(store, kalshi, tmp_path).place(request())["status"] == "resting"


def test_no_answer_to_the_send_holds_the_full_cost(store, tmp_path):
    class Silent(ScriptedKalshi):
        def request(self, method, path, body=None):
            if method == "POST":
                raise OSError("timed out")
            return super().request(method, path, body)
    result = desk_for(store, Silent(), tmp_path).place(request())
    assert result["status"] == "unknown"
    assert store.open_exposure() == pytest.approx(10 * 0.587157 + 0.01, abs=1e-3)


class SilentPost(ScriptedKalshi):
    def request(self, method, path, body=None):
        if method == "POST":
            self.calls.append((method, path, body))
            raise OSError("timed out")
        return super().request(method, path, body)


def test_a_lost_order_found_resting_is_adopted_cancelled_and_read(store, tmp_path):
    clock = {"now": NOW}
    cid = REQUEST_BODY["clientOrderId"]
    kalshi = SilentPost(filled=3, resting=[{"order_id": "k-9", "client_order_id": cid, "ticker": TICKER, "exchange_index": 0}])
    desk = desk_for(store, kalshi, tmp_path, clock=lambda: clock["now"])
    assert desk.place(request())["status"] == "unknown"
    clock["now"] = NOW + timedelta(seconds=30)
    desk.refresh()
    row = order_row(store, cid)
    assert (row["status"], row["kalshi_order_id"], row["exchange_index"], row["count_filled"]) == ("partial", "k-9", 0, 3)
    assert ("DELETE", f"/portfolio/events/orders/k-9?market_ticker={TICKER}&exchange_index=0", None) in kalshi.calls


def test_a_lost_order_listed_nowhere_after_ten_minutes_holds_nothing(store, tmp_path):
    clock = {"now": NOW}
    kalshi = SilentPost(resting=[], listed=[])
    desk = desk_for(store, kalshi, tmp_path, clock=lambda: clock["now"])
    desk.place(request())
    clock["now"] = NOW + timedelta(minutes=5)
    desk.refresh()
    assert store.open_exposure() > 0, "too soon to give up on it"
    clock["now"] = NOW + timedelta(minutes=11)
    desk.refresh()
    assert store.open_exposure() == 0
    assert store.open_orders() == []
    assert any(path.startswith(f"/portfolio/orders?ticker={TICKER}") for _, path, _ in kalshi.calls)


def test_a_lost_order_holds_while_the_listing_is_down(store, tmp_path):
    clock = {"now": NOW + timedelta(minutes=11)}
    kalshi = SilentPost(resting=None, listed=None)
    desk = desk_for(store, kalshi, tmp_path, clock=lambda: clock["now"])
    desk.place(request(decidedAt=clock["now"].isoformat()))
    clock["now"] += timedelta(seconds=30)
    desk.refresh()
    assert store.open_exposure() > 0


def test_an_untracked_bot_order_resting_on_kalshi_is_cancelled_by_the_refresh(store, tmp_path):
    stray = {"order_id": "k-stray", "client_order_id": "ubot-" + "f" * 32, "ticker": OTHER_TICKER, "exchange_index": 2}
    mine = {"order_id": "k-1", "client_order_id": REQUEST_BODY["clientOrderId"], "ticker": TICKER}
    theirs = {"order_id": "k-mm", "client_order_id": "mm-quote-1", "ticker": TICKER}
    kalshi = ScriptedKalshi(filled=10, resting=[stray, mine, theirs])
    desk = desk_for(store, kalshi, tmp_path)
    desk.place(request())
    desk.refresh()
    cancelled = [path for method, path, _ in kalshi.calls if method == "DELETE"]
    assert f"/portfolio/events/orders/k-stray?market_ticker={OTHER_TICKER}&exchange_index=2" in cancelled
    assert not any("k-mm" in path for path in cancelled), "another client's order is not the bot's to cancel"


def test_settlement_is_asked_once_per_five_minutes_per_market(store, tmp_path):
    clock = {"now": NOW}
    kalshi = ScriptedKalshi(filled=10)
    desk = desk_for(store, kalshi, tmp_path, clock=lambda: clock["now"])
    desk.place(request())
    clock["now"] = datetime(2026, 10, 15, 22, 0, tzinfo=timezone.utc)
    for _ in range(3):
        desk.refresh()
        clock["now"] += timedelta(seconds=30)
    assert len([path for method, path, _ in kalshi.calls if path.startswith("/markets/")]) == 1


def test_the_refresh_loop_stops_with_the_service(store, tmp_path):
    stop = threading.Event()
    stop.set()
    desk_for(store, ScriptedKalshi(), tmp_path).run_refresh_loop(stop)


@pytest.mark.parametrize("override, message", [
    ({"clientOrderId": "bot-1"}, "clientOrderId must be a 'ubot-' id"),
    ({"side": "maybe"}, "side must be 'yes' or 'no'"),
    ({"ask": 1.2}, "ask must be a price"),
    ({"allIn": 0.5}, "allIn must be at least the ask"),
    ({"count": 0}, "count must be a whole number"),
    ({"eventId": "123742"}, "eventId must be Unabated's event id"),
    ({"heldIds": "kalshi:x:yes"}, "heldIds must be a list"),
    ({"heldIds": [7]}, "heldIds must be a list"),
    ({"eventStart": None}, "eventStart must be an ISO time"),
    ({"decidedAt": None}, "decidedAt must be an ISO time"),
    ({"decidedAt": "2026-10-14T11:59:50"}, "decidedAt must be an ISO time"),
    ({"kalshiBaseUrl": "kalshi.test"}, "kalshiBaseUrl must be the Kalshi API"),
])
def test_a_bad_request_names_the_field(override, message):
    assert validate_order_request({**REQUEST_BODY, **override}).startswith(message)


def test_a_stale_decision_a_started_game_or_another_kalshi_is_refused_unsent(store, tmp_path):
    kalshi = ScriptedKalshi()
    desk = desk_for(store, kalshi, tmp_path)
    stale = desk.place(request(decidedAt=(NOW - timedelta(seconds=61)).isoformat()))
    assert (stale["status"], stale["detail"]) == ("stale_decision", "decided 61 s ago, over 60 s")
    started = desk.place(request(eventStart=NOW.isoformat()))
    assert started["status"] == "started"
    elsewhere = desk.place(request(kalshiBaseUrl="https://demo-api.kalshi.co/trade-api/v2"))
    assert elsewhere["status"] == "env_mismatch"
    assert elsewhere["detail"] == ("the runner priced https://demo-api.kalshi.co/trade-api/v2, "
                                   "this desk trades https://kalshi.test/trade-api/v2")
    assert kalshi.calls == [] and store.open_orders() == []
    # The same request a minute younger, on this Kalshi, before the game: sent.
    assert desk.place(request())["status"] == "filled"


def test_fills_are_read_across_every_page(store, tmp_path):
    kalshi = ScriptedKalshi(fills_pages=[(4, "page2"), (6, "")])
    desk = desk_for(store, kalshi, tmp_path)
    reply = desk.place(request())
    assert (reply["status"], reply["filled"]) == ("filled", 10)
    assert order_row(store, REQUEST_BODY["clientOrderId"])["count_filled"] == 10
    assert [path for method, path, _ in kalshi.calls if path.startswith("/portfolio/fills")] == [
        "/portfolio/fills?order_id=k-1&limit=1000", "/portfolio/fills?order_id=k-1&limit=1000&cursor=page2"]


def test_a_cut_off_answer_from_kalshi_is_no_answer(tmp_path, monkeypatch):
    import http.client
    import urllib.request
    client = signing_client(tmp_path)

    def cut_off(request, timeout):
        raise http.client.IncompleteRead(b"")
    monkeypatch.setattr(urllib.request, "urlopen", cut_off)
    with pytest.raises(OSError, match="IncompleteRead"):
        client.request("GET", "/markets/X")


def test_open_taker_desk_says_why_it_is_off_and_the_route_repeats_it(store, tmp_path, monkeypatch):
    import json
    import urllib.error
    import urllib.request
    from http.server import ThreadingHTTPServer
    from unabated_ticket.bets_service import config, service
    from unabated_ticket.bets_service.store import BetsStore
    monkeypatch.setattr(config, "KALSHI_TAKER_API_KEY_ID", "")
    assert service.open_taker_desk(store)[0] is None
    assert "no bot Kalshi key" in service.open_taker_desk(store)[1]
    monkeypatch.setattr(config, "KALSHI_TAKER_API_KEY_ID", "key-1")
    monkeypatch.setattr(config, "KALSHI_TAKER_PRIVATE_KEY_PATH", str(tmp_path / "missing.pem"))
    assert service.open_taker_desk(None)[1].startswith("the taker log")
    desk, reason = service.open_taker_desk(store)
    assert desk is None and reason == f"no bot key file at {tmp_path / 'missing.pem'}"
    bets_store = BetsStore(tmp_path / "bets.duckdb", 7)
    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(
        bets_store, 0.0, ["kalshi"], taker_store=store, taker_desk=None, taker_desk_off_reason=reason))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        req = urllib.request.Request(f"http://127.0.0.1:{server.server_address[1]}/taker/order.json",
                                     data=json.dumps(REQUEST_BODY).encode(), method="POST",
                                     headers={"Content-Type": "application/json"})
        with pytest.raises(urllib.error.HTTPError) as refused:
            urllib.request.urlopen(req, timeout=5)
        assert refused.value.code == 503
        assert json.loads(refused.value.read()) == {"status": "off", "error": reason}
    finally:
        server.shutdown()
        server.server_close()
        bets_store.close()


def signing_client(tmp_path):
    """A KalshiOrderClient on a fresh RSA key written to tmp_path."""
    from cryptography.hazmat.primitives import serialization
    from cryptography.hazmat.primitives.asymmetric import rsa
    key = rsa.generate_private_key(public_exponent=65537, key_size=2048)
    key_path = tmp_path / "bot.pem"
    key_path.write_bytes(key.private_bytes(serialization.Encoding.PEM, serialization.PrivateFormat.PKCS8,
                                           serialization.NoEncryption()))
    client = taker_orders.KalshiOrderClient("key-id", str(key_path), "https://kalshi.test/trade-api/v2/")
    client.test_key = key
    return client


def test_the_client_signs_with_kalshis_rsa_pss_scheme(tmp_path):
    from cryptography.hazmat.primitives import hashes
    from cryptography.hazmat.primitives.asymmetric import padding
    import base64
    client = signing_client(tmp_path)
    key = client.test_key
    assert client.base_url == "https://kalshi.test/trade-api/v2"
    signature = base64.b64decode(client._signature("1700000000000", "POST", "/portfolio/events/orders?x=1"))
    key.public_key().verify(signature, b"1700000000000POST/trade-api/v2/portfolio/events/orders",
                            padding.PSS(mgf=padding.MGF1(hashes.SHA256()), salt_length=padding.PSS.DIGEST_LENGTH),
                            hashes.SHA256())


def test_the_order_route_refuses_the_tailnet_name(tmp_path):
    import json
    import urllib.error
    import urllib.request
    from http.server import ThreadingHTTPServer
    from unabated_ticket.bets_service import service
    from unabated_ticket.bets_service.store import BetsStore
    bets_store = BetsStore(tmp_path / "bets.duckdb", 7)
    taker_store = TakerStore(tmp_path / "taker.duckdb")
    kalshi = ScriptedKalshi()
    desk = TakerDesk(kalshi, taker_store, tmp_path / "taker.stop", now=lambda: NOW)
    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(
        bets_store, 0.0, ["kalshi"], extra_hosts=("box.tail.ts.net",), taker_store=taker_store, taker_desk=desk))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{server.server_address[1]}/taker/order.json"
    try:
        def post(host):
            req = urllib.request.Request(url, data=json.dumps(REQUEST_BODY).encode(), method="POST",
                                         headers={"Content-Type": "application/json", "Host": host})
            try:
                with urllib.request.urlopen(req, timeout=5) as response:
                    return response.status, json.loads(response.read())
            except urllib.error.HTTPError as error:
                return error.code, json.loads(error.read())
        assert post("box.tail.ts.net")[0] == 403
        assert kalshi.calls == []
        status, reply = post(f"127.0.0.1:{server.server_address[1]}")
        assert (status, reply["status"]) == (200, "filled")
    finally:
        server.shutdown()
        server.server_close()
