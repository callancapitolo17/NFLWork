"""The bot's order desk (taker_orders.py) against a scripted Kalshi: the stop
file, the $500 ceiling on open exposure, the row written before the send,
limit-at-the-ask bodies, cancel-the-rest, the fills read, the final re-read
and settlement freeing the money. No network."""
from datetime import datetime, timedelta, timezone

import pytest

from unabated_ticket.bets_service import taker_orders
from unabated_ticket.bets_service.taker_orders import TakerDesk, order_body, validate_order_request
from unabated_ticket.bets_service.taker_store import TakerStore

NOW = datetime(2026, 10, 14, 12, 0, tzinfo=timezone.utc)
TICKER = "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36"
REQUEST_BODY = {"clientOrderId": "ubot-" + "a" * 32, "ticker": TICKER, "side": "yes", "ask": 0.57,
                "allIn": 0.587157, "count": 10, "rowKey": "420942308:ms105:si1:tid717:alt-35.5",
                "eventStart": "2026-10-15T19:30:00.000Z"}


class ScriptedKalshi:
    """Answers POST/DELETE/GET by path; records every call."""

    def __init__(self, filled: float = 10, market_status: str = "active", post_status: int = 201):
        self.calls = []
        self.filled = filled
        self.market_status = market_status
        self.post_status = post_status

    def request(self, method, path, body=None):
        self.calls.append((method, path, body))
        if method == "POST":
            return (self.post_status, {"order_id": "k-1"}) if self.post_status == 201 else (self.post_status, {"error": "no"})
        if method == "DELETE":
            return 404, {}
        if path.startswith("/portfolio/fills"):
            fills = [{"count_fp": str(self.filled), "yes_price_dollars": "0.5700", "no_price_dollars": "0.4300",
                      "fee_cost": "0.17"}] if self.filled else []
            return 200, {"fills": fills}
        if path.startswith("/markets/"):
            return 200, {"market": {"status": self.market_status}}
        raise AssertionError(f"unexpected {method} {path}")


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


def test_order_body_buys_yes_as_a_bid_and_no_as_an_ask_at_one_minus():
    assert order_body(TICKER, "yes", 0.57, 10, "ubot-x") == {
        "ticker": TICKER, "side": "bid", "price": "0.5700", "count": "10.00", "time_in_force": "good_till_canceled",
        "self_trade_prevention_type": "taker_at_cross", "client_order_id": "ubot-x"}
    assert order_body(TICKER, "no", 0.43, 3, "ubot-y")["side"] == "ask"
    assert order_body(TICKER, "no", 0.43, 3, "ubot-y")["price"] == "0.5700"


def test_a_filled_order_is_sent_once_cancelled_and_its_fills_recorded(store, tmp_path):
    kalshi = ScriptedKalshi()
    result = desk_for(store, kalshi, tmp_path).place(request())
    assert result["status"] == "filled" and result["filled"] == 10
    assert result["cost"] == pytest.approx(10 * 0.57 + 0.17, abs=1e-3)
    assert [call[0] for call in kalshi.calls] == ["POST", "DELETE", "GET"]
    assert kalshi.calls[0][2]["client_order_id"] == REQUEST_BODY["clientOrderId"]
    # Until the final fills read the whole requested cost is held.
    assert store.open_exposure() == pytest.approx(10 * 0.587157 + 0.01, abs=1e-3)


def test_the_same_client_order_id_is_never_sent_twice(store, tmp_path):
    kalshi = ScriptedKalshi()
    desk = desk_for(store, kalshi, tmp_path)
    desk.place(request())
    assert desk.place(request())["status"] == "duplicate"
    assert [call[0] for call in kalshi.calls].count("POST") == 1


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
    refused = desk.place(request(clientOrderId="ubot-" + "b" * 32))
    assert refused["status"] == "ceiling"
    assert [call[0] for call in kalshi.calls].count("POST") == 1


def test_a_refused_order_holds_nothing(store, tmp_path):
    result = desk_for(store, ScriptedKalshi(post_status=400), tmp_path).place(request())
    assert result["status"] == "rejected"
    assert store.open_exposure() == 0


def test_final_fills_then_settlement_free_the_money(store, tmp_path):
    clock = {"now": NOW}
    kalshi = ScriptedKalshi(filled=4)
    desk = desk_for(store, kalshi, tmp_path, clock=lambda: clock["now"])
    assert desk.place(request())["status"] == "partial"
    # A minute later the fills are final: only what filled is held.
    clock["now"] = NOW + timedelta(seconds=61)
    desk.place(request(clientOrderId="ubot-" + "c" * 32, count=1))
    first_cost = 4 * 0.57 + 0.17
    assert store.open_exposure() == pytest.approx(first_cost + 1 * 0.587157 + 0.01, abs=1e-3)
    # After the game starts and the market settles, both orders are freed.
    kalshi.market_status = "finalized"
    clock["now"] = datetime(2026, 10, 16, tzinfo=timezone.utc)
    desk.place(request(clientOrderId="ubot-" + "d" * 32, count=1))
    assert store.open_exposure() == pytest.approx(1 * 0.587157 + 0.01, abs=1e-3)


def test_no_answer_to_the_send_holds_the_full_cost(store, tmp_path):
    class Silent(ScriptedKalshi):
        def request(self, method, path, body=None):
            if method == "POST":
                raise OSError("timed out")
            return super().request(method, path, body)
    result = desk_for(store, Silent(), tmp_path).place(request())
    assert result["status"] == "unknown"
    assert store.open_exposure() == pytest.approx(10 * 0.587157 + 0.01, abs=1e-3)


@pytest.mark.parametrize("override, message", [
    ({"clientOrderId": "bot-1"}, "clientOrderId must be a 'ubot-' id"),
    ({"side": "maybe"}, "side must be 'yes' or 'no'"),
    ({"ask": 1.2}, "ask must be a price"),
    ({"allIn": 0.5}, "allIn must be at least the ask"),
    ({"count": 0}, "count must be a whole number"),
    ({"eventStart": None}, "eventStart must be an ISO time"),
])
def test_a_bad_request_names_the_field(override, message):
    assert validate_order_request({**REQUEST_BODY, **override}).startswith(message)


def test_the_client_signs_with_kalshis_rsa_pss_scheme(tmp_path):
    from cryptography.hazmat.primitives import hashes, serialization
    from cryptography.hazmat.primitives.asymmetric import padding, rsa
    import base64
    key = rsa.generate_private_key(public_exponent=65537, key_size=2048)
    key_path = tmp_path / "bot.pem"
    key_path.write_bytes(key.private_bytes(serialization.Encoding.PEM, serialization.PrivateFormat.PKCS8,
                                           serialization.NoEncryption()))
    client = taker_orders.KalshiOrderClient("key-id", str(key_path), "https://kalshi.test/trade-api/v2")
    signature = base64.b64decode(client._signature("1700000000000", "POST", "/portfolio/events/orders?x=1"))
    key.public_key().verify(signature, b"1700000000000POST/trade-api/v2/portfolio/events/orders",
                            padding.PSS(mgf=padding.MGF1(hashes.SHA256()), salt_length=padding.PSS.DIGEST_LENGTH),
                            hashes.SHA256())


def test_the_order_route_refuses_the_tailnet_name(tmp_path):
    import json
    import threading
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
        bets_store.close()
        taker_store.close()
