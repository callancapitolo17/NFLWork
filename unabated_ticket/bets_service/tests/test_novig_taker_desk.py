"""novig_taker_desk.py on a fake Novig client: the stop files, freshness of
the held Novig bets, the cross-venue unseen-fill rule, the $500 ceiling in 1c
contracts, insert-before-send, one IOC send never retried, the read to final
and its fills, lost sends looked up by clientId, and settlement."""
from datetime import datetime, timedelta, timezone

import pytest

from unabated_ticket.bets_service import novig_taker_desk as desk_module
from unabated_ticket.bets_service.novig_taker_desk import (NovigTakerDesk, account_problem, held_id_of,
                                                           novig_unseen_fill, validate_order_request)
from unabated_ticket.bets_service.taker_store import TakerStore

NOW = datetime(2026, 10, 18, 16, 0, tzinfo=timezone.utc)
BASE = "https://api.novig.com"
OUTCOME = "3f2504e0-4f89-11d3-9a0c-0305e82c3301"
MARKET = "6a1b2c3d-4e5f-4a6b-8c7d-9e0f1a2b3c4d"
CLIENT_ID = "0b3f5f20-9c1b-4d6e-8f3a-5c2d1e7b4a90"
ORDER_ID = "9d8c7b6a-5f4e-4d3c-2b1a-0f9e8d7c6b5a"


class FakeNovig:
    base_url = BASE

    def __init__(self):
        self.placed = []
        self.place_reply = (201, {"orderId": ORDER_ID, "clientId": CLIENT_ID})
        self.order_reply = (200, {"orderId": ORDER_ID, "status": "FILLED", "qty": 2150, "remaining": 0})
        self.fill_items = [{"fillId": "f1", "orderId": ORDER_ID, "qty": 2150, "cost": "9.99750", "taker": True}]
        self.lookup = None
        self.market_reply = {"marketId": MARKET, "status": "OPEN"}
        self.market_reads = 0
        self.cancelled = []

    def place_ioc_order(self, outcome_id, price, qty, client_id):
        self.placed.append((outcome_id, price, qty, client_id))
        if isinstance(self.place_reply, Exception):
            raise self.place_reply
        return self.place_reply

    def get_order(self, order_id):
        return self.order_reply

    def cancel_order(self, order_id):
        self.cancelled.append(order_id)
        return 200, {}

    def fills(self, order_id=None, outcome_id=None, starts_after_ms=None):
        return list(self.fill_items)

    def find_order_by_client_id(self, outcome_id, client_id):
        return self.lookup

    def market(self, market_id):
        self.market_reads += 1
        return self.market_reply


def fresh_status():
    stamp = (NOW - timedelta(seconds=40)).isoformat()
    return {"novig": {"ok": True, "fetchedAt": stamp}, "novig_bot": {"ok": True, "fetchedAt": stamp}}


@pytest.fixture
def setup(tmp_path):
    store = TakerStore(tmp_path / "taker.duckdb")
    client = FakeNovig()
    clock = {"now": NOW}
    status = {"value": fresh_status()}
    desk = NovigTakerDesk(client, store, (tmp_path / "taker.stop", tmp_path / "novig_taker.stop"),
                          lambda: status["value"], now=lambda: clock["now"], sleep=lambda _s: None)
    yield desk, client, store, clock, status, tmp_path
    store.close()


def request(**changes):
    body = {"clientId": CLIENT_ID, "outcomeId": OUTCOME, "marketId": MARKET, "rowKey": "r1", "price": 0.465, "qty": 2150,
            "eventId": 4242, "heldIds": [], "eventStart": (NOW + timedelta(hours=3)).isoformat(),
            "decidedAt": (NOW - timedelta(seconds=5)).isoformat(), "novigBaseUrl": BASE,
            "line": {"league": "nfl", "betType": "spread", "period": "FG", "side": "home", "points": -3.5,
                     "awayTeam": "Bears", "homeTeam": "Packers", "awayTeamId": "11", "homeTeamId": "12"}}
    body.update(changes)
    parsed = validate_order_request(body)
    assert isinstance(parsed, dict), parsed
    return parsed


def test_a_clean_order_fills_and_keeps_the_board_line(setup):
    desk, client, store, *_ = setup
    result = desk.place(request())
    assert result["status"] == "filled" and result["filled"] == 2150 and result["cost"] == pytest.approx(9.9975)
    assert client.placed == [(OUTCOME, 0.465, 2150, CLIENT_ID)]
    context = store.novig_order_contexts()[ORDER_ID]
    assert (context["league"], context["side"], context["points"], context["home_team_id"]) == ("nfl", "home", -3.5, "12")


@pytest.mark.parametrize("stop_name", ["taker.stop", "novig_taker.stop"])
def test_either_stop_file_stops_it(setup, stop_name):
    desk, client, _store, _clock, _status, tmp_path = setup
    (tmp_path / stop_name).touch()
    assert desk.place(request())["status"] == "stopped" and client.placed == []


def test_refuses_a_stale_decision_a_started_game_and_another_novig(setup):
    desk, client, *_ = setup
    assert desk.place(request(decidedAt=(NOW - timedelta(seconds=90)).isoformat()))["status"] == "stale_decision"
    assert desk.place(request(eventStart=(NOW - timedelta(minutes=1)).isoformat()))["status"] == "started"
    assert desk.place(request(novigBaseUrl="https://api.qa.novig.com"))["status"] == "env_mismatch"
    assert client.placed == []


@pytest.mark.parametrize("change", [{"novig_bot": {"ok": False, "error": "HTTP 500"}},
                                    {"novig": {"ok": True, "fetchedAt": (NOW - timedelta(minutes=6)).isoformat()}},
                                    {"novig_bot": None}])
def test_refuses_while_the_held_novig_bets_are_unknown(setup, change):
    desk, client, _store, _clock, status, _ = setup
    value = fresh_status()
    for name, source in change.items():
        if source is None:
            value.pop(name)
        else:
            value[name] = source
    status["value"] = value
    assert desk.place(request())["status"] == "stale_holdings" and client.placed == []


def test_an_add_waits_until_the_earlier_fill_is_held(setup):
    desk, client, store, *_ = setup
    assert desk.place(request())["status"] == "filled"
    second = request(clientId="1b3f5f20-9c1b-4d6e-8f3a-5c2d1e7b4a91")
    assert desk.place(second)["status"] == "unseen_fill"
    held = request(clientId="1b3f5f20-9c1b-4d6e-8f3a-5c2d1e7b4a91", heldIds=[held_id_of(ORDER_ID)])
    client.place_reply = (201, {"orderId": "8d8c7b6a-5f4e-4d3c-2b1a-0f9e8d7c6b5a"})
    assert desk.place(held)["status"] == "filled"
    # Kalshi's desk asks the same question about Novig orders on its game.
    assert novig_unseen_fill(store, 4242, set()) is not None
    assert novig_unseen_fill(store, 9999, set()) is None


def test_a_kalshi_bot_order_on_the_game_blocks_a_novig_add(setup):
    desk, client, store, *_ = setup
    store.insert_order({"client_order_id": "ubot-1", "placed_at": NOW, "row_key": "k", "event_id": 4242, "ticker": "KX-1",
                        "contract_side": "yes", "ask_price": 0.5, "all_in_price": 0.52, "count_requested": 10,
                        "event_start": NOW + timedelta(hours=3), "open_cost": 5.2, "subaccount": 1})
    assert desk.place(request())["status"] == "unseen_fill" and client.placed == []
    assert desk.place(request(heldIds=["kalshi:KX-1:yes:s1"]))["status"] == "filled"


def test_the_ceiling_cuts_qty_in_one_cent_contracts(setup):
    desk, client, store, *_ = setup
    client.fill_items = [{"fillId": "f1", "qty": 100, "cost": "0.46500", "taker": True}]
    # $500 at 0.465 is 107526 contracts; ask for more.
    desk.place(request(qty=200000))
    assert client.placed[0][2] == 107526
    assert store.novig_open_exposure() == pytest.approx(107526 * 0.465 / 100, abs=1e-6)
    client.placed.clear()
    assert desk.place(request(clientId="2b3f5f20-9c1b-4d6e-8f3a-5c2d1e7b4a92", eventId=7,
                              heldIds=[held_id_of(ORDER_ID)]))["status"] == "ceiling"
    assert client.placed == []


def test_a_refusal_holds_nothing_and_a_lost_answer_holds_the_cost(setup):
    desk, client, store, *_ = setup
    client.place_reply = (451, {"code": "GEOLOCATION_EXPIRED"})
    assert desk.place(request())["status"] == "geo_blocked"
    assert store.novig_open_exposure() == 0
    stored = store._con.execute("SELECT status, open_cost, settled_at IS NOT NULL, novig_order_id, error "
                                "FROM novig_taker_orders WHERE client_id = ?", [CLIENT_ID]).fetchone()
    assert stored[:4] == ("rejected", 0.0, True, None) and stored[4].startswith("HTTP 451 GEOLOCATION_EXPIRED")
    client.place_reply = OSError("timed out")
    lost = request(clientId="3b3f5f20-9c1b-4d6e-8f3a-5c2d1e7b4a93", eventId=8)
    assert desk.place(lost)["status"] == "unknown"
    assert store.novig_open_exposure() == pytest.approx(9.9975)
    assert len(client.placed) == 2  # never retried


def test_a_lost_send_is_found_by_client_id_or_freed(setup):
    desk, client, store, clock, *_ = setup
    client.place_reply = (503, {"code": "UNAVAILABLE"})
    assert desk.place(request())["status"] == "unknown"
    clock["now"] = NOW + timedelta(seconds=30)
    desk.refresh()
    assert store.open_novig_orders()[0]["status"] == "unknown"  # too young to look up
    clock["now"] = NOW + timedelta(seconds=61)
    client.lookup = {"orderId": ORDER_ID, "clientId": CLIENT_ID, "status": "FILLED"}
    desk.refresh()
    order = store.open_novig_orders()[0]
    assert (order["novig_order_id"], order["status"]) == (ORDER_ID, "filled")


def test_a_lost_send_listed_nowhere_never_reached_novig(setup):
    desk, client, store, clock, *_ = setup
    client.place_reply = (503, {})
    desk.place(request())
    clock["now"] = NOW + timedelta(seconds=61)
    desk.refresh()
    assert store.novig_open_exposure() == 0 and store.open_novig_orders() == []


def test_an_order_not_final_yet_is_read_by_the_refresh_then_fills_go_final(setup):
    desk, client, store, clock, *_ = setup
    client.order_reply = (404, {"code": "NOT_FOUND"})
    assert desk.place(request())["status"] == "sent"
    client.order_reply = (200, {"orderId": ORDER_ID, "status": "CANCELED", "qty": 2150, "remaining": 0})
    client.fill_items = [{"fillId": "f1", "qty": 1000, "cost": "4.65000", "taker": True}]
    clock["now"] = NOW + timedelta(seconds=31)
    desk.refresh()
    assert store.open_novig_orders()[0]["status"] == "partial"
    clock["now"] = NOW + timedelta(seconds=92)
    desk.refresh()
    order = store.open_novig_orders()[0]
    assert order["fills_final"] and order["open_cost"] == pytest.approx(4.65)


def test_an_open_ioc_is_cancelled_and_held(setup):
    desk, client, store, *_ = setup
    client.order_reply = (200, {"orderId": ORDER_ID, "status": "OPEN", "remaining": 100})
    assert desk.place(request())["status"] == "resting"
    assert client.cancelled == [ORDER_ID] and store.novig_open_exposure() == pytest.approx(9.9975)


def test_a_started_settled_market_frees_its_cost(setup):
    desk, client, store, clock, *_ = setup
    desk.place(request())
    clock["now"] = NOW + timedelta(hours=4)
    client.market_reply = {"marketId": MARKET, "status": "SETTLED"}
    desk.refresh()
    assert store.open_novig_orders() == [] and store.novig_open_exposure() == 0


def test_two_orders_on_a_settled_market_free_on_one_read(setup):
    desk, client, store, clock, status, _ = setup
    desk.place(request())
    held = [held_id_of(ORDER_ID)]
    second = request(clientId="1b3f5f20-9c1b-4d6e-8f3a-5c2d1e7b4a91", heldIds=held)
    client.place_reply = (201, {"orderId": "8d8c7b6a-5f4e-4d3c-2b1a-0f9e8d7c6b5b"})
    client.fill_items = [{"fillId": "f2", "qty": 100, "cost": "0.46500", "taker": True}]
    assert desk.place(second)["status"] == "partial"
    clock["now"] = NOW + timedelta(hours=4)
    status["value"] = {name: {**source, "fetchedAt": clock["now"].isoformat()} for name, source in fresh_status().items()}
    client.market_reply = {"marketId": MARKET, "status": "SETTLED"}
    client.market_reads = 0
    desk.refresh()
    assert store.open_novig_orders() == [] and client.market_reads == 1


def test_an_order_novig_never_shows_final_is_settled_off_its_listing(setup):
    desk, client, store, clock, *_ = setup
    client.order_reply = (404, {"code": "NOT_FOUND"})
    assert desk.place(request())["status"] == "sent"
    clock["now"] = NOW + timedelta(seconds=30)
    desk.refresh()
    assert store.open_novig_orders()[0]["status"] == "sent"  # too young to look up
    clock["now"] = NOW + timedelta(seconds=61)
    client.lookup = {"orderId": ORDER_ID, "clientId": CLIENT_ID, "status": "CANCELED"}
    client.fill_items = [{"fillId": "f1", "qty": 1000, "cost": "4.65000", "taker": True}]
    desk.refresh()
    assert store.open_novig_orders()[0]["status"] == "partial"
    clock["now"] = NOW + timedelta(seconds=125)
    desk.refresh()
    order = store.open_novig_orders()[0]
    assert order["fills_final"] and order["open_cost"] == pytest.approx(4.65)


def test_a_duplicate_client_id_is_never_sent_twice(setup):
    desk, client, *_ = setup
    desk.place(request())
    assert desk.place(request())["status"] == "duplicate" and len(client.placed) == 1


@pytest.mark.parametrize("change, field", [({"price": 0.467}, "price"), ({"qty": 0}, "qty"), ({"clientId": "x"}, "clientId"),
                                           ({"line": {"league": "nfl"}}, "betType")])
def test_validation_names_the_field(change, field):
    body = {"clientId": CLIENT_ID, "outcomeId": OUTCOME, "marketId": MARKET, "rowKey": "r1", "price": 0.465, "qty": 10,
            "eventId": 1, "heldIds": [], "eventStart": NOW.isoformat(), "decidedAt": NOW.isoformat(), "novigBaseUrl": BASE,
            "line": {"league": "nfl", "betType": "spread", "period": "FG", "side": "home", "points": -3.5}}
    body.update(change)
    assert field in validate_order_request(body)


class ScopeClient:
    def __init__(self, management, balance=500.0):
        self._management, self._balance = management, balance

    def balance(self):
        if isinstance(self._balance, Exception):
            raise self._balance
        return self._balance

    def is_management_key(self):
        return self._management


def test_arms_only_on_a_trading_key_that_reads_its_balance():
    assert account_problem(ScopeClient(False)) == (None, 500.0)
    assert "management key" in account_problem(ScopeClient(True))[0]
    assert account_problem(ScopeClient(None))[0]
    assert "balance" in account_problem(ScopeClient(False, desk_module.NovigError("HTTP 401")))[0]


def _post(port, path, body, host=None):
    import json
    import urllib.error
    import urllib.request
    headers = {"Content-Type": "application/json"}
    if host:
        headers["Host"] = host
    req = urllib.request.Request(f"http://127.0.0.1:{port}{path}", data=json.dumps(body).encode(), method="POST",
                                 headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


def test_the_service_says_why_the_novig_desk_is_off_and_routes_orders_on_loopback_only(tmp_path, monkeypatch):
    import threading
    from http.server import ThreadingHTTPServer
    from unabated_ticket.bets_service import config, service
    from unabated_ticket.bets_service.store import BetsStore
    taker_store = TakerStore(tmp_path / "taker.duckdb")
    monkeypatch.setattr(config, "NOVIG_TAKER_KEY_ID", "")
    desk, reason = service.open_novig_taker_desk(taker_store, dict)
    assert desk is None and "no bot Novig key" in reason
    monkeypatch.setattr(config, "NOVIG_TAKER_KEY_ID", "8f14e45f-ceea-467a-9b1c-3f2a51c8d7e0")
    monkeypatch.setattr(config, "NOVIG_TAKER_PRIVATE_KEY_PATH", str(tmp_path / "missing.pem"))
    assert "no bot Novig key file" in service.open_novig_taker_desk(taker_store, dict)[1]
    from unabated_ticket.bets_service.novig_subaccount_setup import write_new_private_key
    write_new_private_key(tmp_path / "novig_bot.pem")
    monkeypatch.setattr(config, "NOVIG_TAKER_PRIVATE_KEY_PATH", str(tmp_path / "novig_bot.pem"))
    refused = service.open_novig_taker_desk(taker_store, dict, check=lambda _client: ("it is a management key", None))
    assert refused[0] is None and "management key" in refused[1]
    armed, _ = service.open_novig_taker_desk(taker_store, dict, check=lambda _client: (None, 500.0))
    assert isinstance(armed, NovigTakerDesk)
    bets_store = BetsStore(tmp_path / "bets.duckdb", 7)
    assert [source.name for source in service.novig_subaccount_sources(taker_store, bets_store)] == ["novig_bot"]

    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(
        bets_store, 0.0, ["kalshi"], extra_hosts=("box.ts.net",), taker_store=taker_store,
        novig_taker_desk=None, novig_taker_desk_off_reason=reason))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    port = server.server_address[1]
    try:
        assert _post(port, "/taker/novig_order.json", {}) == (503, {"status": "off", "error": reason})
        status, body = _post(port, "/taker/novig_order.json", {}, host="box.ts.net")
        assert status == 403 and "loopback-only" in body["error"]
    finally:
        server.shutdown()
        server.server_close()
        bets_store.close()
        taker_store.close()
