"""The Ticket's Place button on Novig (novig_place.py) against a scripted
order client and a scripted runner quote: the request shape, the refusals
before anything is sent (started, price over 0.1c worse, other outcome held, a
stale Novig read, an earlier Place the stake has not seen at either venue, a
repeat click), the 1c contract count, and what each client answer reads as.
No network."""
import json
import threading
import urllib.error
import urllib.request
from datetime import datetime, timedelta, timezone
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.novig_place import (
    ManualNovigDesk, novig_contracts_for, validate_novig_place_request,
)
from unabated_ticket.bets_service.store import BetsStore

NOW = datetime(2026, 10, 16, 21, 0, tzinfo=timezone.utc)
EVENT_ID = 123742
OUTCOME = "01a0a05e-e864-7b73-9b80-9e502e4e98d5"
OTHER_OUTCOME = "01a0a05e-e864-7b73-9b80-9e64b163cc27"
NOVIG_URL = "https://api.novig.com"
BODY = {"clickId": "n1a2b3c4-0001", "eventId": EVENT_ID, "periodTypeId": 1, "betTypeId": 2, "sideIndex": 1,
        "points": -35.5, "stake": 50, "ticketAllIn": 0.675, "novigReadAt": "2026-10-16T20:59:30Z",
        "kalshiReadAt": "2026-10-16T20:59:40Z"}
QUOTE = {"ok": True, "outcomeId": OUTCOME, "otherOutcomeId": OTHER_OUTCOME, "price": 0.675, "depth": 9000,
         "eventId": EVENT_ID, "eventStart": "2026-10-17T19:30:00.000Z", "label": "Washington St -35.5 full game",
         "novigBaseUrl": NOVIG_URL}


class ScriptedNovig:
    base_url = NOVIG_URL

    def __init__(self, result=None, raises=False):
        self.calls = []
        self.result = result
        self.raises = raises

    def place_ioc_buy(self, outcome_id, price, quantity, client_id):
        self.calls.append((outcome_id, price, quantity, client_id))
        if self.raises:
            raise OSError("timed out")
        if self.result is not None:
            return dict(self.result)
        return {"outcome": "filled", "orderId": "nv-1", "filledQuantity": quantity, "cost": quantity * price / 100}


@pytest.fixture
def store(tmp_path):
    held = BetsStore(tmp_path / "bets.duckdb", 0)
    held.log_source_run("novig", NOW - timedelta(seconds=40), NOW - timedelta(seconds=30), True, None, 3)
    yield held
    held.close()


def desk_with(store, novig, quote=None, now=NOW):
    quotes = []

    def reader(line):
        quotes.append(line)
        return dict(QUOTE if quote is None else quote)

    desk = ManualNovigDesk(novig, store, reader, now=lambda: now)
    desk.quotes = quotes
    return desk


def request(**changes):
    return validate_novig_place_request({**BODY, **changes})


def test_validate_names_what_is_wrong():
    assert validate_novig_place_request({**BODY, "betTypeId": 1}).startswith("betTypeId")
    assert validate_novig_place_request({**BODY, "points": None}).startswith("points")
    assert validate_novig_place_request({**BODY, "ticketAllIn": 1.0}).startswith("ticketAllIn")
    assert validate_novig_place_request({**BODY, "novigReadAt": None}).startswith("novigReadAt")
    assert validate_novig_place_request({**BODY, "kalshiReadAt": "yesterday"}).startswith("kalshiReadAt")
    assert validate_novig_place_request({**BODY, "kalshiReadAt": None})["kalshi_read_at"] is None


def test_a_contract_pays_a_cent_so_50_dollars_at_675_buys_7407():
    assert novig_contracts_for(50, 0.675) == 7407
    assert novig_contracts_for(50, 0.5) == 10000
    assert novig_contracts_for(0.001, 0.5) == 0


def test_a_click_sends_one_ioc_buy_at_the_live_price(store):
    novig = ScriptedNovig()
    desk = desk_with(store, novig)
    result = desk.place(request())
    assert desk.quotes == [{"eventId": EVENT_ID, "periodTypeId": 1, "betTypeId": 2, "sideIndex": 1, "points": -35.5}]
    assert novig.calls == [(OUTCOME, 0.675, 7407, "nman-n1a2b3c4-0001")]
    assert result["status"] == "filled"
    assert result["cost"] == pytest.approx(49.99725)
    assert result["message"] == "Bet $50.00 at 67.5¢ on Novig to win $24.07."


def test_the_same_click_is_never_sent_twice(store):
    novig = ScriptedNovig()
    desk_with(store, novig).place(request())
    again = desk_with(store, novig).place(request(novigReadAt="2026-10-16T21:30:00Z", kalshiReadAt=None))
    assert again["status"] == "refused" and "already sent" in again["message"]
    assert len(novig.calls) == 1


@pytest.mark.parametrize("quote_changes, words", [
    ({"eventStart": "2026-10-16T20:00:00Z"}, "has started"),
    ({"price": 0.6761}, "more than 0.1"),
    ({"ok": False, "reason": "novig_off", "detail": "this runner has no Novig book reader yet"}, "no Novig book reader"),
    ({"novigBaseUrl": "https://api.qa.novig.com"}, "this desk trades"),
    ({"price": None}, "is not a price"),
    ({"otherOutcomeId": None}, "both outcomes"),
])
def test_refusals_send_nothing(store, quote_changes, words):
    novig = ScriptedNovig()
    result = desk_with(store, novig, {**QUOTE, **quote_changes}).place(request())
    assert result["status"] == "refused" and words in result["message"]
    assert novig.calls == []


def test_a_price_within_a_tenth_of_a_cent_worse_or_any_better_places(store):
    novig = ScriptedNovig()
    assert desk_with(store, novig, {**QUOTE, "price": 0.676}).place(request())["status"] == "filled"
    assert desk_with(store, novig, {**QUOTE, "price": 0.6}).place(
        request(clickId="n1a2b3c4-0002", novigReadAt=(NOW + timedelta(seconds=1)).isoformat()))["status"] == "filled"
    assert novig.calls[1][2] == novig_contracts_for(50, 0.6)


def test_the_other_outcome_held_is_refused_but_a_parlay_leg_on_it_is_not(store):
    leg = {"id": "novig:card:0", "venue": "novig", "status": "open", "isParlayLeg": True,
           "venueIds": {"outcomeId": OTHER_OUTCOME}}
    store.upsert_bets([leg], NOW)
    assert desk_with(store, ScriptedNovig()).place(request())["status"] == "filled"
    straight = {"id": "novig:o-1", "venue": "novig", "status": "open", "venueIds": {"outcomeId": OTHER_OUTCOME}}
    store.upsert_bets([straight], NOW)
    novig = ScriptedNovig()
    result = desk_with(store, novig).place(request(clickId="n1a2b3c4-0003", novigReadAt=NOW.isoformat()))
    assert result["status"] == "refused" and "other side" in result["message"]
    assert novig.calls == []


def test_a_failing_or_old_novig_read_is_refused(store):
    novig = ScriptedNovig()
    old = desk_with(store, novig, now=NOW + timedelta(minutes=6)).place(request())
    assert old["status"] == "refused" and "over 5 minutes old" in old["message"]
    store.log_source_run("novig", NOW, NOW, False, "HTTP 401", 0)
    failing = desk_with(store, novig).place(request())
    assert failing["status"] == "refused" and "read is failing" in failing["message"]
    assert novig.calls == []


def test_an_add_waits_until_the_last_place_on_the_game_was_read_at_either_venue(store):
    novig = ScriptedNovig()
    desk_with(store, novig).place(request())
    stale = desk_with(store, novig, now=NOW + timedelta(seconds=20)).place(request(clickId="n1a2b3c4-0009"))
    assert stale["status"] == "refused" and "not on the Bets tab yet" in stale["message"]
    fresh = desk_with(store, novig, now=NOW + timedelta(seconds=90)).place(
        request(clickId="n1a2b3c4-0010", novigReadAt=(NOW + timedelta(seconds=60)).isoformat()))
    assert fresh["status"] == "filled"
    # A Kalshi Place on the game after the Kalshi read the stake knew blocks a Novig Place too.
    store.insert_manual_order({
        "client_order_id": "uman-k1", "placed_at": NOW + timedelta(seconds=100), "venue": "kalshi", "event_id": EVENT_ID,
        "ticker": "KX", "contract_side": "yes", "label": None, "ticket_all_in": 0.5, "ask_price": 0.5,
        "all_in_price": 0.51, "stake": 10, "count_requested": 19, "status": "filled", "updated_at": NOW})
    blocked = desk_with(store, novig, now=NOW + timedelta(seconds=120)).place(
        request(clickId="n1a2b3c4-0011", novigReadAt=(NOW + timedelta(seconds=110)).isoformat()))
    assert blocked["status"] == "refused"
    unseen_by_no_kalshi_read = desk_with(store, novig, now=NOW + timedelta(seconds=120)).place(
        request(clickId="n1a2b3c4-0012", novigReadAt=(NOW + timedelta(seconds=110)).isoformat(), kalshiReadAt=None))
    assert unseen_by_no_kalshi_read["status"] == "refused"


@pytest.mark.parametrize("result, status, words", [
    ({"outcome": "geo_blocked", "detail": "GEOLOCATION_EXPIRED"}, "rejected", "Open the Novig app on your phone"),
    ({"outcome": "rejected", "detail": "insufficient balance"}, "rejected", "insufficient balance"),
    ({"outcome": "unknown", "detail": "lookup found no order"}, "unknown", "Check Novig"),
    ({"outcome": "filled"}, "unknown", "Check Novig"),
    ({"outcome": "unfilled", "orderId": "nv-2", "filledQuantity": 0, "cost": 0}, "unfilled", "price was gone"),
    ({"outcome": "partial", "orderId": "nv-3", "filledQuantity": 4000, "cost": 27.0}, "partial", "Bet $27.00 of $50.00"),
])
def test_what_each_answer_reads_as(store, result, status, words):
    answer = desk_with(store, ScriptedNovig(result=result)).place(request())
    assert answer["status"] == status and words in answer["message"]


def test_no_answer_is_unknown(store):
    answer = desk_with(store, ScriptedNovig(raises=True)).place(request())
    assert answer["status"] == "unknown" and "Check Novig" in answer["message"]


# ---- the route ----
class StubDesk:
    def __init__(self):
        self.requests = []

    def place(self, placed):
        self.requests.append(placed)
        return {"status": "filled", "message": "Bet $50.00 at 67.5¢ on Novig to win $24.07."}


@pytest.fixture
def serve_with(store):
    servers = []

    def start(desk):
        handler_options = {"novig_place_desk": desk} if desk is not None else {}
        server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(
            store, 0.0, ["novig"], extra_hosts=("mac.tail1234.ts.net",), **handler_options))
        threading.Thread(target=server.serve_forever, daemon=True).start()
        servers.append(server)
        return f"http://127.0.0.1:{server.server_address[1]}"

    yield start
    for server in servers:
        server.shutdown()
        server.server_close()


def post_place(url, body):
    sent = urllib.request.Request(f"{url}/place_novig.json", data=json.dumps(body).encode(), method="POST",
                                  headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(sent, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


def test_the_route_places_and_refuses_a_bad_body_or_no_order_client(serve_with):
    desk = StubDesk()
    url = serve_with(desk)
    status, reply = post_place(url, BODY)
    assert status == 200 and reply["ok"] is True and reply["status"] == "filled"
    assert desk.requests == [validate_novig_place_request(BODY)]
    status, reply = post_place(url, {**BODY, "betTypeId": 1})
    assert status == 400 and "betTypeId" in reply["error"]
    status, reply = post_place(serve_with(None), BODY)
    assert status == 503 and "no Novig order client" in reply["error"]
