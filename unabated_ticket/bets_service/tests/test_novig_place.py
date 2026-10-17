"""The Ticket's Place button on Novig (novig_place.py) against a scripted
NovigOrderClient and a scripted runner quote: the request shape, the refusals
before anything is sent (started, price over 0.1c worse or off Novig's grid,
other outcome held, a stale Novig read, an earlier Place the stake has not
seen at either venue, a repeat click), the 1c contract count, what each Novig
answer reads as, and the board line kept for Place's sub-account read. No network."""
import json
import threading
import urllib.error
import urllib.request
from datetime import datetime, timedelta, timezone
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.novig_orders import qty_for_stake
from unabated_ticket.bets_service.novig_place import (PLACE_FINAL_READ_ATTEMPTS, ManualNovigDesk,
                                                      validate_novig_place_request)
from unabated_ticket.bets_service.store import BetsStore

NOW = datetime(2026, 10, 16, 21, 0, tzinfo=timezone.utc)
EVENT_ID = 123742
OUTCOME = "01a0a05e-e864-7b73-9b80-9e502e4e98d5"
OTHER_OUTCOME = "01a0a05e-e864-7b73-9b80-9e64b163cc27"
MARKET = "0b0b0b0b-1111-4222-8333-444444444444"
ORDER_ID = "0c0c0c0c-1111-4222-8333-666666666666"
CLICK_ID = "00000000-0000-4000-8000-000000000001"
NOVIG_URL = "https://api.novig.com"
BODY = {"clickId": CLICK_ID, "eventId": EVENT_ID, "periodTypeId": 1, "betTypeId": 2, "sideIndex": 1,
        "points": -35.5, "stake": 50, "ticketAllIn": 0.675, "novigReadAt": "2026-10-16T20:59:30Z",
        "kalshiReadAt": "2026-10-16T20:59:40Z"}
LINE = {"league": "cfb", "betType": "spread", "period": "full game", "side": "home", "points": -35.5,
        "awayTeam": "Duquesne", "homeTeam": "Washington St", "awayTeamId": "880", "homeTeamId": "717"}
QUOTE = {"ok": True, "outcomeId": OUTCOME, "otherOutcomeId": OTHER_OUTCOME, "marketId": MARKET, "price": 0.675,
         "depth": 9000, "eventId": EVENT_ID, "eventStart": "2026-10-17T19:30:00.000Z",
         "label": "Washington St -35.5 full game", "novigBaseUrl": NOVIG_URL, "line": LINE}


class ScriptedNovig:
    """NovigOrderClient's order calls, scripted. By default the IOC is taken
    (201), reads FILLED and fills every contract asked for."""
    base_url = NOVIG_URL

    def __init__(self, send=None, order_reads=None, fills=None, raises=False, cancel_raises=False):
        self.calls = []
        self.cancel_raises = cancel_raises
        self.send = send
        self.order_reads = list(order_reads or [])
        self.fill_items = fills
        self.raises = raises
        self.cancels = []

    def place_ioc_order(self, outcome_id, price, qty, client_id):
        self.calls.append((outcome_id, price, qty, client_id))
        self.last_qty = qty
        if self.raises:
            raise OSError("timed out")
        return self.send or (201, {"orderId": ORDER_ID})

    def get_order(self, order_id):
        if self.order_reads:
            read = self.order_reads.pop(0)
            if isinstance(read, Exception):
                raise read
            return read
        return 200, {"orderId": order_id, "status": "FILLED"}

    def cancel_order(self, order_id):
        self.cancels.append(order_id)
        if self.cancel_raises:
            raise OSError("reset")
        return 200, {}

    def fills(self, order_id=None):
        if self.fill_items is not None:
            return self.fill_items
        return [{"fillId": "f1", "orderId": order_id, "qty": self.last_qty, "cost": f"{self.last_qty * 0.675 / 100:.5f}"}]


@pytest.fixture
def store(tmp_path):
    held = BetsStore(tmp_path / "bets.duckdb", 0)
    for source in ("novig", "novig_place"):
        held.log_source_run(source, NOW - timedelta(seconds=40), NOW - timedelta(seconds=30), True, None, 3)
    yield held
    held.close()


def desk_with(store, novig, quote=None, now=NOW):
    quotes = []

    def reader(line):
        quotes.append(line)
        return dict(QUOTE if quote is None else quote)

    desk = ManualNovigDesk(novig, store, reader, now=lambda: now, sleep=lambda _sec: None)
    desk.quotes = quotes
    return desk


def request(**changes):
    return validate_novig_place_request({**BODY, **changes})


def test_validate_names_what_is_wrong():
    assert validate_novig_place_request({**BODY, "clickId": "abc12345-0001"}).startswith("clickId")
    assert validate_novig_place_request({**BODY, "betTypeId": 1}).startswith("betTypeId")
    assert validate_novig_place_request({**BODY, "points": None}).startswith("points")
    assert validate_novig_place_request({**BODY, "ticketAllIn": 1.0}).startswith("ticketAllIn")
    assert validate_novig_place_request({**BODY, "novigReadAt": None}).startswith("novigReadAt")
    assert validate_novig_place_request({**BODY, "kalshiReadAt": "yesterday"}).startswith("kalshiReadAt")
    assert validate_novig_place_request({**BODY, "kalshiReadAt": None})["kalshi_read_at"] is None
    assert validate_novig_place_request([BODY]).startswith("body")
    assert validate_novig_place_request({**BODY, "eventId": True}).startswith("eventId")
    assert validate_novig_place_request({**BODY, "stake": 0}).startswith("stake")
    assert validate_novig_place_request({**BODY, "points": float("inf")}).startswith("points")
    assert validate_novig_place_request({**BODY, "ticketAllIn": 0}).startswith("ticketAllIn")
    assert validate_novig_place_request({**BODY, "clickId": CLICK_ID.upper()})["click_id"] == CLICK_ID.upper()


def test_a_click_sends_one_ioc_buy_at_the_live_price_in_1c_contracts(store):
    novig = ScriptedNovig()
    desk = desk_with(store, novig)
    result = desk.place(request())
    assert desk.quotes == [{"eventId": EVENT_ID, "periodTypeId": 1, "betTypeId": 2, "sideIndex": 1, "points": -35.5}]
    assert novig.calls == [(OUTCOME, 0.675, 7407, CLICK_ID)]
    assert result["status"] == "filled"
    assert result["cost"] == pytest.approx(49.99725)
    assert result["message"] == "Bet $50.00 at 67.5¢ on Novig to win $24.07."


def test_the_board_line_is_kept_so_places_subaccount_read_names_the_bet(store):
    desk_with(store, ScriptedNovig()).place(request())
    contexts = store.novig_place_order_contexts()
    assert contexts == {ORDER_ID: {
        "league": "cfb", "bet_type": "spread", "period": "full game", "side": "home", "points": -35.5,
        "away_team": "Duquesne", "home_team": "Washington St", "away_team_id": "880", "home_team_id": "717",
        "event_start": datetime(2026, 10, 17, 19, 30, tzinfo=timezone.utc), "event_id": EVENT_ID, "market_id": MARKET,
        "outcome_id": OUTCOME, "price": 0.675, "settled_at": None}}
    first = store.note_novig_place_settled(ORDER_ID, NOW)
    again = store.note_novig_place_settled(ORDER_ID, NOW + timedelta(hours=1))
    assert first == again == NOW, "only the first SETTLED sighting is kept"
    assert store.novig_place_order_contexts()[ORDER_ID]["settled_at"] == NOW


def test_the_same_click_is_never_sent_twice(store):
    novig = ScriptedNovig()
    desk_with(store, novig).place(request())
    again = desk_with(store, novig).place(request(novigReadAt="2026-10-16T21:30:00Z", kalshiReadAt=None))
    assert again["status"] == "refused" and "already sent" in again["message"]
    assert len(novig.calls) == 1


@pytest.mark.parametrize("quote_changes, words", [
    ({"eventStart": "2026-10-16T20:00:00Z"}, "has started"),
    ({"eventStart": None}, "start time is unknown"),
    ({"eventStart": "2026-10-17T19:30:00"}, "start time is unknown"),
    ({"price": 0.68}, "more than 0.1"),
    ({"ok": False, "reason": "market_mismatch", "detail": "Novig market x: event status \"LIVE\""}, "event status"),
    ({"novigBaseUrl": "https://api.qa.novig.com"}, "this desk trades"),
    ({"price": None}, "is not a price"),
    ({"price": 0.672}, "on its grid"),
    ({"otherOutcomeId": None}, "both its outcomes"),
    ({"marketId": None}, "both its outcomes"),
    ({"line": None}, "line must be the board line"),
    ({"line": {**LINE, "side": "sideways"}}, "line.side must be"),
    ({"line": {**LINE, "league": ""}}, "line.league must be"),
    ({"line": {**LINE, "period": None}}, "line.period must be"),
    ({"line": {**LINE, "points": "x"}}, "line.points must be"),
])
def test_refusals_send_nothing(store, quote_changes, words):
    novig = ScriptedNovig()
    result = desk_with(store, novig, {**QUOTE, **quote_changes}).place(request())
    assert result["status"] == "refused" and words in result["message"]
    assert novig.calls == []


def test_a_stake_that_buys_no_contract_is_refused(store):
    novig = ScriptedNovig()
    result = desk_with(store, novig).place(request(stake=0.005))
    assert result["status"] == "refused" and "buys no contract" in result["message"]
    assert novig.calls == []


def test_a_price_within_a_tenth_of_a_cent_worse_or_any_better_places(store):
    novig = ScriptedNovig()
    assert desk_with(store, novig).place(request(ticketAllIn=0.6745))["status"] == "filled"
    assert desk_with(store, novig, {**QUOTE, "price": 0.6}).place(
        request(clickId="00000000-0000-4000-8000-000000000002", novigReadAt=(NOW + timedelta(seconds=1)).isoformat()))["status"] == "filled"
    assert novig.calls[1][2] == qty_for_stake(50, 0.6)


def test_the_other_outcome_held_is_refused_but_a_parlay_leg_on_it_is_not(store):
    leg = {"id": "novig:card:0", "venue": "novig", "status": "open", "isParlayLeg": True,
           "venueIds": {"outcomeId": OTHER_OUTCOME}}
    store.upsert_bets([leg], NOW)
    assert desk_with(store, ScriptedNovig()).place(request())["status"] == "filled"
    straight = {"id": "novig:o-1", "venue": "novig", "status": "open", "venueIds": {"outcomeId": OTHER_OUTCOME}}
    store.upsert_bets([straight], NOW)
    novig = ScriptedNovig()
    result = desk_with(store, novig).place(request(clickId="00000000-0000-4000-8000-000000000003", novigReadAt=NOW.isoformat()))
    assert result["status"] == "refused" and "other side" in result["message"]
    assert novig.calls == []


def test_a_failing_or_old_novig_read_is_refused(store):
    novig = ScriptedNovig()
    old = desk_with(store, novig, now=NOW + timedelta(minutes=6)).place(request())
    assert old["status"] == "refused" and "over 5 minutes old" in old["message"]
    store.log_source_run("novig_place", NOW, NOW, False, "HTTP 401", 0)
    failing = desk_with(store, novig).place(request())
    assert failing["status"] == "refused" and "novig_place read is failing" in failing["message"]
    # First run: a counted source that has never polled (the bot's, just armed) refuses too.
    never = ManualNovigDesk(novig, store, lambda _line: dict(QUOTE), read_sources=("novig", "novig_bot"), now=lambda: NOW)
    assert "novig_bot read is failing" in never.place(request())["message"]
    assert novig.calls == []


def test_an_add_waits_until_the_last_place_on_the_game_was_read_at_either_venue(store):
    novig = ScriptedNovig()
    desk_with(store, novig).place(request())
    stale = desk_with(store, novig, now=NOW + timedelta(seconds=20)).place(request(clickId="00000000-0000-4000-8000-000000000009"))
    assert stale["status"] == "refused" and "not on the Bets tab yet" in stale["message"]
    fresh = desk_with(store, novig, now=NOW + timedelta(seconds=90)).place(
        request(clickId="00000000-0000-4000-8000-000000000010", novigReadAt=(NOW + timedelta(seconds=60)).isoformat()))
    assert fresh["status"] == "filled"
    # A Kalshi Place on the game after the Kalshi read the stake knew blocks a Novig Place too.
    store.insert_manual_order({
        "client_order_id": "uman-k1", "placed_at": NOW + timedelta(seconds=100), "venue": "kalshi", "event_id": EVENT_ID,
        "ticker": "KX", "contract_side": "yes", "label": None, "ticket_all_in": 0.5, "ask_price": 0.5,
        "all_in_price": 0.51, "stake": 10, "count_requested": 19, "status": "filled", "updated_at": NOW})
    blocked = desk_with(store, novig, now=NOW + timedelta(seconds=120)).place(
        request(clickId="00000000-0000-4000-8000-000000000011", novigReadAt=(NOW + timedelta(seconds=110)).isoformat()))
    assert blocked["status"] == "refused"
    unseen_by_no_kalshi_read = desk_with(store, novig, now=NOW + timedelta(seconds=120)).place(
        request(clickId="00000000-0000-4000-8000-000000000012", novigReadAt=(NOW + timedelta(seconds=110)).isoformat(), kalshiReadAt=None))
    assert unseen_by_no_kalshi_read["status"] == "refused"


@pytest.mark.parametrize("novig, status, words", [
    (ScriptedNovig(send=(451, {"code": "GEOLOCATION_EXPIRED"})), "rejected", "Open the Novig app on your phone"),
    (ScriptedNovig(send=(422, {"code": "INSUFFICIENT_BALANCE"})), "rejected", "INSUFFICIENT_BALANCE"),
    (ScriptedNovig(send=(500, {})), "unknown", "Check Novig"),
    (ScriptedNovig(send=(201, {})), "unknown", "Check Novig"),
    (ScriptedNovig(raises=True), "unknown", "Novig did not answer"),
    (ScriptedNovig(order_reads=[(404, {})] * PLACE_FINAL_READ_ATTEMPTS, fills=[]), "unknown", "has not said what filled"),
    (ScriptedNovig(order_reads=[OSError("reset")] * PLACE_FINAL_READ_ATTEMPTS, fills=[]), "unknown",
     "has not said what filled"),
    # Never read final, but every contract filled: it cannot fill more, so it is filled.
    (ScriptedNovig(order_reads=[(404, {})] * PLACE_FINAL_READ_ATTEMPTS), "filled", "Bet $"),
    # Never read final, some contracts filled: what filled so far, still unknown.
    (ScriptedNovig(order_reads=[(200, {"status": "PENDING"})] * PLACE_FINAL_READ_ATTEMPTS,
                   fills=[{"fillId": "f1", "qty": 4000, "cost": "27.00000"}]), "unknown", "Bet $27.00 at 67.5¢ on Novig so far"),
    (ScriptedNovig(order_reads=[(404, {})] * PLACE_FINAL_READ_ATTEMPTS, fills=[{"fillId": "f1", "qty": 100}]), "unknown",
     "has not said what filled"),
    (ScriptedNovig(order_reads=[(200, {"status": "REJECTED"})], fills=[]), "unfilled", "price was gone"),
    (ScriptedNovig(order_reads=[(200, {"status": "CANCELED"})], fills=[{"fillId": "f1", "qty": 4000, "cost": "27.00000"}]),
     "partial", "Bet $27.00 of $50.00"),
    (ScriptedNovig(fills=[{"fillId": "f1", "qty": 100}]), "unknown", "fills could not be read"),
])
def test_what_each_answer_reads_as(store, novig, status, words):
    answer = desk_with(store, novig).place(request())
    assert answer["status"] == status and words in answer["message"]
    assert len(novig.calls) == 1, "an order is sent once, never retried"


def test_an_order_left_resting_is_cancelled_and_reads_unknown(store):
    novig = ScriptedNovig(order_reads=[(200, {"status": "PENDING"}), (200, {"status": "OPEN"})])
    answer = desk_with(store, novig).place(request())
    assert answer["status"] == "unknown" and "cancel was sent" in answer["message"]
    assert novig.cancels == [ORDER_ID]
    lost_cancel = ScriptedNovig(order_reads=[(200, {"status": "OPEN"})], cancel_raises=True)
    assert desk_with(store, lost_cancel).place(request(clickId="00000000-0000-4000-8000-000000000030", eventId=EVENT_ID + 1))["status"] == "unknown"


def test_a_rejected_send_leaves_the_game_open_but_an_unknown_one_blocks_the_next_place(store):
    desk_with(store, ScriptedNovig(send=(422, {"code": "INSUFFICIENT_BALANCE"}))).place(request())
    later = {"novigReadAt": (NOW - timedelta(seconds=30)).isoformat()}
    assert desk_with(store, ScriptedNovig()).place(
        request(clickId="00000000-0000-4000-8000-000000000020", **later))["status"] == "filled"
    other_game = {"eventId": EVENT_ID + 1, **later}
    desk_with(store, ScriptedNovig(raises=True)).place(request(clickId="00000000-0000-4000-8000-000000000021", **other_game))
    blocked = desk_with(store, ScriptedNovig()).place(request(clickId="00000000-0000-4000-8000-000000000022", **other_game))
    assert blocked["status"] == "refused" and "not on the Bets tab yet" in blocked["message"]
    # A lost send (no Novig order id) blocks the game however late the read: no read can ever have seen it.
    read_later = NOW + timedelta(seconds=90)
    much_later = request(clickId="00000000-0000-4000-8000-000000000023", eventId=EVENT_ID + 1,
                         novigReadAt=(NOW + timedelta(seconds=60)).isoformat())
    assert desk_with(store, ScriptedNovig(), now=read_later).place(much_later)["status"] == "refused"
    # Until the novig_place read finds it on Novig by clientId and names it (or finds it nowhere).
    lost = store.lost_novig_place_orders(read_later)
    assert lost == [{"client_order_id": "nman-00000000-0000-4000-8000-000000000021", "outcome_id": OUTCOME}]
    store.adopt_novig_place_order(lost[0]["client_order_id"], ORDER_ID, read_later)
    assert store.lost_novig_place_orders(read_later) == []
    assert store.novig_place_order_contexts()[ORDER_ID]["event_id"] == EVENT_ID + 1
    assert desk_with(store, ScriptedNovig(), now=read_later).place(much_later)["status"] == "filled"
    # One Novig lists nowhere never reached Novig: rejected, and the game is free.
    desk_with(store, ScriptedNovig(raises=True)).place(request(clickId="00000000-0000-4000-8000-000000000024", eventId=EVENT_ID + 2))
    store.adopt_novig_place_order("nman-00000000-0000-4000-8000-000000000024", None, NOW)
    assert store.manual_orders_maybe_filled_since(EVENT_ID + 2, NOW - timedelta(hours=1), venue="novig") == []


def test_a_novig_unfilled_counts_as_maybe_filled_until_the_next_read(store):
    """Novig's fills feed can lag the order read, so an 'unfilled' may hold contracts."""
    desk_with(store, ScriptedNovig(order_reads=[(200, {"status": "REJECTED"})], fills=[])).place(request())
    soon = desk_with(store, ScriptedNovig(), now=NOW + timedelta(seconds=20)).place(request(clickId="00000000-0000-4000-8000-000000000040"))
    assert soon["status"] == "refused" and "not on the Bets tab yet" in soon["message"]
    read = desk_with(store, ScriptedNovig(), now=NOW + timedelta(seconds=90)).place(
        request(clickId="00000000-0000-4000-8000-000000000041", novigReadAt=(NOW + timedelta(seconds=60)).isoformat()))
    assert read["status"] == "filled"


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
    assert status == 503 and "no Place Novig key" in reply["error"]


def test_the_service_arms_place_only_on_its_own_proven_trading_key_with_both_reads(store, tmp_path, monkeypatch):
    from unabated_ticket.bets_service import config
    from unabated_ticket.bets_service.novig_subaccount_setup import write_new_private_key
    monkeypatch.setattr(config, "NOVIG_PLACE_KEY_ID", "")
    assert "no Place Novig key" in service.open_novig_place_desk(store, ["novig", "novig_place"])[1]
    write_new_private_key(tmp_path / "novig_place.pem")
    monkeypatch.setattr(config, "NOVIG_PLACE_KEY_ID", "8f14e45f-ceea-467a-9b1c-3f2a51c8d7e1")
    monkeypatch.setattr(config, "NOVIG_PLACE_PRIVATE_KEY_PATH", str(tmp_path / "novig_place.pem"))
    trading = lambda _client: (None, 250.0)  # noqa: E731
    assert "does not read novig_place" in service.open_novig_place_desk(store, ["novig"], check=trading)[1]
    refused = service.open_novig_place_desk(store, ["novig", "novig_place"], check=lambda _client: ("management key", None))
    assert refused[0] is None and "management key" in refused[1]
    armed, _ = service.open_novig_place_desk(store, ["novig", "novig_place", "novig_bot"], check=trading)
    assert isinstance(armed, ManualNovigDesk)
    assert armed._read_sources == ("novig", "novig_place", "novig_bot"), "the bot's bets count when they are read"
    sources = service.novig_subaccount_sources(None, store)
    assert [source.name for source in sources] == ["novig_place"]


def test_an_order_never_read_final_stays_unknown_until_the_sub_account_read_confirms_it(store):
    novig = ScriptedNovig(order_reads=[(404, {})] * PLACE_FINAL_READ_ATTEMPTS, fills=[])
    assert desk_with(store, novig).place(request())["status"] == "unknown"
    client_order_id = f"nman-{CLICK_ID}"
    assert store.unconfirmed_novig_place_orders(NOW - timedelta(seconds=1)) == [], "placed after the cutoff: not yet"
    [unconfirmed] = store.unconfirmed_novig_place_orders(NOW)
    assert unconfirmed == {"client_order_id": client_order_id, "novig_order_id": ORDER_ID, "count_requested": novig.last_qty}
    assert store.unconfirmed_novig_place_orders(NOW + timedelta(hours=25)) == [], "re-read for a day, then left"
    store.confirm_novig_place_order(client_order_id, "filled", novig.last_qty, 50.0, NOW + timedelta(minutes=1))
    assert store.unconfirmed_novig_place_orders(NOW + timedelta(minutes=1)) == []


def test_order_reads_stop_at_the_window_however_slow_each_read_is(store):
    """Each read hanging to its 10 s timeout: reads stop once the window has
    passed, so the click answers within the phone's wait."""
    elapsed = [0.0]

    class HangingNovig(ScriptedNovig):
        def get_order(self, order_id):
            elapsed[0] += 10.0
            return 404, {}

    novig = HangingNovig(fills=[])
    desk = ManualNovigDesk(novig, store, lambda _line: dict(QUOTE), now=lambda: NOW, sleep=lambda sec: None,
                           clock=lambda: elapsed[0])
    assert desk.place(request())["status"] == "unknown"
    assert elapsed[0] == 20.0, "two reads: the second ends past the 12 s window, so no third starts"
