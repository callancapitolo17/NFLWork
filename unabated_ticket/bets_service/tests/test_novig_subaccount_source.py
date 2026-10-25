"""sources/novig_subaccount.py: a subaccount's fills -> one record per order,
named by the board line the desk kept, graded off the subaccount's own open
positions and SETTLEMENT ledger, ids apart from Cal's app account ("novig:<order>:bot")."""
from datetime import datetime, timedelta, timezone

import pytest

from unabated_ticket.bets_service.novig_orders import NovigError
from unabated_ticket.bets_service.sources.novig_subaccount import NovigSubaccountSource, normalize_order

NOW = datetime(2026, 10, 18, 16, 0, tzinfo=timezone.utc)
MARKET = "6a1b2c3d-4e5f-4a6b-8c7d-9e0f1a2b3c4d"
OUTCOME = "3f2504e0-4f89-11d3-9a0c-0305e82c3301"
CONTEXT = {"league": "nfl", "bet_type": "spread", "period": "FG", "side": "home", "points": -3.5, "away_team": "Bears",
           "home_team": "Packers", "away_team_id": "11", "home_team_id": "12", "event_start": NOW + timedelta(hours=3),
           "event_id": 4242, "market_id": MARKET, "outcome_id": OUTCOME, "price": 0.465}


class FakeClient:
    """The position stays open and nothing is paid unless a test says otherwise."""

    def __init__(self, fills, positions=None, settlement_rows=()):
        self._fills, self.fill_queries = fills, []
        self.positions = [{"marketId": MARKET}] if positions is None else positions
        self.settlement_rows, self.ledger_reads = list(settlement_rows), 0

    def fills(self, order_id=None, outcome_id=None, starts_after_ms=None):
        self.fill_queries.append(starts_after_ms)
        return self._fills

    def open_positions(self):
        self.ledger_reads += 1
        if isinstance(self.positions, Exception):
            raise self.positions
        return list(self.positions)

    def settlements(self, starts_after_ms):
        return list(self.settlement_rows)


def fill(order_id="o1", qty=1000, cost="4.65000", fill_id="f1", ts=1760800000000):
    return {"fillId": fill_id, "orderId": order_id, "marketId": MARKET, "outcomeId": OUTCOME, "qty": qty, "cost": cost,
            "taker": True, "ts": ts}


def test_an_order_is_one_record_named_by_its_board_line():
    client = FakeClient([fill(), fill(qty=1150, cost="5.34750", fill_id="f2")])
    source = NovigSubaccountSource("bot", client, lambda: {"o1": CONTEXT}, now=lambda: NOW)
    [record] = source.fetch()
    assert source.name == "novig_bot" and record["id"] == "novig:o1:bot" and record["venue"] == "novig"
    assert (record["stake"], record["contracts"], record["toWin"]) == (pytest.approx(9.9975, abs=0.01), 21.5, pytest.approx(11.5, abs=0.01))
    assert record["price"] == 115  # 0.465
    assert (record["league"], record["betType"], record["side"], record["points"], record["status"]) == ("nfl", "spread", "home", -3.5, "open")
    assert record["homeTeamVenue"]["unabatedId"] == "12" and record["venueIds"]["outcomeId"] == OUTCOME
    assert client.fill_queries[0] == int((NOW - timedelta(days=14)).timestamp() * 1000)


def test_no_context_is_recorded_unmatchable():
    source = NovigSubaccountSource("place", FakeClient([fill()]), lambda: {}, now=lambda: NOW)
    [record] = source.fetch()
    assert record["id"] == "novig:o1:place" and record["unmatchable"] and record["status"] == "open"


def test_a_lost_place_order_is_adopted_by_client_id_a_minute_on_or_marked_never_sent():
    class LookupClient(FakeClient):
        def __init__(self, fills, listed):
            super().__init__(fills)
            self.listed, self.lookups = listed, []

        def find_order_by_client_id(self, outcome_id, client_id):
            self.lookups.append((outcome_id, client_id))
            if isinstance(self.listed, Exception):
                raise self.listed
            return self.listed

    adopted = []
    lost = [{"client_order_id": "nman-00000000-0000-4000-8000-000000000001", "outcome_id": OUTCOME}]
    found = LookupClient([], {"orderId": "o9", "status": "FILLED", "clientId": "00000000-0000-4000-8000-000000000001"})
    source = NovigSubaccountSource("place", found, dict, now=lambda: NOW, lost_orders=lambda before: list(lost) if before <= NOW else [],
                                   adopt_order=lambda client_order_id, order_id: adopted.append((client_order_id, order_id)))
    source.fetch()
    assert found.lookups == [(OUTCOME, "00000000-0000-4000-8000-000000000001")]
    assert adopted == [("nman-00000000-0000-4000-8000-000000000001", "o9")]
    nowhere = LookupClient([], None)
    NovigSubaccountSource("place", nowhere, dict, now=lambda: NOW, lost_orders=lambda _before: list(lost),
                          adopt_order=lambda client_order_id, order_id: adopted.append((client_order_id, order_id))).fetch()
    assert adopted[-1] == ("nman-00000000-0000-4000-8000-000000000001", None)
    unread = LookupClient([], NovigError("HTTP 500"))
    NovigSubaccountSource("place", unread, dict, now=lambda: NOW, lost_orders=lambda _before: list(lost),
                          adopt_order=lambda *_args: adopted.append("never")).fetch()
    assert "never" not in adopted, "a failed listing read leaves the order lost for the next poll"


def test_an_unconfirmed_place_order_is_confirmed_once_novig_shows_it_final():
    class OrderClient(FakeClient):
        """The poll's own fills read sees `poll_fills`; the order's fills, read after its final status, `order_fills`."""
        def __init__(self, order_read, order_fills, listed=None):
            super().__init__([])
            self.order_read, self.order_fills, self.listed, self.order_fill_reads = order_read, order_fills, listed, []

        def fills(self, order_id=None, outcome_id=None, starts_after_ms=None):
            if order_id is None:
                return super().fills(starts_after_ms=starts_after_ms)
            self.order_fill_reads.append(order_id)
            return self.order_fills

        def get_order(self, order_id):
            if isinstance(self.order_read, Exception):
                raise self.order_read
            return self.order_read

        def find_order_by_client_id(self, outcome_id, client_id):
            assert (outcome_id, client_id) == (OUTCOME, "00000000-0000-4000-8000-000000000001")
            return self.listed

    def confirmed_by(client, requested=1000):
        confirmed = []
        unconfirmed = [{"client_order_id": "nman-00000000-0000-4000-8000-000000000001", "novig_order_id": "o1",
                        "outcome_id": OUTCOME, "count_requested": requested}]
        NovigSubaccountSource("place", client, lambda: {"o1": CONTEXT}, now=lambda: NOW,
                              unconfirmed_orders=lambda before: list(unconfirmed) if before <= NOW - timedelta(seconds=59) else [],
                              confirm_order=lambda *args: confirmed.append(args[1:])).fetch()
        return confirmed

    filled = OrderClient((200, {"status": "FILLED"}), [fill(qty=1000, cost="4.65000")])
    assert confirmed_by(filled) == [("filled", 1000, 4.65)]
    assert filled.order_fill_reads == ["o1"], "the order's fills are read after its final status"
    assert confirmed_by(OrderClient((200, {"status": "CANCELED"}), [fill(qty=400, cost="1.86000")])) == [("partial", 400, 1.86)]
    assert confirmed_by(OrderClient((200, {"status": "CANCELED"}), [])) == [("unfilled", 0, 0.0)]
    # Novig 404s the order: its listing by clientId decides; another order id or none listed waits.
    assert confirmed_by(OrderClient((404, {}), [fill()], listed={"orderId": "o1", "status": "FILLED"})) == [("filled", 1000, 4.65)]
    assert confirmed_by(OrderClient((404, {}), [fill()], listed={"orderId": "o2", "status": "FILLED"})) == []
    assert confirmed_by(OrderClient((404, {}), [fill()], listed=None)) == []
    # Not final yet, a FILLED order whose fills have not all landed, a failed read: left for the next poll.
    assert confirmed_by(OrderClient((200, {"status": "PENDING"}), [fill()])) == []
    assert confirmed_by(OrderClient((200, {"status": "FILLED"}), [fill(qty=400)])) == []
    assert confirmed_by(OrderClient(OSError("reset"), [fill()])) == []


def settlement(amount, ts=1792303897282, market=MARKET):
    return {"transactionId": "t1", "kind": "SETTLEMENT", "amount": amount, "ref": market, "ts": ts}


def place_source(client, clock, kept):
    def note_settled(order_id, at):
        kept.setdefault(order_id, at)
        return kept[order_id]

    def contexts():
        started = {**CONTEXT, "event_start": NOW - timedelta(hours=4)}
        return {"o1": {**started, "settled_at": kept.get("o1")}, "o2": {**started, "settled_at": kept.get("o2")}}

    return NovigSubaccountSource("place", client, contexts, now=lambda: clock["now"], note_settled=note_settled)


def test_a_settlement_row_grades_each_order_on_its_market_by_what_it_paid():
    # Two orders on one market: 10 and 5 dollars of contracts; the market paid this subaccount $15.
    fills = [fill(), fill(order_id="o2", qty=500, cost="2.32500", fill_id="f2")]
    client = FakeClient(fills, [], [settlement("15.00000")])
    clock, kept = {"now": NOW}, {}
    first, second = place_source(client, clock, kept).fetch()
    expected = normalize_order("place", "o1", [fill()], {**CONTEXT, "event_start": NOW - timedelta(hours=4)},
                               first["sourceFetchedAt"])
    assert first == {**expected, "status": "won", "closedAt": "2026-10-18T06:11:37.282Z", "toWin": 5.35}
    assert (second["status"], second["toWin"], second["closedAt"]) == ("won", 2.68, "2026-10-18T06:11:37.282Z")
    assert kept == {"o1": NOW, "o2": NOW}  # kept, so a restart whose first read fails leaves them as stored


def test_a_fair_value_settlement_pays_part_and_a_refund_is_a_push():
    clock = {"now": NOW}
    partial = place_source(FakeClient([fill()], [], [settlement("6.00000")]), clock, {}).fetch()[0]
    assert (partial["status"], partial["toWin"]) == ("won", 1.35)
    short = place_source(FakeClient([fill()], [], [settlement("2.00000")]), clock, {}).fetch()[0]
    assert (short["status"], short["toWin"]) == ("won", -2.65)  # P&L -2.65, not the full stake
    refund = place_source(FakeClient([fill()], [], [settlement("4.65000")]), clock, {}).fetch()[0]
    assert refund["status"] == "push"


def test_a_closed_position_with_no_payout_is_lost_once_confirmed_at_its_first_sighting():
    client = FakeClient([fill()], [], [])
    clock, kept = {"now": NOW}, {}
    source = place_source(client, clock, kept)
    assert source.fetch()[0]["status"] == "open" and kept == {"o1": NOW}  # first sighting kept
    clock["now"] = NOW + timedelta(minutes=14)
    assert source.fetch() == []  # closed, not yet a loss: the stored row stands
    clock["now"] = NOW + timedelta(minutes=15)
    lost = source.fetch()[0]
    assert (lost["status"], lost["closedAt"]) == ("lost", "2026-10-18T16:00:00.000Z")
    # A restart reads the kept sighting back and grades at once, with the same closedAt.
    again = place_source(FakeClient([fill()], [], []), clock, kept).fetch()[0]
    assert (again["status"], again["closedAt"]) == ("lost", "2026-10-18T16:00:00.000Z")
    # A settlement row landing after the close still wins it.
    late = place_source(FakeClient([fill()], [], [settlement("10.00000")]), clock, kept).fetch()[0]
    assert late["status"] == "won"


def test_a_loss_first_seen_long_after_its_game_is_dated_at_the_game_end_bound():
    # 2026-10-19: the first grading read after a deploy saw Oct 17-18 losses closed days late.
    client = FakeClient([fill()], [], [])
    long_ago = {**CONTEXT, "event_start": NOW - timedelta(days=2), "settled_at": NOW - timedelta(hours=1)}
    source = NovigSubaccountSource("place", client, lambda: {"o1": long_ago}, now=lambda: NOW,
                                   note_settled=lambda order_id, at: at)
    lost = source.fetch()[0]
    assert (lost["status"], lost["closedAt"]) == ("lost", "2026-10-16T20:00:00.000Z")  # start + 4 h


def test_an_open_position_a_read_from_before_the_start_or_a_bot_order_the_desk_has_not_closed_stays_open():
    clock, kept = {"now": NOW}, {}
    held = place_source(FakeClient([fill()], [{"marketId": MARKET}], []), clock, kept).fetch()[0]
    assert (held["status"], kept) == ("open", {})
    # The bot's desk keeps its settle time itself: until it has, the order stays open.
    started = {**CONTEXT, "event_start": NOW - timedelta(hours=4), "settled_at": None}
    bot = NovigSubaccountSource("bot", FakeClient([fill()], [], []), lambda: {"o1": started}, now=lambda: clock["now"])
    assert bot.fetch()[0]["status"] == "open"
    # A positions read taken before the start (cached for LEDGER_RECHECK) says nothing of a later close.
    early = {"now": NOW}
    soon = {**CONTEXT, "event_start": NOW + timedelta(minutes=1)}
    source = NovigSubaccountSource("place", FakeClient([fill()], [], []), lambda: {"o1": soon},
                                   now=lambda: early["now"], note_settled=lambda order_id, at: at)
    source.fetch()
    early["now"] = NOW + timedelta(minutes=1, seconds=30)
    assert source.fetch()[0]["status"] == "open"


def test_the_ledger_is_read_every_two_minutes_and_a_failed_read_keeps_the_last():
    client = FakeClient([fill()], [{"marketId": MARKET}], [])
    clock, kept = {"now": NOW}, {}
    source = place_source(client, clock, kept)
    source.fetch()
    clock["now"] = NOW + timedelta(minutes=1)
    source.fetch()
    assert client.ledger_reads == 1
    client.positions = NovigError("HTTP 429")
    clock["now"] = NOW + timedelta(minutes=2)
    assert source.fetch()[0]["status"] == "open" and kept == {}
    client.positions, client.settlement_rows = [], [settlement("10.00000")]
    clock["now"] = NOW + timedelta(minutes=3)
    assert source.fetch()[0]["status"] == "open"  # a failed read waits LEDGER_RECHECK too
    clock["now"] = NOW + timedelta(minutes=4)
    assert source.fetch()[0]["status"] == "won"


def test_a_restart_whose_first_ledger_read_fails_leaves_graded_bets_as_stored():
    graded = {**CONTEXT, "event_start": NOW - timedelta(hours=4), "settled_at": NOW - timedelta(hours=1)}
    client = FakeClient([fill()], NovigError("HTTP 429"))
    assert NovigSubaccountSource("place", client, lambda: {"o1": graded}, now=lambda: NOW).fetch() == []
    # A started bot order (no settle time kept yet) or one with no context may be stored won: skipped too.
    started = {**CONTEXT, "event_start": NOW - timedelta(hours=4), "settled_at": None}
    assert NovigSubaccountSource("bot", client, lambda: {"o1": started}, now=lambda: NOW).fetch() == []
    assert NovigSubaccountSource("bot", client, lambda: {}, now=lambda: NOW).fetch() == []
    upcoming = {**CONTEXT, "event_start": NOW + timedelta(hours=4)}
    assert NovigSubaccountSource("bot", client, lambda: {"o1": upcoming}, now=lambda: NOW).fetch()[0]["status"] == "open"


def test_an_order_with_no_context_grades_only_off_a_settlement_row():
    unnamed = NovigSubaccountSource("place", FakeClient([fill()], []), lambda: {}, now=lambda: NOW).fetch()[0]
    assert (unnamed["status"], bool(unnamed["unmatchable"])) == ("open", True)
    paid = NovigSubaccountSource("place", FakeClient([fill()], [], [settlement("10.00000")]), lambda: {},
                                 now=lambda: NOW).fetch()[0]
    assert paid["status"] == "won"


def test_an_unreadable_fill_fails_the_run():
    with pytest.raises(NovigError):
        NovigSubaccountSource("bot", FakeClient([{**fill(), "orderId": None}]), lambda: {}, now=lambda: NOW).fetch()


class BalanceClient:
    def __init__(self, reply):
        self._reply = reply

    def balance(self):
        if isinstance(self._reply, Exception):
            raise self._reply
        return self._reply


def test_balance_is_the_subaccounts_own_wallet_with_no_in_bets_split():
    source = NovigSubaccountSource("place", BalanceClient(2446.0684), lambda: {})
    assert source.balance() == {"cash": 2446.0684, "inBets": None, "available": 2446.0684, "creditLimit": None}


def test_a_failed_balance_read_shows_unavailable_under_the_sources_name():
    from unabated_ticket.bets_service.balances import BalanceBook, read_balance_once
    book = BalanceBook()
    read_balance_once(NovigSubaccountSource("bot", BalanceClient(NovigError("HTTP 401")), lambda: {}), book)
    reading = book.snapshot()["novig_bot"]
    assert reading["cash"] is None and "HTTP 401" in reading["error"]
