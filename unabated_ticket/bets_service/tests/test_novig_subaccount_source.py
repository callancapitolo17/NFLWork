"""sources/novig_subaccount.py: a subaccount's fills -> one record per order,
named by the board line the desk kept, graded off the public market once
settled, ids apart from Cal's app account ("novig:<order>:bot")."""
from datetime import datetime, timedelta, timezone

import pytest

from unabated_ticket.bets_service.novig_orders import NovigError
from unabated_ticket.bets_service.sources.novig_subaccount import NovigSubaccountSource, settled_status

NOW = datetime(2026, 10, 18, 16, 0, tzinfo=timezone.utc)
MARKET = "6a1b2c3d-4e5f-4a6b-8c7d-9e0f1a2b3c4d"
OUTCOME = "3f2504e0-4f89-11d3-9a0c-0305e82c3301"
CONTEXT = {"league": "nfl", "bet_type": "spread", "period": "FG", "side": "home", "points": -3.5, "away_team": "Bears",
           "home_team": "Packers", "away_team_id": "11", "home_team_id": "12", "event_start": NOW + timedelta(hours=3),
           "event_id": 4242, "market_id": MARKET, "outcome_id": OUTCOME, "price": 0.465}


class FakeClient:
    def __init__(self, fills, market=None):
        self._fills, self._market, self.market_reads, self.fill_queries = fills, market, 0, []

    def fills(self, order_id=None, outcome_id=None, starts_after_ms=None):
        self.fill_queries.append(starts_after_ms)
        return self._fills

    def market(self, market_id):
        self.market_reads += 1
        return self._market


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
    assert client.market_reads == 0  # before the start nothing can have settled
    assert client.fill_queries[0] == int((NOW - timedelta(days=14)).timestamp() * 1000)


def test_a_settled_market_grades_it_and_is_read_once():
    started = {**CONTEXT, "event_start": NOW - timedelta(hours=4), "settled_at": NOW - timedelta(minutes=2)}
    market = {"marketId": MARKET, "eventId": "ev-1", "status": "SETTLED", "outcomes": [{"outcomeId": OUTCOME, "status": "WIN"}]}
    client = FakeClient([fill()], market)
    source = NovigSubaccountSource("bot", client, lambda: {"o1": started}, now=lambda: NOW)
    assert source.fetch()[0]["status"] == "won"
    source.fetch()
    assert client.market_reads == 1


def test_a_settled_bet_keeps_one_closed_at_poll_after_poll():
    desk_saw_it = NOW - timedelta(hours=1)
    started = {**CONTEXT, "event_start": NOW - timedelta(hours=4), "settled_at": desk_saw_it}
    market = {"marketId": MARKET, "eventId": "ev-1", "status": "SETTLED", "outcomes": [{"outcomeId": OUTCOME, "status": "LOSS"}]}
    clock = {"now": NOW}
    source = NovigSubaccountSource("bot", FakeClient([fill()], market), lambda: {"o1": started}, now=lambda: clock["now"])
    first = source.fetch()[0]
    clock["now"] = NOW + timedelta(minutes=5)
    assert first == {**source.fetch()[0], "sourceFetchedAt": first["sourceFetchedAt"]}
    assert first["closedAt"] == "2026-10-18T15:00:00.000Z"
    # Before the desk has stored its settle time a desk order stays open; with no context the first
    # poll that saw it settled is the closedAt, and it holds.
    unseen = NovigSubaccountSource("bot", FakeClient([fill()], market), lambda: {"o1": {**started, "settled_at": None}},
                                   now=lambda: clock["now"])
    assert (unseen.fetch()[0]["status"], unseen.fetch()[0]["closedAt"]) == ("open", None)
    place = NovigSubaccountSource("place", FakeClient([fill()], market), lambda: {}, now=lambda: clock["now"])
    assert place.fetch()[0]["closedAt"] == "2026-10-18T16:05:00.000Z"
    clock["now"] = NOW + timedelta(minutes=10)
    assert place.fetch()[0]["closedAt"] == "2026-10-18T16:05:00.000Z"


def test_place_keeps_its_first_settled_sighting_through_note_settled():
    """Place has no desk sweep: its read stores the settle time itself, and a
    later poll (or a restart) reads the stored one back through the context."""
    market = {"marketId": MARKET, "eventId": "ev-1", "status": "SETTLED", "outcomes": [{"outcomeId": OUTCOME, "status": "WIN"}]}
    kept = {}

    def note_settled(order_id, at):
        kept.setdefault(order_id, at)
        return kept[order_id]

    contexts = lambda: {"o1": {**CONTEXT, "event_start": NOW - timedelta(hours=4), "settled_at": kept.get("o1")}}  # noqa: E731
    clock = {"now": NOW}
    source = NovigSubaccountSource("place", FakeClient([fill()], market), contexts, now=lambda: clock["now"],
                                   note_settled=note_settled)
    first = source.fetch()[0]
    assert (first["id"], first["status"], first["closedAt"], first["homeTeamVenue"]["unabatedId"]) == (
        "novig:o1:place", "won", "2026-10-18T16:00:00.000Z", "12")
    clock["now"] = NOW + timedelta(minutes=5)
    restarted = NovigSubaccountSource("place", FakeClient([fill()], market), contexts, now=lambda: clock["now"],
                                      note_settled=note_settled)
    assert restarted.fetch()[0]["closedAt"] == "2026-10-18T16:00:00.000Z"


def test_no_context_is_recorded_unmatchable():
    source = NovigSubaccountSource("place", FakeClient([fill()], {"marketId": MARKET, "status": "OPEN"}), lambda: {}, now=lambda: NOW)
    [record] = source.fetch()
    assert record["id"] == "novig:o1:place" and record["unmatchable"] and record["status"] == "open"


def test_an_unreadable_settlement_or_fill_fails_the_run():
    market = {"marketId": MARKET, "status": "SETTLED", "outcomes": [{"outcomeId": OUTCOME, "status": "WEIRD"}]}
    started = {**CONTEXT, "event_start": NOW - timedelta(hours=4), "settled_at": NOW - timedelta(minutes=2)}
    with pytest.raises(NovigError):
        NovigSubaccountSource("bot", FakeClient([fill()], market), lambda: {"o1": started}, now=lambda: NOW).fetch()
    with pytest.raises(NovigError):
        NovigSubaccountSource("bot", FakeClient([{**fill(), "orderId": None}]), lambda: {}, now=lambda: NOW).fetch()


def test_settled_status_reads_keywords_and_fair_value():
    assert settled_status("LOSS", 4.65, 1000) == "lost" and settled_status("PUSH", 4.65, 1000) == "push"
    assert settled_status("TBD", 4.65, 1000) is None
    assert settled_status("0.731", 4.65, 1000) == "won" and settled_status("0.200", 4.65, 1000) == "lost"


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


def test_a_market_no_context_names_is_re_read_every_ten_minutes_not_every_poll():
    market = {"marketId": MARKET, "eventId": "ev-1", "status": "OPEN", "outcomes": [{"outcomeId": OUTCOME, "status": "TBD"}]}
    client = FakeClient([fill()], market)
    clock = {"now": NOW}
    source = NovigSubaccountSource("place", client, dict, now=lambda: clock["now"])
    source.fetch()
    clock["now"] = NOW + timedelta(seconds=30)
    source.fetch()
    assert client.market_reads == 1
    clock["now"] = NOW + timedelta(minutes=11)
    source.fetch()
    assert client.market_reads == 2
