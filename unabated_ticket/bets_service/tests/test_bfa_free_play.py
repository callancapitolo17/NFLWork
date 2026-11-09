"""BFA free play (bfa_free_play.py) on both Place buttons: the balance read, each bet type's
free play rules, the refusals before anything is sent, the UseFreePlay picks, confirmation on a
"(FP)" ticket only, the manual_orders row, and GET /bfa_free_play.json."""
import copy
import json
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import bfa_free_play, bfa_place, bfa_teaser, service
from unabated_ticket.bets_service.bfa_free_play import (
    FreePlayReader, free_play_problem, free_play_rules_of, read_free_play_balance)
from unabated_ticket.bets_service.bfa_teaser import BFATeaserPlacer, validate_place_request
from unabated_ticket.bets_service.store import BetsStore
from unabated_ticket.bets_service.tests import test_bfa_place as straight
from unabated_ticket.bets_service.tests import test_bfa_teaser as teaser
from unabated_ticket.bets_service.tests.test_bfa_place import store  # noqa: F401 — the shared fixture
from unabated_ticket.bets_service.tests.test_bfa_teaser import PLAYER_ID, FakeBFA, FakeResponse

# Free play limits as the account's metadata carried them on 2026-11-09, on every type.
FREE_PLAY_FIELDS = {"freePlay": True, "fpOddsLimit": 800, "fpMaxPayout": 1000.0}


def with_free_play(metadata: dict, **fields) -> dict:
    """The metadata with the free play limits on every betType row (changed by `fields`)."""
    copied = copy.deepcopy(metadata)
    for row in copied["betType"]:
        row.update({**FREE_PLAY_FIELDS, **fields})
    return copied


class FreePlaySession:
    """Wraps a test session: answers BFA's balance and the metadata with free play limits."""

    def __init__(self, inner, metadata: dict, free_play_amount: object = 530.0):
        self.inner = inner
        self.metadata = metadata
        self.free_play_amount = free_play_amount
        self.balance_reads = 0

    @property
    def posts(self):
        return self.inner.posts

    def get(self, url, params=None, headers=None, timeout=None):
        if url == bfa_free_play.BALANCE_URL:
            assert params == {"playerId": PLAYER_ID} and headers["Authorization"] == "Bearer token"
            self.balance_reads += 1
            return FakeResponse(200, {"playerId": 315172, "freePlayAmount": self.free_play_amount,
                                      "availableBalance": 1200.0})
        if url == bfa_teaser.PLAYER_METADATA_URL:
            return FakeResponse(200, self.metadata)
        return self.inner.get(url, params, headers, timeout)

    def post(self, url, params=None, json=None, headers=None, timeout=None):
        return self.inner.post(url, params, json, headers, timeout)


# ---- the rules -------------------------------------------------------------------------

def test_the_balance_is_free_play_amount_and_anything_else_raises():
    assert read_free_play_balance(FreePlaySession(None, {}), {"Authorization": "Bearer token"}, PLAYER_ID) == 530.0
    with pytest.raises(RuntimeError, match="freePlayAmount"):
        read_free_play_balance(FreePlaySession(None, {}, free_play_amount=None), {"Authorization": "Bearer token"},
                               PLAYER_ID)


def test_a_types_rules_need_free_play_on_every_sport_and_take_the_tightest_limits():
    metadata = with_free_play(teaser.METADATA)
    metadata["betType"][-1]["fpMaxPayout"] = 600.0
    assert free_play_rules_of(metadata, 117608, {"NFL", "CFB"}) == {"oddsLimit": 800, "maxPayout": 600.0}
    metadata["betType"][-2]["freePlay"] = False
    assert free_play_rules_of(metadata, 117608, {"NFL", "CFB"}) == "BFA does not offer free play on this bet type"
    assert "not ['NBA']" in free_play_rules_of(metadata, 117608, {"NBA"})
    assert "no free play limits" in free_play_rules_of(with_free_play(teaser.METADATA, fpMaxPayout=None), 117608, {"NFL"})


def test_a_free_play_over_the_balance_the_odds_limit_or_the_payout_cap_refuses():
    rules = {"oddsLimit": 800, "maxPayout": 1000.0}
    assert free_play_problem(rules, 530.0, 100, 90.91, -110) is None
    assert "more than your $530.00 of BFA free play" in free_play_problem(rules, 530.0, 531, 482.73, -110)
    assert "+800 or shorter, this is +900" in free_play_problem(rules, 530.0, 50, 450.0, 900)
    assert free_play_problem(rules, 530.0, 50, 400.0, 800) is None
    assert "win at most $1,000" in free_play_problem(rules, 530.0, 400, 1200.0, None)


# ---- the straight -----------------------------------------------------------------------

def free_play_straight_session(**changes):
    return FreePlaySession(straight.StraightSession(), with_free_play(straight.METADATA), **changes)


GIANTS_FREE_PLAY = {**straight.GIANTS_PLACED, "headerDescription": "STRAIGHT BET (FP)", "riskAmount": 0.0}


def test_a_free_play_straight_sends_use_free_play_and_confirms_on_the_fp_ticket(store):
    session = free_play_straight_session()
    # The cash Giants straight that also shows is not the free play's ticket.
    placer = straight.make_placer(FakeBFA([[], [straight.GIANTS_PLACED, GIANTS_FREE_PLAY]]), session, store)
    result = straight.place(placer, freePlay=True)
    assert result["status"] == "filled" and result["legsMatch"] is True, result
    assert result["message"] == ("Bet $100.00 of free play on New York Giants +2.5 at -110 to win $90.91 "
                                 "(BFA ticket 9001).")
    assert session.balance_reads == 1
    [[wager]] = session.posts
    assert [pick["UseFreePlay"] for pick in wager["Picks"]] == [True]
    assert (wager["RiskWin"], wager["Picks"][0]["RiskOrWinType"], wager["Picks"][0]["Amount"]) == (2, 2, 100.0)
    row = store._con.execute("SELECT free_play, count_filled, cost_filled, status FROM manual_orders").fetchone()
    assert row == (True, 100.0, 0.0, "filled")


def test_a_cash_straight_never_confirms_on_a_free_play_ticket(store):
    session = straight.StraightSession()
    result = straight.place(straight.make_placer(FakeBFA([[], [GIANTS_FREE_PLAY]]), session, store))
    assert result["status"] == "unknown"
    [[wager]] = session.posts
    assert wager["Picks"][0]["UseFreePlay"] is False
    row = store._con.execute("SELECT free_play FROM manual_orders").fetchone()
    assert row == (False,)


@pytest.mark.parametrize("changes, metadata_fields, words", [
    ({"stake": 600}, {}, "more than your $530.00 of BFA free play"),
    ({}, {"freePlay": False}, "does not offer free play"),
    ({}, {"fpMaxPayout": 50.0}, "win at most $50"),
])
def test_a_free_play_straight_bfa_would_refuse_sends_nothing(store, changes, metadata_fields, words):
    session = FreePlaySession(straight.StraightSession(), with_free_play(straight.METADATA, **metadata_fields))
    result = straight.place(straight.make_placer(FakeBFA([[]]), session, store), freePlay=True, **changes)
    assert result["status"] == "refused" and words in result["message"], result
    assert session.posts == [] and straight.manual_rows(store) == []


def test_a_held_free_play_on_the_other_side_blocks_a_place(store):
    cardinals = {**straight.open_straight(8000, straight.GAME_ID,
                                          " [257] ARIZONA CARDINALS -2½-110 [Sport:Football][League:NFL]"),
                 "headerDescription": "STRAIGHT BET (FP)"}
    result = straight.place(straight.make_placer(FakeBFA([[cardinals]]), straight.StraightSession(), store))
    assert "other side of this line" in result["message"]


def test_freeplay_must_be_true_or_false():
    assert "freePlay must be true or false" in bfa_place.validate_bfa_place_request({**straight.REQUEST, "freePlay": 1})
    assert bfa_place.validate_bfa_place_request(straight.REQUEST)["free_play"] is False
    assert "freePlay must be true or false" in validate_place_request({**teaser.REQUEST, "freePlay": "yes"})
    assert validate_place_request(teaser.REQUEST)["freePlay"] is False


# ---- the teaser -------------------------------------------------------------------------

def teaser_placer(bfa, session):
    clock = teaser.Clock()
    return BFATeaserPlacer(bfa, session_factory=lambda: session, clock=clock, sleep=clock.sleep,
                           new_transaction_id=lambda: teaser.TRANSACTION_ID)


def free_play_teaser(ticket):
    placed = teaser.open_teaser(ticket, teaser.PLACED_LEGS, risk=0.0)
    return {**placed, "headerDescription": "4 TEAM TEASERS (FP)"}


def test_a_free_play_teaser_sends_use_free_play_on_every_pick_and_confirms_on_the_fp_ticket():
    session = FreePlaySession(teaser.FakeSession(), with_free_play(teaser.METADATA))
    cash = teaser.open_teaser(356323495, teaser.PLACED_LEGS)
    bfa = FakeBFA([[], [], [cash, free_play_teaser(356323496)]])
    result = teaser_placer(bfa, session).place(validate_place_request({**teaser.REQUEST, "freePlay": True}))
    assert result["status"] == "placed" and result["legsMatch"] is True, result
    assert result["message"] == "Placed · free play · ticket 356323496 · all 4 legs match BFA"
    [[wager]] = session.posts
    expected = copy.deepcopy(teaser.RECORDED_BODY[0])
    for pick in expected["Picks"]:
        pick["UseFreePlay"] = True
    assert wager == expected


def test_a_free_play_teaser_over_the_balance_or_the_payout_cap_sends_nothing():
    over_balance = FreePlaySession(teaser.FakeSession(), with_free_play(teaser.METADATA), free_play_amount=150.0)
    result = teaser_placer(FakeBFA([[]]), over_balance).place(validate_place_request({**teaser.REQUEST, "freePlay": True}))
    assert "more than your $150.00 of BFA free play" in result["message"] and over_balance.posts == []
    # $400 at +300 wins $1,200, over the $1,000 cap (the type's own $500 limit allows it).
    capped = FreePlaySession(teaser.FakeSession(), with_free_play(teaser.METADATA), free_play_amount=900.0)
    result = teaser_placer(FakeBFA([[]]), capped).place(
        validate_place_request({**teaser.REQUEST, "stake": 400, "freePlay": True}))
    assert "win at most $1,000" in result["message"] and capped.posts == []


def test_a_cash_teaser_does_not_read_the_balance_or_confirm_on_a_free_play_ticket():
    session = FreePlaySession(teaser.FakeSession(), with_free_play(teaser.METADATA))
    clock = teaser.Clock()
    placer = BFATeaserPlacer(FakeBFA([[], [free_play_teaser(356323496)]]), session_factory=lambda: session,
                             clock=clock, sleep=clock.sleep, new_transaction_id=lambda: teaser.TRANSACTION_ID)
    result = placer.place(validate_place_request(teaser.REQUEST))
    assert result["status"] == "unconfirmed" and session.balance_reads == 0


def test_a_free_play_tickets_risk_may_read_zero_or_the_stake():
    expected = teaser.expected_for()
    placed = free_play_teaser(356323496)
    assert bfa_teaser.leg_mismatches(placed, expected, 200, free_play=True) == []
    assert bfa_teaser.leg_mismatches({**placed, "riskAmount": 200.0}, expected, 200, free_play=True) == []
    assert bfa_teaser.leg_mismatches(placed, expected, 200) != [], "a cash ticket risking 0 is wrong"


# ---- GET /bfa_free_play.json ------------------------------------------------------------

def test_the_reader_holds_one_read_for_30_seconds_and_says_when_it_fails():
    now = [0.0]
    session = free_play_straight_session()
    reader = FreePlayReader(FakeBFA([[]]), lambda: session, bfa_teaser.read_account_metadata, clock=lambda: now[0])
    assert reader.read() == {"ok": True, "freePlay": 530.0, "maxPayout": 1000.0, "oddsLimit": 800}
    now[0] = 29.0
    reader.read()
    assert session.balance_reads == 1
    now[0] = 31.0
    session.free_play_amount = "oops"
    failed = reader.read()
    assert failed["ok"] is False and "freePlayAmount" in failed["error"]


class StubFreePlay:
    def free_play(self) -> dict:
        return {"ok": True, "freePlay": 530.0, "maxPayout": 1000.0, "oddsLimit": 800}


@pytest.mark.parametrize("placer, status", [(StubFreePlay(), 200), (None, 503)])
def test_the_route_serves_the_free_play_read(tmp_path, placer, status):
    bets_store = BetsStore(tmp_path / "served.duckdb", 7)
    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(
        bets_store, 0.0, ["bfa"], extra_hosts=(), bfa_placer=placer))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        url = f"http://127.0.0.1:{server.server_address[1]}/bfa_free_play.json"
        try:
            with urllib.request.urlopen(url, timeout=5) as response:
                got = (response.status, json.loads(response.read()))
        except urllib.error.HTTPError as error:
            got = (error.code, json.loads(error.read()))
    finally:
        server.shutdown()
        server.server_close()
        bets_store.close()
    assert got[0] == status
    if status == 200:
        assert got[1]["freePlay"] == 530.0
