"""bfa_place.py: the Ticket's request, the account's straight type, the line on BFA's board, the
wager body (the recorded teaser's shape with one pick and the straight type), and the placer's
refusals and confirmation against a fake BFA and a real store."""
import copy
import json
import threading
import urllib.error
import urllib.request
from datetime import datetime, timedelta, timezone
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import bfa_place, bfa_teaser, service
from unabated_ticket.bets_service.bfa_place import (
    BFAStraightPlacer, match_straight, straight_mismatches, straight_type_of, validate_bfa_place_request)
from unabated_ticket.bets_service.store import BetsStore
from unabated_ticket.bets_service.tests.test_bfa_teaser import (
    AGENT_ID, NFL_BOARD, NOW, PLAYER_ID, FakeBFA, FakeResponse, FakeSession)

# The Giants home +2.5 at -110 on the NFL board (game 38431521, rotation 258); Cardinals 257.
GAME_ID = 38431521
START = "2026-10-04T17:00:00Z"
READ_AT = (NOW - timedelta(seconds=30)).isoformat().replace("+00:00", "Z")
STRAIGHT_TYPE = 4633

REQUEST = {"clickId": "click-0001", "eventId": 555, "periodTypeId": 1, "betTypeId": 2, "sideIndex": 1, "points": 2.5,
           "stake": 100, "ticketAllIn": 0.5238095, "bfaReadAt": READ_AT, "kalshiReadAt": READ_AT, "novigReadAt": None,
           "polymarketReadAt": None, "prophetxReadAt": None}


def quote(**changes):
    """The runner's /bfa_quote.json answer for the Giants +2.5."""
    line = {"league": "nfl", "betType": "spread", "period": "FG", "side": "home", "points": 2.5,
            "awayTeam": "Arizona Cardinals", "homeTeam": "New York Giants", "awayTeamId": "1", "homeTeamId": "2"}
    line.update(changes.pop("line", {}))
    answer = {"ok": True, "eventId": 555, "eventStart": START, "label": "New York Giants +2.5", "price": -110,
              "isAlt": False, "rotation": 258, "awayRotation": 257, "homeRotation": 258, "line": line}
    answer.update(changes)
    return answer


METADATA = {"playerId": 315172, "agentId": AGENT_ID, "betType": [
    {"wagerTypeId": STRAIGHT_TYPE, "sportId": "NFL", "numTeams": 1, "wagerTypes": 0, "description": "Straight Bet",
     "maxRisk": 2000.0},
    {"wagerTypeId": 4634, "sportId": "CFB", "numTeams": 1, "wagerTypes": 0, "description": "Straight Bet"},
    {"wagerTypeId": 117608, "sportId": "NFL", "numTeams": 4, "wagerTypes": 2, "description": "4 TEAM TEASERS"},
]}


class StraightSession(FakeSession):
    """FakeSession with the account's straight types in its metadata."""

    def get(self, url, params=None, headers=None, timeout=None):
        if url == bfa_teaser.PLAYER_METADATA_URL:
            self.gets.append(url)
            return FakeResponse(200, METADATA)
        return super().get(url, params, headers, timeout)


def open_straight(ticket, game_id, text, risk=100.0, to_win=90.91):
    return {"idWager": ticket, "ticketNumber": ticket, "headerDescription": "STRAIGHT BET", "riskAmount": risk,
            "winAmount": to_win, "placedDate": "2026-10-03T09:19:08.000",
            "betDetails": [{"idGame": game_id, "detailDescription": text}]}


GIANTS_PLACED = open_straight(9001, GAME_ID, " [258] NEW YORK GIANTS +2½-110 [Sport:Football][League:NFL]")


class Clock:
    def __init__(self):
        self.now = 0.0

    def __call__(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.now += seconds


@pytest.fixture
def store(tmp_path):
    bets_store = BetsStore(tmp_path / "bets.duckdb", 7)
    bets_store.log_source_run("bfa", NOW - timedelta(seconds=31), NOW - timedelta(seconds=30), True, None, 0)
    yield bets_store
    bets_store.close()


def make_placer(bfa, session, store, answer=None, now=NOW):
    clock = Clock()
    reads = []

    def reader(request):
        reads.append(request)
        return answer or quote()
    placer = BFAStraightPlacer(bfa, store, reader, session_factory=lambda: session, now=lambda: now, clock=clock,
                               sleep=clock.sleep, new_transaction_id=lambda: "tx-1")
    placer.reads = reads
    return placer


def place(placer, **changes):
    return placer.place(validate_bfa_place_request({**REQUEST, **changes}))


def manual_rows(store):
    return store._con.execute("SELECT client_order_id, venue, status, stake, ticker, side, error FROM manual_orders").fetchall()


# ---- the request ------------------------------------------------------------------------

def test_a_good_request_validates_and_a_bad_one_names_its_problem():
    request = validate_bfa_place_request(REQUEST)
    assert (request["stake"], request["points"], request["novig_read_at"]) == (100, 2.5, None)
    assert validate_bfa_place_request({**REQUEST, "betTypeId": 1, "points": 7})["points"] is None
    assert "whole dollars" in validate_bfa_place_request({**REQUEST, "stake": 100.5})
    assert "bfaReadAt" in validate_bfa_place_request({**REQUEST, "bfaReadAt": None})
    assert "clickId" in validate_bfa_place_request({**REQUEST, "clickId": "x"})
    assert "points" in validate_bfa_place_request({**REQUEST, "points": None})


# ---- the account and the board ------------------------------------------------------------

def test_the_straight_type_is_the_sports_one_row_with_its_limit():
    assert straight_type_of(METADATA, "NFL") == {"wagerTypeId": STRAIGHT_TYPE, "maxRisk": 2000.0}
    assert straight_type_of(METADATA, "CFB") == {"wagerTypeId": 4634, "maxRisk": None}
    missing = straight_type_of(METADATA, "NBA")
    assert "found none" in missing and "['CFB', 'NFL']" in missing


def test_the_board_line_must_be_the_same_team_side_and_number():
    line = quote()["line"]
    start = datetime(2026, 10, 4, 17, tzinfo=timezone.utc)
    match = match_straight(NFL_BOARD, line, 258, start, NOW, "Giants +2.5")
    assert match["odds"]["price"] == -110 and match["odds"]["dGmId"] == GAME_ID
    assert match_straight(NFL_BOARD, {**line, "points": 3}, 258, start, NOW, "Giants +3") == (
        "Giants +3: BFA has +2.5 now, the Ticket has +3")
    assert "as home, the list has the away team" in match_straight(NFL_BOARD, {**line, "side": "away"}, 258, start, NOW, "x")
    assert match_straight(NFL_BOARD, line, 258, start, start, "x") == "x: the game has started"
    total = match_straight(NFL_BOARD, {**line, "betType": "total", "side": "under", "points": 44.5}, 258, start, NOW, "u")
    assert total["odds"]["side"] == 5
    # The fixture board lists no moneyline market.
    assert "lists no full-game moneyline" in match_straight(NFL_BOARD, {**line, "betType": "moneyline", "points": None},
                                                             258, start, NOW, "x")


# ---- placing ------------------------------------------------------------------------------

def test_a_placed_straight_sends_the_teasers_shape_with_one_pick_and_reads_filled(store):
    session = StraightSession()
    bfa = FakeBFA([[], [], [GIANTS_PLACED]])
    result = place(make_placer(bfa, session, store))
    assert result["status"] == "filled" and result["legsMatch"] is True
    assert result["message"] == "Bet $100.00 on New York Giants +2.5 at -110 to win $90.91 (BFA ticket 9001)."
    assert session.posts == [[{
        "IdTransaction": "tx-1", "IdPlayer": 315172, "FillIdWager": -1, "WagerType": 0, "OpenSpots": 0,
        "Picks": [{"EventId": 3167690, "FixtureId": 18443763, "MarketType": 2, "PeriodNumber": 0, "Side": 1,
                   "Index": 0, "ContestantId": 64424829, "Line": 2.5, "TeaserPoints": 0.0, "Price": -110,
                   "Amount": 100.0, "PointsPurchased": 0.0, "PitcherAction": 0, "RoundRobinCombinations": 0,
                   "UseFreePlay": False, "RiskOrWinType": 2, "IdGameType": 1, "IdSport": "NFL"}],
        "RiskWin": 2, "AcceptChanges": 0, "IdWagerType": STRAIGHT_TYPE, "IsLive": False, "PhoneLine": None,
        "UserId": None}]]
    assert manual_rows(store) == [("bman-click-0001", "bfa", "filled", 100.0, str(GAME_ID), "home", None)]
    assert result["openWagers"] == [GIANTS_PLACED]


def test_a_to_win_a_cent_off_from_the_rounding_still_reads_back_as_sent():
    kansas = open_straight(3, GAME_ID, " [258] NEW YORK GIANTS +2½-110 [Sport:Football][League:NFL]",
                           risk=17.0, to_win=15.46)
    want = {"rotations": {258}, "betType": "spread", "side": "home", "points": 2.5, "price": -110}
    assert straight_mismatches(kansas, want, 17) == []
    # The ticket that booked the wrong size reads as different.
    assert straight_mismatches({**kansas, "riskAmount": 18.7, "winAmount": 17.0}, want, 17) != []


def test_the_same_click_twice_sends_once(store):
    session = StraightSession()
    placer = make_placer(FakeBFA([[], [GIANTS_PLACED]]), session, store)
    assert place(placer)["status"] == "filled"
    again = place(placer)
    assert again["status"] == "refused" and len(session.posts) == 1


@pytest.mark.parametrize("answer, words", [
    (quote(isAlt=True), "is an alt line"),
    (quote(line={"period": "1H"}), "is a 1H line"),
    (quote(line={"league": "soccer"}), "does not read BFA's SOCCER board"),
    (quote(rotation=None), "no rotation number"),
    (quote(eventStart="2026-10-03T16:00:00Z"), "the game has started"),
    ({"ok": False, "reason": "no_line", "detail": "the board holds 0 BFA lines"}, "the board holds 0 BFA lines"),
])
def test_a_line_bfa_place_cannot_bet_sends_nothing(store, answer, words):
    session = StraightSession()
    result = place(make_placer(FakeBFA([[]]), session, store, answer=answer))
    assert result["status"] == "refused" and words in result["message"]
    assert session.posts == [] and manual_rows(store) == []


def test_a_worse_price_the_other_side_held_a_limit_and_a_stale_read_refuse(store):
    session = StraightSession()
    worse = place(make_placer(FakeBFA([[]]), session, store), ticketAllIn=0.51)
    assert "more than 0.1¢ worse" in worse["message"]
    cardinals = open_straight(8000, GAME_ID, " [257] ARIZONA CARDINALS -2½-110 [Sport:Football][League:NFL]")
    held = place(make_placer(FakeBFA([[cardinals]]), session, store))
    assert "other side of this line" in held["message"]
    over = place(make_placer(FakeBFA([[]]), session, store), stake=2500)
    assert "over BFA's $2000 limit" in over["message"]
    stale = place(make_placer(FakeBFA([[]]), session, store, now=NOW + timedelta(minutes=10)))
    assert "over 5 minutes old" in stale["message"]
    assert session.posts == [] and manual_rows(store) == []


def test_an_earlier_place_after_the_bfa_read_refuses_the_next(store):
    session = StraightSession()
    placer = make_placer(FakeBFA([[], [GIANTS_PLACED]]), session, store)
    assert place(placer)["status"] == "filled"
    # The panel's stake is still sized on the read from before that bet.
    second = place(placer, clickId="click-0002")
    assert second["status"] == "refused" and "last Place on this game" in second["message"]
    assert len(session.posts) == 1


def test_a_bfa_straight_counts_against_a_place_at_any_venue_until_the_bfa_read_after_it(store):
    from unabated_ticket.bets_service.manual_orders import unseen_manual_fill, validate_place_request
    assert place(make_placer(FakeBFA([[], [GIANTS_PLACED]]), StraightSession(), store))["status"] == "filled"
    kalshi = validate_place_request({**REQUEST, "stake": 10.0, "ticketAllIn": 0.5})
    assert kalshi["bfa_read_at"] == datetime.fromisoformat(READ_AT.replace("Z", "+00:00"))
    reads = [kalshi["kalshi_read_at"], None, None, None]
    assert unseen_manual_fill(store, 555, *reads, kalshi["bfa_read_at"]) is True
    assert unseen_manual_fill(store, 555, *reads, NOW + timedelta(seconds=1)) is False
    assert unseen_manual_fill(store, 556, *reads, kalshi["bfa_read_at"]) is False, "another game"


def test_a_4xx_is_rejected_and_no_ticket_is_unknown(store):
    refused = place(make_placer(FakeBFA([[]]), StraightSession(wager_status=400, wager_reply={"error": "limit"}), store))
    assert refused["status"] == "rejected" and "HTTP 400" in refused["message"]
    lost = place(make_placer(FakeBFA([[]]), StraightSession(wager_raises=True), store), clickId="click-0002")
    assert lost["status"] == "unknown" and "Check BFA's open bets" in lost["message"]
    assert [row[2] for row in manual_rows(store)] == ["rejected", "unknown"]


def test_a_placed_ticket_that_differs_is_filled_and_says_so(store):
    other_price = open_straight(9002, GAME_ID, " [258] NEW YORK GIANTS +2½-115 [Sport:Football][League:NFL]", to_win=86.96)
    result = place(make_placer(FakeBFA([[], [other_price]]), StraightSession(), store))
    assert result["status"] == "filled" and result["legsMatch"] is False
    assert "BFA has [258] NEW YORK GIANTS +2½-115" in result["message"]
    assert "risks 100.0 to win 86.96" in result["message"]


def test_a_total_and_a_moneyline_read_back_as_sent():
    total = open_straight(1, GAME_ID, " [257] TOTAL u44½-110 \r(ARIZONA CARDINALS vrs NEW YORK GIANTS) [Sport:Football]")
    want = {"rotations": {257, 258}, "betType": "total", "side": "under", "points": 44.5, "price": -110}
    assert straight_mismatches(total, want, 100) == []
    moneyline = open_straight(2, GAME_ID, " [258] NEW YORK GIANTS +120 [Sport:Football][League:NFL]", to_win=120.0)
    want = {"rotations": {258}, "betType": "moneyline", "side": "home", "points": None, "price": 120}
    assert straight_mismatches(moneyline, want, 100) == []
    assert straight_mismatches(moneyline, {**want, "side": "away"}, 100) != []


def test_a_free_play_or_a_parlay_is_not_the_new_straight():
    free_play = {**GIANTS_PLACED, "headerDescription": "STRAIGHT BET (FP)"}
    parlay = {**GIANTS_PLACED, "headerDescription": "PARLAY (2 TEAMS)"}
    assert bfa_place.new_straights_on_game([free_play, parlay, GIANTS_PLACED], GAME_ID, set()) == [GIANTS_PLACED]
    assert bfa_place.new_straights_on_game([GIANTS_PLACED], GAME_ID, {9001}) == []


class StubPlacer:
    def __init__(self, result: dict):
        self.result = result
        self.requests = []

    def place(self, request: dict) -> dict:
        self.requests.append(request)
        return copy.deepcopy(self.result)


@pytest.fixture
def serve_with(tmp_path):
    servers = []

    def start(placer):
        bets_store = BetsStore(tmp_path / "served.duckdb", 7)
        server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(
            bets_store, 0.0, ["bfa"], extra_hosts=(), bfa_placer=placer))
        threading.Thread(target=server.serve_forever, daemon=True).start()
        servers.append((server, bets_store))
        return f"http://127.0.0.1:{server.server_address[1]}", bets_store

    yield start
    for server, bets_store in servers:
        server.shutdown()
        server.server_close()
        bets_store.close()


def post_place(url: str, body: object) -> tuple[int, dict]:
    request = urllib.request.Request(f"{url}/place_bfa.json", data=json.dumps(body).encode(), method="POST",
                                     headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


def test_the_route_places_and_stores_the_open_bets_it_read(serve_with):
    placer = StubPlacer({"status": "filled", "message": "Bet $100.00 ...", "ticket": {"ticketNumber": 9001},
                         "openWagers": [GIANTS_PLACED]})
    url, bets_store = serve_with(placer)
    status, reply = post_place(url, REQUEST)
    assert status == 200
    assert reply == {"ok": True, "status": "filled", "message": "Bet $100.00 ...", "ticket": {"ticketNumber": 9001}}
    assert placer.requests == [validate_bfa_place_request(REQUEST)]
    stored = {record["id"] for record in bets_store.load_bets(30, datetime.now(timezone.utc))}
    assert any(bet_id.startswith("bfa:9001") for bet_id in stored)


def test_the_route_refuses_a_bad_body_and_answers_503_without_a_bfa_account(serve_with):
    url, _ = serve_with(StubPlacer({"status": "filled", "message": ""}))
    status, reply = post_place(url, {**REQUEST, "stake": 10.5})
    assert status == 400 and "whole dollars" in reply["error"]
    url, _ = serve_with(None)
    status, reply = post_place(url, REQUEST)
    assert status == 503 and "BFA_USERNAME" in reply["error"]
