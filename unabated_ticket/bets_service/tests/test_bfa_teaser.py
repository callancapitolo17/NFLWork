"""bfa_teaser.py: the panel's request, the account's teaser type, each leg on BFA's board, the
wager body against the one Cal placed by hand on 2026-10-03, and the placer's never-twice rules
against a fake BFA."""
import copy
import json
import threading
import urllib.error
import urllib.request
from datetime import datetime, timezone
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import bfa_teaser, service
from unabated_ticket.bets_service.bfa_teaser import (
    BFATeaserPlacer, match_leg, pick_of, teaser_type_of, ticket_key, validate_place_request, wager_body)
from unabated_ticket.bets_service.store import BetsStore

PLAYER_ID = "315172"
AGENT_ID = 28080
TRANSACTION_ID = "ccfa8e18-e4ed-432b-9354-af2f134507a6"
# 2026-10-03 09:19 Pacific, when the recorded ticket went in.
NOW = datetime(2026, 10, 3, 16, 19, tzinfo=timezone.utc)


def odds(contestant_id, line, price, side, sport, game_id):
    return {"index": 0, "contestantId": contestant_id, "line": line, "price": price, "side": side, "status": 0,
            "dSpId": sport, "dGtId": 1, "dGmId": game_id, "dGmTn": None}


def board_game(event_id, fixture_id, away, home, start, away_spread, home_price_away_price, total, sport, game_id):
    """One game of BFA's popular list, the live shape: away contestant side 2, home side 1."""
    (away_id, away_name, away_rot), (home_id, home_name, home_rot) = away, home
    away_price, home_price = home_price_away_price
    return {
        "id": event_id, "name": f"{away_name} vs {home_name}", "type": 1, "isLive": False,
        "fixtures": [{"id": fixture_id, "isMain": True, "isLive": False, "date": start, "contestants": [
            {"id": away_id, "side": 2, "name": away_name, "rotNum": away_rot},
            {"id": home_id, "side": 1, "name": home_name, "rotNum": home_rot}]}],
        "markets": [
            {"id": event_id * 10 + 2, "fixtureId": fixture_id, "periodNumber": 0, "type": 2, "status": 1, "odds": [
                odds(away_id, away_spread, away_price, 2, sport, game_id),
                odds(home_id, -away_spread, home_price, 1, sport, game_id)]},
            {"id": event_id * 10 + 3, "fixtureId": fixture_id, "periodNumber": 0, "type": 3, "status": 1, "odds": [
                odds(None, total, -110, 4, sport, game_id), odds(None, total, -110, 5, sport, game_id)]},
        ],
    }


NFL_BOARD = [
    board_game(3167690, 18443763, (64424828, "Arizona Cardinals", 257), (64424829, "New York Giants", 258),
               "2026-10-04T17:00:00Z", -2.5, (-110, -110), 44.5, "NFL", 38431521),
    board_game(3167698, 18443730, (64424689, "Denver Broncos", 271), (64424690, "San Francisco 49ers", 272),
               "2026-10-04T20:25:00Z", 2.5, (-105, -115), 41.5, "NFL", 38431485),
    board_game(3167701, 18443731, (64424691, "Los Angeles Chargers", 275), (64424692, "Seattle Seahawks", 276),
               "2026-10-04T20:25:00Z", 7, (-105, -115), 47.5, "NFL", 38431490),
]
CFB_BOARD = [
    board_game(3179927, 18548061, (64817450, "Old Dominion", 143), (64817451, "Georgia State", 144),
               "2026-10-03T19:30:00Z", 2.5, (-110, -110), 52.5, "CFB", 38502300),
]

# The ticket as the panel sends it: Buckeye's numbers BEFORE the 6 points.
REQUEST = {"stake": 200, "legs": [
    {"league": "nfl", "betType": "spread", "side": "home", "rotation": 258, "points": 2.5,
     "eventStart": "2026-10-04T17:00:00Z", "label": "Giants +8.5"},
    {"league": "nfl", "betType": "spread", "side": "away", "rotation": 271, "points": 2.5,
     "eventStart": "2026-10-04T20:25:00Z", "label": "Broncos +8.5"},
    {"league": "nfl", "betType": "spread", "side": "home", "rotation": 276, "points": -7,
     "eventStart": "2026-10-04T20:25:00Z", "label": "Seahawks -1"},
    {"league": "cfb", "betType": "spread", "side": "away", "rotation": 143, "points": 2.5,
     "eventStart": "2026-10-03T19:30:00Z", "label": "Old Dominion +8.5"},
]}


def recorded_pick(event_id, fixture_id, side, contestant_id, line, price, sport):
    return {"EventId": event_id, "FixtureId": fixture_id, "MarketType": 2, "PeriodNumber": 0, "Side": side,
            "Index": 0, "ContestantId": contestant_id, "Line": line, "TeaserPoints": 6.0, "Price": price,
            "Amount": 200.0, "PointsPurchased": 0.0, "PitcherAction": 0, "RoundRobinCombinations": 0,
            "UseFreePlay": False, "RiskOrWinType": 1, "IdGameType": 1, "IdSport": sport}


# What the site POSTed for ticket 356323496 (decoded from the recording, 2026-10-03).
RECORDED_BODY = [{
    "IdTransaction": TRANSACTION_ID, "IdPlayer": 315172, "FillIdWager": -1, "WagerType": 2, "OpenSpots": 0,
    "Picks": [
        recorded_pick(3167690, 18443763, 1, 64424829, 2.5, -110, "NFL"),
        recorded_pick(3167698, 18443730, 2, 64424689, 2.5, -105, "NFL"),
        recorded_pick(3167701, 18443731, 1, 64424692, -7.0, -115, "NFL"),
        recorded_pick(3179927, 18548061, 2, 64817450, 2.5, -110, "CFB"),
    ],
    "RiskWin": 1, "AcceptChanges": 0, "IdWagerType": 117608, "IsLive": False, "PhoneLine": None, "UserId": None,
}]


def type_row(type_id, sport, teams, description, points, payout, max_risk=500.0):
    side_field, total_field = {"NFL": ("nflSide", "nflTotal"), "CFB": ("cfbSide", "cfbTotal"),
                               "CBB": ("cbbSide", "cbbTotal"), "NBA": ("nbaSide", "nbaTotal")}[sport]
    return {"wagerTypeId": type_id, "sportId": sport, "numTeams": teams, "wagerTypes": 2, "description": description,
            "maxRisk": max_risk, "payOuts": payout, "nflSide": 0, "nflTotal": 0, "cfbSide": 0, "cfbTotal": 0,
            "cbbSide": 0, "cbbTotal": 0, "nbaSide": 0, "nbaTotal": 0, side_field: points, total_field: points}


# The account's teaser rows as its metadata lists them (2026-10-03), straight bets around them.
METADATA = {"playerId": 315172, "agentId": AGENT_ID, "betType": [
    {"wagerTypeId": 4633, "sportId": "NFL", "numTeams": 1, "wagerTypes": 0, "description": "Straight Bet"},
    type_row(117606, "NFL", 3, "3 TEAM TEASERS", 6, 180.0),
    type_row(117607, "NFL", 4, "4 TEAM SWEETHEART TEASERS", 13, -130.0),
    type_row(117607, "CFB", 4, "4 TEAM SWEETHEART TEASERS", 13, -130.0),
    type_row(117608, "CBB", 4, "4 TEAM TEASERS", 4, 300.0),
    type_row(117608, "CFB", 4, "4 TEAM TEASERS", 6, 300.0),
    type_row(117608, "NFL", 4, "4 TEAM TEASERS", 6, 300.0),
]}


# ---- the panel's request ---------------------------------------------------------------

def test_a_good_request_validates_to_its_clean_shape():
    request = validate_place_request(REQUEST)
    assert request["stake"] == 200 and len(request["legs"]) == 4
    assert request["legs"][2] == {"league": "nfl", "betType": "spread", "side": "home", "rotation": 276,
                                  "points": -7.0, "eventStart": "2026-10-04T20:25:00Z", "label": "Seahawks -1"}


@pytest.mark.parametrize("change, error", [
    (lambda body: body.update(stake=200.5), "stake must be a whole number"),
    (lambda body: body.update(stake=0), "stake must be a whole number"),
    (lambda body: body["legs"].pop(), "legs must be a list of 4 legs, got 3"),
    (lambda body: body["legs"][0].update(league="nba"), "legs[0].league"),
    (lambda body: body["legs"][1].update(side="over"), "legs[1].side must be one of ['away', 'home'] on a spread"),
    (lambda body: body["legs"][1].update(betType="total", side="away"), "legs[1].side must be one of ['over', 'under']"),
    (lambda body: body["legs"][2].update(rotation="276"), "legs[2].rotation"),
    (lambda body: body["legs"][2].update(points=None), "legs[2].points"),
    (lambda body: body["legs"][3].update(eventStart="tomorrow"), "legs[3].eventStart"),
    (lambda body: body["legs"][3].update(eventStart=None), "legs[3].eventStart"),
    (lambda body: body["legs"][3].update(label=""), "legs[3].label"),
])
def test_a_bad_request_names_the_first_problem(change, error):
    body = copy.deepcopy(REQUEST)
    change(body)
    assert error in validate_place_request(body)


def test_ticket_key_ignores_the_order_of_the_legs():
    request = validate_place_request(REQUEST)
    assert ticket_key(request["legs"]) == ticket_key(list(reversed(request["legs"])))


# ---- the account's teaser type ---------------------------------------------------------

def test_the_four_team_six_point_type_covering_both_sports_is_the_accounts_own():
    assert teaser_type_of(METADATA, {"NFL", "CFB"}) == {"wagerTypeId": 117608, "maxRisk": 500.0}
    assert teaser_type_of(METADATA, {"NFL"}) == {"wagerTypeId": 117608, "maxRisk": 500.0}


def test_a_type_paying_other_than_300_or_missing_refuses():
    paying_less = {"betType": [type_row(117608, "NFL", 4, "4 TEAM TEASERS", 6, 250.0)]}
    assert "pays [250.0], the Teasers tab prices +300" in teaser_type_of(paying_less, {"NFL"})
    assert "found none" in teaser_type_of({"betType": []}, {"NFL"})
    assert "no betType list" in teaser_type_of({"playerId": 1}, {"NFL"})


# ---- legs on BFA's board -----------------------------------------------------------------

def test_the_recorded_ticket_builds_the_body_the_site_sent():
    request = validate_place_request(REQUEST)
    boards = {"nfl": NFL_BOARD, "cfb": CFB_BOARD}
    matches = [match_leg(boards[leg["league"]], leg, NOW) for leg in request["legs"]]
    picks = [pick_of(match, request["stake"]) for match in matches]
    assert wager_body(TRANSACTION_ID, PLAYER_ID, 117608, picks) == RECORDED_BODY


def test_a_total_leg_is_its_side_4_or_5_row_with_no_contestant():
    leg = {**REQUEST["legs"][0], "betType": "total", "side": "under", "points": 44.5, "label": "Under 50.5"}
    pick = pick_of(match_leg(NFL_BOARD, leg, NOW), 200)
    assert (pick["MarketType"], pick["Side"], pick["ContestantId"], pick["Line"]) == (3, 5, None, 44.5)


def test_a_moved_number_refuses_and_says_both_numbers():
    leg = {**REQUEST["legs"][0], "points": 3.0}
    assert match_leg(NFL_BOARD, leg, NOW) == "Giants +8.5: BFA has +2.5 now, the list has +3"
    total = {**leg, "betType": "total", "side": "over", "points": 45.0, "label": "Over 39"}
    assert match_leg(NFL_BOARD, total, NOW) == "Over 39: BFA has 44.5 now, the list has 45"


def test_no_game_a_started_game_a_closed_market_and_a_start_far_off_refuse():
    leg = REQUEST["legs"][3]
    assert "lists no game with rotation 999" in match_leg(CFB_BOARD, {**leg, "rotation": 999}, NOW)
    after_kickoff = datetime(2026, 10, 3, 19, 31, tzinfo=timezone.utc)
    assert match_leg(CFB_BOARD, leg, after_kickoff) == "Old Dominion +8.5: the game has started"
    closed = copy.deepcopy(CFB_BOARD)
    closed[0]["markets"][0]["status"] = 2
    assert match_leg(closed, leg, NOW) == "Old Dominion +8.5: BFA has the spread closed"
    next_week = {**leg, "eventStart": "2026-10-10T19:30:00Z"}
    assert "lists no game with rotation 143" in match_leg(CFB_BOARD, next_week, NOW)


def test_an_odds_row_without_a_game_id_and_a_non_game_entry_refuse():
    no_game_id = copy.deepcopy(CFB_BOARD)
    for odds_row in no_game_id[0]["markets"][0]["odds"]:
        del odds_row["dGmId"]
    leg = REQUEST["legs"][3]
    assert match_leg(no_game_id, leg, NOW) == "Old Dominion +8.5: BFA's odds row carries no price, sport, game type or game id"
    not_a_game = copy.deepcopy(CFB_BOARD)
    not_a_game[0]["type"] = 2
    assert "lists no game with rotation 143" in match_leg(not_a_game, leg, NOW)


# ---- the placer against a fake BFA -------------------------------------------------------

class FakeResponse:
    def __init__(self, status_code: int, body: object):
        self.status_code = status_code
        self._body = body
        self.text = str(body)

    def json(self) -> object:
        return self._body


class FakeBFA:
    """The BFA source's two placer-facing methods, its open list scripted per read."""

    def __init__(self, open_reads: list[list[dict]]):
        self.open_reads = open_reads
        self.open_calls = 0

    def authorized(self) -> tuple[dict, str]:
        return {"Authorization": "Bearer token"}, PLAYER_ID

    def open_wagers(self) -> list[dict]:
        reply = self.open_reads[min(self.open_calls, len(self.open_reads) - 1)]
        self.open_calls += 1
        return reply


class FakeSession:
    def __init__(self, wager_reply: object = None, wager_status: int = 200, wager_raises: bool = False):
        self.wager_reply = {"state": 0} if wager_reply is None else wager_reply
        self.wager_status = wager_status
        self.wager_raises = wager_raises
        self.posts: list[dict] = []
        self.gets: list[str] = []

    def get(self, url, params=None, headers=None, timeout=None):
        self.gets.append(url)
        if url == bfa_teaser.PLAYER_METADATA_URL:
            assert params == {"playerId": PLAYER_ID}
            return FakeResponse(200, METADATA)
        if url == bfa_teaser.BOARD_URL.format(slug="nfl"):
            assert params["agentId"] == AGENT_ID and params["fixtureType"] == 1
            return FakeResponse(200, {"games": NFL_BOARD})
        if url == bfa_teaser.BOARD_URL.format(slug="ncaa_f"):
            return FakeResponse(200, {"games": CFB_BOARD})
        raise AssertionError(f"unexpected GET {url}")

    def post(self, url, params=None, json=None, headers=None, timeout=None):
        assert url == bfa_teaser.WAGER_URL and params == {"playerId": PLAYER_ID}
        assert headers["Authorization"] == "Bearer token"
        self.posts.append(json)
        if self.wager_raises:
            raise TimeoutError("read timed out")
        return FakeResponse(self.wager_status, self.wager_reply)


def open_teaser(ticket, legs):
    return {"idWager": ticket, "ticketNumber": ticket, "headerDescription": "4 TEAM TEASERS", "riskAmount": 200.0,
            "winAmount": 600.0, "placedDate": "2026-10-03T09:19:08.000",
            "betDetails": [{"idGame": game_id, "detailDescription": f" [{rotation}] TEAM +8½-110 (B+6)"}
                           for game_id, rotation in legs]}


PLACED_LEGS = [(38431521, 258), (38431485, 271), (38431490, 276), (38502300, 143)]
OLDER_TICKET = open_teaser(356205894, [(38431480, 102), (38431485, 271), (38431517, 261), (38502187, 201)])


class Clock:
    def __init__(self):
        self.now = NOW.timestamp()

    def __call__(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.now += seconds


def make_placer(bfa: FakeBFA, session: FakeSession, clock: Clock) -> BFATeaserPlacer:
    return BFATeaserPlacer(bfa, session_factory=lambda: session, clock=clock, sleep=clock.sleep,
                           new_transaction_id=lambda: TRANSACTION_ID)


def test_placing_posts_the_recorded_body_once_and_confirms_off_the_open_bets():
    placed = open_teaser(356323496, PLACED_LEGS)
    bfa = FakeBFA([[OLDER_TICKET], [OLDER_TICKET], [OLDER_TICKET, placed]])
    session, clock = FakeSession(), Clock()
    result = make_placer(bfa, session, clock).place(validate_place_request(REQUEST))
    assert session.posts == [RECORDED_BODY]
    assert result["status"] == "placed" and result["message"] == "Placed · ticket 356323496"
    assert result["ticket"] == {"ticketNumber": 356323496, "risk": 200.0, "toWin": 600.0,
                                "placedDate": "2026-10-03T09:19:08.000"}
    assert result["openWagers"] == [OLDER_TICKET, placed]


def test_the_same_ticket_never_goes_out_twice():
    placed = open_teaser(356323496, PLACED_LEGS)
    bfa = FakeBFA([[OLDER_TICKET], [OLDER_TICKET, placed]])
    session, clock = FakeSession(), Clock()
    placer = make_placer(bfa, session, clock)
    assert placer.place(validate_place_request(REQUEST))["status"] == "placed"
    again = placer.place(validate_place_request(REQUEST))
    assert again["status"] == "refused" and "was placed" in again["message"] and "356323496" in again["message"]
    assert len(session.posts) == 1


def test_a_ticket_already_open_at_bfa_is_refused_before_any_post():
    bfa = FakeBFA([[open_teaser(356323496, PLACED_LEGS)]])
    session = FakeSession()
    result = make_placer(bfa, session, Clock()).place(validate_place_request(REQUEST))
    assert result == {"status": "refused", "message": "Not placed: this ticket is already open at BFA (ticket 356323496)"}
    assert session.posts == []


def test_a_moved_number_refuses_before_any_post():
    body = copy.deepcopy(REQUEST)
    body["legs"][1]["points"] = 3.0
    session = FakeSession()
    result = make_placer(FakeBFA([[]]), session, Clock()).place(validate_place_request(body))
    assert result == {"status": "refused", "message": "Not placed: Broncos +8.5: BFA has +2.5 now, the list has +3"}
    assert session.posts == []


def test_a_stake_over_the_types_limit_refuses():
    body = {**copy.deepcopy(REQUEST), "stake": 600}
    session = FakeSession()
    result = make_placer(FakeBFA([[]]), session, Clock()).place(validate_place_request(body))
    assert result["status"] == "refused" and "over BFA's $500 limit" in result["message"]
    assert session.posts == []


def test_a_4xx_is_refused_and_frees_the_legs_for_another_try():
    bfa = FakeBFA([[]])
    session, clock = FakeSession(wager_reply={"message": "Line changed"}, wager_status=400), Clock()
    placer = make_placer(bfa, session, clock)
    first = placer.place(validate_place_request(REQUEST))
    assert first["status"] == "refused" and "HTTP 400" in first["message"] and "Line changed" in first["message"]
    second = placer.place(validate_place_request(REQUEST))
    assert second["status"] == "refused" and len(session.posts) == 2  # it went out again: nothing was booked


def test_a_lost_reply_is_unconfirmed_and_holds_the_legs():
    bfa = FakeBFA([[OLDER_TICKET]])
    session, clock = FakeSession(wager_raises=True), Clock()
    placer = make_placer(bfa, session, clock)
    result = placer.place(validate_place_request(REQUEST))
    assert result["status"] == "unconfirmed"
    assert "the request to BFA failed (TimeoutError)" in result["message"] and "Check BFA's open bets" in result["message"]
    again = placer.place(validate_place_request(REQUEST))
    assert again["status"] == "refused" and "no ticket has shown yet" in again["message"]
    assert len(session.posts) == 1
    clock.now += bfa_teaser.RECENT_TICKET_SEC
    assert placer.place(validate_place_request(REQUEST))["status"] == "unconfirmed"  # the hold ends
    assert len(session.posts) == 2


def test_the_metadata_is_read_once_an_hour():
    placed = open_teaser(356323496, PLACED_LEGS)
    bfa = FakeBFA([[OLDER_TICKET], [OLDER_TICKET, placed]])
    session, clock = FakeSession(), Clock()
    placer = make_placer(bfa, session, clock)
    placer.place(validate_place_request(REQUEST))
    body = copy.deepcopy(REQUEST)
    body["legs"][0]["points"] = 3.0  # refused at the board, after the metadata
    placer.place(validate_place_request(body))
    assert session.gets.count(bfa_teaser.PLAYER_METADATA_URL) == 1


# ---- the service route -------------------------------------------------------------------

class StubPlacer:
    def __init__(self, result: dict):
        self.result = result
        self.requests: list[dict] = []

    def place(self, request: dict) -> dict:
        self.requests.append(request)
        return copy.deepcopy(self.result)


@pytest.fixture
def serve_with(tmp_path):
    servers = []

    def start(placer):
        store = BetsStore(tmp_path / "bets.duckdb", 7)
        server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(store, 0.0, ["bfa"], placer))
        threading.Thread(target=server.serve_forever, daemon=True).start()
        servers.append((server, store))
        return f"http://127.0.0.1:{server.server_address[1]}", store

    yield start
    for server, store in servers:
        server.shutdown()
        server.server_close()
        store.close()


def post_place(url: str, body: object) -> tuple[int, dict]:
    request = urllib.request.Request(f"{url}/place_teaser.json", data=json.dumps(body).encode(), method="POST",
                                     headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


def test_the_route_places_and_stores_the_open_bets_it_read(serve_with):
    placed = open_teaser(356323496, PLACED_LEGS)
    placer = StubPlacer({"status": "placed", "message": "Placed · ticket 356323496",
                         "ticket": {"ticketNumber": 356323496}, "openWagers": [placed]})
    url, store = serve_with(placer)
    status, reply = post_place(url, REQUEST)
    assert status == 200
    assert reply == {"ok": True, "status": "placed", "message": "Placed · ticket 356323496",
                     "ticket": {"ticketNumber": 356323496}}
    assert placer.requests == [validate_place_request(REQUEST)]
    stored = {record["id"] for record in store.load_bets(30, datetime.now(timezone.utc))}
    assert {f"bfa:356323496:leg{index}" for index in range(4)} <= stored


def test_the_route_refuses_a_bad_body_and_answers_503_without_a_bfa_account(serve_with):
    url, _store = serve_with(StubPlacer({"status": "placed", "message": ""}))
    status, reply = post_place(url, {**REQUEST, "stake": "200"})
    assert status == 400 and "stake" in reply["error"]
    url, _store = serve_with(None)
    status, reply = post_place(url, REQUEST)
    assert status == 503 and "BFA_USERNAME" in reply["error"]


def test_an_error_before_the_post_reads_not_placed(serve_with):
    class Broken:
        def place(self, request):
            raise RuntimeError("BFA NFL board: HTTP 503")

    url, _store = serve_with(Broken())
    status, reply = post_place(url, REQUEST)
    assert status == 200
    assert reply == {"ok": False, "status": "refused", "message": "Not placed: RuntimeError: BFA NFL board: HTTP 503"}


def test_a_5xx_or_an_odd_reply_is_unconfirmed_and_holds_the_legs():
    session, clock = FakeSession(wager_reply={"error": "busy"}, wager_status=503), Clock()
    placer = make_placer(FakeBFA([[OLDER_TICKET]]), session, clock)
    result = placer.place(validate_place_request(REQUEST))
    assert result["status"] == "unconfirmed" and "BFA answered HTTP 503" in result["message"]
    assert placer.place(validate_place_request(REQUEST))["status"] == "refused"
    assert len(session.posts) == 1


def test_a_failure_after_the_post_is_unconfirmed_never_not_placed():
    class BrokenAfterPost(FakeBFA):
        def open_wagers(self):
            if self.open_calls:
                raise AssertionError("never reached: the confirm read swallows its own errors")
            self.open_calls += 1
            return [OLDER_TICKET]

    session, clock = FakeSession(), Clock()
    placer = make_placer(BrokenAfterPost([]), session, clock)
    placer._confirm = lambda *args: {}["boom"]  # a bug inside the confirmation itself
    result = placer.place(validate_place_request(REQUEST))
    assert result["status"] == "unconfirmed" and "reading the open bets failed (KeyError)" in result["message"]
    assert placer.place(validate_place_request(REQUEST))["status"] == "refused"  # held
    assert len(session.posts) == 1


def test_the_hold_covers_the_same_sides_at_a_moved_number():
    session, clock = FakeSession(wager_raises=True), Clock()
    placer = make_placer(FakeBFA([[OLDER_TICKET]]), session, clock)
    assert placer.place(validate_place_request(REQUEST))["status"] == "unconfirmed"
    moved = copy.deepcopy(REQUEST)
    moved["legs"][0]["points"] = 3.0
    again = placer.place(validate_place_request(moved))
    assert again["status"] == "refused" and "no ticket has shown yet" in again["message"]
    assert len(session.posts) == 1


def test_a_second_ticket_while_one_is_placing_is_refused():
    session = FakeSession()
    placer = make_placer(FakeBFA([[]]), session, Clock())
    placer._placing.acquire()
    result = placer.place(validate_place_request(REQUEST))
    assert result == {"status": "refused", "message": "Not placed: another ticket is being placed at BFA right now"}
    assert session.posts == [] and session.gets == []
