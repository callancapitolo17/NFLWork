"""DraftKings parser on tests/fixtures/bets/draftkings_bets.json (the push body the
extension sends; its first three open bets are Cal's My Bets capture of 2026-10-10, the
rest hand-written for shapes not yet seen): both lists and their merge, the push
validation, the closed-by-absence fallback, and POST /draftkings.json over a real HTTP
server."""
import json
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

import pytest

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.sources import draftkings
from unabated_ticket.bets_service.sources.draftkings import (
    american_of, closed_by_absence, market_of, moment_of, normalize_draftkings, validate_push)
from unabated_ticket.bets_service.store import BetsStore

FIXTURE_PATH = Path(__file__).parents[2] / "tests" / "fixtures" / "bets" / "draftkings_bets.json"
CONTRACT_KEYS = {"id", "source", "venue", "league", "eventStart", "eventDate", "awayTeam", "homeTeam",
                 "awayKey", "homeKey", "rotation", "betType", "period", "side", "points", "price", "stake",
                 "toWin", "contracts", "placedAt", "status", "closedAt", "isParlayLeg", "parlayId",
                 "legIndex", "legCount", "approx", "unmatchable", "sourceFetchedAt", "raw"}
# 13 open bets, two of them 2-leg parlays, one also in the settled list (settled wins): 14.
OPEN_RECORDS = 14
# 11 settled bets, one a 2-leg parlay, one skipped (settlementStatus None): 11.
SETTLED_RECORDS = 11
ALL_RECORDS = OPEN_RECORDS + SETTLED_RECORDS
NOW_ISO = "2026-10-10T19:05:00Z"


@pytest.fixture
def push() -> dict:
    checked = validate_push(json.loads(FIXTURE_PATH.read_text()))
    assert isinstance(checked, dict), checked
    return checked


@pytest.fixture
def records(push) -> list[dict]:
    found, _skipped = normalize_draftkings(push)
    return found


def by_id(records: list[dict], record_id: str) -> dict:
    found = [record for record in records if record["id"] == record_id]
    assert len(found) == 1, f"record {record_id} missing"
    return found[0]


def test_the_captured_alternate_spreads_read_as_open_away_and_home_favourites(records):
    ndsu = by_id(records, "draftkings:639272551505388534")
    assert (ndsu["league"], ndsu["period"], ndsu["betType"], ndsu["side"], ndsu["points"], ndsu["price"]) == (
        "cfb", "FG", "spread", "away", -14.5, 566)
    assert (ndsu["awayTeam"], ndsu["homeTeam"]) == ("North Dakota State", "UNLV")
    assert (ndsu["stake"], ndsu["toWin"], ndsu["status"], ndsu["pnl"]) == (75, 424.5, "open", None)
    assert (ndsu["placedAt"], ndsu["eventStart"], ndsu["eventDate"]) == (
        "2026-10-10T18:52:30Z", "2026-10-10T23:00:00Z", "2026-10-10")
    assert ndsu["closedAt"] is None and ndsu["unmatchable"] is None
    temple = by_id(records, "draftkings:639272551036408222")
    assert (temple["side"], temple["points"], temple["price"], temple["homeTeam"]) == ("home", -17.5, 485, "Temple")
    jmu = by_id(records, "draftkings:639272549945262625")
    assert (jmu["side"], jmu["points"], jmu["toWin"]) == ("away", -22.5, 1221.75)


def test_a_total_a_moneyline_and_a_first_half_spread(records):
    total = by_id(records, "draftkings:700000000000000001")
    assert (total["league"], total["betType"], total["side"], total["points"], total["price"]) == (
        "nfl", "total", "over", 44.5, -110)
    moneyline = by_id(records, "draftkings:700000000000000002")
    assert (moneyline["betType"], moneyline["side"], moneyline["points"], moneyline["price"]) == (
        "moneyline", "away", None, 135)
    first_half = by_id(records, "draftkings:700000000000000003")
    assert (first_half["period"], first_half["betType"], first_half["side"], first_half["points"]) == (
        "1H", "spread", "home", -3)


def test_a_parlay_is_one_record_per_leg_on_the_tickets_stake(records):
    legs = [by_id(records, f"draftkings:700000000000000004:leg{index}") for index in range(2)]
    assert [(leg["league"], leg["betType"], leg["side"]) for leg in legs] == [("nfl", "total", "under"),
                                                                            ("cfb", "spread", "home")]
    for index, leg in enumerate(legs):
        assert (leg["isParlayLeg"], leg["parlayId"], leg["legIndex"], leg["legCount"]) == (
            True, "draftkings:700000000000000004", index, 2)
        assert (leg["stake"], leg["toWin"], leg["raw"]["parlayPrice"]) == (20, 52.8, 264)


@pytest.mark.parametrize("record_id, reason", [
    ("draftkings:700000000000000005", draftkings.REASON_ROUND_ROBIN),
    ("draftkings:700000000000000006", "market 'Team Total Points' is not a game line"),
    ("draftkings:700000000000000007", "league id not supported ('40253')"),
    ("draftkings:700000000000000008:leg0", draftkings.REASON_SGP_GROUP),
])
def test_unseen_shapes_fail_closed_with_the_reason(records, record_id, reason):
    record = by_id(records, record_id)
    assert (record["unmatchable"], record["betType"], record["league"]) == (reason, "other", None)


def test_the_other_leg_of_a_parlay_with_an_sgp_group_still_reads(records):
    assert by_id(records, "draftkings:700000000000000008:leg1")["unmatchable"] is None


def test_a_settled_bet_replaces_its_open_copy_with_the_venues_result(records):
    won = by_id(records, "draftkings:700000000000000009")
    assert (won["status"], won["toWin"], won["pnl"], won["closedAt"], won["raw"]["read"]) == (
        "won", 50, 50, "2026-10-04T23:40:00Z", "settled")
    assert (won["side"], won["points"]) == ("away", 3.5)


@pytest.mark.parametrize("native_id, status, pnl", [
    ("700000000000000010", "lost", -110),
    ("700000000000000011", "push", 0),
    ("700000000000000012", "void", 0),
    ("700000000000000013", "closed", 30.5),  # cashed out: the venue's P&L, no won / lost
    ("700000000000000014", "closed", 50),    # half won: the venue's P&L, never forced
])
def test_each_settlement_status_books_the_venues_result(records, native_id, status, pnl):
    record = by_id(records, f"draftkings:{native_id}")
    assert (record["status"], record["pnl"]) == (status, pnl)
    assert record["closedAt"] is not None


def test_a_free_bet_stakes_nothing_so_its_loss_costs_nothing_and_its_win_is_all_profit(records):
    lost = by_id(records, "draftkings:700000000000000015")
    assert (lost["stake"], lost["pnl"], lost["status"]) == (0, 0, "lost")
    won = by_id(records, "draftkings:700000000000000016")
    assert (won["stake"], won["toWin"], won["pnl"], won["status"]) == (0, 50, 50, "won")


def test_a_risk_free_bet_is_cals_own_stake_so_its_loss_costs_it(records):
    risk_free = by_id(records, "draftkings:700000000000000019")
    assert (risk_free["stake"], risk_free["toWin"], risk_free["pnl"], risk_free["status"]) == (100, 90.91, -100, "lost")
    assert risk_free["raw"]["bonusType"] == "RiskFreeBet"


def test_a_bet_on_the_open_list_stays_open_whatever_its_settlement_status_says(records):
    partial = by_id(records, "draftkings:700000000000000020")
    assert (partial["status"], partial["closedAt"], partial["pnl"]) == ("open", None, None)
    assert partial["raw"]["settlementStatus"] == "PartialCashOut"


def test_a_settled_status_the_table_does_not_name_is_skipped_with_the_reason(push):
    _records, skipped = normalize_draftkings(push)
    assert skipped == ["bet 700000000000000017: settlementStatus 'None' in the settled list"]


def test_every_record_carries_the_contract_keys_and_ids_are_unique(records):
    assert len(records) == ALL_RECORDS
    assert len({record["id"] for record in records}) == ALL_RECORDS
    for record in records:
        assert CONTRACT_KEYS <= set(record), record["id"]
        assert (record["venue"], record["source"]) == ("draftkings", "draftkings_extension")


@pytest.mark.parametrize("text, american", [
    ("+566", 566), ("−110", -110), ("-110", -110), ("EVEN", 100), ("+50", None), ("", None), (None, None)])
def test_american_of(text, american):
    assert american_of(text) == american


@pytest.mark.parametrize("name, parsed", [
    ("Spread", ("spread", "FG")), ("Spread Alternate", ("spread", "FG")), ("1st Half Spread", ("spread", "1H")),
    ("Total - 1st Half", ("total", "1H")), ("Run Line", ("spread", "FG")), ("Moneyline", ("moneyline", "FG")),
    ("1st 5 Innings Total Runs", ("total", "F5")), ("2nd Quarter Moneyline", ("moneyline", "2Q")),
])
def test_market_of_reads_the_period_and_the_market(name, parsed):
    assert market_of(name) == parsed


def test_market_of_names_what_it_cannot_read():
    assert market_of("Player Passing Yards") == "market 'Player Passing Yards' is not a game line"


def test_moment_of_never_guesses_a_zone():
    assert moment_of("2026-10-10T18:52:30.706Z").isoformat() == "2026-10-10T18:52:30.706000+00:00"
    assert moment_of("2026-10-10T18:52:30") is None
    assert moment_of(None) is None


def test_validate_push_accepts_a_complete_push_and_an_error(push):
    assert set(push) == {"fetchedAt", "open", "settled", "events"}
    assert validate_push({"error": " not logged in "}) == {"error": "not logged in"}


@pytest.mark.parametrize("body, message", [
    ([], "body must be an object, got list"),
    ({"error": ""}, "error must be a non-empty string"),
    ({"open": [], "settled": [], "events": {}}, "fetchedAt must be an ISO time, got None"),
    ({"fetchedAt": NOW_ISO, "settled": [], "events": {}}, "open must be a list, got NoneType"),
    ({"fetchedAt": NOW_ISO, "open": [], "events": {}}, "settled must be a list, got NoneType"),
    ({"fetchedAt": NOW_ISO, "open": [{"stake": 1}], "settled": [], "events": {}},
     "open[0] carries no betId; keys seen: ['stake']"),
    ({"fetchedAt": NOW_ISO, "open": [], "settled": [], "events": []}, "events must be an object of event objects, each with a participants list of objects"),
    ({"fetchedAt": NOW_ISO, "open": [{"betId": "1", "selections": ["x"]}], "settled": [], "events": {}},
     "open[0].selections must be a list of objects with a participants list"),
    ({"fetchedAt": NOW_ISO, "open": [], "settled": [], "events": {"1": {"participants": ["x"]}}},
     "events must be an object of event objects, each with a participants list of objects"),
])
def test_validate_push_says_what_is_wrong(body, message):
    assert validate_push(body) == message


def test_closed_by_absence_closes_only_the_open_draftkings_records_the_push_lacks():
    stored = [
        {"id": "draftkings:1", "venue": "draftkings", "status": "open", "raw": {"nativeId": "1"}},
        {"id": "draftkings:2", "venue": "draftkings", "status": "open", "raw": {"nativeId": "2"}},
        {"id": "draftkings:3", "venue": "draftkings", "status": "won", "raw": {"nativeId": "3"}},
        {"id": "bet105:prematch:4", "venue": "bet105", "status": "open", "raw": {}},
    ]
    closed = closed_by_absence(stored, {"draftkings:1"}, NOW_ISO)
    assert [record["id"] for record in closed] == ["draftkings:2"]
    assert (closed[0]["status"], closed[0]["closedAt"], closed[0]["raw"]["closedReason"]) == (
        "closed", NOW_ISO, draftkings.REASON_CLOSED_BY_ABSENCE)


@pytest.fixture
def store(tmp_path):
    bets_store = BetsStore(tmp_path / "bets.duckdb", 7)
    yield bets_store
    bets_store.close()


@pytest.fixture
def http_server(store):
    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(store, 0.0, ["kalshi", draftkings.VENUE]))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()
    server.server_close()


def request_json(method: str, url: str, body: object = None, content_type: str = "application/json") -> tuple[int, dict]:
    data = None if body is None else (body if isinstance(body, bytes) else json.dumps(body).encode())
    headers = {"Content-Type": content_type} if data is not None else {}
    request = urllib.request.Request(url, data=data, method=method, headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


def draftkings_bets(payload: dict) -> dict[str, dict]:
    return {bet["id"]: bet for bet in payload["bets"] if bet["venue"] == "draftkings"}


def test_http_push_stores_settles_closes_what_neither_list_carries_and_records_an_error_run(http_server, push):
    status, reply = request_json("POST", f"{http_server}/draftkings.json", push)
    assert (status, reply) == (200, {"ok": True, "count": ALL_RECORDS, "closed": 0, "skipped": 1})
    status, payload = request_json("GET", f"{http_server}/bets.json?days=3650")
    source = payload["sources"]["draftkings"]
    assert (source["ok"], source["error"], source["count"]) == (True, None, ALL_RECORDS)
    assert len(draftkings_bets(payload)) == ALL_RECORDS

    # The NDSU bet leaves the open list and the settled list does not carry it: closed, no result.
    push["open"] = [bet for bet in push["open"] if bet["betId"] != "639272551505388534"]
    status, reply = request_json("POST", f"{http_server}/draftkings.json", push)
    assert (status, reply) == (200, {"ok": True, "count": ALL_RECORDS - 1, "closed": 1, "skipped": 1})
    bets = draftkings_bets(request_json("GET", f"{http_server}/bets.json?days=3650")[1])
    gone = bets["draftkings:639272551505388534"]
    assert (gone["status"], gone["raw"]["closedReason"]) == ("closed", draftkings.REASON_CLOSED_BY_ABSENCE)
    assert bets["draftkings:639272551036408222"]["status"] == "open"

    status, reply = request_json("POST", f"{http_server}/draftkings.json", {"error": draftkings_error()})
    assert (status, reply) == (200, {"ok": True, "recorded": "error"})
    source = request_json("GET", f"{http_server}/bets.json?days=3650")[1]["sources"]["draftkings"]
    assert (source["ok"], source["error"]) == (False, draftkings_error())


def draftkings_error() -> str:
    return "not logged in at sportsbook.draftkings.com"


def test_http_push_refuses_bad_bodies_and_other_verbs(http_server):
    status, reply = request_json("POST", f"{http_server}/draftkings.json", {"open": []})
    assert (status, reply["error"]) == (400, "fetchedAt must be an ISO time, got None")
    status, _reply = request_json("POST", f"{http_server}/draftkings.json", b'{"error": "x"}', "text/plain")
    assert status == 415
    status, _reply = request_json("GET", f"{http_server}/draftkings.json")
    assert status == 404


def test_the_pushed_source_is_listed_before_its_first_push(http_server):
    status, payload = request_json("GET", f"{http_server}/bets.json")
    assert status == 200 and payload["sources"]["draftkings"] == service.NO_POLL_YET
