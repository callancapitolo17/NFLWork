"""Bet105 parser on tests/fixtures/bets/bet105_history.json (the push body the
extension sends), the push validation, the closed-by-absence rule, and
POST /bet105.json over a real HTTP server."""
import json
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

import pytest

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.sources import bet105
from unabated_ticket.bets_service.sources.bet105 import (
    american_from_decimal, closed_by_absence, normalize_bet105, normalize_group, status_of, validate_push)
from unabated_ticket.bets_service.store import BetsStore

FIXTURE_PATH = Path(__file__).parents[2] / "tests" / "fixtures" / "bets" / "bet105_history.json"
CONTRACT_KEYS = {"id", "source", "venue", "league", "eventStart", "eventDate", "awayTeam", "homeTeam",
                 "awayKey", "homeKey", "rotation", "betType", "period", "side", "points", "price", "stake",
                 "toWin", "contracts", "placedAt", "status", "closedAt", "isParlayLeg", "parlayId",
                 "legIndex", "legCount", "approx", "unmatchable", "sourceFetchedAt", "raw"}
# 14 prematch groups (one a 2-leg parlay) + 1 live group.
FIXTURE_RECORDS = 16
NOW_ISO = "2026-09-29T05:00:00Z"


@pytest.fixture
def push() -> dict:
    return json.loads(FIXTURE_PATH.read_text())


@pytest.fixture
def records(push) -> list[dict]:
    return normalize_bet105(push["feeds"], push["fetchedAt"])


def by_id(records: list[dict], record_id: str) -> dict:
    found = [record for record in records if record["id"] == record_id]
    assert len(found) == 1, f"record {record_id} missing"
    return found[0]


# ---- parser ---------------------------------------------------------------------------

def test_the_captured_first_half_total_reads_as_an_open_over(records, push):
    record = by_id(records, "bet105:prematch:91000001")
    assert (record["source"], record["venue"], record["status"], record["closedAt"]) == ("bet105_extension", "bet105", "open", None)
    assert (record["league"], record["betType"], record["period"]) == ("nfl", "total", "1H")
    assert (record["side"], record["points"], record["price"]) == ("over", 24, -105)
    assert (record["awayTeam"], record["homeTeam"]) == ("Jacksonville Jaguars", "Cincinnati Bengals")
    assert (record["awayKey"], record["homeKey"], record["rotation"]) == (None, None, None)  # the panel resolves keys
    assert (record["stake"], record["toWin"], record["contracts"]) == (105, 100, None)
    assert record["placedAt"] == "2026-09-28T02:00:31Z"
    assert (record["eventStart"], record["eventDate"]) == ("2026-10-04T17:00:00Z", "2026-10-04")  # 1 PM ET Sunday
    assert (record["approx"], record["unmatchable"]) == ([], None)
    assert (record["isParlayLeg"], record["parlayId"], record["legIndex"], record["legCount"]) == (False, None, None, None)
    assert record["sourceFetchedAt"] == push["fetchedAt"]
    assert (record["raw"]["nativeId"], record["raw"]["feed"], record["raw"]["betGroupId"]) == ("91000001", "prematch", 91000001)
    assert (record["raw"]["marketId"], record["raw"]["key"], record["raw"]["subKey"], record["raw"]["periodId"]) == ("5", "24", "1", "h1")


def test_a_monday_night_kickoff_dates_on_its_eastern_day(records):
    record = by_id(records, "bet105:prematch:91000002")
    assert (record["eventStart"], record["eventDate"]) == ("2026-10-06T00:15:00Z", "2026-10-05")  # 8:15 PM ET Monday
    assert (record["side"], record["points"], record["price"]) == ("over", 22.5, -117)
    assert (record["stake"], record["toWin"]) == (117, 100)


def test_the_third_capture_and_the_under(records):
    record = by_id(records, "bet105:prematch:91000003")
    assert (record["side"], record["points"], record["price"], record["period"]) == ("over", 23, -110, "1H")
    assert (record["awayTeam"], record["homeTeam"]) == ("Dallas Cowboys", "Houston Texans")
    under = by_id(records, "bet105:prematch:91000004")
    assert (under["betType"], under["side"], under["points"], under["price"], under["period"]) == ("total", "under", 47.5, -110, "FG")


def test_a_spread_reads_the_key_as_the_away_number_and_negates_it_for_home(records):
    away = by_id(records, "bet105:prematch:91000005")
    assert (away["betType"], away["side"], away["points"], away["price"]) == ("spread", "away", -3.5, -110)
    assert (away["awayTeam"], away["homeTeam"]) == ("Los Angeles Rams", "San Francisco 49ers")
    home = by_id(records, "bet105:prematch:91000006")
    assert (home["betType"], home["side"], home["points"]) == ("spread", "home", 3.5)


def test_a_moneyline_has_no_points(records):
    record = by_id(records, "bet105:prematch:91000007")
    assert (record["betType"], record["side"], record["points"], record["price"]) == ("moneyline", "home", None, 135)
    assert (record["stake"], record["toWin"]) == (100, 135)


@pytest.mark.parametrize("record_id, reason, raw_key, raw_value", [
    ("bet105:prematch:91000008", "market 7 (Team Total 1) is not a game line", "marketId", "7"),
    ("bet105:prematch:91000009", "unknown total side code (3)", "subKey", "3"),
    ("bet105:prematch:91000010", "unknown period code ('x9')", "periodId", "x9"),
    ("bet105:prematch:91000011", "league not supported (leagueName 'NCAA Football')", "leagueName", "NCAA Football"),
    ("bet105:prematch:91000012", "styled market (marketStyleId 616, market 5) is not a game line", "marketStyleId", 616),
])
def test_unobserved_shapes_fail_closed_with_the_reason_and_keep_their_codes(records, record_id, reason, raw_key, raw_value):
    record = by_id(records, record_id)
    assert record["unmatchable"] == reason
    assert (record["league"], record["betType"], record["side"], record["points"], record["price"]) == (None, "other", None, None, None)
    assert record["status"] == "open" and record["stake"] == 110
    assert record["raw"][raw_key] == raw_value


def test_a_group_without_legs_fails_closed(records):
    record = by_id(records, "bet105:prematch:91000014")
    assert record["unmatchable"] == bet105.REASON_NO_LEGS
    assert record["status"] == "open" and record["betType"] == "other"


def test_a_parlay_is_one_record_per_leg_on_the_tickets_stake_and_price(records):
    first = by_id(records, "bet105:prematch:91000013:leg0")
    second = by_id(records, "bet105:prematch:91000013:leg1")
    for record in (first, second):
        assert (record["isParlayLeg"], record["parlayId"], record["legCount"]) == (True, "bet105:prematch:91000013", 2)
        assert (record["stake"], record["toWin"], record["raw"]["parlayPrice"]) == (100, 264, 264)
    assert (first["legIndex"], first["betType"], first["side"], first["points"], first["price"]) == (0, "total", "over", 44.5, -110)
    assert first["unmatchable"] is None
    assert (second["legIndex"], second["unmatchable"]) == (1, "market 30 (Outright Winner) is not a game line")


def test_the_live_feed_prefixes_the_id(records):
    record = by_id(records, "bet105:live:92000001")
    assert record["raw"]["feed"] == "live"
    assert (record["side"], record["points"], record["price"], record["period"]) == ("over", 21.5, -115, "1H")


def test_every_record_carries_the_contract_keys_and_ids_are_unique(records):
    assert len(records) == FIXTURE_RECORDS
    assert len({record["id"] for record in records}) == FIXTURE_RECORDS
    for record in records:
        assert set(record) == CONTRACT_KEYS, record["id"]


def test_a_group_without_a_bet_group_id_raises():
    with pytest.raises(RuntimeError, match="betGroupId"):
        normalize_group({"risk": 100, "componentBets": []}, "prematch", NOW_ISO)


def test_the_older_dict_shaped_bet_groups_are_read_too(push):
    feeds = {"prematch": {"91000001": push["feeds"]["prematch"][0]}, "live": {}}
    [record] = normalize_bet105(feeds, NOW_ISO)
    assert record["id"] == "bet105:prematch:91000001"


@pytest.mark.parametrize("decimal_odds, american", [
    (1.9523799419403076, -105), (1.854699969291687, -117), (1.9090900421142578, -110), (1.87, -115),
    (2.0, 100), (2.35, 135), (3.64, 264), (1.0, None), (0.5, None), (None, None), (True, None), ("2.0", None),
])
def test_american_from_decimal(decimal_odds, american):
    assert american_from_decimal(decimal_odds) == american


@pytest.mark.parametrize("state, status", [(0, "open"), (-1, "void"), ("0", "open"), (7, "unknown"), (None, "unknown"), (0.5, "unknown")])
def test_status_of(state, status):
    assert status_of(state) == status


# ---- the push -------------------------------------------------------------------------

def test_validate_push_accepts_a_complete_push_and_an_error(push):
    clean = validate_push(push)
    assert isinstance(clean, dict) and clean["fetchedAt"] == push["fetchedAt"]
    assert [len(clean["feeds"][feed]) for feed in bet105.FEEDS] == [14, 1]
    assert validate_push({"error": "  not logged in at app.bet105.ag  "}) == {"error": "not logged in at app.bet105.ag"}
    dict_shaped = {"fetchedAt": NOW_ISO, "feeds": {"prematch": {"1": {"betGroupId": 1}}, "live": []}}
    assert validate_push(dict_shaped) == {"fetchedAt": NOW_ISO, "feeds": {"prematch": [{"betGroupId": 1}], "live": []}}


@pytest.mark.parametrize("body, message", [
    ([], "body must be an object, got list"),
    ({"error": ""}, "error must be a non-empty string"),
    ({"feeds": {}}, "fetchedAt must be an ISO time, got None"),
    ({"fetchedAt": "yesterday", "feeds": {}}, "fetchedAt must be an ISO time, got 'yesterday'"),
    ({"fetchedAt": NOW_ISO, "feeds": {"prematch": []}}, "feeds must be an object with exactly ['prematch', 'live'], got ['prematch']"),
    ({"fetchedAt": NOW_ISO, "feeds": []}, "feeds must be an object with exactly ['prematch', 'live'], got list"),
    ({"fetchedAt": NOW_ISO, "feeds": {"prematch": "x", "live": []}}, "feeds.prematch must be a list, got str"),
    ({"fetchedAt": NOW_ISO, "feeds": {"prematch": [1], "live": []}}, "feeds.prematch[0] must be an object, got int"),
    ({"fetchedAt": NOW_ISO, "feeds": {"prematch": [], "live": [{"risk": 1}]}}, "feeds.live[0] carries no betGroupId; keys seen: ['risk']"),
    ({"fetchedAt": NOW_ISO, "feeds": {"prematch": [{"betGroupId": i} for i in range(1, 2002)], "live": []}},
     "at most 2000 bet groups per push, got 2001"),
])
def test_validate_push_says_what_is_wrong(body, message):
    assert validate_push(body) == message


def test_closed_by_absence_closes_only_the_open_bet105_records_the_push_lacks():
    stored = [
        {"id": "bet105:prematch:1", "venue": "bet105", "status": "open", "closedAt": None, "raw": {"nativeId": "1"}},
        {"id": "bet105:prematch:2", "venue": "bet105", "status": "open", "closedAt": None, "raw": {"nativeId": "2"}},
        {"id": "bet105:prematch:3", "venue": "bet105", "status": "closed", "closedAt": "2026-09-20T00:00:00Z", "raw": {}},
        {"id": "kalshi:x:yes", "venue": "kalshi", "status": "open", "closedAt": None, "raw": {}},
    ]
    [closed] = closed_by_absence(stored, {"bet105:prematch:1"}, NOW_ISO)
    assert (closed["id"], closed["status"], closed["closedAt"], closed["sourceFetchedAt"]) == ("bet105:prematch:2", "closed", NOW_ISO, NOW_ISO)
    assert closed["raw"] == {"nativeId": "2", "closedReason": bet105.REASON_CLOSED_BY_ABSENCE}
    assert stored[1]["status"] == "open"  # the input is not mutated


# ---- HTTP -----------------------------------------------------------------------------

@pytest.fixture
def store(tmp_path):
    bets_store = BetsStore(tmp_path / "bets.duckdb", 7)
    yield bets_store
    bets_store.close()


@pytest.fixture
def http_server(store):
    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(store, 0.0, ["kalshi", bet105.VENUE]))
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


def bet105_bets(payload: dict) -> dict[str, dict]:
    return {bet["id"]: bet for bet in payload["bets"] if bet["venue"] == "bet105"}


def test_http_push_stores_the_records_closes_the_absent_and_records_an_error_run(http_server, push):
    status, reply = request_json("POST", f"{http_server}/bet105.json", push)
    assert (status, reply) == (200, {"ok": True, "count": FIXTURE_RECORDS, "closed": 0})
    status, payload = request_json("GET", f"{http_server}/bets.json")
    assert status == 200
    source = payload["sources"]["bet105"]
    assert (source["ok"], source["error"], source["count"]) == (True, None, FIXTURE_RECORDS)
    assert source["fetchedAt"] is not None
    bets = bet105_bets(payload)
    assert len(bets) == FIXTURE_RECORDS and bets["bet105:prematch:91000002"]["status"] == "open"

    # The next complete push no longer lists group 91000002: it is closed, with no result.
    push["feeds"]["prematch"] = [group for group in push["feeds"]["prematch"] if group["betGroupId"] != 91000002]
    status, reply = request_json("POST", f"{http_server}/bet105.json", push)
    assert (status, reply) == (200, {"ok": True, "count": FIXTURE_RECORDS - 1, "closed": 1})
    bets = bet105_bets(request_json("GET", f"{http_server}/bets.json")[1])
    closed = bets["bet105:prematch:91000002"]
    assert (closed["status"], closed["raw"]["closedReason"]) == ("closed", bet105.REASON_CLOSED_BY_ABSENCE)
    assert closed["closedAt"] is not None and bets["bet105:prematch:91000001"]["status"] == "open"

    # An error push is a failed run: the panel shows the text; the records stand.
    status, reply = request_json("POST", f"{http_server}/bet105.json", {"error": "not logged in at app.bet105.ag"})
    assert (status, reply) == (200, {"ok": True, "recorded": "error"})
    status, payload = request_json("GET", f"{http_server}/bets.json")
    source = payload["sources"]["bet105"]
    assert (source["ok"], source["error"], source["count"]) == (False, "not logged in at app.bet105.ag", FIXTURE_RECORDS - 1)
    assert len(bet105_bets(payload)) == FIXTURE_RECORDS


def test_http_push_refuses_bad_bodies_and_other_verbs(http_server):
    status, reply = request_json("POST", f"{http_server}/bet105.json", {"feeds": {}})
    assert (status, reply["error"]) == (400, "fetchedAt must be an ISO time, got None")
    status, reply = request_json("POST", f"{http_server}/bet105.json", b'{"error": "x"}', "text/plain")
    assert status == 415
    status, _reply = request_json("GET", f"{http_server}/bet105.json")
    assert status == 404
    status, _reply = request_json("DELETE", f"{http_server}/bet105.json")
    assert status == 404


def test_the_pushed_source_is_listed_before_its_first_push(http_server):
    status, payload = request_json("GET", f"{http_server}/bets.json")
    assert status == 200 and payload["sources"]["bet105"] == service.NO_POLL_YET
