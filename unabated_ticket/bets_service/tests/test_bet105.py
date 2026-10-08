"""Bet105 parser on tests/fixtures/bets/bet105_history.json (the push body the
extension sends): the open read, the settled read and their merge, the push
validation, the closed-by-absence fallback, and POST /bet105.json over a real
HTTP server."""
import json
import threading
import urllib.error
import urllib.request
from datetime import datetime
from http.server import ThreadingHTTPServer
from pathlib import Path
from zoneinfo import ZoneInfo

import pytest

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.sources import bet105
from unabated_ticket.bets_service.sources.bet105 import (
    american_from_decimal, american_from_text, closed_by_absence, feeds_by_native_id, merge_settled, moment_of_iso,
    normalize_bet105, normalize_group, normalize_settled, settled_status_of, status_of, validate_push)
from unabated_ticket.bets_service.store import BetsStore

FIXTURE_PATH = Path(__file__).parents[2] / "tests" / "fixtures" / "bets" / "bet105_history.json"
CONTRACT_KEYS = {"id", "source", "venue", "league", "eventStart", "eventDate", "awayTeam", "homeTeam",
                 "awayKey", "homeKey", "rotation", "betType", "period", "side", "points", "price", "stake",
                 "toWin", "contracts", "placedAt", "status", "closedAt", "isParlayLeg", "parlayId",
                 "legIndex", "legCount", "approx", "unmatchable", "sourceFetchedAt", "raw"}
# 14 prematch groups (one a 2-leg parlay) + 1 live group.
FIXTURE_RECORDS = 16
# 22 settled-list wagers: 19 graded (one a 3-leg parlay) read, one Pending left to the
# open read, two skipped (an unknown productCode, an unknown wagerStatus).
SETTLED_RECORDS = 21
SETTLED_SKIPPED = 2
# The settled read carries open groups 91000001-91000003, graded since.
MERGED_RECORDS = FIXTURE_RECORDS - 3 + SETTLED_RECORDS
NOW_ISO = "2026-09-29T05:00:00Z"
PACIFIC = ZoneInfo("America/Los_Angeles")


@pytest.fixture
def push() -> dict:
    return json.loads(FIXTURE_PATH.read_text())


@pytest.fixture
def records(push) -> list[dict]:
    return normalize_bet105(push["feeds"], push["fetchedAt"])


@pytest.fixture
def settled(push) -> list[dict]:
    records, _skipped = normalize_settled(push["settled"], push["fetchedAt"], {})
    return records


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
    ("bet105:prematch:91000011", "league not supported ('NCAA Football')", "leagueName", "NCAA Football"),
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


# ---- the settled read -----------------------------------------------------------------

def tracker_pnl(record: dict) -> float | None:
    """What the Bet Tracker books (server/tracker/trackerstats.js pnlOf)."""
    if record["status"] == "lost":
        return -record["stake"]
    if record["status"] in ("push", "void"):
        return 0
    if record["status"] == "won":
        return record["toWin"]
    return None


def test_a_graded_open_bet_settles_on_the_same_id_with_the_venues_result_and_settle_time(settled):
    record = by_id(settled, "bet105:prematch:91000002")
    assert (record["status"], record["stake"], record["toWin"]) == ("won", 117, 100)
    assert record["closedAt"] == "2026-10-06T01:56:19Z"  # gradeTime, never the push's clock
    assert record["placedAt"] == "2026-09-28T02:01:09Z"
    # The tracker buckets by the Pacific day it settled: Monday night, not UTC Tuesday.
    settle_day = datetime.fromisoformat(record["closedAt"]).astimezone(PACIFIC).date().isoformat()
    assert settle_day == "2026-10-05"
    assert (record["league"], record["betType"], record["period"], record["side"], record["points"], record["price"]) == (
        "nfl", "total", "1H", "over", 22.5, -117)
    assert (record["awayTeam"], record["homeTeam"], record["eventStart"], record["eventDate"]) == (
        "Atlanta Falcons", "New Orleans Saints", "2026-10-06T00:15:00Z", "2026-10-05")
    assert (record["raw"]["read"], record["raw"]["wagerStatus"], record["raw"]["result"], record["raw"]["legGrade"]) == (
        "settled", "Win", 100, "W")
    assert (record["raw"]["nativeId"], record["raw"]["feed"], record["raw"]["wagerId"]) == ("91000002", "prematch", 81000002)
    lost = by_id(settled, "bet105:prematch:91000001")
    assert (lost["status"], lost["closedAt"], tracker_pnl(lost)) == ("lost", "2026-10-04T18:30:46Z", -105)


def test_every_graded_wager_books_exactly_the_venues_result(settled):
    for record in settled:
        assert record["pnl"] == record["raw"]["result"], record["id"]  # the venue's own P&L, counted first
        if record["status"] == "closed":
            continue  # a cash-out: only its pnl says what it made
        assert tracker_pnl(record) == record["raw"]["result"], record["id"]


def test_a_win_pays_what_the_venue_paid_not_the_quoted_to_win(settled):
    record = by_id(settled, "bet105:prematch:93000014")
    assert (record["toWin"], record["raw"]["toWin"], record["raw"]["result"]) == (10.4, 10.401, 10.4)


def test_a_free_play_stakes_nothing_so_its_loss_costs_nothing(settled):
    record = by_id(settled, "bet105:prematch:93000007")
    assert (record["status"], record["stake"], record["raw"]["risk"], record["raw"]["isFreePlay"]) == ("lost", 0, 112, True)
    assert tracker_pnl(record) == 0


@pytest.mark.parametrize("record_id, bet_type, side, points, price", [
    ("bet105:prematch:93000001", "spread", "home", -2, -103),    # Merrimack -2: figure 2 is the away number
    ("bet105:prematch:93000004", "spread", "home", 9, -104),     # Maryland +9: figure -9
    ("bet105:prematch:93000011", "spread", "away", -18.5, -117), # San Diego -18.5
    ("bet105:prematch:93000003", "moneyline", "home", None, 192),
    ("bet105:prematch:93000005", "total", "under", 12.5, 116),
    ("bet105:prematch:93000002", "total", "over", 31, 100),      # older leg: odds "+100" as text
])
def test_settled_legs_read_their_side_and_line_the_open_reads_way(settled, record_id, bet_type, side, points, price):
    record = by_id(settled, record_id)
    assert (record["betType"], record["side"], record["points"], record["price"], record["unmatchable"]) == (
        bet_type, side, points, price, None)


@pytest.mark.parametrize("record_id, league, period", [
    ("bet105:prematch:93000005", "nfl", "4Q"), ("bet105:prematch:93000006", "nfl", "1Q"),
    ("bet105:prematch:93000002", "cfb", "1H"), ("bet105:prematch:93000011", "cfb", "FG"),
    ("bet105:prematch:93000001", "cbb", "1H"), ("bet105:prematch:93000012", "cbb", "1H"),
    ("bet105:prematch:93000013", "nba", "FG"), ("bet105:prematch:93000015", "mlb", "FG"),
])
def test_the_leagues_and_periods_seen_on_bets(settled, league, period, record_id):
    record = by_id(settled, record_id)
    assert (record["league"], record["period"]) == (league, period)


def test_a_settled_wager_that_is_not_a_game_line_still_carries_its_result(settled):
    team_total = by_id(settled, "bet105:prematch:93000009")
    assert team_total["unmatchable"] == "market 8 (Team Total 2) is not a game line"
    assert (team_total["status"], tracker_pnl(team_total)) == ("won", 100)
    soccer = by_id(settled, "bet105:prematch:93000014")
    assert soccer["unmatchable"] == "league not supported ('UEFA Champions League')"


def test_a_settled_parlay_is_one_record_per_leg_on_the_tickets_stake_and_status(settled):
    legs = [by_id(settled, f"bet105:prematch:93000008:leg{index}") for index in range(3)]
    for index, record in enumerate(legs):
        assert (record["isParlayLeg"], record["parlayId"], record["legIndex"], record["legCount"]) == (
            True, "bet105:prematch:93000008", index, 3)
        assert (record["status"], record["stake"], record["toWin"], record["raw"]["parlayPrice"]) == ("lost", 0, 0, 623)
    assert [record["raw"]["legGrade"] for record in legs] == ["W", "L", "W"]
    assert [record["side"] for record in legs] == ["over", "under", "over"]


def test_void_and_cash_out_come_from_the_sites_own_status_rules(settled):
    void = by_id(settled, "bet105:prematch:94000001")
    assert (void["status"], void["closedAt"], tracker_pnl(void)) == ("void", "2026-09-28T17:00:00Z", 0)
    cashed_out = by_id(settled, "bet105:prematch:94000002")
    assert (cashed_out["status"], cashed_out["closedAt"], cashed_out["pnl"]) == ("closed", "2026-09-28T01:30:00Z", 40)


def test_pending_wagers_are_left_to_the_open_read_and_unreadable_ones_name_the_reason(push):
    records, skipped = normalize_settled(push["settled"], push["fetchedAt"], {})
    assert len(records) == SETTLED_RECORDS
    assert not [record for record in records if record["raw"]["nativeId"] == "93000010"]  # the March Pending wager
    assert skipped == [
        "ticket 94000003: unknown productCode 'Unseen' and no open record to take the feed from",
        "ticket 94000004: unknown wagerStatus 'Unseen'",
    ]
    assert len(skipped) == SETTLED_SKIPPED


def test_a_wager_settles_under_the_feed_its_tickets_record_sits_under(push):
    # 94000003's productCode names no feed; held under the live feed, it settles there.
    records, skipped = normalize_settled(push["settled"], push["fetchedAt"], {"94000003": "live", "91000002": "live"})
    assert by_id(records, "bet105:live:94000003")["status"] == "won"
    assert by_id(records, "bet105:live:91000002")["status"] == "won"  # the record's feed beats the productCode's
    assert skipped == ["ticket 94000004: unknown wagerStatus 'Unseen'"]


def test_feeds_by_native_id_reads_bet105_records_and_the_later_one_wins():
    records = [
        {"venue": "bet105", "raw": {"nativeId": "1", "feed": "prematch"}},
        {"venue": "bet105", "raw": {"nativeId": "1", "feed": "live"}},
        {"venue": "bet105", "raw": {"nativeId": "2", "feed": "prematch"}},
        {"venue": "kalshi", "raw": {"nativeId": "3", "feed": "prematch"}},
        {"venue": "bet105", "raw": {"nativeId": "4"}},
    ]
    assert feeds_by_native_id(records) == {"1": "live", "2": "prematch"}


def test_settled_records_carry_the_contract_keys_and_pnl_and_ids_are_unique(settled):
    assert len({record["id"] for record in settled}) == SETTLED_RECORDS
    for record in settled:
        assert set(record) == CONTRACT_KEYS | {"pnl"}, record["id"]


def test_merge_settled_lets_a_graded_ticket_replace_its_open_records(records, settled):
    merged = merge_settled(records, settled)
    assert len(merged) == MERGED_RECORDS and len({record["id"] for record in merged}) == MERGED_RECORDS
    assert by_id(merged, "bet105:prematch:91000002")["status"] == "won"
    assert by_id(merged, "bet105:prematch:91000004")["status"] == "open"
    # The ticket decides, not the feed: an open record under the other feed gives way too.
    live_open = {**by_id(records, "bet105:prematch:91000001"), "id": "bet105:live:91000001"}
    live_open["raw"] = {**live_open["raw"], "feed": "live"}
    assert "bet105:live:91000001" not in {record["id"] for record in merge_settled([live_open], settled)}


@pytest.mark.parametrize("wager, status", [
    ({"wagerStatus": "Win"}, "won"), ({"wagerStatus": "Loss"}, "lost"), ({"wagerStatus": "Push"}, "push"),
    ({"wagerStatus": "Cancel"}, "void"), ({"wagerStatus": "NO_ACTION"}, "void"),
    ({"wagerStatus": "Loss", "isCashout": True}, "closed"), ({"wagerStatus": "Pending"}, None),
    ({"wagerStatus": "Graded"}, None), ({}, None),
])
def test_settled_status_of(wager, status):
    assert settled_status_of(wager) == status


@pytest.mark.parametrize("text, american", [
    ("-117", -117), ("+192", 192), ("100", 100), ("+100", 100), ("-99", None), ("EV", None), ("", None), (None, None), (-117, None),
])
def test_american_from_text(text, american):
    assert american_from_text(text) == american


@pytest.mark.parametrize("text, iso", [
    ("2026-10-06T01:56:19Z", "2026-10-06T01:56:19+00:00"), ("2026-10-06 00:15:00+00:00", "2026-10-06T00:15:00+00:00"),
    ("2026-10-06 00:15:00", None), ("soon", None), ("", None), (None, None),
])
def test_moment_of_iso_never_guesses_a_zone(text, iso):
    moment = moment_of_iso(text)
    assert (None if moment is None else moment.isoformat()) == iso


# ---- the push -------------------------------------------------------------------------

def test_validate_push_accepts_a_complete_push_and_an_error(push):
    clean = validate_push(push)
    assert isinstance(clean, dict) and clean["fetchedAt"] == push["fetchedAt"]
    assert [len(clean["feeds"][feed]) for feed in bet105.FEEDS] == [14, 1]
    assert len(clean["settled"]) == len(push["settled"])
    assert validate_push({"error": "  not logged in at app.bet105.ag  "}) == {"error": "not logged in at app.bet105.ag"}
    dict_shaped = {"fetchedAt": NOW_ISO, "feeds": {"prematch": {"1": {"betGroupId": 1}}, "live": []}, "settled": []}
    assert validate_push(dict_shaped) == {"fetchedAt": NOW_ISO, "feeds": {"prematch": [{"betGroupId": 1}], "live": []}, "settled": []}


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
    ({"fetchedAt": NOW_ISO, "feeds": {"prematch": [], "live": []}}, "settled must be a list, got NoneType"),
    ({"fetchedAt": NOW_ISO, "feeds": {"prematch": [], "live": []}, "settled": {}}, "settled must be a list, got dict"),
    ({"fetchedAt": NOW_ISO, "feeds": {"prematch": [], "live": []}, "settled": ["x"]}, "settled[0] must be an object, got str"),
    ({"fetchedAt": NOW_ISO, "feeds": {"prematch": [], "live": []}, "settled": [{"wagerStatus": "Win"}]},
     "settled[0] carries no ticketNumber; keys seen: ['wagerStatus']"),
    ({"fetchedAt": NOW_ISO, "feeds": {"prematch": [], "live": []}, "settled": [{"ticketNumber": " "}]},
     "settled[0] carries no ticketNumber; keys seen: ['ticketNumber']"),
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


def test_http_push_settles_from_the_venue_closes_only_what_neither_read_lists_and_records_an_error_run(http_server, push):
    settled_list = push["settled"]
    # Before the bets were graded: every open group is open, the wager list shows none of them.
    push["settled"] = []
    status, reply = request_json("POST", f"{http_server}/bet105.json", push)
    assert (status, reply) == (200, {"ok": True, "count": FIXTURE_RECORDS, "settled": 0, "closed": 0, "skipped": 0})
    status, payload = request_json("GET", f"{http_server}/bets.json")
    assert status == 200
    source = payload["sources"]["bet105"]
    assert (source["ok"], source["error"], source["count"]) == (True, None, FIXTURE_RECORDS)
    assert source["fetchedAt"] is not None
    bets = bet105_bets(payload)
    assert len(bets) == FIXTURE_RECORDS and bets["bet105:prematch:91000002"]["status"] == "open"

    # Graded: 91000001-91000003 leave the open list and the wager list carries their results.
    # 91000004 leaves the open list too, but no wager carries it: the closed-by-absence fallback.
    # The live group 92000001 is graded under a productCode that names no feed: it settles
    # under the live feed its stored record sits under, not closed with no result.
    graded = {91000001, 91000002, 91000003, 91000004}
    push["feeds"]["prematch"] = [group for group in push["feeds"]["prematch"] if group["betGroupId"] not in graded]
    push["feeds"]["live"] = []
    live_wager = {**settled_list[1], "ticketNumber": "92000001", "productCode": "Unseen"}
    push["settled"] = settled_list + [live_wager]
    status, reply = request_json("POST", f"{http_server}/bet105.json", push)
    count = FIXTURE_RECORDS - len(graded) - 1 + SETTLED_RECORDS + 1
    assert (status, reply) == (200, {"ok": True, "count": count, "settled": SETTLED_RECORDS + 1, "closed": 1,
                                     "skipped": SETTLED_SKIPPED})
    bets = bet105_bets(request_json("GET", f"{http_server}/bets.json?days=3650")[1])
    won = bets["bet105:prematch:91000002"]
    assert (won["status"], won["toWin"], won["closedAt"]) == ("won", 100, "2026-10-06T01:56:19Z")
    assert "closedReason" not in won["raw"]
    assert (bets["bet105:prematch:91000001"]["status"], bets["bet105:prematch:91000001"]["closedAt"]) == (
        "lost", "2026-10-04T18:30:46Z")
    fallback = bets["bet105:prematch:91000004"]
    assert (fallback["status"], fallback["raw"]["closedReason"]) == ("closed", bet105.REASON_CLOSED_BY_ABSENCE)
    assert fallback["closedAt"] is not None and bets["bet105:prematch:91000005"]["status"] == "open"
    live = bets["bet105:live:92000001"]
    assert (live["status"], live["pnl"], live["closedAt"]) == ("won", 100, "2026-10-06T01:56:19Z")
    assert "bet105:prematch:92000001" not in bets
    assert len(bets) == FIXTURE_RECORDS - 3 + SETTLED_RECORDS  # history the store never saw open is kept too

    # An error push is a failed run: the panel shows the text; the records stand.
    status, reply = request_json("POST", f"{http_server}/bet105.json", {"error": "not logged in at app.bet105.ag"})
    assert (status, reply) == (200, {"ok": True, "recorded": "error"})
    status, payload = request_json("GET", f"{http_server}/bets.json?days=3650")
    source = payload["sources"]["bet105"]
    assert (source["ok"], source["error"], source["count"]) == (False, "not logged in at app.bet105.ag", count)
    assert len(bet105_bets(payload)) == FIXTURE_RECORDS - 3 + SETTLED_RECORDS


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
