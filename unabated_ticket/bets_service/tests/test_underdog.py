"""Underdog (2026-11-07; pushed by the extension since 2026-11-09): the prediction
slips' records (sources/underdog.py) and one POST /underdog.json push. The slip is Cal's
$1 49ers moneyline (trimmed), open on 2026-11-07 and settled won on 2026-11-09 with
payout "1.23" = its contract_quantity. No network."""
import json
import threading
import urllib.error
import urllib.request
from datetime import datetime, timezone
from http.server import ThreadingHTTPServer
from pathlib import Path

import pytest

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.store import BetsStore
from unabated_ticket.bets_service.sources.underdog import (
    MAX_SETTLED_PAGES, closed_by_absence, normalize_underdog, records_of_push, validate_push)

FETCHED_AT = "2026-11-07T06:20:00Z"
SF_ID, LV_ID, CIN_ID, ATL_ID = "t-sf", "t-lv", "t-cin", "t-atl"
TEAMS = [{"id": SF_ID, "name": "San Francisco 49ers"}, {"id": LV_ID, "name": "Las Vegas Raiders"},
         {"id": CIN_ID, "name": "Cincinnati Bengals"}, {"id": ATL_ID, "name": "Atlanta Falcons"}]
GAMES = [{"id": 179088, "sport_id": "NFL", "away_team_id": LV_ID, "home_team_id": SF_ID,
          "scheduled_at": "2026-11-08T21:25:00Z"},
         {"id": 179090, "sport_id": "NFL", "away_team_id": CIN_ID, "home_team_id": ATL_ID,
          "scheduled_at": "2026-11-08T18:00:00Z"},
         {"id": 501, "sport_id": "PGA", "away_team_id": None, "home_team_id": None, "scheduled_at": None}]
OVER_UNDERS = [{"id": "ou-ml", "title": "LV @ SF Moneyline", "appearance_stat": {"stat": "moneyline"}},
               {"id": "ou-sp", "title": "CIN @ ATL Spread", "appearance_stat": {"stat": "spread"}},
               {"id": "ou-tot", "title": "CIN @ ATL Total Points O/U", "appearance_stat": {"stat": "points"}},
               {"id": "ou-1h", "title": "CIN @ ATL 1H Spread", "appearance_stat": {"stat": "spread"}},
               {"id": "ou-prop", "title": "Joe Burrow Passing Yards", "appearance_stat": {"stat": "passing_yds"}}]
LINES = [{"id": "l-ml", "over_under_id": "ou-ml", "provider_id": "udx", "stat_value": None},
         {"id": "l-sp", "over_under_id": "ou-sp", "provider_id": "udx", "stat_value": "3.5"},
         {"id": "l-tot", "over_under_id": "ou-tot", "provider_id": "udx", "stat_value": "52.5"},
         {"id": "l-1h", "over_under_id": "ou-1h", "provider_id": "udx", "stat_value": "1.5"},
         {"id": "l-prop", "over_under_id": "ou-prop", "provider_id": "udx", "stat_value": "250.5"}]
OPTIONS = [{"id": "o-sf", "over_under_line_id": "l-ml", "choice": "home", "choice_display": "49ers to win",
            "american_price": "-435", "decimal_price": "1.23"},
           {"id": "o-cin", "over_under_line_id": "l-sp", "choice": "away", "choice_display": "CIN -3.5",
            "american_price": "+101", "decimal_price": "2.01"},
           {"id": "o-atl", "over_under_line_id": "l-sp", "choice": "home", "choice_display": "ATL +3.5",
            "american_price": "-132", "decimal_price": "1.76"},
           {"id": "o-under", "over_under_line_id": "l-tot", "choice": "lower", "choice_display": "Lower",
            "american_price": "-112", "decimal_price": "1.9"},
           {"id": "o-1h", "over_under_line_id": "l-1h", "choice": "home", "choice_display": "ATL +1.5"},
           {"id": "o-prop", "over_under_line_id": "l-prop", "choice": "higher", "choice_display": "Higher"}]


def slip(slip_id, option_ids, match_ids, fee="1.0", contracts="1.23", status="active", payout="0.0",
         product="prediction_markets", created_at="2026-11-07T06:09:41Z", payout_at=None):
    return {"id": slip_id, "product": product, "type": "Prediction", "status": status, "fee": fee,
            "payout": payout, "created_at": created_at, "payout_at": payout_at, "display_result": status,
            "prediction": {"contract_quantity": contracts, "open_contract_quantity": contracts,
                           "trade_cost": "0.994", "total_fees": "0.01"},
            "selection_groups": [{"match_id": match_id, "match_type": "Game",
                                  "selections": [{"option_id": option_id, "payout_multiplier": "1.23",
                                                  "result": "pending"}]}
                                 for option_id, match_id in zip(option_ids, match_ids)]}


def page(slips, next_page=None, number=1):
    return {"data": {"entry_slips": slips, "games": GAMES, "teams": TEAMS, "over_unders": OVER_UNDERS,
                     "over_under_lines": LINES, "over_under_options": OPTIONS},
            "meta": {"count": len(slips), "items": len(slips), "next": next_page, "page": number}}


def read(active_slips, settled_slips=(), stored=None):
    """normalize_underdog on one active page and one settled page."""
    return normalize_underdog(page(list(active_slips)), [page(list(settled_slips))], stored or {}, FETCHED_AT)


def by_id(records):
    return {record["id"]: record for record in records}


def test_cals_live_moneyline_reads_as_the_49ers_at_minus_435():
    [record] = read([slip("s1", ["o-sf"], [179088])])
    assert record["id"] == "underdog:s1" and record["venue"] == "underdog" and record["status"] == "open"
    assert (record["league"], record["awayTeam"], record["homeTeam"]) == ("nfl", "Las Vegas Raiders", "San Francisco 49ers")
    assert (record["betType"], record["period"], record["side"], record["points"]) == ("moneyline", "FG", "home", None)
    assert (record["stake"], record["toWin"], record["price"]) == (1, 0.23, -435)
    assert record["eventStart"] == "2026-11-08T21:25:00Z" and record["eventDate"] == "2026-11-08"
    assert record["placedAt"] == "2026-11-07T06:09:41Z" and record["closedAt"] is None
    assert record["unmatchable"] is None


def test_spread_takes_its_sign_from_the_pick_and_a_total_reads_under():
    records = by_id(read([
        slip("fav", ["o-cin"], [179090], fee="10.0", contracts="20.1"),
        slip("dog", ["o-atl"], [179090], fee="10.0", contracts="17.6"),
        slip("und", ["o-under"], [179090], fee="10.0", contracts="19.0")]))
    assert (records["underdog:fav"]["side"], records["underdog:fav"]["points"]) == ("away", -3.5)
    assert (records["underdog:dog"]["side"], records["underdog:dog"]["points"]) == ("home", 3.5)
    assert (records["underdog:und"]["betType"], records["underdog:und"]["side"], records["underdog:und"]["points"]) == \
        ("total", "under", 52.5)
    assert records["underdog:fav"]["price"] == 101 and records["underdog:und"]["price"] == -111


@pytest.mark.parametrize("paid, status, pnl", [("0.0", "lost", None), ("1.23", "won", None), ("1.0", "void", None),
                                               ("0.6", "closed", -0.4), (None, "unknown", None)])
def test_a_settled_slip_grades_off_what_it_paid(paid, status, pnl):
    [record] = read([], [slip("s", ["o-sf"], [179088], status="settled", payout=paid,
                              payout_at="2026-11-09T00:40:00Z")])
    assert (record["status"], record.get("pnl"), record["closedAt"]) == (status, pnl, "2026-11-09T00:40:00Z")
    assert (record["stake"], record["toWin"], record["side"]) == (1, 0.23, "home")


def test_with_no_payout_time_a_settled_slip_keeps_its_stored_close_and_never_the_poll_clock_twice():
    settled = slip("s", ["o-sf"], [179088], status="settled", payout="0.0")
    [first] = read([], [settled])
    assert first["closedAt"] == FETCHED_AT
    [again] = normalize_underdog(page([]), [page([settled])], {"underdog:s": first}, "2026-11-07T07:00:00Z")
    assert again["closedAt"] == FETCHED_AT


def test_anything_on_the_active_list_is_open_whatever_its_status_word():
    [record] = read([slip("s", ["o-sf"], [179088], status="live")])
    assert (record["status"], record["closedAt"]) == ("open", None)


def test_a_settled_answer_without_the_game_keeps_the_stored_game_fields():
    [stored] = read([slip("s", ["o-sf"], [179088])])
    settled = page([slip("s", ["o-sf"], [179088], status="settled", payout="1.23")])
    settled["data"]["games"] = []
    [record] = normalize_underdog(page([]), [settled], {"underdog:s": stored}, FETCHED_AT)
    assert (record["status"], record["homeTeam"], record["side"], record["unmatchable"]) == \
        ("won", "San Francisco 49ers", "home", None)
    assert not record["raw"].get("parseFailed")


def test_pickem_is_skipped_and_a_slip_on_both_lists_keeps_its_settled_copy():
    records = read([slip("s1", ["o-sf"], [179088])],
                   [slip("s1", ["o-sf"], [179088], status="settled", payout="1.23"),
                    slip("pk", ["o-sf"], [179088], product="fantasy")])
    assert [(record["id"], record["status"]) for record in records] == [("underdog:s1", "won")]


def test_what_does_not_place_is_unmatchable_and_only_an_unreadable_game_line_is_a_code_fix():
    records = by_id(read([
        slip("combo", ["o-sf", "o-cin"], [179088, 179090], contracts="2.47"),
        slip("half", ["o-1h"], [179090]),
        slip("prop", ["o-prop"], [179090]),
        slip("gone", ["o-missing"], [179090])]))
    assert records["underdog:combo"]["unmatchable"].startswith("combo of 2 picks")
    assert records["underdog:combo"]["raw"]["sportId"] == "NFL"
    assert "names a period" in records["underdog:half"]["unmatchable"]
    assert "not a game line" in records["underdog:prop"]["unmatchable"]
    for quiet in ("combo", "half", "prop"):
        assert records[f"underdog:{quiet}"]["betType"] == "other"
        assert not records[f"underdog:{quiet}"]["raw"].get("parseFailed")
    assert records["underdog:gone"]["raw"]["parseFailed"] is True


def test_a_spread_whose_pick_disagrees_with_the_line_fails_closed():
    options = OPTIONS + [{"id": "o-bad", "over_under_line_id": "l-sp", "choice": "away", "choice_display": "CIN -4.5"}]
    bad = page([slip("bad", ["o-bad"], [179090])])
    bad["data"]["over_under_options"] = options
    [record] = normalize_underdog(bad, [], {}, FETCHED_AT)
    assert record["raw"]["parseFailed"] is True and record["side"] is None


def test_an_answer_without_a_slip_list_raises():
    for answer in ({"error": "nope"}, {"data": {}}):
        with pytest.raises(ValueError):
            normalize_underdog(answer, [], {}, FETCHED_AT)


def test_a_stored_open_slip_neither_list_returns_closes_with_no_result_once_the_walk_passed_it():
    stored = {"underdog:a": {"id": "underdog:a", "status": "open", "placedAt": "2026-11-06T00:00:00Z", "raw": {}},
              "underdog:b": {"id": "underdog:b", "status": "open", "placedAt": "2026-10-01T00:00:00Z", "raw": {}},
              "underdog:c": {"id": "underdog:c", "status": "won", "placedAt": "2026-11-06T00:00:00Z", "raw": {}}}
    walked_to = datetime(2026, 11, 1, tzinfo=timezone.utc).timestamp() * 1000
    [closed] = closed_by_absence(stored, set(), walked_to, FETCHED_AT)
    assert (closed["id"], closed["status"], closed["closedAt"], closed["raw"]["closedByAbsence"]) == \
        ("underdog:a", "closed", FETCHED_AT, True)
    assert closed_by_absence(stored, {"underdog:a"}, walked_to, FETCHED_AT) == []
    assert {record["id"] for record in closed_by_absence(stored, set(), -1, FETCHED_AT)} == {"underdog:a", "underdog:b"}


# ---- one push (POST /underdog.json) ---------------------------------------------------

def push(active_slips=(), settled_pages=None, complete=False):
    return {"fetchedAt": FETCHED_AT, "active": page(list(active_slips)),
            "settled": settled_pages if settled_pages is not None else [page([])], "complete": complete}


def test_validate_push_takes_a_whole_read_or_an_error_and_names_what_is_wrong():
    body = push([slip("a", ["o-sf"], [179088])])
    assert validate_push(body) == body
    assert validate_push({"error": "  no tab  "}) == {"error": "no tab"}
    assert validate_push({"error": ""}) == "error must be a non-empty string"
    assert validate_push([]) == "body must be an object, got list"
    assert validate_push({**body, "fetchedAt": "soon"}).startswith("fetchedAt must be an ISO time")
    assert validate_push({**body, "active": {"data": {}}}) == "active.data.entry_slips must be a list of objects"
    assert validate_push({**body, "active": None}) == "active must be an answer object with data"
    assert validate_push({**body, "settled": []}) == f"settled must be a list of 1 to {MAX_SETTLED_PAGES} pages"
    assert validate_push({**body, "settled": [page([])] * (MAX_SETTLED_PAGES + 1)}).startswith("settled must be")
    assert validate_push({**body, "settled": [{"data": {"entry_slips": [], "games": {}}}]}) == \
        "settled[0].data.games must be a list, got dict"
    assert validate_push({k: v for k, v in body.items() if k != "complete"}) == "complete must be true or false, got None"
    many = [slip(f"s{n}", ["o-sf"], [179088]) for n in range(501)]
    assert validate_push(push(many)) == "at most 500 slips per push, got 501"


def test_the_page_cap_is_the_extensions_own():
    extension = (Path(__file__).parents[2] / "extension" / "underdog.js").read_text()
    assert f"const MAX_SETTLED_PAGES = {MAX_SETTLED_PAGES};" in extension


def test_cals_49ers_slip_settles_won_off_the_push():
    won = slip("7ed9", ["o-sf"], [179088], status="settled", payout="1.23", payout_at="2026-11-09T00:25:08Z")
    records, closed = records_of_push(push(settled_pages=[page([won], next_page=2)]), {})
    [record] = records
    assert (record["id"], record["status"], record["closedAt"], record["toWin"]) == \
        ("underdog:7ed9", "won", "2026-11-09T00:25:08Z", 0.23)
    assert closed == []


def test_a_push_closes_held_slips_the_pages_passed_and_a_complete_one_closes_every_unseen_one():
    stored = {"underdog:new": {"id": "underdog:new", "status": "open", "placedAt": "2026-11-07T00:00:00Z", "raw": {}},
              "underdog:old": {"id": "underdog:old", "status": "open", "placedAt": "2026-10-01T00:00:00Z", "raw": {}},
              "underdog:seen": {"id": "underdog:seen", "status": "open", "placedAt": "2026-11-07T00:00:00Z", "raw": {}}}
    oldest = slip("x", ["o-sf"], [179088], status="settled", payout="0.0", created_at="2026-11-01T00:00:00Z")
    body = push([slip("seen", ["o-sf"], [179088])], [page([oldest], next_page=2)])
    _, closed = records_of_push(body, stored)
    assert [record["id"] for record in closed] == ["underdog:new"]
    _, closed = records_of_push({**body, "complete": True}, stored)
    assert {record["id"] for record in closed} == {"underdog:new", "underdog:old"}
    _, closed = records_of_push(push([slip("seen", ["o-sf"], [179088])]), stored)
    assert closed == [], "pages with no slips have passed nothing"


# ---- POST /underdog.json over a real HTTP server --------------------------------------

@pytest.fixture
def http_server(tmp_path):
    store = BetsStore(tmp_path / "bets.duckdb", 3650)
    server = ThreadingHTTPServer(("127.0.0.1", 0), service.make_handler(store, 0.0, ["underdog"]))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()
    server.server_close()
    store.close()


def request_json(method, url, body=None):
    data = None if body is None else json.dumps(body).encode()
    request = urllib.request.Request(url, data=data, method=method,
                                     headers={"Content-Type": "application/json"} if data else {})
    try:
        with urllib.request.urlopen(request, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


def test_http_push_stores_the_slip_then_settles_it_and_records_an_error_run(http_server):
    opened = slip("7ed9", ["o-sf"], [179088], created_at="2026-11-07T06:09:41Z")
    assert request_json("POST", f"{http_server}/underdog.json", push([opened])) == \
        (200, {"ok": True, "count": 1, "closed": 0})
    won = slip("7ed9", ["o-sf"], [179088], status="settled", payout="1.23", payout_at="2026-11-09T00:25:08Z")
    assert request_json("POST", f"{http_server}/underdog.json", push([], [page([won])], complete=True)) == \
        (200, {"ok": True, "count": 1, "closed": 0})
    payload = request_json("GET", f"{http_server}/bets.json?days=3650")[1]
    [bet] = [bet for bet in payload["bets"] if bet["venue"] == "underdog"]
    assert (bet["id"], bet["status"], bet["awayTeam"], bet["homeTeam"]) == \
        ("underdog:7ed9", "won", "Las Vegas Raiders", "San Francisco 49ers")
    # A second slip goes open, then a complete push lists it nowhere: closed, no result.
    assert request_json("POST", f"{http_server}/underdog.json", push([slip("b2", ["o-sf"], [179088])])) == \
        (200, {"ok": True, "count": 1, "closed": 0})
    assert request_json("POST", f"{http_server}/underdog.json", push([], [page([won])], complete=True)) == \
        (200, {"ok": True, "count": 1, "closed": 1})
    bets = {bet["id"]: bet for bet in request_json("GET", f"{http_server}/bets.json?days=3650")[1]["bets"]}
    assert (bets["underdog:b2"]["status"], bets["underdog:b2"]["raw"]["closedByAbsence"]) == ("closed", True)
    assert bets["underdog:7ed9"]["status"] == "won"
    assert request_json("POST", f"{http_server}/underdog.json", {"error": "no app.underdogsports.com tab"}) == \
        (200, {"ok": True, "recorded": "error"})
    source = request_json("GET", f"{http_server}/bets.json?days=3650")[1]["sources"]["underdog"]
    assert (source["ok"], source["error"]) == (False, "no app.underdogsports.com tab")
    status, reply = request_json("POST", f"{http_server}/underdog.json", {"active": {}})
    assert (status, reply["error"]) == (400, "fetchedAt must be an ISO time, got None")
