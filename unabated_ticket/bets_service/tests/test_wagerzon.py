"""Wagerzon parser on tests/fixtures/bets/wagerzon_history.json, and WagerzonSource's
in-memory session and never-partial rules against a fake session."""
import json
from pathlib import Path

import pytest

from unabated_ticket.bets_service.sources import wagerzon
from unabated_ticket.bets_service.sources.wagerzon import (
    WagerzonSource, normalize_open_bets, normalize_wager, normalize_wagerzon, open_tickets_of, parse_leg, period_for,
    wagers_of_history)

FIXTURE_PATH = Path(__file__).parents[2] / "tests" / "fixtures" / "bets" / "wagerzon_history.json"
FETCHED_AT = "2026-09-23T15:00:00Z"
CONTRACT_KEYS = {"id", "source", "venue", "league", "eventStart", "eventDate", "awayTeam", "homeTeam",
                 "awayKey", "homeKey", "rotation", "betType", "period", "side", "points", "price", "stake",
                 "toWin", "contracts", "placedAt", "status", "closedAt", "isParlayLeg", "parlayId",
                 "legIndex", "legCount", "approx", "unmatchable", "sourceFetchedAt", "raw"}


@pytest.fixture
def history() -> dict:
    return json.loads(FIXTURE_PATH.read_text())


@pytest.fixture
def wagers(history) -> list[dict]:
    return [wager for week in history["weeks"] for wager in wagers_of_history(week)]


@pytest.fixture
def records(wagers) -> list[dict]:
    return normalize_wagerzon(wagers, FETCHED_AT)


def by_id(records: list[dict], record_id: str) -> dict:
    found = [record for record in records if record["id"] == record_id]
    assert len(found) == 1, f"record {record_id} missing"
    return found[0]


# ---- parser ---------------------------------------------------------------------------

def test_pending_mlb_total_carries_its_start_from_the_leg(records):
    record = by_id(records, "wagerzon:300000001")
    assert record["status"] == "open" and record["closedAt"] is None
    assert (record["league"], record["betType"], record["period"]) == ("mlb", "total", "FG")
    assert (record["rotation"], record["side"], record["points"], record["price"]) == (967, "over", 7.5, -120)
    assert (record["awayTeam"], record["homeTeam"]) == ("CHI CUBS", "TB RAYS")
    assert (record["awayKey"], record["homeKey"]) == (None, None)  # the panel resolves keys
    assert (record["stake"], record["toWin"]) == (120, 100)
    assert record["placedAt"] == "2026-09-23T15:02:00Z"  # 11:02 AM Eastern
    assert (record["eventStart"], record["eventDate"]) == ("2026-09-24T23:10:00Z", "2026-09-24")  # 07:10 PM Eastern
    assert record["approx"] == []  # a start is served and a total names both teams
    assert record["unmatchable"] is None
    assert (record["source"], record["venue"], record["sourceFetchedAt"]) == ("wagerzon_api", "wagerzon", FETCHED_AT)
    assert (record["raw"]["nativeId"], record["raw"]["sportCode"], record["raw"]["legResult"]) == ("300000001", "MLB", "")
    assert record["raw"]["gameNumber"] is None


def test_first_five_run_line_is_placed_by_rotation_parity_and_closes_at_its_start(records):
    record = by_id(records, "wagerzon:300000002")
    assert (record["status"], record["league"], record["betType"], record["period"]) == ("won", "mlb", "spread", "F5")
    assert (record["rotation"], record["side"], record["points"], record["price"]) == (969, "away", -0.5, 115)
    assert (record["awayTeam"], record["homeTeam"]) == ("ARI DBACKS", None)  # 969 is odd = away
    assert record["approx"] == ["side_from_rotation_parity"]
    assert (record["eventStart"], record["closedAt"]) == ("2026-09-20T20:10:00Z", "2026-09-20T20:10:00Z")


def test_doubleheader_total_keeps_the_game_number(records):
    record = by_id(records, "wagerzon:300000003")
    assert (record["status"], record["betType"], record["side"], record["points"], record["price"], record["period"]) == \
        ("lost", "total", "over", 5, 100, "F5")
    assert (record["awayTeam"], record["homeTeam"], record["raw"]["gameNumber"]) == ("HOU ASTROS", "BAL ORIOLES", 2)
    assert record["eventStart"] == "2026-09-19T22:35:00Z"


def test_pending_parlay_is_one_open_record_per_leg(records):
    spread, total = [record for record in records if record["parlayId"] == "wagerzon:300000004"]
    assert [spread["id"], total["id"]] == ["wagerzon:300000004:leg0", "wagerzon:300000004:leg1"]
    assert all(record["isParlayLeg"] and record["legCount"] == 2 and record["status"] == "open" for record in (spread, total))
    assert (spread["betType"], spread["side"], spread["points"], spread["price"], spread["awayTeam"]) == \
        ("spread", "away", -1.5, 151, "BOS RED SOX")
    assert (total["betType"], total["side"], total["points"], total["price"], total["period"]) == \
        ("total", "under", 7.5, -105, "FG")
    assert (total["awayTeam"], total["homeTeam"], total["raw"]["gameNumber"]) == ("COL ROCKIES", "NY METS", 1)
    assert all(record["stake"] == 50 and record["toWin"] == 214 and record["raw"]["parlayPrice"] == 428
               for record in (spread, total))
    assert all(record["eventStart"] == "2026-09-24T22:40:00Z" and record["closedAt"] is None for record in (spread, total))
    assert not any(record["id"] == "wagerzon:300000004" for record in records)


def test_nfl_spread_cfb_moneyline_nhl_push_and_cancel(records):
    jaguars = by_id(records, "wagerzon:300000005")
    assert (jaguars["status"], jaguars["league"], jaguars["side"], jaguars["points"], jaguars["price"]) == \
        ("open", "nfl", "away", 2.5, -110)
    assert (jaguars["eventStart"], jaguars["eventDate"]) == ("2026-09-28T17:00:00Z", "2026-09-28")
    moneyline = by_id(records, "wagerzon:300000006")
    assert (moneyline["status"], moneyline["league"], moneyline["betType"], moneyline["side"], moneyline["points"],
            moneyline["price"], moneyline["homeTeam"]) == ("won", "cfb", "moneyline", "home", None, 168, "TEXAS TECH")
    assert moneyline["closedAt"] == "2026-09-19T23:30:00Z"
    push = by_id(records, "wagerzon:300000008")
    assert (push["status"], push["league"], push["side"], push["points"], push["price"]) == ("push", "nhl", "home", -1.5, 180)
    void = by_id(records, "wagerzon:300000010")
    assert (void["status"], void["side"], void["awayTeam"]) == ("void", "away", "DENVER BRONCOS")


def test_props_unsupported_sports_and_a_postponed_leg_fail_closed_with_a_reason(records):
    special = by_id(records, "wagerzon:219337542")
    assert special["unmatchable"] == "not a game market (IdSport PROP)"
    assert (special["status"], special["stake"], special["toWin"], special["betType"]) == ("lost", 250, 4838, "other")
    assert (special["placedAt"], special["closedAt"]) == ("2026-09-20T20:09:00Z", "2026-09-20T20:09:00Z")
    assert special["raw"]["headerDesc"] == "NFL WEEK 2 - SPECIALS"
    assert by_id(records, "wagerzon:219292648")["unmatchable"] == "not a game market (IdSport RBL)"
    assert by_id(records, "wagerzon:300000009")["unmatchable"] == "league not supported (IdSport TNS)"
    postponed, mariners = [record for record in records if record["parlayId"] == "wagerzon:300000007"]
    assert postponed["unmatchable"] == "postponed (( NYM vs COL Has Been Postponed. NO Action ))"
    assert (mariners["unmatchable"], mariners["betType"], mariners["awayTeam"], mariners["status"]) == \
        (None, "total", "SEA MARINERS", "push")


def test_transfer_rows_are_not_bets(wagers, records):
    assert any(wager["WagerOrTrans"] == "TRAN" for wager in wagers)
    assert not any(record["raw"]["nativeId"] == "0" for record in records)


def test_every_fixture_record_carries_the_contract_keys(records):
    # 15 rows: a transfer dropped, two 2-leg parlays, 12 single bets.
    assert len(records) == 16
    for record in records:
        assert set(record) == CONTRACT_KEYS, record["id"]
        assert record["id"].startswith("wagerzon:")


@pytest.mark.parametrize("text, expected", [
    ("[967] TOTAL o7½-120 (CHI CUBS vrs TB RAYS)<BR>( JAMESON TAILLON - R / SHANE MCCLANAHAN - L )",
     ("total", "over", 7.5, -120, None, "CHI CUBS", "TB RAYS", None)),
    ("[969] 1H ARI DBACKS -½+115<BR>( ZAC GALLEN - R / SHOHEI OHTANI - R )",
     ("spread", "away", -0.5, 115, "1H", "ARI DBACKS", None, None)),
    ("[915] TOTAL o5EV (1H HOU ASTROS GM#2 vrs 1H BAL ORIOLES GM#2)",
     ("total", "over", 5, 100, "1H", "HOU ASTROS", "BAL ORIOLES", 2)),
    ("[970] 1H BAL ORIOLES GM#1 -½+120", ("spread", "home", -0.5, 120, "1H", None, "BAL ORIOLES", 1)),
    ("[212] TEXAS TECH +168", ("moneyline", "home", None, 168, None, None, "TEXAS TECH", None)),
    ("[42] BOS BRUINS -1½+180", ("spread", "home", -1.5, 180, None, None, "BOS BRUINS", None)),
    ("[968] TB RAYS +120 (CHI CUBS vrs TB RAYS)", ("moneyline", "home", None, 120, None, "CHI CUBS", "TB RAYS", None)),
    ("[464] HOUSTON TEXANS pk-105", ("spread", "home", 0, -105, None, None, "HOUSTON TEXANS", None)),
])
def test_parse_leg_grammar(text, expected):
    leg = parse_leg(text)
    assert isinstance(leg, dict), leg
    assert (leg["betType"], leg["side"], leg["points"], leg["price"], leg["period"],
            leg["awayTeam"], leg["homeTeam"], leg["gameNumber"]) == expected


@pytest.mark.parametrize("text, reason", [
    ("( NYM vs COL Has Been Postponed. NO Action )", "postponed"),
    ("[1] TOTAL o8-110", "unrecognised selection"),
    ("[1] TOTAL o8-110 (SEA at SD)", "no 'AWAY vrs HOME' bracket"),
    ("[1] TOTAL o8-110 (1H SEA MARINERS vrs 2H SD PADRES)", "periods differ (1H vs 2H)"),
    ("[1] TOTAL o8-110 (SEA MARINERS GM#1 vrs SD PADRES GM#2)", "doubleheader games differ (1 vs 2)"),
    ("[1] 3H SEA MARINERS -1½+120", "unknown period (3H)"),
    ("[1] TB RAYS +120 (CHI CUBS vrs NY YANKEES)", "team TB RAYS is not in its bracket"),
    ("TB RAYS +120", "no rotation to place the side of TB RAYS"),
    ("Mark Andrews (BAL) 1+ Touchdowns <br /> New Orleans Saints vs Baltimore Ravens", "unrecognised selection"),
])
def test_parse_leg_fails_closed_with_the_reason(text, reason):
    result = parse_leg(text)
    assert isinstance(result, str) and result.startswith(reason), result


def test_period_for_maps_mlb_first_half_to_first_five():
    assert period_for("mlb", "1H") == "F5"
    assert period_for("nfl", "1H") == "1H"
    assert period_for("mlb", None) == "FG"


def test_wager_without_an_id_fails_the_poll_loudly():
    wager = {"WagerOrTrans": "WAGER", "IdWager": 0, "Result": "", "details": [], "RiskAmount": "1", "WinAmount": "1"}
    with pytest.raises(RuntimeError, match="carries no IdWager"):
        normalize_wager(wager, FETCHED_AT)


def test_wager_without_legs_is_unmatchable_not_dropped():
    wager = {"WagerOrTrans": "WAGER", "IdWager": 5, "Result": "", "details": [], "RiskAmount": "10", "WinAmount": "9",
             "PlacedDate": "09/23/2026", "PlacedTime": "10:00 AM"}
    [record] = normalize_wager(wager, FETCHED_AT)
    assert (record["id"], record["unmatchable"], record["status"]) == ("wagerzon:5", "wager carries no legs", "open")


# ---- open bets ------------------------------------------------------------------------

@pytest.fixture
def open_records(history) -> list[dict]:
    return normalize_open_bets(open_tickets_of(history["openBets"]), FETCHED_AT)


def test_open_tickets_of_groups_rows_by_ticket_and_refuses_unknown_shapes(history):
    tickets = open_tickets_of(history["openBets"])
    assert [(rows[0]["TicketNumber"], len(rows)) for rows in tickets] == \
        [(400000001, 1), (400000002, 2), (400000003, 1), (300000005, 1)]
    with pytest.raises(RuntimeError, match=r"expected \{result: \[...\]\}, got \['rows'\]"):
        open_tickets_of({"rows": []})
    with pytest.raises(RuntimeError, match=r"unknown shape \(keys \['TicketNumber', 'legs'\]\)"):
        open_tickets_of({"result": [{"TicketNumber": 1, "legs": []}]})
    assert open_tickets_of({"result": []}) == []


def test_open_mlb_total_reads_its_start_from_the_helper_row(open_records):
    record = by_id(open_records, "wagerzon:400000001")
    assert (record["status"], record["closedAt"], record["unmatchable"]) == ("open", None, None)
    assert (record["league"], record["betType"], record["period"], record["side"], record["points"], record["price"]) == \
        ("mlb", "total", "FG", "over", 7.5, -120)
    assert (record["awayTeam"], record["homeTeam"], record["rotation"]) == ("CHI CUBS", "TB RAYS", 967)
    assert (record["eventStart"], record["eventDate"], record["placedAt"]) == \
        ("2026-09-24T23:10:00Z", "2026-09-24", "2026-09-23T15:02:00Z")
    assert (record["stake"], record["toWin"], record["approx"]) == (120, 100, [])
    assert record["raw"]["openBet"] is True and record["raw"]["rotationNumbers"] == "967"


def test_open_parlay_rows_of_one_ticket_are_its_legs(open_records):
    spread, total = [record for record in open_records if record["parlayId"] == "wagerzon:400000002"]
    assert [spread["id"], total["id"]] == ["wagerzon:400000002:leg0", "wagerzon:400000002:leg1"]
    assert all(record["status"] == "open" and record["legCount"] == 2 and record["stake"] == 50
               and record["raw"]["parlayPrice"] == 428 for record in (spread, total))
    assert (spread["betType"], spread["side"], spread["points"], spread["price"]) == ("spread", "away", -1.5, 151)
    assert (total["betType"], total["side"], total["points"], total["raw"]["gameNumber"]) == ("total", "under", 7.5, 1)
    assert by_id(open_records, "wagerzon:400000003")["unmatchable"] == "not a game market (IdSport RBL)"
    assert len(open_records) == 5


# ---- network half ---------------------------------------------------------------------

class FakeResponse:
    def __init__(self, status_code: int, body: object = None, text: str = "", content_type: str = "text/html",
                 headers: dict | None = None):
        self.status_code = status_code
        self._body = body
        self.text = text
        self.headers = {"Content-Type": content_type, **(headers or {})}

    def json(self) -> object:
        if self._body is None:
            raise ValueError("not JSON")
        return self._body


LOGIN_PAGE = ('<form><input type="hidden" name="__VIEWSTATE" id="__VIEWSTATE" value="vs-1" />'
              '<input type="hidden" name="__VIEWSTATEGENERATOR" id="__VIEWSTATEGENERATOR" value="gen" />'
              '<input type="hidden" name="__EVENTVALIDATION" id="__EVENTVALIDATION" value="ev" /></form>')
EMPTY_WEEK = {"details": [], "ErrorMsg": ""}


class FakeSession:
    """Answers the login page, the form POST and the two helpers; a logged-out
    session gets the site's redirect to the root instead of JSON."""

    def __init__(self, weeks: list[dict], open_rows: list[dict], login_takes: bool = True):
        self.weeks = weeks
        self.open_rows = open_rows
        self.login_takes = login_takes
        self.logged_in = False
        self.logins = 0
        self.posted_fields: list[dict] = []
        self.history_calls: list[tuple[int, dict]] = []
        self.open_calls = 0
        self.error_week: int | None = None

    def get(self, url: str, params=None, headers=None, timeout=None, allow_redirects=True) -> FakeResponse:
        if url == wagerzon.BASE_URL:
            return FakeResponse(200, text=LOGIN_PAGE)
        assert headers == wagerzon.XHR_HEADERS and allow_redirects is False
        if not self.logged_in:
            return FakeResponse(302, headers={"Location": "/"})
        if url == wagerzon.HISTORY_URL:
            week = params["week"]
            self.history_calls.append((week, dict(headers)))
            if week == self.error_week:
                return FakeResponse(200, {"result": {"details": [], "ErrorMsg": "Session expired"}}, content_type="application/json")
            body = self.weeks[week] if week < len(self.weeks) else EMPTY_WEEK
            return FakeResponse(200, {"result": body}, content_type="application/json; charset=utf-8")
        if url == wagerzon.OPEN_BETS_URL:
            self.open_calls += 1
            return FakeResponse(200, {"result": self.open_rows}, content_type="application/json; charset=utf-8")
        raise AssertionError(f"unexpected GET {url}")

    def post(self, url: str, data=None, timeout=None) -> FakeResponse:
        assert url == wagerzon.BASE_URL
        self.logins += 1
        self.posted_fields.append(dict(data))
        self.logged_in = self.login_takes
        return FakeResponse(200, text="<html>welcome</html>")


def make_source(session: FakeSession) -> WagerzonSource:
    return WagerzonSource(username="user", password="secret", history_weeks=6, poll_sec=300,
                          session_factory=lambda: session)


def test_fetch_logs_in_once_pulls_six_weeks_and_the_open_list_which_wins_on_a_shared_ticket(history):
    session = FakeSession(history["weeks"], history["openBets"]["result"])
    source = make_source(session)
    records = source.fetch()
    assert session.logins == 1
    assert session.posted_fields[0]["__VIEWSTATE"] == "vs-1" and session.posted_fields[0]["Account"] == "user"
    assert session.posted_fields[0]["Password"] == "secret" and session.posted_fields[0]["BtnSubmit"] == ""
    assert [week for week, _ in session.history_calls] == [0, 1, 2, 3, 4, 5] and session.open_calls == 1
    assert len(records) == 20  # 16 from the weeks + 5 open, the pending Jaguars spread listed by both
    jaguars = by_id(records, "wagerzon:300000005")
    assert jaguars["status"] == "open" and jaguars["raw"]["openBet"] is True
    assert by_id(records, "wagerzon:400000001")["status"] == "open"


def test_session_is_kept_across_polls_and_a_dead_one_is_replaced_by_one_login(history):
    session = FakeSession(history["weeks"], [])
    source = make_source(session)
    source.fetch()
    source.fetch()
    assert session.logins == 1
    session.logged_in = False  # the site dropped the session between polls
    records = source.fetch()
    assert session.logins == 2 and len(records) == 16


def test_a_login_that_does_not_take_raises_after_one_retry(history):
    session = FakeSession(history["weeks"], [], login_takes=False)
    with pytest.raises(RuntimeError, match="returned no JSON after a fresh login"):
        make_source(session).fetch()
    assert session.logins == 2


def test_an_error_message_in_a_week_fails_the_poll(history):
    session = FakeSession(history["weeks"], [])
    session.error_week = 3
    with pytest.raises(RuntimeError, match="HistoryHelper week 3: Session expired"):
        make_source(session).fetch()


def test_open_bets_rows_of_unknown_shape_fail_the_poll(history):
    session = FakeSession(history["weeks"], [{"WagerId": 1, "Legs": []}])
    with pytest.raises(RuntimeError, match=r"unknown shape \(keys \['Legs', 'WagerId'\]\)"):
        make_source(session).fetch()
    session = FakeSession(history["weeks"], [{"TicketNumber": 0, "DetailDescription": "x"}])
    with pytest.raises(RuntimeError, match="carries no TicketNumber"):
        make_source(session).fetch()


def test_missing_credentials_name_the_env_keys(history):
    source = WagerzonSource(username="", password="", history_weeks=6, poll_sec=300,
                            session_factory=lambda: FakeSession(history["weeks"], []))
    with pytest.raises(RuntimeError, match="WAGERZONC_USERNAME and WAGERZONC_PASSWORD"):
        source.fetch()


def test_source_if_configured_needs_both_credentials(monkeypatch):
    monkeypatch.setattr(wagerzon.config, "WAGERZON_USERNAME", "user")
    monkeypatch.setattr(wagerzon.config, "WAGERZON_PASSWORD", None)
    assert wagerzon.source_if_configured() is None
    monkeypatch.setattr(wagerzon.config, "WAGERZON_PASSWORD", "secret")
    source = wagerzon.source_if_configured()
    assert isinstance(source, WagerzonSource) and (source.name, source.poll_sec) == ("wagerzon", 300.0)
