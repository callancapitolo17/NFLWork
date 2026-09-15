"""BetOnline parser on tests/fixtures/bets/betonline_history.json, and
BetOnlineSource's token-rotation and never-partial rules against a fake session."""
import json
from pathlib import Path

import pytest

from unabated_ticket.bets_service.sources import betonline
from unabated_ticket.bets_service.sources.betonline import (
    BetOnlineSource, normalize_betonline, normalize_row, parse_leg)

FIXTURE_PATH = Path(__file__).parents[2] / "tests" / "fixtures" / "bets" / "betonline_history.json"
FETCHED_AT = "2026-09-11T21:00:00Z"


@pytest.fixture
def report() -> dict:
    return json.loads(FIXTURE_PATH.read_text())


@pytest.fixture
def records(report) -> list[dict]:
    return normalize_betonline(report["Data"], FETCHED_AT)


def by_id(records: list[dict], record_id: str) -> dict:
    found = [record for record in records if record["id"] == record_id]
    assert len(found) == 1, f"record {record_id} missing"
    return found[0]


# ---- parser ---------------------------------------------------------------------------

def test_pending_nfl_total_names_both_teams_and_stays_open(records):
    record = by_id(records, "betonline:995909271-1")
    assert record["status"] == "open" and record["closedAt"] is None
    assert (record["league"], record["betType"], record["period"]) == ("nfl", "total", "FG")
    assert (record["rotation"], record["side"], record["points"], record["price"]) == (273, "over", 39.5, -108)
    assert (record["awayTeam"], record["homeTeam"]) == ("Seattle Seahawks", "Arizona Cardinals")
    assert (record["awayKey"], record["homeKey"]) == (None, None)  # the panel resolves keys
    assert (record["stake"], record["toWin"]) == (540, 500)
    assert record["placedAt"] == "2026-09-14T00:19:21Z"  # naive report time read as UTC-8
    assert record["eventStart"] is None and record["eventDate"] is None
    assert record["approx"] == ["game_date_unknown"]  # no parity guess on a total
    assert record["unmatchable"] is None
    assert (record["source"], record["venue"], record["sourceFetchedAt"]) == ("betonline_api", "betonline", FETCHED_AT)
    assert record["raw"]["ticketNumber"] == 995909271 and record["raw"]["wagerNumber"] == 1


def test_college_first_half_total_is_cfb_with_the_home_rotation(records):
    record = by_id(records, "betonline:995445752-1")
    assert (record["league"], record["period"], record["side"], record["points"], record["price"]) == \
        ("cfb", "1H", "under", 26.5, 105)
    assert (record["awayTeam"], record["homeTeam"], record["rotation"]) == ("Tulsa", "Sam Houston St", 374)
    assert (record["status"], record["closedAt"]) == ("won", "2026-09-12T19:33:46Z")  # placed time: a lower bound


def test_one_team_spread_is_placed_by_rotation_parity(records):
    record = by_id(records, "betonline:993434702-1")
    assert (record["betType"], record["side"], record["points"], record["price"]) == ("spread", "home", 15, 104)
    assert (record["awayTeam"], record["homeTeam"]) == (None, "San Diego")  # 307122 is even = home
    assert record["approx"] == ["game_date_unknown", "side_from_rotation_parity"]
    assert record["status"] == "lost"
    negative = by_id(records, "betonline:996000003-1")
    assert (negative["side"], negative["points"], negative["homeTeam"]) == ("home", -3, "Green Bay Packers")


def test_same_game_parlay_is_one_record_per_leg(records):
    legs = [record for record in records if record["parlayId"] == "betonline:996000001-1"]
    assert [record["id"] for record in legs] == ["betonline:996000001-1:leg0", "betonline:996000001-1:leg1"]
    assert all(record["isParlayLeg"] and record["legCount"] == 2 for record in legs)
    first, second = legs
    assert (first["betType"], first["side"], first["points"], first["price"], first["period"]) == \
        ("spread", "away", 3.5, -110, "1H")
    assert first["awayTeam"] == "Chicago Bears" and first["league"] == "nfl"
    assert (second["betType"], second["side"], second["points"], second["price"], second["period"]) == \
        ("total", "under", 44.5, -115, "FG")
    assert all(record["stake"] == 50 and record["toWin"] == 130 for record in legs)
    assert all(record["raw"]["parlayPrice"] == 260 for record in legs)
    assert not any(record["id"] == "betonline:996000001-1" for record in records)


def test_moneyline_push_and_cancelled_statuses(records):
    push = by_id(records, "betonline:996000002-1")
    assert (push["betType"], push["points"], push["price"], push["status"], push["side"]) == \
        ("moneyline", None, 140, "push", "away")
    assert by_id(records, "betonline:996000003-1")["status"] == "void"


def test_unknown_period_fails_closed_with_a_reason(records):
    record = by_id(records, "betonline:996000004-1")
    assert record["unmatchable"] == "unknown period (1ST PERIOD)"
    assert record["betType"] == "other" and record["league"] is None


def test_basketball_quarter_and_baseball_first_five(records):
    quarter = by_id(records, "betonline:996000005-1")
    assert (quarter["league"], quarter["period"], quarter["points"]) == ("nba", "1Q", -4.5)
    first_five = by_id(records, "betonline:996000006-1")
    assert (first_five["league"], first_five["period"], first_five["points"], first_five["price"]) == \
        ("mlb", "F5", -0.5, -120)
    assert first_five["awayTeam"] == "New York Yankees"


def test_every_fixture_record_carries_the_contract_keys(records):
    expected = {"id", "source", "venue", "league", "eventStart", "eventDate", "awayTeam", "homeTeam",
                "awayKey", "homeKey", "rotation", "betType", "period", "side", "points", "price", "stake",
                "toWin", "contracts", "placedAt", "status", "closedAt", "isParlayLeg", "parlayId",
                "legIndex", "legCount", "approx", "unmatchable", "sourceFetchedAt", "raw"}
    assert len(records) == 15  # 14 rows, the parlay is two
    for record in records:
        assert set(record) == expected, record["id"]


@pytest.mark.parametrize("text, expected", [
    ("465 Chicago Bears +3½ -110 for GAME", ("spread", "away", 3.5, -110, "FG")),
    ("466 Carolina Panthers -3 -110 for 2ND HALF", ("spread", "home", -3, -110, "2H")),
    ("466 Carolina Panthers pk -105 for GAME", ("spread", "home", 0, -105, "FG")),
    ("465 Chicago Bears/Carolina Panthers over 21½ -110 for 1ST QUARTER", ("total", "over", 21.5, -110, "1Q")),
    ("465 Chicago Bears +145 for game ", ("moneyline", "away", None, 145, "FG")),
])
def test_parse_leg_grammar(text, expected):
    leg = parse_leg(text)
    assert isinstance(leg, dict), leg
    assert (leg["betType"], leg["side"], leg["points"], leg["price"], leg["period"]) == expected


def test_leg_without_the_for_period_tail_does_not_parse():
    assert parse_leg("465 Chicago Bears +3.5 -110").startswith("no '<rotation> <selection> for <period>'")


def test_wager_type_disagreeing_with_the_selection_is_unmatchable():
    row = {"Id": "1-1", "Date": "2026-09-11T18:00:00", "Description": "FOOTBALL - 465 Chicago Bears +3.5 -110 for GAME",
           "WagerType": "Total", "WagerStatus": "Pending", "Risk": 10, "ToWin": 9}
    [record] = normalize_row(row, FETCHED_AT)
    assert record["unmatchable"] == "wager type Total but the selection reads as spread"


def test_sport_outside_the_scanner_is_unmatchable():
    row = {"Id": "1-1", "Date": "2026-09-11T18:00:00", "Description": "Desktop - TENNIS - 1 Alcaraz -200 for MATCH",
           "WagerType": "Money Line", "WagerStatus": "Pending", "Risk": 10, "ToWin": 5}
    [record] = normalize_row(row, FETCHED_AT)
    assert record["unmatchable"] == "unknown sport prefix (TENNIS)"


def test_parlay_with_an_unparseable_leg_is_one_unmatchable_record():
    row = {"Id": "2-1", "Date": "2026-09-11T18:00:00",
           "Description": "FOOTBALL - FOOTBALL - 465 Chicago Bears +3.5 -110 for GAME, Justin Fields 200+ passing yards -120",
           "WagerType": "Same Game Parlay", "WagerStatus": "Pending", "Risk": 10, "ToWin": 30}
    [record] = normalize_row(row, FETCHED_AT)
    assert record["id"] == "betonline:2-1" and not record["isParlayLeg"]
    assert record["unmatchable"].startswith("parlay legs not parsed")


def test_row_without_a_native_id_fails_the_poll_loudly():
    row = {"Date": "2026-09-11T18:00:00", "Description": "FOOTBALL - 465 Chicago Bears +3.5 -110 for GAME",
           "WagerType": "Spread", "WagerStatus": "Pending", "Risk": 10, "ToWin": 9}
    with pytest.raises(RuntimeError, match="none of the id fields"):
        normalize_betonline([row], FETCHED_AT)


# ---- network half ---------------------------------------------------------------------

class FakeResponse:
    def __init__(self, status_code: int, body: object):
        self.status_code = status_code
        self._body = body

    def json(self) -> object:
        return self._body


class FakeSession:
    """Answers the token and report POSTs; records every call."""

    def __init__(self, report: dict, token_status: int = 200, report_status: int = 200,
                 expires_in: int = 300, rotate: bool = True):
        self.report = report
        self.token_status = token_status
        self.report_status = report_status
        self.expires_in = expires_in
        self.rotate = rotate
        self.calls: list[tuple[str, dict]] = []
        self.refresh_tokens_seen: list[str] = []

    def post(self, url: str, data=None, json=None, headers=None, timeout=None) -> FakeResponse:
        self.calls.append((url, data or json))
        if url == betonline.TOKEN_URL:
            self.refresh_tokens_seen.append(data["refresh_token"])
            if self.token_status != 200:
                return FakeResponse(self.token_status, {"error": "invalid_grant"})
            body = {"access_token": f"access-{len(self.refresh_tokens_seen)}", "expires_in": self.expires_in}
            if self.rotate:
                body["refresh_token"] = f"rotated-{len(self.refresh_tokens_seen)}"
            return FakeResponse(200, body)
        if url == betonline.BET_HISTORY_URL:
            if self.report_status != 200:
                return FakeResponse(self.report_status, "unavailable")
            assert headers["Authorization"].startswith("Bearer access-")
            start = json["StartPosition"]
            page = self.report["Data"][start:start + json["TotalPerPage"]]
            return FakeResponse(200, {"TotalRows": len(self.report["Data"]), "Data": page})
        raise AssertionError(f"unexpected url {url}")


def write_cookie_file(path: Path, refresh_token: str = "initial") -> None:
    path.write_text(json.dumps([
        {"name": "cf_clearance", "value": "cf", "domain": ".betonline.ag", "path": "/"},
        {"name": "krefresh", "value": refresh_token, "domain": "www.betonline.ag", "path": "/"},
    ]))


def make_source(tmp_path: Path, session: FakeSession, clock):
    cookies_path = tmp_path / "recon_betonline_cookies.json"
    write_cookie_file(cookies_path)
    return cookies_path, BetOnlineSource(cookies_path=cookies_path, poll_sec=300, history_days=31,
                                         session_factory=lambda cookies: session, clock=clock)


def test_fetch_normalises_the_paged_report_and_rotates_the_refresh_token(tmp_path, report):
    session = FakeSession(report)
    now = [1_800_000_000.0]
    cookies_path, source = make_source(tmp_path, session, clock=lambda: now[0])
    records = source.fetch()
    assert len(records) == 15 and records[0]["id"] == "betonline:995909271-1"
    assert session.refresh_tokens_seen == ["initial"]
    saved = json.loads(cookies_path.read_text())
    assert [cookie["value"] for cookie in saved if cookie["name"] == "krefresh"] == ["rotated-1"]
    assert not (tmp_path / "recon_betonline_cookies.json.tmp").exists()
    report_calls = [body for url, body in session.calls if url == betonline.BET_HISTORY_URL]
    assert report_calls[0]["StartDate"].endswith("T00:00:00.000Z") and report_calls[0]["StartPosition"] == 0


def test_access_token_is_refreshed_only_within_60s_of_expiry(tmp_path, report):
    session = FakeSession(report, expires_in=300)
    now = [1_800_000_000.0]
    _, source = make_source(tmp_path, session, clock=lambda: now[0])
    source.fetch()
    now[0] += 200  # 100 s of life left: no refresh
    source.fetch()
    assert session.refresh_tokens_seen == ["initial"]
    now[0] += 50  # 50 s left, inside the 60 s margin: refresh, using the ROTATED token from the file
    source.fetch()
    assert session.refresh_tokens_seen == ["initial", "rotated-1"]


def test_refresh_re_reads_a_rotation_made_by_the_launchagent(tmp_path, report):
    session = FakeSession(report, expires_in=10)
    now = [1_800_000_000.0]
    cookies_path, source = make_source(tmp_path, session, clock=lambda: now[0])
    source.fetch()
    write_cookie_file(cookies_path, refresh_token="agent-rotated")  # another process rotated it
    now[0] += 60
    source.fetch()
    assert session.refresh_tokens_seen[-1] == "agent-rotated"


def test_dead_refresh_token_raises_and_leaves_the_file_alone(tmp_path, report):
    session = FakeSession(report, token_status=400)
    cookies_path, source = make_source(tmp_path, session, clock=lambda: 0.0)
    with pytest.raises(RuntimeError, match="HTTP 400, invalid_grant"):
        source.fetch()
    assert "initial" in cookies_path.read_text()


def test_failed_report_page_raises_instead_of_returning_a_partial_list(tmp_path, report):
    session = FakeSession(report, report_status=503)
    _, source = make_source(tmp_path, session, clock=lambda: 0.0)
    with pytest.raises(RuntimeError, match="report page 0 failed: HTTP 503"):
        source.fetch()


def test_missing_cookie_file_names_the_recon_script(tmp_path, report):
    source = BetOnlineSource(cookies_path=tmp_path / "missing.json", poll_sec=300, history_days=31,
                             session_factory=lambda cookies: FakeSession(report))
    with pytest.raises(RuntimeError, match="recon_betonline.py"):
        source.fetch()
