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

def test_pending_first_half_spread_is_kept_open_with_the_side_number(records):
    record = by_id(records, "betonline:900001")
    assert record["status"] == "open" and record["closedAt"] is None
    assert (record["league"], record["betType"], record["period"]) == ("nfl", "spread", "1H")
    assert (record["rotation"], record["side"], record["points"], record["price"]) == (465, "away", 3.5, -110)
    assert (record["awayTeam"], record["homeTeam"]) == ("Chicago Bears", None)
    assert (record["awayKey"], record["homeKey"]) == (None, None)  # the panel resolves keys
    assert (record["stake"], record["toWin"]) == (110, 100)
    assert record["placedAt"] == "2026-09-11T18:12:33Z"
    assert record["eventStart"] is None and record["eventDate"] is None
    assert record["approx"] == ["game_date_unknown", "side_from_rotation_parity"]
    assert record["unmatchable"] is None
    assert (record["source"], record["venue"], record["sourceFetchedAt"]) == ("betonline_api", "betonline", FETCHED_AT)
    assert record["raw"]["description"].startswith("Desktop - NFL - 465")


def test_settled_negative_spread_is_won_with_placed_time_as_closed_at(records):
    record = by_id(records, "betonline:900002")
    assert (record["status"], record["closedAt"]) == ("won", "2026-09-07T15:40:02Z")
    assert (record["side"], record["points"], record["price"], record["period"]) == ("home", -7, -105, "FG")
    assert (record["awayTeam"], record["homeTeam"]) == (None, "Carolina Panthers")


def test_total_reads_over_under_and_the_half_glyph(records):
    record = by_id(records, "betonline:900003")
    assert (record["betType"], record["side"], record["points"], record["price"]) == ("total", "over", 44.5, -110)
    assert record["approx"] == ["game_date_unknown"]  # no parity guess on a total
    assert record["homeTeam"] == "Carolina Panthers"  # rotation 466 is the home team


def test_same_game_parlay_is_one_record_per_leg(records):
    legs = [record for record in records if record["parlayId"] == "betonline:900004"]
    assert [record["id"] for record in legs] == ["betonline:900004:leg0", "betonline:900004:leg1"]
    assert all(record["isParlayLeg"] and record["legCount"] == 2 for record in legs)
    assert [record["legIndex"] for record in legs] == [0, 1]
    first, second = legs
    assert (first["betType"], first["side"], first["points"], first["price"]) == ("spread", "away", 3.5, -110)
    assert (second["betType"], second["side"], second["points"], second["price"]) == ("total", "under", 44.5, -115)
    assert all(record["stake"] == 50 and record["toWin"] == 130 for record in legs)
    assert all(record["raw"]["parlayPrice"] == 260 for record in legs)
    assert not any(record["id"] == "betonline:900004" for record in records)


def test_moneyline_lost_push_and_cancelled_statuses(records):
    lost = by_id(records, "betonline:900005")
    assert (lost["betType"], lost["points"], lost["price"], lost["status"]) == ("moneyline", None, 140, "lost")
    assert by_id(records, "betonline:900007")["status"] == "push"
    assert by_id(records, "betonline:900008")["status"] == "void"


def test_unknown_period_sport_only_prefix_fail_closed_with_a_reason(records):
    period = by_id(records, "betonline:900006")
    assert period["unmatchable"] == "unknown period (1st Period)"
    assert period["betType"] == "other" and period["league"] is None
    sport = by_id(records, "betonline:900009")
    assert sport["unmatchable"] == "league not determined from sport prefix (Football)"


def test_mlb_first_five_period(records):
    record = by_id(records, "betonline:900010")
    assert (record["league"], record["period"], record["points"], record["price"]) == ("mlb", "F5", -0.5, -120)
    assert record["awayTeam"] == "New York Yankees"


def test_every_fixture_record_carries_the_contract_keys(records):
    expected = {"id", "source", "venue", "league", "eventStart", "eventDate", "awayTeam", "homeTeam",
                "awayKey", "homeKey", "rotation", "betType", "period", "side", "points", "price", "stake",
                "toWin", "contracts", "placedAt", "status", "closedAt", "isParlayLeg", "parlayId",
                "legIndex", "legCount", "approx", "unmatchable", "sourceFetchedAt", "raw"}
    assert len(records) == 11
    for record in records:
        assert set(record) == expected, record["id"]


@pytest.mark.parametrize("text, expected", [
    ("465 Chicago Bears +3½ -110", ("spread", "away", 3.5, -110, "FG")),
    ("466 Carolina Panthers -3 -110 - 2nd Half", ("spread", "home", -3, -110, "2H")),
    ("466 Carolina Panthers pk -105", ("spread", "home", 0, -105, "FG")),
    ("465 Chicago Bears over 21½ -110 - 1st Quarter", ("total", "over", 21.5, -110, "1Q")),
    ("465 Chicago Bears +145 for 100.00", ("moneyline", "away", None, 145, "FG")),
])
def test_parse_leg_grammar(text, expected):
    leg = parse_leg(text)
    assert isinstance(leg, dict), leg
    assert (leg["betType"], leg["side"], leg["points"], leg["price"], leg["period"]) == expected


def test_wager_type_disagreeing_with_the_selection_is_unmatchable():
    row = {"Id": 1, "Date": "2026-09-11T18:00:00Z", "Description": "NFL - 465 Chicago Bears +3.5 -110",
           "WagerType": "Total", "WagerStatus": "Pending", "Risk": 10, "ToWin": 9}
    [record] = normalize_row(row, FETCHED_AT)
    assert record["unmatchable"] == "wager type Total but the selection reads as spread"


def test_parlay_with_an_unparseable_leg_is_one_unmatchable_record():
    row = {"Id": 2, "Date": "2026-09-11T18:00:00Z",
           "Description": "NFL - NFL - 465 Chicago Bears +3.5 -110, Justin Fields 200+ passing yards -120",
           "WagerType": "Same Game Parlay", "WagerStatus": "Pending", "Risk": 10, "ToWin": 30}
    [record] = normalize_row(row, FETCHED_AT)
    assert record["id"] == "betonline:2" and not record["isParlayLeg"]
    assert record["unmatchable"].startswith("parlay legs not parsed")


def test_row_without_a_native_id_fails_the_poll_loudly():
    row = {"Date": "2026-09-11T18:00:00Z", "Description": "NFL - 465 Chicago Bears +3.5 -110",
           "WagerType": "Spread", "WagerStatus": "Pending", "Risk": 10, "ToWin": 9}
    with pytest.raises(RuntimeError, match="none of the id fields"):
        normalize_betonline([row], FETCHED_AT)


def test_ticket_number_is_accepted_as_the_native_id():
    row = {"TicketNumber": "AB12", "Date": "2026-09-11T18:00:00Z",
           "Description": "NFL - 465 Chicago Bears +3.5 -110", "WagerType": "Spread",
           "WagerStatus": "Pending", "Risk": 10, "ToWin": 9}
    [record] = normalize_row(row, FETCHED_AT)
    assert record["id"] == "betonline:AB12"


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
    assert len(records) == 11 and records[0]["id"] == "betonline:900001"
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
