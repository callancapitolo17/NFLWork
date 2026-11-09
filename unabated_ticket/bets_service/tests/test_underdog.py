"""Underdog (2026-11-07): the prediction slips' records (sources/underdog.py), the
source's settled-page walk, and the auth's refresh and rotation
(sources/underdog_auth.py). The active slip is Cal's live $1 49ers moneyline
(trimmed); the settled grades are built on the same grammar, since no prediction
slip had settled when this was written. No network."""
import json
from datetime import datetime, timezone

import pytest

from unabated_ticket.bets_service.sources import underdog_auth
from unabated_ticket.bets_service.sources.underdog import UnderdogSource, closed_by_absence, normalize_underdog
from unabated_ticket.bets_service.sources.underdog_auth import UnderdogAuth, UnderdogAuthError

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


# ---- the network half, on a scripted session ----------------------------------------

class Reply:
    def __init__(self, status, payload):
        self.status_code, self._payload = status, payload
        self.text = json.dumps(payload)

    def json(self):
        return self._payload


class ScriptedSession:
    """Answers each API GET from `replies` by (path, page)."""

    def __init__(self, replies):
        self.replies, self.calls = replies, []

    def get(self, url, params=None, headers=None, timeout=None):
        key = (url.split(".com", 1)[1], (params or {}).get("page"))
        self.calls.append((key, headers.get("Authorization")))
        return self.replies[key].pop(0) if isinstance(self.replies[key], list) else self.replies[key]


class FixedAuth:
    def __init__(self):
        self.forgotten = 0

    def token(self):
        return underdog_auth.AccessToken(f"tok{self.forgotten}", 1e12, "auth0|cal")

    def forget(self):
        self.forgotten += 1


NOW = datetime(2026, 11, 7, 6, 20, tzinfo=timezone.utc).timestamp()
OLD = "2026-08-01T00:00:00Z"


def source(replies, stored=None):
    session = ScriptedSession(replies)
    return UnderdogSource(FixedAuth(), lambda: dict(stored or {}), session=session, clock=lambda: NOW), session


ACTIVE = ("/v10/user/active_entry_slips", None)


def settled(number):
    return ("/v10/user/settled_entry_slips", number)


def test_the_first_poll_reads_back_the_history_window_and_later_polls_read_one_settled_page():
    recent, old = "2026-11-01T00:00:00Z", OLD
    underdog, session = source({
        ACTIVE: Reply(200, page([slip("s1", ["o-sf"], [179088])])),
        settled(1): [Reply(200, page([slip("r", ["o-sf"], [179088], status="settled", created_at=recent)], next_page=2)),
                     Reply(200, page([slip("r", ["o-sf"], [179088], status="settled", created_at=recent)], next_page=2))],
        settled(2): Reply(200, page([slip("o", ["o-sf"], [179088], status="settled", created_at=old)], next_page=3))})
    assert {record["id"] for record in underdog.fetch()} == {"underdog:s1", "underdog:r", "underdog:o"}
    assert [call[0] for call in session.calls] == [ACTIVE, settled(1), settled(2)]
    assert session.calls[0][1] == "Bearer tok0"
    session.calls.clear()
    underdog.fetch()
    assert [call[0] for call in session.calls] == [ACTIVE, settled(1)]


def test_a_slip_still_held_open_is_chased_down_the_settled_pages():
    stored = {"underdog:late": {"id": "underdog:late", "status": "open", "placedAt": OLD, "raw": {}}}
    underdog, session = source({
        ACTIVE: Reply(200, page([])),
        settled(1): Reply(200, page([slip("x", ["o-sf"], [179088], status="settled", created_at="2026-11-05T00:00:00Z")],
                                    next_page=2)),
        settled(2): Reply(200, page([slip("late", ["o-sf"], [179088], status="settled", payout="1.23", created_at=OLD)],
                                    next_page=3)),
    }, stored)
    underdog._read_history = True
    records = by_id(underdog.fetch())
    assert records["underdog:late"]["status"] == "won"
    assert [call[0] for call in session.calls] == [ACTIVE, settled(1), settled(2)]


def test_a_stored_open_slip_the_walk_passes_without_seeing_is_closed():
    stored = {"underdog:gone": {"id": "underdog:gone", "status": "open", "placedAt": "2026-11-06T00:00:00Z", "raw": {}}}
    underdog, session = source({
        ACTIVE: Reply(200, page([])),
        settled(1): Reply(200, page([slip("x", ["o-sf"], [179088], status="settled", created_at=OLD)], next_page=2))},
        stored)
    underdog._read_history = True
    records = by_id(underdog.fetch())
    assert records["underdog:gone"]["status"] == "closed"
    assert [call[0] for call in session.calls] == [ACTIVE, settled(1)]


def test_a_refused_token_is_refreshed_once_and_a_failed_read_raises():
    underdog, session = source({ACTIVE: [Reply(401, {}), Reply(200, page([]))], settled(1): Reply(503, {"error": "down"})})
    with pytest.raises(RuntimeError, match="HTTP 503"):
        underdog.fetch()
    assert [call[1] for call in session.calls[:2]] == ["Bearer tok0", "Bearer tok1"]


# ---- auth -----------------------------------------------------------------------------

def jwt(claims):
    import base64
    body = base64.urlsafe_b64encode(json.dumps(claims).encode()).rstrip(b"=").decode()
    return f"h.{body}.s"


def test_refresh_saves_a_rotated_token_and_caches_until_near_expiry(tmp_path):
    path = tmp_path / "underdog_token.json"
    underdog_auth.save_token_file(path, "r1", "auth0|cal")
    sent = []

    def request(body):
        sent.append(body)
        return {"access_token": jwt({"sub": "auth0|cal"}), "refresh_token": f"r{len(sent) + 1}", "expires_in": 600}

    clock = [1000.0]
    auth = UnderdogAuth(path, request=request, clock=lambda: clock[0])
    auth.token()
    auth.token()
    assert len(sent) == 1 and sent[0]["grant_type"] == "refresh_token" and sent[0]["refresh_token"] == "r1"
    assert json.loads(path.read_text())["refresh_token"] == "r2"
    clock[0] += 600 - underdog_auth.REFRESH_MARGIN_SEC + 1
    auth.token()
    assert sent[1]["refresh_token"] == "r2"


def test_connect_uses_the_password_grant_and_stores_only_the_refresh_token(tmp_path):
    path = tmp_path / "underdog_token.json"
    sent = []

    def request(body):
        sent.append(body)
        return {"access_token": jwt({"sub": "auth0|cal"}), "refresh_token": "r1"}

    assert underdog_auth.connect(path, "cal@example.com", "hunter2", request=request) == "auth0|cal"
    assert sent[0]["grant_type"] == "password" and sent[0]["audience"] == "https://api.underdogfantasy.com"
    stored = path.read_text()
    assert "hunter2" not in stored and json.loads(stored)["refresh_token"] == "r1"


def test_no_token_file_is_a_clear_error(tmp_path):
    with pytest.raises(UnderdogAuthError, match="underdog_auth connect"):
        UnderdogAuth(tmp_path / "missing.json").token()
