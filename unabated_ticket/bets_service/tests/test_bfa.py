"""BFA parser on tests/fixtures/bets/bfa_history.json, and BFASource's in-memory
Keycloak session and never-partial rules against a fake session."""
import base64
import json
from datetime import datetime, timezone
from pathlib import Path

import pytest

from unabated_ticket.bets_service.sources import bfa
from unabated_ticket.bets_service.sources.bfa import (
    BFASource, game_rotation_of, merge_open_over_history, normalize_bfa, normalize_open_bets, normalize_wager,
    parse_leg)
from unabated_ticket.bets_service.store import BetsStore

FIXTURE_PATH = Path(__file__).parents[2] / "tests" / "fixtures" / "bets" / "bfa_history.json"
FETCHED_AT = "2026-09-23T05:00:00Z"
# 2026-09-22 21:46:36 on the account's Pacific clock.
CLOCK = datetime(2026, 9, 23, 4, 46, 36, tzinfo=timezone.utc).timestamp()
CONTRACT_KEYS = {"id", "source", "venue", "league", "eventStart", "eventDate", "awayTeam", "homeTeam",
                 "awayKey", "homeKey", "rotation", "betType", "period", "side", "points", "price", "stake",
                 "toWin", "contracts", "placedAt", "status", "closedAt", "isParlayLeg", "parlayId",
                 "legIndex", "legCount", "approx", "unmatchable", "sourceFetchedAt", "raw"}


@pytest.fixture
def history() -> dict:
    return json.loads(FIXTURE_PATH.read_text())


@pytest.fixture
def records(history) -> list[dict]:
    return normalize_bfa(history["wagers"], FETCHED_AT)


def by_id(records: list[dict], record_id: str) -> dict:
    found = [record for record in records if record["id"] == record_id]
    assert len(found) == 1, f"record {record_id} missing"
    return found[0]


# ---- parser ---------------------------------------------------------------------------

def test_pending_nfl_first_half_total_carries_its_kickoff_from_settled_date(records):
    record = by_id(records, "bfa:990000001")
    assert record["status"] == "open" and record["closedAt"] is None
    assert (record["league"], record["betType"], record["period"]) == ("nfl", "total", "1H")
    assert (record["rotation"], record["side"], record["points"], record["price"]) == (455, "over", 24.5, -110)
    assert (record["awayTeam"], record["homeTeam"]) == ("JACKSONVILLE JAGUARS", "CINCINNATI BENGALS")
    assert (record["awayKey"], record["homeKey"]) == (None, None)  # the panel resolves keys
    assert (record["stake"], record["toWin"]) == (220, 200)
    assert record["placedAt"] == "2026-09-27T14:12:03Z"  # 07:12 Pacific
    assert (record["eventStart"], record["eventDate"]) == ("2026-09-27T17:00:00Z", "2026-09-27")  # 10:00 Pacific = 1 PM ET
    assert record["approx"] == ["event_start_from_settled_date"]  # no parity guess on a total
    assert record["unmatchable"] is None
    assert (record["source"], record["venue"], record["sourceFetchedAt"]) == ("bfa_api", "bfa", FETCHED_AT)
    assert record["raw"]["nativeId"] == "990000001" and record["raw"]["type"] == "STRAIGHT BET"


def test_a_dotnet_default_settled_date_means_no_start_and_the_dateless_window(records):
    record = by_id(records, "bfa:990000002")
    assert record["status"] == "open"  # an empty result is a pending bet
    assert (record["eventStart"], record["eventDate"]) == (None, None)
    assert record["approx"] == ["game_date_unknown", "side_from_rotation_parity"]
    assert (record["betType"], record["side"], record["points"], record["price"]) == ("spread", "home", 7, -115)
    assert (record["awayTeam"], record["homeTeam"]) == (None, "MIAMI DOLPHINS")  # 458 is even = home


def test_college_bets_are_league_unknown_and_never_guessed(records):
    for record_id in ("bfa:354669433", "bfa:354693777", "bfa:354693779", "bfa:354693780", "bfa:354715025",
                      "bfa:355252008", "bfa:990000003"):
        record = by_id(records, record_id)
        assert record["unmatchable"] == bfa.REASON_LEAGUE_UNKNOWN, record_id
        assert record["league"] is None and record["betType"] == "other"
    lost = by_id(records, "bfa:354669433")
    assert (lost["status"], lost["stake"], lost["toWin"]) == ("lost", 200, 200)
    assert lost["placedAt"] == "2026-09-12T02:10:33Z"  # 19:10 Pacific the day before
    assert lost["closedAt"] == "2026-09-12T21:14:20Z"  # lastModification, the grading time
    free_play = by_id(records, "bfa:354693777")
    assert (free_play["status"], free_play["stake"], free_play["toWin"]) == ("won", 0, 200)


def test_nfl_teaser_is_one_record_per_leg_on_the_tickets_status(records):
    legs = [record for record in records if record["parlayId"] == "bfa:354852863"]
    assert [record["id"] for record in legs] == [f"bfa:354852863:leg{index}" for index in range(4)]
    assert all(record["isParlayLeg"] and record["legCount"] == 4 and record["status"] == "lost" for record in legs)
    assert all(record["league"] == "nfl" and record["betType"] == "spread" and record["period"] == "FG"
               for record in legs)
    assert [(record["rotation"], record["side"], record["points"], record["price"]) for record in legs] == [
        (477, "away", 8, -110), (481, "away", 8.5, -110), (456, "home", -2.5, -110), (470, "home", -1, -110)]
    assert [record["awayTeam"] or record["homeTeam"] for record in legs] == [
        "GREEN BAY PACKERS", "DENVER BRONCOS", "JACKSONVILLE JAGUARS", "DETROIT LIONS"]
    assert all(record["approx"] == ["game_date_unknown", "side_from_rotation_parity"] for record in legs)
    assert all(record["eventStart"] is None and record["closedAt"] == "2026-09-13T23:44:45Z" for record in legs)
    assert all(record["stake"] == 200 and record["toWin"] == 600 and record["raw"]["parlayPrice"] == 300
               for record in legs)
    assert [record["raw"]["legResult"] for record in legs] == [None, "PENDING", "WIN", "PUSH"]
    assert legs[1]["raw"]["legDescription"] == "[481] DENVER BRONCOS +8½-110 (B+6)"
    assert not any(record["id"] == "bfa:354852863" for record in records)


def test_college_parlay_legs_are_each_league_unknown(records):
    legs = [record for record in records if record["parlayId"] == "bfa:355159147"]
    assert [record["id"] for record in legs] == ["bfa:355159147:leg0", "bfa:355159147:leg1"]
    assert all(record["unmatchable"] == bfa.REASON_LEAGUE_UNKNOWN and record["isParlayLeg"] for record in legs)
    assert all(record["raw"]["parlayPrice"] == 255 for record in legs)  # $200 to win $509.45


def test_pending_nfl_parlay_legs_are_open_with_no_start(records):
    spread, total = [record for record in records if record["parlayId"] == "bfa:990000008"]
    assert (spread["status"], spread["betType"], spread["side"], spread["points"], spread["awayTeam"]) == \
        ("open", "spread", "away", -3, "GREEN BAY PACKERS")
    assert (total["betType"], total["side"], total["points"], total["period"]) == ("total", "under", 44.5, "FG")
    assert (total["awayTeam"], total["homeTeam"]) == ("DALLAS COWBOYS", "NEW YORK GIANTS")
    assert total["approx"] == ["game_date_unknown"]
    assert all(record["eventStart"] is None and record["closedAt"] is None for record in (spread, total))
    assert all(record["raw"]["parlayPrice"] == 260 for record in (spread, total))


def test_moneyline_push_and_cancelled_pick_em(records):
    push = by_id(records, "bfa:990000004")
    assert (push["betType"], push["points"], push["price"], push["status"], push["side"]) == \
        ("moneyline", None, 140, "push", "away")
    assert (push["eventStart"], push["closedAt"]) == ("2026-09-20T20:25:00Z", "2026-09-20T23:20:00Z")
    void = by_id(records, "bfa:990000005")
    assert (void["betType"], void["points"], void["price"], void["status"], void["side"], void["homeTeam"]) == \
        ("spread", 0, -105, "void", "home", "HOUSTON TEXANS")


def test_baseball_total_reads_past_the_pitchers_bracket(records):
    record = by_id(records, "bfa:990000007")
    assert (record["league"], record["betType"], record["side"], record["points"], record["price"]) == \
        ("mlb", "total", "over", 7.5, -120)
    assert (record["awayTeam"], record["homeTeam"], record["period"], record["status"]) == \
        ("CHI CUBS", "TB RAYS", "FG", "open")
    assert (record["eventStart"], record["eventDate"]) == ("2026-09-21T23:10:00Z", "2026-09-21")


def test_team_total_short_teaser_and_prop_fail_closed_with_a_reason(records):
    assert by_id(records, "bfa:990000006")["unmatchable"] == "team total"
    assert by_id(records, "bfa:990000010")["unmatchable"].startswith("unrecognised selection (GIANTS SUPERFECTA")


def test_short_teaser_is_still_one_unmatchable_record_per_declared_leg(records):
    # The open list stored it as :leg0..:leg3, so those are the ids that must settle — never a lone bfa:990000009.
    legs = [record for record in records if record["parlayId"] == "bfa:990000009"]
    assert [record["id"] for record in legs] == [f"bfa:990000009:leg{index}" for index in range(4)]
    assert not any(record["id"] == "bfa:990000009" for record in records)
    assert all(record["unmatchable"] == "4 TEAM TEASERS names 4 legs but carries 3" for record in legs)
    assert all(record["isParlayLeg"] and record["legCount"] == 4 and record["status"] == "open"
               and record["betType"] == "other" and record["league"] is None for record in legs)
    assert all(record["stake"] == 100 and record["toWin"] == 300 and record["raw"]["parlayPrice"] == 300
               for record in legs)
    assert [record["raw"]["legDescription"] for record in legs] == [
        "[471] DETROIT LIONS +9-110 (B+6)", "[473] LOS ANGELES RAMS +8½-110 (B+6)",
        "[475] KANSAS CITY CHIEFS +2-110 (B+6)", None]
    assert [record["raw"]["legResult"] for record in legs] == [None, "PENDING", "PENDING", None]


def test_every_fixture_record_carries_the_contract_keys(records):
    # 18 wagers: two 4-leg teasers (one short of its legs), two 2-leg parlays, 14 straight bets.
    assert len(records) == 26
    for record in records:
        assert set(record) == CONTRACT_KEYS, record["id"]
        assert record["id"].startswith("bfa:")


# ---- open bets ------------------------------------------------------------------------

@pytest.fixture
def open_records(history) -> list[dict]:
    return normalize_open_bets(history["openBets"], FETCHED_AT)


def test_open_college_bets_carry_their_league_and_start_from_the_open_list(open_records):
    total = by_id(open_records, "bfa:343243731")
    assert (total["status"], total["closedAt"], total["unmatchable"]) == ("open", None, None)
    assert (total["league"], total["betType"], total["period"], total["side"], total["points"], total["price"]) == \
        ("cbb", "total", "1H", "under", 68.5, 110)
    assert (total["rotation"], total["awayTeam"], total["homeTeam"]) == (674, "NEW MEXICO", "NEVADA")  # 1674 as written
    assert total["raw"]["rotationAsWritten"] == 1674
    # No clock assertion here: the February rows read 7 h ahead of today's clock (module docstring).
    assert total["approx"] == [] and (total["stake"], total["toWin"]) == (150, 165)
    assert total["raw"]["openBet"] is True and total["raw"]["idSport"] == "CBB"
    spread = by_id(open_records, "bfa:343243900")
    assert (spread["league"], spread["betType"], spread["period"], spread["side"], spread["points"], spread["price"]) == \
        ("cbb", "spread", "1H", "home", -2.5, -137)  # 1670 is even = home
    assert (spread["awayTeam"], spread["homeTeam"], spread["approx"]) == (None, "UCLA", ["side_from_rotation_parity"])
    assert (spread["rotation"], spread["raw"]["rotationAsWritten"]) == (670, 1670)


def test_live_open_bet_with_the_two_bracket_suffix_reads_on_the_pacific_clock(open_records):
    spread = by_id(open_records, "bfa:355820026")
    assert spread["unmatchable"] is None
    assert (spread["league"], spread["betType"], spread["period"], spread["side"], spread["points"], spread["price"]) == \
        ("cfb", "spread", "FG", "away", 35.5, -110)  # 309011 is odd = away
    assert (spread["rotation"], spread["awayTeam"], spread["homeTeam"]) == (309011, "EASTERN ILLINOIS", None)
    # 16:00 PDT = the board's 23:00 UTC kickoff; placed 12:42 PDT, first seen on the 12:43 PT poll
    assert (spread["eventStart"], spread["eventDate"]) == ("2026-09-26T23:00:00Z", "2026-09-26")
    assert spread["placedAt"] == "2026-09-26T19:42:10Z"
    assert spread["approx"] == ["side_from_rotation_parity"] and (spread["stake"], spread["toWin"]) == (220, 200)


@pytest.mark.parametrize("rotation, period, expected", [
    (1340, "1H", 340), (1306551, "1H", 306551), (1670, "1H", 670),
    (455, "1H", 455),      # three digits: nothing to strip
    (1340, "FG", 1340), (340, "FG", 340), (2340, "1H", 2340),
])
def test_game_rotation_of_strips_the_first_half_prefix(rotation, period, expected):
    assert game_rotation_of(rotation, period) == expected


def test_parse_leg_serves_the_game_rotation_and_keeps_the_written_one():
    leg = parse_leg("[1340] TOTAL u24EV \r(ARIZONA 1H vrs BYU 1H)")
    assert (leg["rotation"], leg["rotationAsWritten"], leg["period"]) == (340, 1340, "1H")
    full_game = parse_leg("[308945] ILLINOIS ST -3-110")
    assert (full_game["rotation"], full_game["rotationAsWritten"]) == (308945, 308945)


def test_open_parlay_prop_and_first_five(open_records):
    spread, total = [record for record in open_records if record["parlayId"] == "bfa:990100001"]
    assert [spread["id"], total["id"]] == ["bfa:990100001:leg0", "bfa:990100001:leg1"]
    assert all(record["status"] == "open" and record["league"] == "nfl" and record["legCount"] == 2
               and record["raw"]["parlayPrice"] == 260 for record in (spread, total))
    assert (spread["side"], spread["points"], spread["eventStart"], spread["eventDate"]) == \
        ("away", -3, "2026-09-27T20:00:00Z", "2026-09-27")  # 13:00 PDT
    assert (total["side"], total["points"], total["awayTeam"], total["homeTeam"]) == \
        ("under", 44.5, "DALLAS COWBOYS", "NEW YORK GIANTS")
    assert by_id(open_records, "bfa:990100002")["unmatchable"] == "not a game market (idSport PROP)"
    first_five = by_id(open_records, "bfa:990100003")
    assert (first_five["league"], first_five["period"], first_five["side"], first_five["points"]) == ("mlb", "F5", "over", 4.5)


def test_every_open_record_carries_the_contract_keys(open_records):
    assert len(open_records) == 7  # 6 wagers, the parlay is two
    for record in open_records:
        assert set(record) == CONTRACT_KEYS, record["id"]


def test_open_records_replace_the_historys_pending_copy():
    pending = {"id": 343243731, "type": "STRAIGHT BET", "description": "[1674] TOTAL u68\u00bd+110 \r(NEW MEXICO 1H vrs NEVADA 1H)",
               "placedDate": "2026-02-24T15:33:05", "result": "PENDING", "risk": 150.0, "win": 165.0}
    open_wager = {"idWager": 343243731, "headerDescription": "STRAIGHT BET", "riskAmount": 150.0, "winAmount": 165.0,
                  "placedDate": "2026-02-24T22:33:05.13",
                  "betDetails": [{"idSport": "CBB", "gameDateTime": "2026-02-25T03:00:00",
                                  "detailDescription": "CBB - Game <br> [1674] TOTAL u68\u00bd+110 \r(NEW MEXICO 1H vrs NEVADA 1H) [Sport:Basketball, League:NCAA]"}]}
    merged = merge_open_over_history(normalize_bfa([pending], FETCHED_AT), normalize_open_bets([open_wager], FETCHED_AT))
    [record] = merged
    assert (record["league"], record["unmatchable"], record["status"]) == ("cbb", None, "open")


# ---- a ticket the history cannot read -------------------------------------------------

TEASER_LEGS = ["[471] DETROIT LIONS +9-110 (B+6)", "[473] LOS ANGELES RAMS +8½-110 (B+6)",
               "[475] KANSAS CITY CHIEFS +2-110 (B+6)", "[479] NEW ORLEANS SAINTS +7½-110 (B+6)"]


def open_list_ticket(native_id: int, header: str) -> dict:
    """GetPlayerOpenBets' copy of a 4-leg NFL ticket: one betDetails row per leg."""
    return {"idWager": native_id, "headerDescription": header, "riskAmount": 100.0, "winAmount": 300.0,
            "placedDate": "2026-09-27T08:05:00",
            "betDetails": [{"idSport": "NFL", "gameDateTime": "2026-09-27T13:00:00", "idGame": 36810201 + index,
                            "detailDescription": f"NFL - Teasers <br> {text} [Sport:Football][League:NFL]"}
                           for index, text in enumerate(TEASER_LEGS)]}


def history_ticket(native_id: int, wager_type: str, texts: list[str]) -> dict:
    """GetPlayerHistory's pending copy: the first text is the ticket's description, the rest its picks."""
    return {"id": native_id, "type": wager_type, "description": texts[0], "placedDate": "2026-09-27T08:05:00",
            "lastModification": "2026-09-27T08:05:00", "settledDate": "2026-09-27T08:05:00",
            "result": "PENDING", "risk": 100.0, "win": 300.0, "amount": 0.0,
            "picks": [{"description": text, "result": "PENDING"} for text in texts[1:]]}


def settled(wager: dict) -> dict:
    """The same history wager once the ticket lost: graded 16:40 Pacific."""
    return {**wager, "result": "LOSE", "lastModification": "2026-09-27T16:40:00", "amount": -100.0,
            "picks": [{**pick, "result": "LOSE"} for pick in wager.get("picks") or []]}


def unreadable_history_copy(case: str, history: dict) -> tuple[dict, str]:
    """(the history's pending copy, the reason it cannot be read) for a 4-leg ticket."""
    if case == "short of its legs":  # the fixture's own: 4 TEAM TEASERS carrying 3
        [wager] = [wager for wager in history["wagers"] if wager["id"] == 990000009]
        return wager, "4 TEAM TEASERS names 4 legs but carries 3"
    if case == "a leg that does not parse":
        texts = TEASER_LEGS[:3] + ["NEW ORLEANS SAINTS +7½-110 (B+6)"]
        return (history_ticket(990000009, "4 TEAM TEASERS", texts),
                "leg not parsed: no '[rotation]' prefix (NEW ORLEANS SAINTS +7½-110)")
    return (history_ticket(990000009, "IF WIN OR PUSH (4 TEAMS)", TEASER_LEGS),
            "unrecognised multi-leg type (IF WIN OR PUSH (4 TEAMS))")


@pytest.mark.parametrize("case", ["short of its legs", "a leg that does not parse", "an unrecognised type"])
def test_a_ticket_the_history_cannot_read_still_settles_every_open_list_leg(history, tmp_path, case):
    pending_copy, reason = unreadable_history_copy(case, history)
    open_copy = open_list_ticket(990000009, pending_copy["type"])
    leg_ids = [f"bfa:990000009:leg{index}" for index in range(4)]
    store = BetsStore(tmp_path / "bets.duckdb", source_runs_retention_days=7)
    try:
        # While open, the open list's legs replace the history's pending copies one for one: no extra record.
        while_open = merge_open_over_history(normalize_bfa([pending_copy], FETCHED_AT),
                                             normalize_open_bets([open_copy], FETCHED_AT))
        assert [record["id"] for record in while_open] == leg_ids
        assert all(record["unmatchable"] is None and record["league"] == "nfl" for record in while_open)
        store.upsert_bets(while_open, datetime(2026, 9, 27, 15, 10, tzinfo=timezone.utc))
        # Graded: it has left the open list, and the history's copy is the only one.
        after = merge_open_over_history(normalize_bfa([settled(pending_copy)], FETCHED_AT), normalize_open_bets([], FETCHED_AT))
        store.upsert_bets(after, datetime(2026, 9, 28, 0, 0, tzinfo=timezone.utc))
        stored = store.load_bets(30, datetime(2026, 9, 28, 0, 5, tzinfo=timezone.utc))
    finally:
        store.close()
    assert sorted(record["id"] for record in stored) == leg_ids  # never a lone bfa:990000009
    assert all(record["status"] == "lost" and record["closedAt"] == "2026-09-27T23:40:00Z" for record in stored)
    assert all(record["unmatchable"] == reason and record["legCount"] == 4 for record in stored)


def test_a_type_declaring_no_legs_is_one_unmatchable_record_never_dropped():
    wager = settled(history_ticket(990000011, "PARLAY (0 TEAMS)", TEASER_LEGS[:2]))
    [record] = normalize_wager(wager, FETCHED_AT)
    assert (record["id"], record["status"], record["isParlayLeg"]) == ("bfa:990000011", "lost", False)
    assert record["unmatchable"] == "PARLAY (0 TEAMS) names 0 legs but carries 2"


def open_straight(id_sport: str, description: str) -> dict:
    """A one-leg GetPlayerOpenBets wager in the 2026-09-26 live shape."""
    return {"idWager": 355820026, "headerDescription": "STRAIGHT BET", "riskAmount": 220.0, "winAmount": 200.0,
            "placedDate": "2026-09-26T12:42:10",
            "betDetails": [{"idSport": id_sport, "gameDateTime": "2026-09-26T16:00:00", "detailDescription": description}]}


def test_an_open_game_leg_that_does_not_parse_is_marked_a_parse_failure():
    # A suffix BFA has not sent yet: SPORT_SUFFIX_RE strips only [Sport:…] / [League:…] brackets.
    wager = open_straight("CFB", "College Football FCS <br> [309011] EASTERN ILLINOIS +35½-110 "
                                 "[Sport:Football][League:NCAA][Region:US]")
    [record] = normalize_open_bets([wager], FETCHED_AT)
    assert record["unmatchable"].startswith("unrecognised selection (EASTERN ILLINOIS +35½-110")
    assert record["raw"]["parseFailed"] is True
    assert (record["status"], record["league"], record["betType"], record["raw"]["idSport"]) == \
        ("open", None, "other", "CFB")


@pytest.mark.parametrize("id_sport, description, reason", [
    ("CFB", "College Football <br> [1117] TOTAL o17½-115 \r(KENT STATE 1H TEAM PTS vrs OHIO STATE 1H) "
            "[Sport:Football][League:NCAA]", "team total"),
    ("PROP", "NFL Specials <br> [777045] GIANTS SUPERFECTA (SCR 1ST, 1Q, 1H & GM +7½) +1935 "
             "[Sport:Football][League:NFL]", "not a game market (idSport PROP)"),
    ("TNS", "Tennis <br> [5001] DJOKOVIC -150 [Sport:Tennis][League:ATP]", "league not supported (idSport TNS)"),
    # A code the table does not list, its league read off a team nickname: a guess, never marked.
    ("FB", "Pro Football <br> [101] PITTSBURGH STEELERS SPRING SPECIAL [Sport:Football][League:NFL]",
     "unrecognised selection"),
])
def test_deliberate_exclusions_and_guessed_leagues_are_not_marked(id_sport, description, reason):
    [record] = normalize_open_bets([open_straight(id_sport, description)], FETCHED_AT)
    assert record["unmatchable"].startswith(reason), record["unmatchable"]
    assert "parseFailed" not in record["raw"]


def test_a_history_wager_that_does_not_parse_is_not_marked():
    # The history names no league code, so its parse failure is never flagged.
    wager = {"id": 1, "type": "STRAIGHT BET", "description": "[1] SOMETHING NEW", "result": "PENDING",
             "placedDate": "2026-09-26T12:42:10", "risk": 110.0, "win": 100.0}
    [record] = normalize_wager(wager, FETCHED_AT)
    assert record["unmatchable"].startswith("unrecognised selection") and "parseFailed" not in record["raw"]


@pytest.mark.parametrize("text, expected", [
    ("[1340] TOTAL u24EV \r(ARIZONA 1H vrs BYU 1H)", ("total", "under", 24, 100, "1H", "ARIZONA", "BYU")),
    ("[1117] TOTAL O35-110 \r(KENT STATE 1H VRS OHIO STATE 1H)", ("total", "over", 35, -110, "1H", "KENT STATE", "OHIO STATE")),
    ("[308915] TOTAL o47-110 \r(HOLY CROSS vrs MIAMI OHIO)", ("total", "over", 47, -110, "FG", "HOLY CROSS", "MIAMI OHIO")),
    ("[308945] ILLINOIS ST -3-110", ("spread", "away", -3, -110, "FG", "ILLINOIS ST", None)),
    ("[308959] WESTERN ILLINOIS +42-110", ("spread", "away", 42, -110, "FG", "WESTERN ILLINOIS", None)),
    ("[462] PITTSBURGH STEELERS -½-110 (B+6)", ("spread", "home", -0.5, -110, "FG", None, "PITTSBURGH STEELERS")),
    ("[1306564] LAMAR 1H -2-110", ("spread", "home", -2, -110, "1H", None, "LAMAR")),
    ("[1306551] GRAMBLING 1H +168", ("moneyline", "away", None, 168, "1H", "GRAMBLING", None)),
    ("[464] HOUSTON TEXANS pk-105", ("spread", "home", 0, -105, "FG", None, "HOUSTON TEXANS")),
    ("[968] TB RAYS +120 (CHI CUBS vrs TB RAYS)", ("moneyline", "home", None, 120, "FG", "CHI CUBS", "TB RAYS")),
    ("CBB - Alternative Lines <br> [1670] UCLA 1H -2\u00bd-137 [Sport:Basketball, League:NCAA]",
     ("spread", "home", -2.5, -137, "1H", None, "UCLA")),
    ("College Football FCS <br> [309011] EASTERN ILLINOIS +35\u00bd-110 [Sport:Football][League:NCAA]",
     ("spread", "away", 35.5, -110, "FG", "EASTERN ILLINOIS", None)),
    ("NFL - Game <br> [453] TOTAL u44\u00bd-110 \r(DALLAS COWBOYS vrs NEW YORK GIANTS) [League:NFL][Sport:Football]",
     ("total", "under", 44.5, -110, "FG", "DALLAS COWBOYS", "NEW YORK GIANTS")),
])
def test_parse_leg_grammar(text, expected):
    leg = parse_leg(text)
    assert isinstance(leg, dict), leg
    assert (leg["betType"], leg["side"], leg["points"], leg["price"], leg["period"],
            leg["awayTeam"], leg["homeTeam"]) == expected


@pytest.mark.parametrize("text, reason", [
    ("TOTAL o47-110 (HOLY CROSS vrs MIAMI OHIO)", "no '[rotation]' prefix"),
    ("[1] TOTAL o47-110", "total names no teams"),
    ("[1] TOTAL o47-110 (HOLY CROSS at MIAMI OHIO)", "no 'AWAY vrs HOME' bracket"),
    ("[1] TOTAL o47-110 (HOLY CROSS 1H vrs MIAMI OHIO 2H)", "periods differ (1H vs 2H)"),
    ("[1] TOTAL o47-110 (HOLY CROSS 3H vrs MIAMI OHIO 3H)", "unknown period (3H)"),
    ("[1] TB RAYS +120 (CHI CUBS vrs NY YANKEES)", "team TB RAYS is not in its bracket"),
])
def test_parse_leg_fails_closed_with_the_reason(text, reason):
    result = parse_leg(text)
    assert isinstance(result, str) and result.startswith(reason), result


def test_a_player_prop_reads_as_a_moneyline_on_nobody_and_is_league_unknown():
    wager = {"id": 2, "type": "STRAIGHT BET", "description": "[1] Justin Fields 200+ passing yards -120",
             "placedDate": "2026-09-20T07:00:00", "result": "PENDING", "risk": 120.0, "win": 100.0}
    [record] = normalize_wager(wager, FETCHED_AT)
    assert record["unmatchable"] == bfa.REASON_LEAGUE_UNKNOWN and record["betType"] == "other"


def test_wager_without_an_id_fails_the_poll_loudly():
    wager = {"type": "STRAIGHT BET", "description": "[455] TOTAL o24½-110 (JACKSONVILLE JAGUARS vrs CINCINNATI BENGALS)",
             "placedDate": "2026-09-27T07:12:03.5", "result": "PENDING", "risk": 220.0, "win": 200.0}
    with pytest.raises(RuntimeError, match="carries no id"):
        normalize_bfa([wager], FETCHED_AT)


def test_a_settled_bet_without_a_grading_time_closes_at_its_placed_time():
    wager = {"id": 1, "type": "STRAIGHT BET", "description": "[455] TOTAL o24½-110 (JACKSONVILLE JAGUARS vrs CINCINNATI BENGALS)",
             "placedDate": "2026-09-20T07:00:00", "result": "WIN", "risk": 110.0, "win": 100.0}
    [record] = normalize_wager(wager, FETCHED_AT)
    assert (record["status"], record["closedAt"], record["placedAt"]) == ("won", "2026-09-20T14:00:00Z", "2026-09-20T14:00:00Z")
    assert record["eventStart"] is None and record["approx"] == ["game_date_unknown"]


# ---- network half ---------------------------------------------------------------------

def jwt_with(player_id: str) -> str:
    payload = base64.urlsafe_b64encode(json.dumps({"player_id": player_id}).encode()).rstrip(b"=").decode()
    return f"header.{payload}.signature"


class FakeResponse:
    def __init__(self, status_code: int, body: object = None, text: str = "", headers: dict | None = None):
        self.status_code = status_code
        self._body = body
        self.text = text
        self.headers = headers or {}

    def json(self) -> object:
        return self._body


class FakeSession:
    """Answers the Keycloak login, the token grants and the paged history; records every call."""

    LOGIN_ACTION = "https://auth.bfagaming.com/realms/players_realm/login-actions/authenticate?session_code=abc&execution=def"

    def __init__(self, wagers: list[dict], login_ok: bool = True, history_status: int = 200, expires_in: int = 300,
                 open_wagers: list[dict] | None = None):
        self.wagers = wagers
        self.open_wagers = open_wagers or []
        self.login_ok = login_ok
        self.history_status = history_status
        self.expires_in = expires_in
        self.refresh_ok = True
        self.logins = 0
        self.tokens_minted = 0
        self.refresh_tokens_seen: list[str] = []
        self.history_calls: list[dict] = []
        self.passwords_seen: list[str] = []

    def _tokens(self) -> dict:
        self.tokens_minted += 1
        return {"access_token": jwt_with("777"), "expires_in": self.expires_in,
                "refresh_token": f"refresh-{self.tokens_minted}"}

    def get(self, url: str, params=None, headers=None, timeout=None) -> FakeResponse:
        if url == bfa.AUTH_URL:
            assert params["code_challenge_method"] == "S256" and params["client_id"] == "bfagaming"
            return FakeResponse(200, text=f'<form action="{self.LOGIN_ACTION.replace("&", "&amp;")}" method="post">')
        if url == bfa.OPEN_BETS_URL:
            assert headers["Authorization"] == f"Bearer {jwt_with('777')}" and params == {"playerId": "777"}
            return FakeResponse(200, self.open_wagers)
        if url == bfa.HISTORY_URL:
            assert headers["Authorization"] == f"Bearer {jwt_with('777')}" and params["playerId"] == "777"
            self.history_calls.append(params)
            if self.history_status != 200:
                return FakeResponse(self.history_status, "unavailable")
            start = params["page"] * params["recordsByPage"]
            page = self.wagers[start:start + params["recordsByPage"]]
            # totalRecords counts transactions too, as live (18 against 16 wagers).
            return FakeResponse(200, {"totalRecords": len(self.wagers) + 2, "wagers": page})
        raise AssertionError(f"unexpected GET {url}")

    def post(self, url: str, data=None, headers=None, allow_redirects=True, timeout=None) -> FakeResponse:
        if url == self.LOGIN_ACTION:
            assert allow_redirects is False
            self.logins += 1
            self.passwords_seen.append(data["password"])
            if not self.login_ok:
                return FakeResponse(200, text="<span>Invalid username or password.</span>")
            return FakeResponse(302, headers={"Location": f"https://bfagaming.com/?session_state=s&code=code-{self.logins}"})
        if url == bfa.TOKEN_URL and data["grant_type"] == "authorization_code":
            assert data["code"] == f"code-{self.logins}" and data["code_verifier"]
            return FakeResponse(200, self._tokens())
        if url == bfa.TOKEN_URL and data["grant_type"] == "refresh_token":
            self.refresh_tokens_seen.append(data["refresh_token"])
            if not self.refresh_ok:
                return FakeResponse(400, {"error": "invalid_grant"})
            return FakeResponse(200, self._tokens())
        raise AssertionError(f"unexpected POST {url} {data}")


def make_source(session: FakeSession, clock) -> BFASource:
    return BFASource(username="user", password="secret", history_days=31, poll_sec=300,
                     session_factory=lambda: session, clock=clock)


def test_fetch_logs_in_once_reads_the_open_list_and_the_paged_history(history, monkeypatch):
    monkeypatch.setattr(bfa, "RECORDS_PER_PAGE", 5)
    session = FakeSession(history["wagers"], open_wagers=history["openBets"])
    source = make_source(session, clock=lambda: CLOCK)
    records = source.fetch()
    assert len(records) == 33 and records[0]["id"] == "bfa:354669433"  # 26 history + 7 open
    assert sum(1 for record in records if record["status"] == "open") == 16  # 9 pending in the weeks + 7 open
    assert session.logins == 1 and session.passwords_seen == ["secret"]
    # 18 wagers in pages of 5, then the empty page that ends a total padded by transactions.
    assert [call["page"] for call in session.history_calls] == [0, 1, 2, 3, 4]
    assert (session.history_calls[0]["startDate"], session.history_calls[0]["endDate"]) == ("2026-08-22", "2026-09-23")


def test_access_token_is_reused_then_refreshed_within_60s_of_expiry_and_a_refused_refresh_logs_in_again(history):
    session = FakeSession(history["wagers"], expires_in=300)
    now = [CLOCK]
    source = make_source(session, clock=lambda: now[0])
    source.fetch()
    now[0] += 200  # 100 s of life left: no refresh, no login
    source.fetch()
    assert (session.logins, session.refresh_tokens_seen) == (1, [])
    now[0] += 50  # 50 s left, inside the 60 s margin: refresh with the in-memory token
    source.fetch()
    assert (session.logins, session.refresh_tokens_seen) == (1, ["refresh-1"])
    session.refresh_ok = False
    now[0] += 300  # expired; Keycloak refuses the refresh: a new password login, nothing to fix on disk
    source.fetch()
    assert (session.logins, session.refresh_tokens_seen) == (2, ["refresh-1", "refresh-2"])


def test_invalid_credentials_raise(history):
    session = FakeSession(history["wagers"], login_ok=False)
    with pytest.raises(RuntimeError, match=r"BFA login failed: HTTP 200 \(invalid credentials\)"):
        make_source(session, clock=lambda: CLOCK).fetch()


def test_failed_history_page_raises_instead_of_returning_a_partial_list(history):
    session = FakeSession(history["wagers"], history_status=503)
    with pytest.raises(RuntimeError, match="history page 0 failed: HTTP 503"):
        make_source(session, clock=lambda: CLOCK).fetch()


def test_missing_credentials_name_the_env_keys(history):
    source = BFASource(username="", password="", history_days=31, poll_sec=300,
                       session_factory=lambda: FakeSession(history["wagers"]))
    with pytest.raises(RuntimeError, match="BFA_USERNAME and BFA_PASSWORD"):
        source.fetch()


def test_source_if_configured_needs_both_credentials(monkeypatch):
    monkeypatch.setattr(bfa.config, "BFA_USERNAME", "user")
    monkeypatch.setattr(bfa.config, "BFA_PASSWORD", None)
    assert bfa.source_if_configured() is None
    monkeypatch.setattr(bfa.config, "BFA_PASSWORD", "secret")
    source = bfa.source_if_configured()
    assert isinstance(source, BFASource) and (source.name, source.poll_sec) == ("bfa", 300.0)
