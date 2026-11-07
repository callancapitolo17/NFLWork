"""Unabated's stored fair ladders (issue #264): ladders.duckdb::ladder_readings
(its own file apart from bets.duckdb) and POST /ladders.json over a real HTTP
server — the first opener stays, the newest pregame close wins, and a bad row
refuses the whole request."""
import json
import threading
import urllib.error
import urllib.request
from datetime import datetime, timezone
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.ladder_store import LadderStore
from unabated_ticket.bets_service.store import BetsStore

SPREAD_CLOSE = {
    "reading": "close", "eventId": 125807, "leagueId": 1, "periodTypeId": 1, "betTypeId": 2,
    "eventName": "Bears Chicago CHI @ Panthers Carolina CAR", "awayTeamId": 6, "homeTeamId": 5,
    "eventStart": "2026-09-13T17:00:00.000Z", "observedAt": "2026-09-13T16:00:00.000Z",
    "mainPoints0": -3, "mainFair0": -103, "mainPoints1": 3, "mainFair1": 103,
    "rungs": [[0, -3.5, 110, 9, 9], [0, -3, -103, 12, 11], [1, 3, 103, 12, 12]],
}
SPREAD_OPENER = {**SPREAD_CLOSE, "reading": "opener"}
MONEYLINE_CLOSE = {**SPREAD_CLOSE, "betTypeId": 1, "mainPoints0": None, "mainFair0": -156, "mainPoints1": None,
                   "mainFair1": 156, "rungs": [[0, None, -156, 9, 9], [1, None, 156, 8, 8]]}
RECEIVED_AT = datetime(2026, 9, 13, 16, 0, 5, tzinfo=timezone.utc)


@pytest.fixture
def stores(tmp_path):
    bets_store = BetsStore(tmp_path / "bets.duckdb", 7)
    ladder_store = LadderStore(tmp_path / "ladders.duckdb")
    yield bets_store, ladder_store
    bets_store.close()
    ladder_store.close()


def serve(handler) -> tuple[ThreadingHTTPServer, str]:
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, f"http://127.0.0.1:{server.server_address[1]}"


def post(url: str, payload: object) -> tuple[int, dict]:
    request = urllib.request.Request(url, data=json.dumps(payload).encode(), method="POST",
                                     headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


def stored(ladder_store: LadderStore) -> list[tuple]:
    return ladder_store._con.execute(
        "SELECT reading, bet_type_id, CAST(observed_at AS VARCHAR), main_fair_0, rung_count, rungs_json "
        "FROM ladder_readings ORDER BY reading, bet_type_id").fetchall()


def save(ladder_store: LadderStore, rows: list[dict]) -> int:
    validated = service.validate_ladder_rows({"rows": rows})
    assert not isinstance(validated, str), validated
    return ladder_store.save_readings(validated, RECEIVED_AT)


def test_the_first_opener_stays_and_the_newest_close_wins(stores):
    _, ladder_store = stores
    assert save(ladder_store, [SPREAD_OPENER, SPREAD_CLOSE, MONEYLINE_CLOSE]) == 3
    later = {"observedAt": "2026-09-13T16:50:00.000Z", "mainFair0": -110}
    assert save(ladder_store, [{**SPREAD_OPENER, **later}, {**SPREAD_CLOSE, **later}]) == 1
    # An older close (a retry that landed late) never replaces a newer one.
    assert save(ladder_store, [SPREAD_CLOSE]) == 0
    rows = stored(ladder_store)
    assert [row[:5] for row in rows] == [
        ("close", 1, "2026-09-13 16:00:00+00", -156, 2),
        ("close", 2, "2026-09-13 16:50:00+00", -110, 3),
        ("opener", 2, "2026-09-13 16:00:00+00", -103, 3),
    ]
    assert json.loads(rows[0][5]) == [[0, None, -156, 9, 9], [1, None, 156, 8, 8]]


@pytest.mark.parametrize("bad_row, message", [
    ({**SPREAD_CLOSE, "reading": "kickoff"}, "rows[0].reading must be one of"),
    ({**SPREAD_CLOSE, "leagueId": 5}, "rows[0].leagueId must be one of"),
    ({**SPREAD_CLOSE, "periodTypeId": 2}, "rows[0].periodTypeId must be one of"),
    ({**SPREAD_CLOSE, "eventId": "125807"}, "rows[0].eventId must be an integer"),
    ({**SPREAD_CLOSE, "observedAt": "2026-09-13T17:00:00.000Z"}, "rows[0].observedAt must be before eventStart"),
    ({**SPREAD_CLOSE, "mainFair0": 50}, "rows[0].mainFair0 must be a whole American price or null"),
    ({**SPREAD_CLOSE, "rungs": []}, "rows[0].rungs must be a list of 1 to"),
    ({**SPREAD_CLOSE, "rungs": [[2, -3.5, 110, 1, 1]]}, "rows[0].rungs[0][0] (side) must be 0 or 1"),
    ({**SPREAD_CLOSE, "rungs": [[0, None, 110, 1, 1]]}, "rows[0].rungs[0][1] (points) must be a number"),
    ({**MONEYLINE_CLOSE, "rungs": [[0, 0.5, 110, 1, 1]]}, "rows[0].rungs[0][1] (points) must be null on a moneyline"),
    ({**SPREAD_CLOSE, "rungs": [[0, -3.5, 110.5, 1, 1]]}, "rows[0].rungs[0][2] (fair) must be a whole American price"),
    ({**SPREAD_CLOSE, "rungs": [[0, -3.5, 110, 1, 2]]}, "rows[0].rungs[0][3:5] (books, agreeing)"),
])
def test_a_bad_row_refuses_the_whole_request(bad_row, message):
    assert service.validate_ladder_rows({"rows": [bad_row, SPREAD_CLOSE]}).startswith(message)


def test_post_saves_the_rows(stores):
    bets_store, ladder_store = stores
    server, base = serve(service.make_handler(bets_store, 0.0, ["kalshi"], ladder_store=ladder_store))
    try:
        assert post(f"{base}/ladders.json", {"rows": [SPREAD_OPENER, SPREAD_CLOSE]}) == (200, {"ok": True, "saved": 2})
        status, body = post(f"{base}/ladders.json", {"rows": [{**SPREAD_CLOSE, "reading": "nope"}]})
        assert status == 400 and body["error"].startswith("rows[0].reading must be one of")
    finally:
        server.shutdown()
        server.server_close()
    assert len(stored(ladder_store)) == 2


def test_without_a_ladder_store_the_route_is_503(stores):
    bets_store, _ = stores
    server, base = serve(service.make_handler(bets_store, 0.0, ["kalshi"]))
    try:
        assert post(f"{base}/ladders.json", {"rows": [SPREAD_CLOSE]}) == (
            503, {"error": "the ladder store is not open on this service"})
    finally:
        server.shutdown()
        server.server_close()
