"""The automated taker's decision log (automated_taker_plan.md, T1):
taker.duckdb::taker_decisions (append, its own file apart from bets.duckdb)
and POST /taker_decisions.json over a real HTTP server, with the body's
validation refusing the whole request on its first bad row."""
import json
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.store import BetsStore
from unabated_ticket.bets_service.taker_store import TakerStore

WOULD_BET = {
    "rowKey": "420942308:ms105:si1:tid717:alt-35.5", "league": "CFB", "eventId": 123742,
    "eventStart": "2026-09-19T19:30:00.000Z", "decidedAt": "2026-09-15T19:00:00.000Z", "reason": "bet",
    "detail": None, "ticker": "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36", "side": "yes", "unabatedPrice": -143,
    "liveAllIn": 0.587157, "liveEdge": 0.05, "contracts": 12, "stake": 7.5,
}


@pytest.fixture
def stores(tmp_path):
    bets_store = BetsStore(tmp_path / "bets.duckdb", 7)
    taker_store = TakerStore(tmp_path / "taker.duckdb")
    yield bets_store, taker_store
    bets_store.close()
    taker_store.close()


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


def test_post_appends_every_row_with_its_whole_decision(stores, tmp_path):
    bets_store, taker_store = stores
    server, base = serve(service.make_handler(bets_store, 0.0, ["kalshi"], taker_store=taker_store))
    try:
        assert post(f"{base}/taker_decisions.json", {"rows": [WOULD_BET, {**WOULD_BET, "reason": "price_moved"}]}) == (
            200, {"ok": True, "logged": 2})
        assert post(f"{base}/taker_decisions.json", {"rows": [WOULD_BET]})[0] == 200
    finally:
        server.shutdown()
        server.server_close()
    assert taker_store.count_decisions() == 3
    row = taker_store._con.execute(
        "SELECT reason, ticker, contract_side, unabated_price, contracts, event_id, decision FROM taker_decisions "
        "ORDER BY received_at LIMIT 1").fetchone()
    assert row[:6] == ("bet", WOULD_BET["ticker"], "yes", -143, 12, 123742)
    assert json.loads(row[6])["stake"] == 7.5


@pytest.mark.parametrize("bad_row, message", [
    ({**WOULD_BET, "reason": "placed"}, "rows[0].reason must be one of"),
    ({**WOULD_BET, "decidedAt": "yesterday"}, "rows[0].decidedAt must be an ISO time"),
    ({**WOULD_BET, "rowKey": ""}, "rows[0].rowKey must be a non-empty string"),
    ({**WOULD_BET, "contracts": 1.5}, "rows[0].contracts must be an integer or null"),
    ({**WOULD_BET, "liveEdge": "big"}, "rows[0].liveEdge must be a number or null"),
])
def test_a_bad_row_refuses_the_whole_request(bad_row, message):
    assert service.validate_taker_decisions({"rows": [bad_row, WOULD_BET]}).startswith(message)


def test_without_a_taker_store_the_route_is_503(stores):
    bets_store, _ = stores
    server, base = serve(service.make_handler(bets_store, 0.0, ["kalshi"]))
    try:
        status, body = post(f"{base}/taker_decisions.json", {"rows": [WOULD_BET]})
    finally:
        server.shutdown()
        server.server_close()
    assert status == 503 and "taker log" in body["error"]


def test_an_unopenable_taker_log_is_none_not_a_crash(tmp_path):
    (tmp_path / "taker.duckdb").mkdir()  # a directory where the file should be
    assert service.open_taker_store(tmp_path / "taker.duckdb") is None
